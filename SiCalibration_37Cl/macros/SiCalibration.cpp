#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "MulginCalibration.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <map>
#include <vector>

struct RunFitInfo {
  Double_t mu = 0;
  Double_t mu_err = 0;
  Double_t sigma = 0;
  Double_t sigma_err = 0;
  Double_t reduced_chi2 = 0;
};

void SiCalibration() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");

  TFile *fit_file = IO::OpenForReading(project_root + "/root_files/" +
                                       Constants::FIT_RESULTS_FILE);
  if (!fit_file || fit_file->IsZombie()) {
    std::cerr << "ERROR: cannot open " << Constants::FIT_RESULTS_FILE
              << " (run SiFits first)" << std::endl;
    return;
  }

  TTree *fit_tree = static_cast<TTree *>(fit_file->Get("FitResults"));
  if (!fit_tree) {
    std::cerr << "ERROR: FitResults tree not found" << std::endl;
    fit_file->Close();
    delete fit_file;
    return;
  }

  Int_t b_run;
  Double_t b_mu, b_mu_err, b_sigma, b_sigma_err, b_chi2;
  fit_tree->SetBranchAddress("RunNumber", &b_run);
  fit_tree->SetBranchAddress("Mu", &b_mu);
  fit_tree->SetBranchAddress("MuError", &b_mu_err);
  fit_tree->SetBranchAddress("Sigma", &b_sigma);
  fit_tree->SetBranchAddress("SigmaError", &b_sigma_err);
  fit_tree->SetBranchAddress("ReducedChi2", &b_chi2);

  std::map<Int_t, RunFitInfo> fits;
  Long64_t n = fit_tree->GetEntries();
  for (Long64_t i = 0; i < n; i++) {
    fit_tree->GetEntry(i);
    RunFitInfo info;
    info.mu = b_mu;
    info.mu_err = b_mu_err;
    info.sigma = b_sigma;
    info.sigma_err = b_sigma_err;
    info.reduced_chi2 = b_chi2;
    fits[b_run] = info;
  }

  const Bool_t INCLUDE_PPAC = Constants::INCLUDE_PPAC;
  const TString hypo_tag = INCLUDE_PPAC ? "PPAC_IN" : "PPAC_OUT";
  TString config_label =
      Form("Mulgin calibration, PPAC %s", INCLUDE_PPAC ? "IN" : "OUT");
  std::cout << "Running calibration: " << config_label << std::endl;

  // ── Load Mulgin fit parameters for the selected PPAC hypothesis ──────────
  TString mulgin_path =
      project_root + "/root_files/" + Constants::MULGIN_CALIBRATION_FILE;
  TFile *mulgin_file = IO::OpenForReading(mulgin_path);
  if (!mulgin_file || mulgin_file->IsZombie()) {
    std::cerr << "ERROR: cannot open " << mulgin_path
              << " (run MulginCalibration first)" << std::endl;
    fit_file->Close();
    delete fit_file;
    return;
  }
  TString mulgin_tree_name = TString("MulginFit_") + hypo_tag;
  TTree *mulgin_tree = static_cast<TTree *>(mulgin_file->Get(mulgin_tree_name));
  if (!mulgin_tree || mulgin_tree->GetEntries() < 1) {
    std::cerr << "ERROR: tree " << mulgin_tree_name << " missing or empty in "
              << mulgin_path << std::endl;
    mulgin_file->Close();
    delete mulgin_file;
    fit_file->Close();
    delete fit_file;
    return;
  }
  Double_t B_A = 0, B_A_err = 0, C_A = 0, C_A_err = 0;
  Double_t B_B = 0, B_B_err = 0, C_B = 0, C_B_err = 0;
  Double_t gamma = 0, gamma_err = 0;
  mulgin_tree->SetBranchAddress("B_A", &B_A);
  mulgin_tree->SetBranchAddress("B_A_err", &B_A_err);
  mulgin_tree->SetBranchAddress("C_A", &C_A);
  mulgin_tree->SetBranchAddress("C_A_err", &C_A_err);
  mulgin_tree->SetBranchAddress("B_B", &B_B);
  mulgin_tree->SetBranchAddress("B_B_err", &B_B_err);
  mulgin_tree->SetBranchAddress("C_B", &C_B);
  mulgin_tree->SetBranchAddress("C_B_err", &C_B_err);
  mulgin_tree->SetBranchAddress("Gamma", &gamma);
  mulgin_tree->SetBranchAddress("Gamma_err", &gamma_err);
  mulgin_tree->GetEntry(0);

  std::printf("  Mulgin params (%s):\n", hypo_tag.Data());
  std::printf("    Group A: B = %.4f +/- %.4f   C = %.4f +/- %.4f\n", B_A,
              B_A_err, C_A, C_A_err);
  std::printf("    Group B: B = %.4f +/- %.4f   C = %.4f +/- %.4f\n", B_B,
              B_B_err, C_B, C_B_err);
  std::printf("    gamma  = %.6f +/- %.6f\n", gamma, gamma_err);

  // ── Invert: ADC channel → energy (MeV), per gain group ───────────────────
  // mu_err is propagated through dE/d(adc); B/C/gamma errors are not folded
  // in here (they are correlated and dominated by the calibration constraint
  // set, not by individual mu measurements).
  std::map<Int_t, Double_t> energy_mev, energy_err_mev;
  for (size_t i = 0; i < Constants::RUN_NUMBERS.size(); i++) {
    Int_t run = Constants::RUN_NUMBERS[i];
    if (fits.count(run) == 0)
      continue;
    Bool_t group_A = MulginIsGroupA(run);
    Double_t B = group_A ? B_A : B_B;
    Double_t C = group_A ? C_A : C_B;
    Double_t adc = fits[run].mu;
    energy_mev[run] = MulginInverseEnergy(adc, B, C, gamma);
    energy_err_mev[run] =
        MulginInverseEnergyError(adc, B, C, gamma, fits[run].mu_err);
  }

  mulgin_file->Close();
  delete mulgin_file;

  // dE relative to beam energy (only meaningful for pressure-scan runs).
  std::map<Int_t, Double_t> dE_mev, dE_err_mev;
  for (size_t i = 0; i < Constants::RUN_NUMBERS.size(); i++) {
    Int_t run = Constants::RUN_NUMBERS[i];
    if (Constants::GAS_PRESSURE_TORR.count(run) == 0)
      continue;
    if (energy_mev.count(run) == 0)
      continue;
    Double_t dE = Constants::BEAM_ENERGY_MEV - energy_mev[run];
    Double_t dE_err =
        TMath::Sqrt(TMath::Power(Constants::BEAM_ENERGY_ERR_MEV, 2) +
                    TMath::Power(energy_err_mev[run], 2));
    dE_mev[run] = dE;
    dE_err_mev[run] = dE_err;
  }

  // Annotated per-run plots: histogram + fit overlay + E/dE/mu/sigma legend.
  Bool_t logy = kTRUE;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(logy);

  for (size_t i = 0; i < Constants::RUN_NUMBERS.size(); i++) {
    Int_t run = Constants::RUN_NUMBERS[i];
    if (fits.count(run) == 0)
      continue;

    TString hist_name = Form("Run_%d", run);
    TH1F *hist = static_cast<TH1F *>(fit_file->Get(hist_name));
    if (!hist) {
      std::cerr << "WARNING: histogram " << hist_name << " missing"
                << std::endl;
      continue;
    }

    TF1 *fit_func =
        static_cast<TF1 *>(fit_file->Get(Form("FitFunc_Run_%d", run)));

    canvas->cd();
    PlottingUtils::ConfigureAndDrawHistogram(hist, kBlue, Form("Run %d", run));
    hist->SetStats(0);

    if (fit_func) {
      fit_func->SetLineColor(kRed);
      fit_func->SetLineWidth(2);
      fit_func->Draw("SAME");
    }

    TLegend *legend = run <= 28
                          ? PlottingUtils::AddLegend(0.2, 0.475, 0.65, 0.88)
                          : PlottingUtils::AddLegend(0.6, 0.9, 0.65, 0.88);
    if (fit_func)
      legend->AddEntry(fit_func, "Fit", "l");
    legend->AddEntry(
        (TObject *)0,
        Form("E = %.2f #pm %.2f MeV", energy_mev[run], energy_err_mev[run]),
        "");
    if (dE_mev.find(run) != dE_mev.end()) {
      legend->AddEntry(
          (TObject *)0,
          Form("#DeltaE = %.2f #pm %.2f MeV", dE_mev[run], dE_err_mev[run]),
          "");
    }
    legend->AddEntry(
        (TObject *)0,
        Form("#mu = %.1f #pm %.1f ADC", fits[run].mu, fits[run].mu_err), "");
    legend->AddEntry((TObject *)0,
                     Form("#sigma = %.1f #pm %.1f ADC", fits[run].sigma,
                          fits[run].sigma_err),
                     "");
    legend->SetMargin(0.05);
    legend->Draw();

    canvas->Update();
    PlottingUtils::SaveFigure(canvas, Form("Run_%d", run), "",
                              PlotSaveOptions::kLOG);
  }

  // dE vs pressure (Torr).
  std::vector<Double_t> pressures, dE_vals, dE_errs, p_errs;
  for (size_t i = 0; i < Constants::RUN_NUMBERS.size(); i++) {
    Int_t run = Constants::RUN_NUMBERS[i];
    if (Constants::GAS_PRESSURE_TORR.count(run) == 0)
      continue;
    if (dE_mev.count(run) == 0)
      continue;
    pressures.push_back(Constants::GAS_PRESSURE_TORR.at(run));
    dE_vals.push_back(dE_mev[run]);
    dE_errs.push_back(dE_err_mev[run]);
    p_errs.push_back(0);
  }

  TGraphErrors *dE_vs_pressure =
      new TGraphErrors(pressures.size(), pressures.data(), dE_vals.data(),
                       p_errs.data(), dE_errs.data());
  dE_vs_pressure->SetTitle(
      "Energy Loss vs. Gas Pressure;Pressure (Torr);#DeltaE (MeV)");
  dE_vs_pressure->SetMarkerStyle(20);
  dE_vs_pressure->SetMarkerSize(1.2);
  dE_vs_pressure->SetMarkerColor(kBlue);
  dE_vs_pressure->SetLineColor(kBlue);

  TCanvas *dE_canvas = PlottingUtils::GetConfiguredCanvas();
  dE_vs_pressure->Draw("AP");
  dE_vs_pressure->GetXaxis()->SetTitleSize(0.05);
  dE_vs_pressure->GetYaxis()->SetTitleSize(0.05);
  dE_vs_pressure->GetXaxis()->SetLabelSize(0.045);
  dE_vs_pressure->GetYaxis()->SetLabelSize(0.045);
  dE_vs_pressure->GetXaxis()->SetTitleOffset(1.2);
  dE_vs_pressure->GetYaxis()->SetTitleOffset(1.3);
  dE_canvas->Update();
  PlottingUtils::SaveFigure(dE_canvas, "dE_vs_Pressure", "",
                            PlotSaveOptions::kLINEAR);

  // Resolution (% FWHM) vs run.
  std::vector<Double_t> run_x, res_y, res_yerr, run_xerr;
  for (size_t i = 0; i < Constants::RUN_NUMBERS.size(); i++) {
    Int_t run = Constants::RUN_NUMBERS[i];
    if (fits.count(run) == 0)
      continue;
    Double_t fwhm = 2.355 * fits[run].sigma;
    Double_t fwhm_err = 2.355 * fits[run].sigma_err;
    Double_t res = (fwhm / fits[run].mu) * 100.0;
    Double_t res_err =
        100.0 * TMath::Sqrt(TMath::Power(fwhm_err / fits[run].mu, 2) +
                            TMath::Power(fwhm * fits[run].mu_err /
                                             TMath::Power(fits[run].mu, 2),
                                         2));
    run_x.push_back(run);
    res_y.push_back(res);
    res_yerr.push_back(res_err);
    run_xerr.push_back(0);
  }

  TGraphErrors *resolution_vs_run =
      new TGraphErrors(run_x.size(), run_x.data(), res_y.data(),
                       run_xerr.data(), res_yerr.data());
  resolution_vs_run->SetTitle(
      "Energy Resolution vs. Run Number;Run Number;Resolution (% FWHM)");
  resolution_vs_run->SetMarkerStyle(20);
  resolution_vs_run->SetMarkerSize(1.2);
  resolution_vs_run->SetMarkerColor(kRed);
  resolution_vs_run->SetLineColor(kRed);

  TCanvas *res_canvas = PlottingUtils::GetConfiguredCanvas();
  resolution_vs_run->Draw("AP");
  resolution_vs_run->GetXaxis()->SetTitleSize(0.05);
  resolution_vs_run->GetYaxis()->SetTitleSize(0.05);
  resolution_vs_run->GetXaxis()->SetLabelSize(0.045);
  resolution_vs_run->GetYaxis()->SetLabelSize(0.045);
  resolution_vs_run->GetXaxis()->SetTitleOffset(1.2);
  resolution_vs_run->GetYaxis()->SetTitleOffset(1.3);
  resolution_vs_run->GetXaxis()->SetNdivisions(516);
  res_canvas->Update();
  PlottingUtils::SaveFigure(res_canvas, "Resolution_vs_Run", "",
                            PlotSaveOptions::kLINEAR);

  // Persist derived results.
  TFile *out_file = IO::OpenForWriting(Constants::CALIBRATION_RESULTS_FILE);

  TTree *out_tree =
      new TTree("CalibrationResults", "Silicon Detector Calibration Results");
  Int_t run_number;
  Double_t gas_pressure_torr;
  Double_t delta_E_MeV;
  Double_t delta_E_err_MeV;
  Double_t fwhm_percent;
  Double_t fwhm_percent_err;
  Double_t centroid_ch;
  Double_t centroid_err_ch;
  Double_t sigma_ch;
  Double_t sigma_err_ch;
  out_tree->Branch("RunNumber", &run_number, "RunNumber/I");
  out_tree->Branch("GasPressure", &gas_pressure_torr, "GasPressure/D");
  out_tree->Branch("DeltaE", &delta_E_MeV, "DeltaE/D");
  out_tree->Branch("DeltaE_Error", &delta_E_err_MeV, "DeltaE_Error/D");
  out_tree->Branch("FWHM_Percent", &fwhm_percent, "FWHM_Percent/D");
  out_tree->Branch("FWHM_Percent_Error", &fwhm_percent_err,
                   "FWHM_Percent_Error/D");
  out_tree->Branch("Centroid", &centroid_ch, "Centroid/D");
  out_tree->Branch("Centroid_Error", &centroid_err_ch, "Centroid_Error/D");
  out_tree->Branch("Sigma", &sigma_ch, "Sigma/D");
  out_tree->Branch("Sigma_Error", &sigma_err_ch, "Sigma_Error/D");

  for (size_t i = 0; i < Constants::RUN_NUMBERS.size(); i++) {
    Int_t run = Constants::RUN_NUMBERS[i];
    if (Constants::GAS_PRESSURE_TORR.count(run) == 0)
      continue;
    if (fits.count(run) == 0)
      continue;
    run_number = run;
    gas_pressure_torr = Constants::GAS_PRESSURE_TORR.at(run);
    delta_E_MeV = dE_mev[run];
    delta_E_err_MeV = dE_err_mev[run];
    Double_t fwhm = 2.355 * fits[run].sigma;
    Double_t fwhm_err = 2.355 * fits[run].sigma_err;
    fwhm_percent = (fwhm / fits[run].mu) * 100.0;
    fwhm_percent_err =
        100.0 * TMath::Sqrt(TMath::Power(fwhm_err / fits[run].mu, 2) +
                            TMath::Power(fwhm * fits[run].mu_err /
                                             TMath::Power(fits[run].mu, 2),
                                         2));
    centroid_ch = fits[run].mu;
    centroid_err_ch = fits[run].mu_err;
    sigma_ch = fits[run].sigma;
    sigma_err_ch = fits[run].sigma_err;
    out_tree->Fill();
  }

  out_file->cd();
  out_tree->Write("CalibrationResults", TObject::kOverwrite);
  dE_vs_pressure->Write("DeltaE_vs_Pressure", TObject::kOverwrite);
  resolution_vs_run->Write("Resolution_vs_Run", TObject::kOverwrite);
  out_file->Close();

  delete out_file;
  delete res_canvas;
  delete dE_canvas;
  delete canvas;
  fit_file->Close();
  delete fit_file;
}
