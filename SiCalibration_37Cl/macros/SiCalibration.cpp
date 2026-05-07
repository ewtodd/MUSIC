#include "FittingUtils.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <vector>

void SiCalibration() {

  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                TString(gSystem->pwd()) + "/plots",
                                TString(gSystem->pwd()) + "/root_files");

  std::vector<std::string> filepaths = {
      "/home/e-work/LabData/MUSIC/37Cl/SiCalibration/run_21/ROOT/DataR_run_21.root",
      "/home/e-work/LabData/MUSIC/37Cl/SiCalibration/run_22/ROOT/DataR_run_22.root",
      "/home/e-work/LabData/MUSIC/37Cl/SiCalibration/run_23/ROOT/DataR_run_23.root",
      "/home/e-work/LabData/MUSIC/37Cl/SiCalibration/run_24/ROOT/DataR_run_24.root",
      "/home/e-work/LabData/MUSIC/37Cl/SiCalibration/run_25/ROOT/DataR_run_25.root",
      "/home/e-work/LabData/MUSIC/37Cl/SiCalibration/run_26/ROOT/DataR_run_26.root",
      "/home/e-work/LabData/MUSIC/37Cl/SiCalibration/run_27/ROOT/DataR_run_27.root",
      "/home/e-work/LabData/MUSIC/37Cl/SiCalibration/run_28/ROOT/DataR_run_28.root",
      "/home/e-work/LabData/MUSIC/37Cl/SiCalibration/run_29/ROOT/DataR_run_29.root",
      "/home/e-work/LabData/MUSIC/37Cl/SiCalibration/run_30/ROOT/DataR_run_30.root",
      "/home/e-work/LabData/MUSIC/37Cl/SiCalibration/run_31/ROOT/DataR_run_31.root",
      "/home/e-work/LabData/MUSIC/37Cl/SiCalibration/run_32/ROOT/DataR_run_32.root",
      "/home/e-work/LabData/MUSIC/37Cl/SiCalibration/run_33/ROOT/DataR_run_33.root",
      "/home/e-work/LabData/MUSIC/37Cl/SiCalibration/run_34/ROOT/DataR_run_34.root",
      "/home/e-work/LabData/MUSIC/37Cl/SiCalibration/run_35/ROOT/DataR_run_35.root",
      "/home/e-work/LabData/MUSIC/37Cl/SiCalibration/run_36/ROOT/DataR_run_36.root"};

  std::map<int, double> mu_guesses;
  mu_guesses[21] = 15300;
  mu_guesses[22] = 10000;
  mu_guesses[23] = 9500;
  mu_guesses[24] = 13200;
  mu_guesses[25] = 10150;
  mu_guesses[26] = 9200;
  mu_guesses[27] = 8500;
  mu_guesses[28] = 7800;
  mu_guesses[29] = 6800;
  mu_guesses[30] = 6000;
  mu_guesses[31] = 5100;
  mu_guesses[32] = 4200;
  mu_guesses[33] = 3400;
  mu_guesses[34] = 2520;
  mu_guesses[35] = 1760;
  mu_guesses[36] = 1500;

  std::vector<Double_t> mu_values;
  std::vector<Double_t> sigma_values;
  std::vector<Double_t> mu_errors;
  std::vector<Double_t> sigma_errors;
  std::vector<Double_t> energy_values;
  std::vector<Double_t> energy_errors;

  // First pass: get mu values for all runs via FittingUtils::FitSinglePeak.
  // Plots are written by FittingUtils to plots/fits/<peak>_Run_<run>.{png,pdf}.
  Double_t mu_21 = 0, mu_22 = 0;
  Double_t mu_21_err = 0, mu_22_err = 0;

  for (int run = 21; run <= 36; run++) {
    int idx = run - 21;

    TFile *file = IO::OpenForReading(filepaths[idx].c_str());
    if (!file || file->IsZombie()) {
      std::cerr << "Error opening file for run " << run << std::endl;
      continue;
    }

    TTree *tree = static_cast<TTree *>(file->Get("Data_R"));
    if (!tree) {
      std::cerr << "Error getting tree for run " << run << std::endl;
      file->Close();
      continue;
    }

    UShort_t channel, energy;
    tree->SetBranchAddress("Channel", &channel);
    tree->SetBranchAddress("Energy", &energy);

    TH1F *hist =
        new TH1F(Form("Run_%d", run), Form("Run %d ; ADC Channel; Counts", run),
                 600, 0, 16384);

    Int_t number_samples = tree->GetEntries();
    for (int i = 0; i < number_samples; i++) {
      tree->GetEntry(i);
      if (channel == 11) {
        hist->Fill(energy);
      }
    }

    Double_t mu_guess = mu_guesses[run];
    Double_t fit_min = TMath::Max(0.0, mu_guess - 2000);
    Double_t fit_max = TMath::Min(16384.0, mu_guess + 2000);

    FittingUtils fitter(hist, fit_min, fit_max,
                        /*use_flat_background=*/kTRUE);

    FitResult fit_result = fitter.FitSinglePeak(Form("Run_%d", run), "alpha");

    if (!fit_result.peaks.empty() && fit_result.peaks[0].mu > 0) {
      mu_values.push_back(fit_result.peaks[0].mu);
      sigma_values.push_back(fit_result.peaks[0].sigma);
      mu_errors.push_back(fit_result.peaks[0].mu_error);
      sigma_errors.push_back(fit_result.peaks[0].sigma_error);
    } else {
      std::cerr << "WARNING: FitSinglePeak failed for run " << run << std::endl;
      mu_values.push_back(mu_guess);
      sigma_values.push_back(200);
      mu_errors.push_back(0);
      sigma_errors.push_back(0);
    }

    if (run == 21) {
      mu_21 = mu_values.back();
      mu_21_err = mu_errors.back();
    }
    if (run == 22) {
      mu_22 = mu_values.back();
      mu_22_err = mu_errors.back();
    }

    delete hist;
    delete file;
  }

  Double_t cal_x[2] = {mu_21, mu_22};
  Double_t cal_y[2] = {91.87, 71.35};

  TGraph *cal_graph = new TGraph(2, cal_x, cal_y);
  cal_graph->SetTitle(";ADC Channel;TOF Energy (MeV)");
  cal_graph->SetMarkerStyle(20);
  cal_graph->SetMarkerSize(1.5);
  cal_graph->SetMarkerColor(kBlue);

  TF1 *cal_fit = new TF1("cal_fit", "[0] + [1]*x", 0, 16384);
  cal_fit->SetParameter(0, 0);
  cal_graph->Fit("cal_fit", "QRL+");
  cal_fit->SetRange(0, 16384);

  Double_t p0 = cal_fit->GetParameter(0);
  Double_t p1 = cal_fit->GetParameter(1);
  Double_t p0_err = cal_fit->GetParError(0);
  Double_t p1_err = cal_fit->GetParError(1);

  TCanvas *cal_canvas = PlottingUtils::GetConfiguredCanvas();

  cal_graph->Draw("AP");
  cal_graph->GetXaxis()->SetLimits(-1, 16384);
  cal_fit->SetRange(0, 16384);
  cal_canvas->Modified();
  cal_canvas->Update();
  cal_graph->SetMinimum(0);
  cal_graph->SetMaximum(100);
  cal_graph->GetXaxis()->SetTitleSize(0.05);
  cal_graph->GetYaxis()->SetTitleSize(0.05);
  cal_graph->GetXaxis()->SetLabelSize(0.045);
  cal_graph->GetYaxis()->SetLabelSize(0.045);

  TLegend *cal_legend = PlottingUtils::AddLegend(0.15, 0.60, 0.65, 0.9);
  cal_legend->AddEntry(cal_graph, "Calibration points", "p");
  cal_legend->AddEntry(cal_fit, Form("Fit: E = p_{0} + p_{1}x + p_{2}x^2"),
                       "l");
  cal_legend->AddEntry((TObject *)0, Form("p_{0} = %.2f #pm %.3f", p0, p0_err),
                       "");
  cal_legend->AddEntry((TObject *)0, Form("p_{1} = %.5f #pm %.6f", p1, p1_err),
                       "");
  cal_legend->Draw();

  cal_canvas->Update();
  PlottingUtils::SaveFigure(cal_canvas, "Energy_Calibration", "",
                            PlotSaveOptions::kLINEAR);
  for (size_t i = 0; i < mu_values.size(); i++) {
    Double_t adc = mu_values[i];
    Double_t energy_MeV = p0 + p1 * adc;

    Double_t energy_err_MeV =
        TMath::Sqrt(TMath::Power(p0_err, 2) + TMath::Power(adc * p1_err, 2) +
                    TMath::Power(p1 * mu_errors[i], 2));
    energy_values.push_back(energy_MeV);
    energy_errors.push_back(energy_err_MeV);
  }

  std::map<int, double> gas_pressure;
  gas_pressure[25] = 0;
  gas_pressure[26] = 20;
  gas_pressure[27] = 40;
  gas_pressure[28] = 60;
  gas_pressure[29] = 80;
  gas_pressure[30] = 100;
  gas_pressure[31] = 120;
  gas_pressure[32] = 140;
  gas_pressure[33] = 160;
  gas_pressure[34] = 180;
  gas_pressure[35] = 200;
  gas_pressure[36] = 220;

  Double_t beam_energy = 92.00;    // MeV
  Double_t beam_energy_err = 0.05; // MeV
  std::vector<Double_t> dE_values;
  std::vector<Double_t> dE_errors;

  for (size_t i = 0; i < energy_values.size(); i++) {
    Double_t dE = beam_energy - energy_values[i];
    Double_t dE_err = TMath::Sqrt(TMath::Power(beam_energy_err, 2) +
                                  TMath::Power(energy_errors[i], 2));
    dE_values.push_back(dE);
    dE_errors.push_back(dE_err);
  }

  Bool_t logy = kTRUE;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(logy);

  for (int run = 21; run <= 36; run++) {
    int idx = run - 21;

    TFile *file = IO::OpenForReading(filepaths[idx].c_str());
    if (!file || file->IsZombie())
      continue;

    TTree *tree = static_cast<TTree *>(file->Get("Data_R"));
    if (!tree) {
      file->Close();
      continue;
    }

    UShort_t channel, energy;
    tree->SetBranchAddress("Channel", &channel);
    tree->SetBranchAddress("Energy", &energy);

    TH1F *hist =
        new TH1F(Form("Run_%d_annot", run),
                 Form("Run %d ; ADC Channel; Counts", run), 600, 0, 16384);

    Int_t number_samples = tree->GetEntries();
    for (int i = 0; i < number_samples; i++) {
      tree->GetEntry(i);
      if (channel == 11) {
        hist->Fill(energy);
      }
    }

    canvas->cd();
    PlottingUtils::ConfigureAndDrawHistogram(hist, kBlue, Form("Run %d", run));
    hist->SetStats(0);

    TLegend *legend = run <= 28
                          ? PlottingUtils::AddLegend(0.2, 0.475, 0.70, 0.88)
                          : PlottingUtils::AddLegend(0.625, 0.9, 0.70, 0.88);
    legend->AddEntry(
        (TObject *)0,
        Form("E = %.2f #pm %.2f MeV", energy_values[idx], energy_errors[idx]),
        "");
    if (run >= 25) {
      legend->AddEntry(
          (TObject *)0,
          Form("#DeltaE = %.2f #pm %.2f MeV", dE_values[idx], dE_errors[idx]),
          "");
    }
    legend->AddEntry(
        (TObject *)0,
        Form("#mu = %.1f #pm %.1f ADC", mu_values[idx], mu_errors[idx]), "");
    legend->AddEntry((TObject *)0,
                     Form("#sigma = %.1f #pm %.1f ADC", sigma_values[idx],
                          sigma_errors[idx]),
                     "");
    legend->SetMargin(0.05);
    legend->Draw();

    canvas->Update();
    PlottingUtils::SaveFigure(canvas, Form("Run_%d", run), "",
                              PlotSaveOptions::kLOG);
    delete hist;
    delete file;
  }
  std::vector<Double_t> pressure_vec;
  std::vector<Double_t> dE_for_plot;
  std::vector<Double_t> dE_err_for_plot;
  std::vector<Double_t> pressure_err_vec;

  for (int run = 25; run <= 36; run++) {
    int idx = run - 21;
    pressure_vec.push_back(gas_pressure[run]);
    dE_for_plot.push_back(dE_values[idx]);
    dE_err_for_plot.push_back(dE_errors[idx]);
    pressure_err_vec.push_back(0);
  }

  TGraphErrors *dE_vs_pressure =
      new TGraphErrors(pressure_vec.size(), &pressure_vec[0], &dE_for_plot[0],
                       &pressure_err_vec[0], &dE_err_for_plot[0]);

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

  std::vector<Double_t> run_numbers;
  std::vector<Double_t> resolution_percent;
  std::vector<Double_t> resolution_errors;
  std::vector<Double_t> run_err_vec;

  for (int run = 21; run <= 36; run++) {
    int idx = run - 21;

    Double_t fwhm = 2.355 * sigma_values[idx];
    Double_t fwhm_err = 2.355 * sigma_errors[idx];

    Double_t resolution = (fwhm / mu_values[idx]) * 100.0;

    Double_t res_err =
        100.0 * TMath::Sqrt(TMath::Power(fwhm_err / mu_values[idx], 2) +
                            TMath::Power(fwhm * mu_errors[idx] /
                                             TMath::Power(mu_values[idx], 2),
                                         2));

    run_numbers.push_back(run);
    resolution_percent.push_back(resolution);
    resolution_errors.push_back(res_err);
    run_err_vec.push_back(0);
  }

  TGraphErrors *resolution_vs_run = new TGraphErrors(
      run_numbers.size(), &run_numbers[0], &resolution_percent[0],
      &run_err_vec[0], &resolution_errors[0]);

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

  TFile *outFile = IO::OpenForWriting("SiCalibration_Results.root");

  TTree *results =
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

  results->Branch("RunNumber", &run_number, "RunNumber/I");
  results->Branch("GasPressure", &gas_pressure_torr, "GasPressure/D");
  results->Branch("DeltaE", &delta_E_MeV, "DeltaE/D");
  results->Branch("DeltaE_Error", &delta_E_err_MeV, "DeltaE_Error/D");
  results->Branch("FWHM_Percent", &fwhm_percent, "FWHM_Percent/D");
  results->Branch("FWHM_Percent_Error", &fwhm_percent_err,
                  "FWHM_Percent_Error/D");
  results->Branch("Centroid", &centroid_ch, "Centroid/D");
  results->Branch("Centroid_Error", &centroid_err_ch, "Centroid_Error/D");
  results->Branch("Sigma", &sigma_ch, "Sigma/D");
  results->Branch("Sigma_Error", &sigma_err_ch, "Sigma_Error/D");

  for (Int_t run = 25; run <= 36; run++) {
    Int_t idx = run - 21;
    run_number = run;
    gas_pressure_torr = gas_pressure[run];
    delta_E_MeV = dE_values[idx];
    delta_E_err_MeV = dE_errors[idx];
    Double_t fwhm = 2.355 * sigma_values[idx];
    Double_t fwhm_err = 2.355 * sigma_errors[idx];
    fwhm_percent = (fwhm / mu_values[idx]) * 100.0;
    fwhm_percent_err =
        100.0 * TMath::Sqrt(TMath::Power(fwhm_err / mu_values[idx], 2) +
                            TMath::Power(fwhm * mu_errors[idx] /
                                             TMath::Power(mu_values[idx], 2),
                                         2));

    centroid_ch = mu_values[idx];
    centroid_err_ch = mu_errors[idx];
    sigma_ch = sigma_values[idx];
    sigma_err_ch = sigma_errors[idx];

    results->Fill();
  }

  results->Write();

  cal_graph->Write("EnergyCalibration");
  dE_vs_pressure->Write("DeltaE_vs_Pressure");
  resolution_vs_run->Write("Resolution_vs_Run");

  outFile->Close();

  delete outFile;
  delete res_canvas;
  delete dE_canvas;
  delete canvas;
  delete cal_canvas;
}
