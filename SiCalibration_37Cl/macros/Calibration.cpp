#include "Calibration.hpp"
#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "StoppingPowerLISE.hpp"
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <cstdio>
#include <iostream>

struct RunConfig {
  Int_t run;
  Double_t tof_MeV;
  Bool_t front_Ti; // 0.9 mg/cm^2 Ti entrance window
  Bool_t exit_Ti;  // 1.3 mg/cm^2 Ti exit window
  Double_t adc_centroid;
  Double_t adc_sigma;
  Int_t group; // 0 = A, 1 = B
};

const Int_t kNRuns = 5;
const Int_t kNPar = 4; // B_A, C_A, B_B, C_B

static Double_t g_E_true[kNRuns];
static Double_t g_P_meas[kNRuns];
static Double_t g_P_sig[kNRuns];
static Int_t g_group[kNRuns];

// par[0] = B_A, par[1] = C_A, par[2] = B_B, par[3] = C_B.
static void MulginFcn(Int_t & /*npar*/, Double_t * /*gin*/, Double_t &f,
                      Double_t *par, Int_t /*iflag*/) {
  Double_t B_A = par[0];
  Double_t C_A = par[1];
  Double_t B_B = par[2];
  Double_t C_B = par[3];

  Double_t chi2 = 0.0;
  for (Int_t i = 0; i < kNRuns; i++) {
    Double_t B = (g_group[i] == 0) ? B_A : B_B;
    Double_t C = (g_group[i] == 0) ? C_A : C_B;
    Double_t P_pred = MulginPredictChannel(g_E_true[i], B, C);
    Double_t r = (g_P_meas[i] - P_pred) / g_P_sig[i];
    chi2 += r * r;
  }
  f = chi2;
}

// Propagate the TOF beam energy through the active materials for one run,
// returning the energy at the Si front face (MeV).
static Double_t ComputeETrue(const RunConfig &cfg, Bool_t use_ppac,
                             const TString &mylar_file,
                             const TString &ti_file) {
  const Int_t model_index = 1; // Ziegler
  const Double_t A_Cl37 = 37.0;
  const Double_t ppac_layer_micron = 1.5;
  const Int_t ppac_n_layers = 4;
  const Double_t ti_front_mg_cm2 = 0.9;
  const Double_t ti_exit_mg_cm2 = 1.3;

  Double_t E = cfg.tof_MeV;

  if (use_ppac) {
    for (Int_t k = 0; k < ppac_n_layers; k++) {
      Double_t dedx =
          GetStoppingPowerFromLISE(mylar_file, model_index, E / A_Cl37);
      E -= dedx * ppac_layer_micron;
    }
  }

  if (cfg.front_Ti) {
    Double_t dedx = GetStoppingPowerFromLISE(ti_file, model_index, E / A_Cl37);
    E -= dedx * ti_front_mg_cm2;
  }

  if (cfg.exit_Ti) {
    Double_t dedx = GetStoppingPowerFromLISE(ti_file, model_index, E / A_Cl37);
    E -= dedx * ti_exit_mg_cm2;
  }

  return E;
}

static void RunOneHypothesis(Bool_t use_ppac, const RunConfig *runs,
                             const TString &mylar_file, const TString &ti_file,
                             TFile *out_file) {
  const TString hypo_tag = use_ppac ? "PPAC_IN" : "PPAC_OUT";
  const TString hypo_pretty = use_ppac ? "PPAC IN" : "PPAC OUT";

  for (Int_t i = 0; i < kNRuns; i++) {
    g_E_true[i] = ComputeETrue(runs[i], use_ppac, mylar_file, ti_file);
    g_P_meas[i] = runs[i].adc_centroid;
    g_P_sig[i] = runs[i].adc_sigma;
    g_group[i] = runs[i].group;
  }

  Double_t dEA = g_E_true[0] - g_E_true[1];
  Double_t slopeA = (dEA != 0.0) ? (g_P_meas[0] - g_P_meas[1]) / dEA : 160.0;
  Double_t intA = g_P_meas[0] - slopeA * g_E_true[0];
  Double_t dEB = g_E_true[3] - g_E_true[4];
  Double_t slopeB = (dEB != 0.0) ? (g_P_meas[3] - g_P_meas[4]) / dEB : 160.0;
  Double_t intB = g_P_meas[3] - slopeB * g_E_true[3];

  TMinuit *minuit = new TMinuit(kNPar);
  minuit->SetFCN(MulginFcn);

  Int_t ierflg = 0;
  Double_t arglist[10];

  arglist[0] = 1.0;
  minuit->mnexcm("SET ERR", arglist, 1, ierflg);
  arglist[0] = -1.0;
  minuit->mnexcm("SET PRI", arglist, 1, ierflg);
  arglist[0] = 0.0;
  minuit->mnexcm("SET NOW", arglist, 0, ierflg);

  minuit->mnparm(0, "B_A", slopeA, 1.0, 0.0, 0.0, ierflg);
  minuit->mnparm(1, "C_A", intA, 50.0, 0.0, 0.0, ierflg);
  minuit->mnparm(2, "B_B", slopeB, 1.0, 0.0, 0.0, ierflg);
  minuit->mnparm(3, "C_B", intB, 50.0, 0.0, 0.0, ierflg);

  arglist[0] = 10000.0;
  arglist[1] = 0.01;
  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  minuit->mnexcm("HESSE", arglist, 1, ierflg);

  Double_t par[kNPar];
  Double_t err[kNPar];
  for (Int_t i = 0; i < kNPar; i++) {
    minuit->GetParameter(i, par[i], err[i]);
  }

  Double_t fmin = 0.0;
  Double_t fedm = 0.0;
  Double_t errdef = 0.0;
  Int_t npari = 0;
  Int_t nparx = 0;
  Int_t istat = 0;
  minuit->mnstat(fmin, fedm, errdef, npari, nparx, istat);
  Int_t ndf = kNRuns - kNPar;

  const char *pname[kNPar] = {"B_A", "C_A", "B_B", "C_B"};
  std::printf("Mulgin fit, %s:\n", hypo_pretty.Data());
  for (Int_t i = 0; i < kNPar; i++) {
    std::printf("  %-5s = %14.6f  +/-  %12.6f\n", pname[i], par[i], err[i]);
  }
  if (ndf > 0) {
    std::printf("  chi2/ndf = %.4f / %d = %.4f\n", fmin, ndf, fmin / ndf);
  } else {
    std::printf("  chi2 = %.4g (ndf = 0)\n", fmin);
  }

  std::printf("  %-3s  %-11s  %-10s  %-10s  %-7s\n", "Run", "E_true(MeV)",
              "P_meas", "P_pred", "resid/s");
  for (Int_t i = 0; i < kNRuns; i++) {
    Double_t B = (g_group[i] == 0) ? par[0] : par[2];
    Double_t C = (g_group[i] == 0) ? par[1] : par[3];
    Double_t P_pred = MulginPredictChannel(g_E_true[i], B, C);
    Double_t resid_sigma = (g_P_meas[i] - P_pred) / g_P_sig[i];
    std::printf("  %3d  %11.4f  %10.2f  %10.2f  %+7.3f\n", runs[i].run,
                g_E_true[i], g_P_meas[i], P_pred, resid_sigma);
  }

  Double_t E_25_on_A = MulginInverseEnergy(g_P_meas[4], par[0], par[1]);
  std::printf("  run 25 on group-A gain: inferred E = %.4f MeV (true E = "
              "%.4f)\n",
              E_25_on_A, g_E_true[4]);

  Double_t emat[kNPar * kNPar];
  minuit->mnemat(emat, kNPar);
  Double_t den_A = TMath::Sqrt(emat[0 * kNPar + 0] * emat[1 * kNPar + 1]);
  Double_t den_B = TMath::Sqrt(emat[2 * kNPar + 2] * emat[3 * kNPar + 3]);
  Double_t rho_BA_CA = (den_A > 0.0) ? emat[0 * kNPar + 1] / den_A : 0.0;
  Double_t rho_BB_CB = (den_B > 0.0) ? emat[2 * kNPar + 3] / den_B : 0.0;
  std::printf("  corr(B_A, C_A) = %+.4f   corr(B_B, C_B) = %+.4f\n", rho_BA_CA,
              rho_BB_CB);
  std::printf("  C_A - C_B = %+.4f  (err %.4f)\n\n", par[1] - par[3],
              TMath::Sqrt(err[1] * err[1] + err[3] * err[3]));

  out_file->cd();
  const TString tree_name = TString("MulginFit_") + hypo_tag;
  TTree *fit_tree = new TTree(tree_name, tree_name);
  Bool_t include_ppac = use_ppac;
  Double_t B_A = par[0], B_A_err = err[0];
  Double_t C_A = par[1], C_A_err = err[1];
  Double_t B_B = par[2], B_B_err = err[2];
  Double_t C_B = par[3], C_B_err = err[3];
  Double_t chi2 = fmin;
  Int_t ndf_out = ndf;
  fit_tree->Branch("IncludePPAC", &include_ppac, "IncludePPAC/O");
  fit_tree->Branch("B_A", &B_A, "B_A/D");
  fit_tree->Branch("B_A_err", &B_A_err, "B_A_err/D");
  fit_tree->Branch("C_A", &C_A, "C_A/D");
  fit_tree->Branch("C_A_err", &C_A_err, "C_A_err/D");
  fit_tree->Branch("B_B", &B_B, "B_B/D");
  fit_tree->Branch("B_B_err", &B_B_err, "B_B_err/D");
  fit_tree->Branch("C_B", &C_B, "C_B/D");
  fit_tree->Branch("C_B_err", &C_B_err, "C_B_err/D");
  fit_tree->Branch("Chi2", &chi2, "Chi2/D");
  fit_tree->Branch("Ndf", &ndf_out, "Ndf/I");
  fit_tree->Fill();
  fit_tree->Write(tree_name, TObject::kOverwrite);

  Int_t nA = 0;
  Int_t nB = 0;
  for (Int_t i = 0; i < kNRuns; i++) {
    if (g_group[i] == 0)
      nA++;
    else
      nB++;
  }
  TGraphErrors *gA = new TGraphErrors(nA);
  TGraphErrors *gB = new TGraphErrors(nB);
  Int_t iA = 0;
  Int_t iB = 0;
  Double_t xmin = 1.0e9;
  Double_t xmax = -1.0e9;
  Double_t ymax_data = 0.0;
  for (Int_t i = 0; i < kNRuns; i++) {
    if (g_E_true[i] < xmin)
      xmin = g_E_true[i];
    if (g_E_true[i] > xmax)
      xmax = g_E_true[i];
    if (g_P_meas[i] > ymax_data)
      ymax_data = g_P_meas[i];
    if (g_group[i] == 0) {
      gA->SetPoint(iA, g_E_true[i], g_P_meas[i]);
      gA->SetPointError(iA, 0.0, g_P_sig[i]);
      iA++;
    } else {
      gB->SetPoint(iB, g_E_true[i], g_P_meas[i]);
      gB->SetPointError(iB, 0.0, g_P_sig[i]);
      iB++;
    }
  }
  Double_t xlo = xmin - 5.0;
  if (xlo < 0.0)
    xlo = 0.0;
  Double_t xhi = xmax + 5.0;

  TString title = ";Deposited Energy [MeV];ADC Channel";
  PlottingUtils::ConfigureGraph(gA, kBlue + 1, title);
  PlottingUtils::ConfigureGraph(gB, kRed + 1, title);
  gB->SetMarkerStyle(21);

  TCanvas *cal_canvas = PlottingUtils::GetConfiguredCanvas();
  cal_canvas->SetRightMargin(0.28);
  gA->GetXaxis()->SetLimits(xlo, xhi);
  gA->SetMinimum(0.0);
  gA->SetMaximum(1.15 * ymax_data);
  gA->GetYaxis()->SetTitleOffset(1.3);
  gA->Draw("AP");
  gB->Draw("P SAME");

  TString formula = "[0]*(x - 0.55*x/(1.0 + 0.37568*x)) + [1]";
  TF1 *fA = new TF1(Form("fA_%s", hypo_tag.Data()), formula, xlo, xhi);
  fA->SetParameter(0, par[0]);
  fA->SetParameter(1, par[1]);
  fA->SetLineColor(kBlue + 1);
  fA->SetLineWidth(PlottingUtils::GetLineWidth());
  fA->SetNpx(500);
  fA->Draw("SAME");

  TF1 *fB = new TF1(Form("fB_%s", hypo_tag.Data()), formula, xlo, xhi);
  fB->SetParameter(0, par[2]);
  fB->SetParameter(1, par[3]);
  fB->SetLineColor(kRed + 1);
  fB->SetLineWidth(PlottingUtils::GetLineWidth());
  fB->SetNpx(500);
  fB->Draw("SAME");

  TLegend *leg = PlottingUtils::AddLegend(0.735, 0.985, 0.12, 0.88);
  leg->SetTextFont(42);
  leg->SetTextSize(0.025);
  leg->SetMargin(0.1);
  leg->AddEntry((TObject *)0, Form("%s", hypo_pretty.Data()), "");
  leg->AddEntry(gA, "Group A (runs 21, 22)", "lep");
  leg->AddEntry(gB, "Group B (runs 23, 24, 25)", "lep");
  for (Int_t i = 0; i < kNPar; i++) {
    leg->AddEntry((TObject *)0,
                  Form("%s = %.4f #pm %.4f", pname[i], par[i], err[i]), "");
  }
  Double_t gain_X = (par[2] != 0.0) ? par[0] / par[2] : 0.0;
  leg->AddEntry((TObject *)0, Form("Gain X = %.3f ", gain_X), "");
  leg->Draw();

  cal_canvas->Update();
  PlottingUtils::SaveFigure(cal_canvas, "MulginCalibration_" + hypo_tag, "",
                            PlotSaveOptions::kLINEAR);

  delete minuit;
}

void Calibration() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");

  const TString lise_dir = project_root + "/macros/lise";
  const TString mylar_file = lise_dir + "/37Cl_in_Mylar.lise";
  const TString ti_file = lise_dir + "/37Cl_in_Ti_NOT_MICRON.lise";

  RunConfig runs[kNRuns];
  // run, TOF_MeV, front_Ti, exit_Ti, adc_centroid, adc_sigma, group(0=A,1=B)
  runs[0] = {21, 91.87, kFALSE, kFALSE, 15332.8, 48.8, 0};
  runs[1] = {22, 71.35, kFALSE, kFALSE, 10061.9, 42.4, 0};
  runs[2] = {23, 71.35, kTRUE, kFALSE, 9604.6, 120.4, 1};
  runs[3] = {24, 92.00, kTRUE, kFALSE, 13344.3, 119.4, 1};
  runs[4] = {25, 92.00, kTRUE, kTRUE, 10189.3, 140.8, 1};

  TFile *out_file = IO::OpenForWriting(Constants::MULGIN_CALIBRATION_FILE);

  RunOneHypothesis(kFALSE, runs, mylar_file, ti_file, out_file);
  RunOneHypothesis(kTRUE, runs, mylar_file, ti_file, out_file);

  out_file->Close();
  delete out_file;
}
