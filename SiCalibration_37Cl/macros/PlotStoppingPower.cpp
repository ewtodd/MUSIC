#include "InitUtils.hpp"
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

void PlotStoppingPower() {
  InitUtils::SetROOTPreferences();
  TString input_filepath = "root_files/SiCalibration_Results.root";
  TFile *inFile = new TFile(input_filepath, "READ");
  TTree *calibration_results =
      static_cast<TTree *>(inFile->Get("CalibrationResults"));

  TTree *results = (TTree *)inFile->Get("StoppingPower");

  Int_t run_number;
  Double_t gas_pressure_torr;
  Double_t delta_E_experimental;
  Double_t deltaE_model0, deltaE_model1, deltaE_model2, deltaE_model3;
  Double_t deltaE_model4, deltaE_model5, deltaE_model6;

  calibration_results->SetBranchAddress("RunNumber", &run_number);
  calibration_results->SetBranchAddress("GasPressure", &gas_pressure_torr);
  results->SetBranchAddress("DeltaE", &delta_E_experimental);
  results->SetBranchAddress("DeltaE_Model0_Hubert", &deltaE_model0);
  results->SetBranchAddress("DeltaE_Model1_Ziegler", &deltaE_model1);
  results->SetBranchAddress("DeltaE_Model2_ATIMA12LS", &deltaE_model2);
  results->SetBranchAddress("DeltaE_Model3_ATIMA12NoLS", &deltaE_model3);
  results->SetBranchAddress("DeltaE_Model4_ATIMA14Weick", &deltaE_model4);
  results->SetBranchAddress("DeltaE_Model5_Electrical", &deltaE_model5);
  results->SetBranchAddress("DeltaE_Model6_Nuclear", &deltaE_model6);

  Int_t colors[7] = {kRed,  kBlue,   kGreen + 2,  kMagenta,
                     kCyan, kOrange, kViolet + 10};
  Int_t markerStyles[7] = {24, 25, 26, 27, 28, 30, 32};

  std::vector<Double_t> pressure_vec, exp_dE_vec, model0_vec, model1_vec,
      model2_vec;
  std::vector<Double_t> model3_vec, model4_vec, model5_vec, model6_vec;

  Int_t nentries = results->GetEntries();
  Int_t debug_nentries = calibration_results->GetEntries();
  if (nentries != debug_nentries) {
    std::cout << "Number of entries in two trees differs. Quitting..."
              << std::endl;
    return;
  };
  for (Int_t i = 0; i < nentries; ++i) {

    calibration_results->GetEntry(i);
    results->GetEntry(i);

    pressure_vec.push_back(gas_pressure_torr);
    exp_dE_vec.push_back(delta_E_experimental);
    model0_vec.push_back(deltaE_model0);
    model1_vec.push_back(deltaE_model1);
    model2_vec.push_back(deltaE_model2);
    model3_vec.push_back(deltaE_model3);
    model4_vec.push_back(deltaE_model4);
    model5_vec.push_back(deltaE_model5);
    model6_vec.push_back(deltaE_model6);
  }

  std::vector<Double_t> pressure_err(pressure_vec.size(), 0.0);
  std::vector<Double_t> dE_err(pressure_vec.size(), 0.5);

  // Create error vectors for model data (10% of the value)
  std::vector<Double_t> model0_err, model1_err, model2_err, model3_err;
  std::vector<Double_t> model4_err, model5_err, model6_err;

  for (size_t i = 0; i < pressure_vec.size(); ++i) {
    model0_err.push_back(model0_vec[i] * 0.10);
    model1_err.push_back(model1_vec[i] * 0.10);
    model2_err.push_back(model2_vec[i] * 0.10);
    model3_err.push_back(model3_vec[i] * 0.10);
    model4_err.push_back(model4_vec[i] * 0.10);
    model5_err.push_back(model5_vec[i] * 0.10);
    model6_err.push_back(model6_vec[i] * 0.10);
  }

  TMultiGraph *mg = new TMultiGraph();

  TGraphErrors *grExp =
      new TGraphErrors(pressure_vec.size(), &pressure_vec[0], &exp_dE_vec[0],
                       &pressure_err[0], &dE_err[0]);
  grExp->SetName("Experimental");
  grExp->SetMarkerStyle(20);
  grExp->SetMarkerSize(1);
  grExp->SetMarkerColor(kBlack);
  grExp->SetLineColor(kBlack);
  grExp->SetLineWidth(2);
  mg->Add(grExp);

  // Model 0 - now with 10% error bars
  TGraphErrors *gr0 =
      new TGraphErrors(pressure_vec.size(), &pressure_vec[0], &model0_vec[0],
                       &pressure_err[0], &model0_err[0]);
  gr0->SetName("Hubert (He-base)");
  gr0->SetMarkerStyle(markerStyles[0]);
  gr0->SetMarkerSize(1.3);
  gr0->SetMarkerColor(colors[0]);
  gr0->SetLineColor(colors[0]);
  mg->Add(gr0);

  // Model 1 - now with 10% error bars
  TGraphErrors *gr1 =
      new TGraphErrors(pressure_vec.size(), &pressure_vec[0], &model1_vec[0],
                       &pressure_err[0], &model1_err[0]);
  gr1->SetName("Ziegler (H-base)");
  gr1->SetLineWidth(2);
  gr1->SetMarkerStyle(markerStyles[1]);
  gr1->SetMarkerSize(2);
  gr1->SetMarkerColor(colors[1]);
  gr1->SetLineColor(colors[1]);
  mg->Add(gr1);

  // Model 2 (Recommended) - now with 10% error bars
  TGraphErrors *gr2 =
      new TGraphErrors(pressure_vec.size(), &pressure_vec[0], &model2_vec[0],
                       &pressure_err[0], &model2_err[0]);
  gr2->SetName("ATIMA 1.2 LS-theory");
  gr2->SetMarkerStyle(markerStyles[2]);
  gr2->SetMarkerSize(1.3);
  gr2->SetMarkerColor(colors[2]);
  gr2->SetLineColor(colors[2]);
  mg->Add(gr2);

  // Model 3 - now with 10% error bars
  TGraphErrors *gr3 =
      new TGraphErrors(pressure_vec.size(), &pressure_vec[0], &model3_vec[0],
                       &pressure_err[0], &model3_err[0]);
  gr3->SetName("ATIMA 1.2 no LS");
  gr3->SetMarkerStyle(markerStyles[3]);
  gr3->SetMarkerSize(1.3);
  gr3->SetMarkerColor(colors[3]);
  gr3->SetLineColor(colors[3]);
  mg->Add(gr3);

  // Model 4 - now with 10% error bars
  TGraphErrors *gr4 =
      new TGraphErrors(pressure_vec.size(), &pressure_vec[0], &model4_vec[0],
                       &pressure_err[0], &model4_err[0]);
  gr4->SetName("ATIMA 1.4 Weick");
  gr4->SetMarkerStyle(markerStyles[4]);
  gr4->SetMarkerSize(1.3);
  gr4->SetMarkerColor(colors[4]);
  gr4->SetLineColor(colors[4]);
  mg->Add(gr4);

  // Model 5 - now with 10% error bars
  TGraphErrors *gr5 =
      new TGraphErrors(pressure_vec.size(), &pressure_vec[0], &model5_vec[0],
                       &pressure_err[0], &model5_err[0]);
  gr5->SetName("Electrical component");
  gr5->SetMarkerStyle(markerStyles[5]);
  gr5->SetMarkerSize(1.3);
  gr5->SetMarkerColor(colors[5]);
  gr5->SetLineColor(colors[5]);
  mg->Add(gr5);

  // Model 6 - now with 10% error bars
  TGraphErrors *gr6 =
      new TGraphErrors(pressure_vec.size(), &pressure_vec[0], &model6_vec[0],
                       &pressure_err[0], &model6_err[0]);
  gr6->SetName("Nuclear component");
  gr6->SetMarkerStyle(markerStyles[6]);
  gr6->SetMarkerSize(1.3);
  gr6->SetMarkerColor(colors[6]);
  gr6->SetLineColor(colors[6]);
  mg->Add(gr6);

  gStyle->SetPadRightMargin(0.4);
  TCanvas *canvas =
      new TCanvas("dE_comparison", "Energy Loss vs Pressure", 1800, 900);
  canvas->SetGridx(1);
  canvas->SetGridy(1);
  canvas->SetTicks(1, 1);

  mg->Draw("ALP");
  mg->SetTitle(";Gas Pressure (Torr);#DeltaE (MeV)");

  TLegend *legend = new TLegend(0.70, 0.15, 0.99, 0.95);
  legend->SetTextFont(42);
  legend->SetTextSize(0.032);
  legend->SetBorderSize(1);
  legend->SetFillStyle(1001);
  legend->SetFillColorAlpha(kWhite, 0.85);

  legend->AddEntry(grExp, "Si Data", "lp");
  legend->AddEntry(gr0, "Hubert (He-base)", "lp");
  legend->AddEntry(gr1, "Ziegler (H-base)", "lp");
  legend->AddEntry(gr2, "ATIMA 1.2 LS-theory", "lp");
  legend->AddEntry(gr3, "ATIMA 1.2 no LS-correction", "lp");
  legend->AddEntry(gr4, "ATIMA 1.4 Weick", "lp");
  legend->AddEntry(gr5, "Electrical component", "lp");
  legend->AddEntry(gr6, "Nuclear component", "lp");

  legend->Draw();

  canvas->Update();

  canvas->SaveAs("plots/DeltaE_vs_Pressure_Comparison.png");

  std::cout << "Plot saved!" << std::endl;

  inFile->Close();
}
