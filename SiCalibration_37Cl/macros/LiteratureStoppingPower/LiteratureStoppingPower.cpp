#include "InitUtils.hpp"
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

void ReadLISEFile(const TString &filename, Double_t exp_min, Double_t exp_max,
                  std::vector<std::vector<Double_t>> &energies,
                  std::vector<std::vector<Double_t>> &stopping_powers) {

  std::ifstream file(filename.Data());
  if (!file.is_open()) {
    std::cerr << "ERROR: Could not open file " << filename << std::endl;
    return;
  }

  std::string line;
  for (int i = 0; i < 3; ++i) {
    std::getline(file, line);
  }

  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '!')
      continue;

    std::istringstream iss(line);
    std::vector<std::string> tokens;
    std::string token;

    while (std::getline(iss, token, '\t')) {
      tokens.push_back(token);
    }

    if (tokens.size() < 10)
      continue;

    try {
      for (Int_t model = 0; model < 5; ++model) {
        Double_t energy = std::stod(tokens[2 * model]);
        Double_t stopping = std::stod(tokens[2 * model + 1]);

        if (energy >= exp_min && energy <= exp_max) {
          energies[model].push_back(energy);
          stopping_powers[model].push_back(stopping);
        }
      }
    } catch (const std::exception &e) {
      continue;
    }
  }

  file.close();
}

void LiteratureStoppingPower() {
  InitUtils::SetROOTPreferences();

  Double_t exp_E_A[] = {2.624, 2.550, 2.450, 2.339, 2.207, 2.134, 2.038, 1.955,
                        1.843, 1.717, 1.624, 1.446, 1.414, 1.352, 1.230, 1.091,
                        0.927, 0.766, 0.763, 0.547, 0.495, 0.474, 0.415, 0.391,
                        0.356, 0.332, 0.301, 0.296, 0.251, 0.247, 0.212, 0.202,
                        0.183, 0.157, 0.134, 0.114};

  Double_t exp_stopping[] = {
      23.1, 25.8, 26.2, 25.7, 27.1, 28.0, 26.9, 28.2, 29.7, 29.5, 31.7, 30.9,
      31.1, 31.5, 32.3, 32.4, 33.5, 32.9, 32.7, 33.6, 30.1, 30.5, 28.5, 27.9,
      26.3, 26.2, 24.9, 25.2, 22.6, 21.9, 21.5, 19.5, 19.2, 17.1, 16.0, 13.2};

  Int_t n_exp = sizeof(exp_E_A) / sizeof(exp_E_A[0]);

  // Find energy range
  Double_t exp_min = *std::min_element(exp_E_A, exp_E_A + n_exp);
  Double_t exp_max = *std::max_element(exp_E_A, exp_E_A + n_exp);

  TString lise_filename = "37Cl_in_He4_230Torr_293K_NOT_MICRON.lise";
  std::vector<std::vector<Double_t>> energies(5);
  std::vector<std::vector<Double_t>> stopping_powers(5);

  ReadLISEFile(lise_filename, exp_min, exp_max, energies, stopping_powers);
  TMultiGraph *mg = new TMultiGraph();

  TGraph *grExp = new TGraph(n_exp, exp_E_A, exp_stopping);
  grExp->SetName("Experimental, Pierce/Blann 1968");
  grExp->SetMarkerStyle(20);
  grExp->SetMarkerSize(1.8);
  grExp->SetMarkerColor(kBlack);
  grExp->SetLineColor(kBlack);
  grExp->SetLineWidth(3);
  mg->Add(grExp);

  Int_t colors[] = {kBlue, kRed, kGreen + 2, kOrange, kMagenta};
  Int_t lineStyles[] = {1, 2, 3, 9, 1};
  const char *model_names[] = {"Hubert et al (He-base)",
                               "Ziegler et al (H-base)", "ATIMA 1.2 LS-theory",
                               "ATIMA 1.2 no LS-correction", "ATIMA 1.4 Weick"};

  TGraph *grModels[5];
  for (Int_t i = 0; i < 5; ++i) {
    if (energies[i].size() > 0) {
      grModels[i] = new TGraph(energies[i].size(), &energies[i][0],
                               &stopping_powers[i][0]);
      grModels[i]->SetName(model_names[i]);
      grModels[i]->SetLineColor(colors[i]);
      grModels[i]->SetLineWidth(2);
      grModels[i]->SetLineStyle(lineStyles[i]);
      mg->Add(grModels[i]);
    }
  }

  TCanvas *canvas =
      new TCanvas("stopping_power", "Stopping Power Comparison", 1800, 900);
  canvas->SetGridx(1);
  canvas->SetGridy(1);
  canvas->SetTicks(1, 1);

  mg->Draw("AL");
  mg->SetTitle(";E/A (MeV/u);dE/dx (MeV/mg/cm^{2})");

  TLegend *legend = new TLegend(0.62, 0.15, 0.95, 0.92);
  legend->SetTextFont(42);
  legend->SetTextSize(0.032);
  legend->SetBorderSize(1);
  legend->SetFillStyle(1001);
  legend->SetFillColorAlpha(kWhite, 0.85);

  legend->AddEntry(grExp, "Experimental, Pierce/Blann 1968", "lp");
  for (Int_t i = 0; i < 5; ++i) {
    if (energies[i].size() > 0) {
      legend->AddEntry(grModels[i], model_names[i], "l");
    }
  }

  legend->Draw();

  canvas->Update();

  // Save plot
  if (gSystem->AccessPathName("plots")) {
    gSystem->mkdir("plots", kTRUE);
  }
  canvas->SaveAs("plots/StoppingPower_vs_Energy_Comparison.png");
}
