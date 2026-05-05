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

Double_t GetStoppingPowerFromLISE(const TString &filename, Int_t modelIndex,
                                  Double_t energy_MeV_per_u) {

  std::ifstream file(filename.Data());
  if (!file.is_open()) {
    std::cerr << "ERROR: Could not open file " << filename << std::endl;
    return -1.0;
  }

  std::string line;
  for (int i = 0; i < 3; ++i) {
    std::getline(file, line);
  }

  std::vector<Double_t> energies;
  std::vector<std::vector<Double_t>> dedx_values;

  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '!')
      continue;

    std::istringstream iss(line);
    std::vector<std::string> tokens;
    std::string token;

    while (std::getline(iss, token, '\t')) {
      tokens.push_back(token);
    }

    if (tokens.size() < 2)
      continue;

    try {
      Double_t E = std::stod(tokens[0]);

      if (energies.empty() || energies.back() != E) {
        energies.push_back(E);
        std::vector<Double_t> dedx_for_models(7);

        for (int m = 0; m < 7; ++m) {
          int col_index = 1 + m * 2;
          if (col_index < tokens.size()) {
            dedx_for_models[m] = std::stod(tokens[col_index]);
          }
        }

        dedx_values.push_back(dedx_for_models);
      }
    } catch (const std::exception &e) {
      continue;
    }
  }

  file.close();

  if (energies.size() < 2) {
    std::cerr << "ERROR: Not enough data points in " << filename << std::endl;
    return -1.0;
  }

  Int_t idx = 0;
  for (idx = 0; idx < energies.size() - 1; ++idx) {
    if (energies[idx] <= energy_MeV_per_u &&
        energy_MeV_per_u <= energies[idx + 1]) {
      break;
    }
  }

  if (idx >= energies.size() - 1) {
    idx = energies.size() - 2;
  }

  Double_t E1 = energies[idx];
  Double_t E2 = energies[idx + 1];
  Double_t dedx1 = dedx_values[idx][modelIndex];
  Double_t dedx2 = dedx_values[idx + 1][modelIndex];

  Double_t dedx_interpolated =
      dedx1 + (energy_MeV_per_u - E1) / (E2 - E1) * (dedx2 - dedx1);

  return dedx_interpolated;
}

void CalcStoppingPower() {

  TString input_filepath = "root_files/SiCalibration_Results.root";
  TFile *inFile = new TFile(input_filepath, "UPDATE");

  if (!inFile || inFile->IsZombie()) {
    std::cerr << "ERROR: Could not read ROOT file." << std::endl;
    return;
  }

  TTree *calibration_results =
      static_cast<TTree *>(inFile->Get("CalibrationResults"));
  TTree *results =
      new TTree("StoppingPower", "Comparison of Stopping Power Tables");

  Int_t run_number;
  Double_t gas_pressure_torr;
  Double_t delta_E_experimental;
  Double_t copy_delta_E_experimental;
  Double_t deltaE_model0, deltaE_model1, deltaE_model2, deltaE_model3;
  Double_t deltaE_model4, deltaE_model5, deltaE_model6;

  calibration_results->SetBranchAddress("RunNumber", &run_number);
  calibration_results->SetBranchAddress("GasPressure", &gas_pressure_torr);
  calibration_results->SetBranchAddress("DeltaE", &delta_E_experimental);

  results->Branch("DeltaE", &copy_delta_E_experimental);
  results->Branch("DeltaE_Model0_Hubert", &deltaE_model0, "DeltaE_Model0/D");
  results->Branch("DeltaE_Model1_Ziegler", &deltaE_model1, "DeltaE_Model1/D");
  results->Branch("DeltaE_Model2_ATIMA12LS", &deltaE_model2, "DeltaE_Model2/D");
  results->Branch("DeltaE_Model3_ATIMA12NoLS", &deltaE_model3,
                  "DeltaE_Model3/D");
  results->Branch("DeltaE_Model4_ATIMA14Weick", &deltaE_model4,
                  "DeltaE_Model4/D");
  results->Branch("DeltaE_Model5_Electrical", &deltaE_model5,
                  "DeltaE_Model5/D");
  results->Branch("DeltaE_Model6_Nuclear", &deltaE_model6, "DeltaE_Model6/D");

  const Double_t detector_length_cm = 35;
  const Double_t detector_length_micron = detector_length_cm * 10000;
  const Int_t A_Cl37 = 37;
  const Double_t beam_energy_TOF = 92;
  const Double_t beam_energy_per_u = beam_energy_TOF / A_Cl37;

  Double_t deltaE_models[7] = {0};

  const Double_t titanium_entrance = 0.9; // mg / cm2
  const Double_t titanium_exit = 1.3;     // mg / cm2
  const TString ti_lise_filename = "lise/37Cl_in_Ti_NOT_MICRON.lise";

  Double_t deltaE_entrance, deltaE_gas, deltaE_exit;
  Double_t beam_energy_per_u_gas, beam_energy_per_u_exit;

  Int_t nentries = calibration_results->GetEntries();
  for (Int_t i = 0; i < nentries; i++) {
    calibration_results->GetEntry(i);
    std::cout << "Run number: " << run_number
              << " with gas pressure: " << gas_pressure_torr << std::endl;

    TString lise_filename =
        Form("lise/37Cl_in_He4_%03.0fTorr_293K.lise", gas_pressure_torr);

    for (Int_t model = 0; model < 7; ++model) {
      Double_t entrance_dedx_MeV_per_mg_cm2 =
          GetStoppingPowerFromLISE(ti_lise_filename, model, beam_energy_per_u);
      deltaE_entrance = entrance_dedx_MeV_per_mg_cm2 * titanium_entrance;
      beam_energy_per_u_gas = (beam_energy_TOF - deltaE_entrance) / A_Cl37;
      if (gas_pressure_torr == 0) {
        deltaE_gas = 0;
      } else {
        Double_t dedx_MeV_per_micron = GetStoppingPowerFromLISE(
            lise_filename, model, beam_energy_per_u_gas);
        Double_t current_energy = beam_energy_TOF - deltaE_entrance;
        deltaE_gas = 0;

        const Double_t segment_micron = 1000;

        for (Double_t distance = 0; distance < detector_length_micron;
             distance += segment_micron) {

          Double_t segment_length =
              std::min(segment_micron, detector_length_micron - distance);
          Double_t current_energy_per_u = current_energy / A_Cl37;
          Double_t dedx = GetStoppingPowerFromLISE(lise_filename, model,
                                                   current_energy_per_u);
          Double_t dE_segment = dedx * segment_length;
          deltaE_gas += dE_segment;
          current_energy -= dE_segment;
        }
      }

      beam_energy_per_u_exit =
          (beam_energy_TOF - deltaE_entrance - deltaE_gas) / A_Cl37;
      Double_t exit_dedx_MeV_per_mg_cm2 = GetStoppingPowerFromLISE(
          ti_lise_filename, model, beam_energy_per_u_exit);
      deltaE_exit = exit_dedx_MeV_per_mg_cm2 * titanium_exit;
      deltaE_models[model] = deltaE_entrance + deltaE_gas + deltaE_exit;
    }

    copy_delta_E_experimental = delta_E_experimental;
    deltaE_model0 = deltaE_models[0];
    deltaE_model1 = deltaE_models[1];
    deltaE_model2 = deltaE_models[2];
    deltaE_model3 = deltaE_models[3];
    deltaE_model4 = deltaE_models[4];
    deltaE_model5 = deltaE_models[5];
    deltaE_model6 = deltaE_models[6];

    results->Fill();
  }

  results->Write();
  inFile->Close();
  delete inFile;

  std::cout << "Stopping power comparison complete!" << std::endl;
}
