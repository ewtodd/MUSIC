#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "StoppingPowerLISE.hpp"
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
#include <iostream>
#include <vector>

void CalcStoppingPower() {

  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");

  TString input_filepath =
      project_root + "/root_files/" + Constants::CALIBRATION_RESULTS_FILE;
  TFile *inFile = IO::OpenForWriting(input_filepath, "UPDATE");

  if (!inFile || inFile->IsZombie()) {
    std::cerr << "ERROR: Could not read ROOT file." << std::endl;
    return;
  }

  const Bool_t INCLUDE_PPAC = Constants::INCLUDE_PPAC;

  TTree *calibration_results =
      static_cast<TTree *>(inFile->Get("CalibrationResults"));
  TTree *results =
      new TTree("StoppingPower", "Comparison of Stopping Power Tables");

  Int_t run_number;
  Double_t gas_pressure_torr;
  Double_t delta_E_experimental;
  Double_t copy_delta_E_experimental;
  Double_t deltaE_model1, deltaE_model2, deltaE_model3;
  Double_t deltaE_model4, deltaE_model5, deltaE_model6;

  calibration_results->SetBranchAddress("RunNumber", &run_number);
  calibration_results->SetBranchAddress("GasPressure", &gas_pressure_torr);
  calibration_results->SetBranchAddress("DeltaE", &delta_E_experimental);

  results->Branch("DeltaE", &copy_delta_E_experimental);
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

  const Double_t ppac_layer_micron = 1.5; // um, single mylar foil
  const Int_t ppac_n_layers = 4;          // four foils per PPAC
  const TString mylar_lise_filename = "lise/37Cl_in_Mylar.lise";

  Double_t deltaE_ppac, deltaE_entrance, deltaE_gas, deltaE_exit;
  Double_t beam_energy_per_u_entrance, beam_energy_per_u_gas,
      beam_energy_per_u_exit;

  Int_t nentries = calibration_results->GetEntries();
  const Double_t segment_micron = 100;

  for (Int_t i = 0; i < nentries; i++) {
    calibration_results->GetEntry(i);
    std::cout << "Run number: " << run_number
              << " with gas pressure: " << gas_pressure_torr << std::endl;

    TString lise_filename =
        Form("lise/37Cl_in_He4_%03.0fTorr_293K.lise", gas_pressure_torr);

    for (Int_t model = 1; model < 7; model++) {
      Double_t current_energy = beam_energy_TOF;
      deltaE_ppac = 0;
      for (Int_t layer = 0; layer < ppac_n_layers; layer++) {
        Double_t dedx = GetStoppingPowerFromLISE(mylar_lise_filename, model,
                                                 current_energy / A_Cl37);
        Double_t dE = dedx * ppac_layer_micron;
        deltaE_ppac += dE;
        current_energy -= dE;
      }
      beam_energy_per_u_entrance = (beam_energy_TOF - deltaE_ppac) / A_Cl37;

      if (!INCLUDE_PPAC) {
        beam_energy_per_u_entrance = beam_energy_per_u;
        deltaE_ppac = 0;
      }

      Double_t entrance_dedx_MeV_per_mg_cm2 = GetStoppingPowerFromLISE(
          ti_lise_filename, model, beam_energy_per_u_entrance);
      deltaE_entrance = entrance_dedx_MeV_per_mg_cm2 * titanium_entrance;
      beam_energy_per_u_gas =
          (beam_energy_TOF - deltaE_ppac - deltaE_entrance) / A_Cl37;
      if (gas_pressure_torr == 0) {
        deltaE_gas = 0;
      } else {
        Double_t current_energy =
            beam_energy_TOF - deltaE_entrance - deltaE_ppac;
        deltaE_gas = 0;

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
          (beam_energy_TOF - deltaE_ppac - deltaE_entrance - deltaE_gas) /
          A_Cl37;
      Double_t exit_dedx_MeV_per_mg_cm2 = GetStoppingPowerFromLISE(
          ti_lise_filename, model, beam_energy_per_u_exit);
      deltaE_exit = exit_dedx_MeV_per_mg_cm2 * titanium_exit;
      // deltaE_ppac is 0 when INCLUDE_PPAC=kFALSE, so this is always consistent
      // with the experimental dE = BEAM_ENERGY_TOF - E_Si.
      deltaE_models[model] =
          deltaE_ppac + deltaE_entrance + deltaE_gas + deltaE_exit;
    }

    copy_delta_E_experimental = delta_E_experimental;
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
