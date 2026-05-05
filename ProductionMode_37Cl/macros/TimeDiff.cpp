#include "Constants.hpp"
#include "InitUtils.hpp"
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <ostream>
#include <vector>

void CheckTimeDiff(std::vector<TString> input_filenames,
                   std::vector<TString> output_names,
                   Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  Int_t n_files = input_filenames.size();
  for (Int_t i = 0; i < n_files; i++) {
    TString input_filename = input_filenames[i];
    TString input_filepath = "root_files/" + input_filename + ".root";
    TFile *input_file = new TFile(input_filepath, "READ");

    if (!input_file || input_file->IsZombie()) {
      std::cerr << "Error opening file: " << input_filepath << std::endl;
      continue;
    }

    TTree *input_tree = static_cast<TTree *>(input_file->Get("Data_R"));
    if (!input_tree) {
      std::cerr << "Error getting tree from: " << input_filepath << std::endl;
      input_file->Close();
      continue;
    }

    UShort_t board, channel, energy;
    UInt_t flags;
    ULong64_t timestamp;
    input_tree->SetBranchAddress("Board", &board);
    input_tree->SetBranchAddress("Channel", &channel);
    input_tree->SetBranchAddress("Energy", &energy);
    input_tree->SetBranchAddress("ShiftedTimestamp", &timestamp);
    input_tree->SetBranchAddress("Flags", &flags);

    std::cout << "Loading baskets into memory..." << std::endl;
    input_tree->LoadBaskets();

    TString output_name = output_names[i];
    TString output_filepath = "root_files/" + output_name + ".root";
    TFile *output_file = new TFile(output_filepath, "RECREATE");

    Long64_t gridTimediff;
    Int_t grid;

    TTree *output_tree = new TTree("event", "MUSIC events");
    output_tree->Branch("Grid", &grid, "Grid/I");
    output_tree->Branch("GridTimediff", &gridTimediff, "GridTimediff/l");
    TH1I *gridTimediffHist = new TH1I(PlottingUtils::GetRandomName(),
                                      "; Time [us]; Counts", 100, 0, 250);

    ULong64_t event_time_zero = 0;
    Bool_t in_event = kFALSE;

    grid = 0;

    Long64_t n_entries = input_tree->GetEntries();
    std::cout << "Building events from " << n_entries << " entries..."
              << std::endl;

    Long64_t lastGrid = 0;
    for (Long64_t j = 0; j < n_entries; j++) {
      input_tree->GetEntry(j);
      TString map_name = Constants::channelMap.at({board, channel});

      if (map_name == "Grid") {
        gridTimediff = timestamp - lastGrid;
        lastGrid = timestamp;
        grid = energy;
        gridTimediffHist->Fill(gridTimediff * 1e-6);
        output_tree->Fill();
      }

      if (j % 10000000 == 0) {
        std::cout << "  Progress: " << j << "/" << n_entries << std::endl;
      }
    }

    TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDrawHistogram(gridTimediffHist, kRed);
    PlottingUtils::SaveFigure(canvas, "gridtimediff");

    output_file->cd();
    output_tree->Write("event", TObject::kOverwrite);
    gridTimediffHist->Write("gridTimediffHist", TObject::kOverwrite);
    output_file->Close();
    input_file->Close();
    std::cout << "Events saved to: " << output_filepath << std::endl;
  }
}

void TimeDiff() {
  Bool_t reprocess_initial = kTRUE;

  std::vector<TString> filenames, output_names;

  for (Int_t i = 0; i < Constants::N_FILES; i++) {
    filenames.push_back(Form("Sorted_Run37_%d", i));
    std::cout << "Processing file: " << std::endl;
    std::cout << filenames[i] << std::endl;
    output_names.push_back(Form("TimeDiff_Run37_%d", i));
  }

  InitUtils::SetROOTPreferences();
  CheckTimeDiff(filenames, output_names, reprocess_initial);
}
