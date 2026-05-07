#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PipelineMutex.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1I.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <cstddef>
#include <iostream>
#include <mutex>
#include <ostream>
#include <vector>

void CheckTimeDiff(std::vector<TString> input_filenames,
                   std::vector<TString> output_names,
                   std::vector<TString> file_labels,
                   Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  Int_t n_files = input_filenames.size();
  for (Int_t i = 0; i < n_files; i++) {
    TString input_filename = input_filenames[i];
    TString input_filepath = input_filename + ".root";
    TString file_label = file_labels[i];
    TFile *input_file = IO::OpenForReading(input_filepath);

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
    ULong64_t timestamp;
    input_tree->SetBranchAddress("Board", &board);
    input_tree->SetBranchAddress("Channel", &channel);
    input_tree->SetBranchAddress("Energy", &energy);
    input_tree->SetBranchAddress("ShiftedTimestamp", &timestamp);

    std::cout << "Loading baskets into memory..." << std::endl;
    input_tree->LoadBaskets();

    TString output_name = output_names[i];
    TString output_filepath = output_name + ".root";
    TFile *output_file = IO::OpenForWriting(output_filepath);

    ULong64_t gridTimediff;
    Int_t grid;

    TTree *output_tree = new TTree("event", "MUSIC events");
    output_tree->Branch("Grid", &grid, "Grid/I");
    output_tree->Branch("GridTimediff", &gridTimediff, "GridTimediff/l");
    TH1I *gridTimediffHist = new TH1I(PlottingUtils::GetRandomName(),
                                      "; Time [us]; Counts", 100, 0, 300);

    Long64_t n_entries = input_tree->GetEntries();
    std::cout << "Walking " << n_entries << " entries..." << std::endl;

    ULong64_t lastGrid = 0;
    Bool_t have_lastGrid = kFALSE;
    for (Long64_t j = 0; j < n_entries; j++) {
      input_tree->GetEntry(j);
      TString map_name = Constants::channelMap.at({board, channel});

      if (map_name == "Grid") {
        if (have_lastGrid) {
          gridTimediff = timestamp - lastGrid;
          grid = energy;
          gridTimediffHist->Fill(gridTimediff * 1e-6);
          output_tree->Fill();
        }
        lastGrid = timestamp;
        have_lastGrid = kTRUE;
      }

      if (j % 10000000 == 0) {
        std::cout << "  Progress: " << j << "/" << n_entries << std::endl;
      }
    }

    {
      std::lock_guard<std::mutex> lock(g_plot_mutex);
      TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDrawHistogram(gridTimediffHist, kRed);
      PlottingUtils::SaveFigure(canvas, "gridtimediff",
                                "timediff/" + file_label,
                                PlotSaveOptions::kLOG);
      delete canvas;
    }

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

  std::vector<TString> filenames, output_names, file_labels;

  std::vector<FileSpec> specs = BuildFileSpecs();
  for (std::size_t k = 0; k < specs.size(); k++) {
    filenames.push_back(SortedName(specs[k]));
    std::cout << "Processing file: " << std::endl;
    std::cout << filenames.back() << std::endl;
    output_names.push_back(TimeDiffName(specs[k]));
    file_labels.push_back(FileLabel(specs[k]));
  }

  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                TString(gSystem->pwd()) + "/plots",
                                TString(gSystem->pwd()) + "/root_files");
  CheckTimeDiff(filenames, output_names, file_labels, reprocess_initial);
}
