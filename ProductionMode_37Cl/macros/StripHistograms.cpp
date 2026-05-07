#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PipelineMutex.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <cstddef>
#include <iostream>
#include <ostream>
#include <vector>

void BuildStripHistograms(std::vector<TString> input_output_filenames,
                          std::vector<TString> file_labels,
                          Bool_t save_plots = kTRUE,
                          Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  const Int_t LR_BINS = 300;
  const Double_t ADC_MAX = 16384.0;
  const Int_t TOTAL_BINS = 300;

  Int_t n_files = input_output_filenames.size();
  for (Int_t i = 0; i < n_files; i++) {
    TString input_output_filename = input_output_filenames[i];
    TString input_output_filepath = input_output_filename + ".root";
    TString file_label = file_labels[i];
    TFile *input_output_file =
        IO::OpenForWriting(input_output_filepath, "UPDATE");

    if (!input_output_file || input_output_file->IsZombie()) {
      std::cerr << "Error opening file: " << input_output_filepath << std::endl;
      continue;
    }

    TTree *input_output_tree =
        static_cast<TTree *>(input_output_file->Get("event"));
    if (!input_output_tree) {
      std::cerr << "Error getting tree from: " << input_output_filepath
                << std::endl;
      input_output_file->Close();
      continue;
    }

    Int_t leftdE[18], rightdE[18], totaldE[18];
    input_output_tree->SetBranchAddress("LeftdE", leftdE);
    input_output_tree->SetBranchAddress("RightdE", rightdE);
    input_output_tree->SetBranchAddress("TotaldE", totaldE);

    std::cout << "Loading baskets into memory..." << std::endl;
    input_output_tree->LoadBaskets();

    Long64_t n_entries = input_output_tree->GetEntries();
    std::cout << "Building strip histograms from " << n_entries << " entries..."
              << std::endl;

    TDirectory *hist_dir =
        input_output_file->mkdir("strip_histograms", "strip_histograms", kTRUE);

    TH1D *h_strip0 =
        new TH1D("h_Strip0", "Strip0;dE [ADC];Counts", TOTAL_BINS, 0, ADC_MAX);
    TH1D *h_strip17 = new TH1D("h_Strip17", "Strip17;dE [ADC];Counts",
                               TOTAL_BINS, 0, ADC_MAX);

    TH1D *h_left[17];
    TH1D *h_right[17];
    for (Int_t s = 1; s <= 16; s++) {
      h_left[s] = new TH1D(Form("h_L%d", s), Form("L%d;dE [ADC];Counts", s),
                           LR_BINS, 0, ADC_MAX);
      h_right[s] = new TH1D(Form("h_R%d", s), Form("R%d;dE [ADC];Counts", s),
                            LR_BINS, 0, ADC_MAX);
    }

    for (Long64_t j = 0; j < n_entries; j++) {
      input_output_tree->GetEntry(j);

      Double_t event_total = 0.0;
      for (Int_t s = 0; s <= 17; s++) {
        event_total += totaldE[s];
      }

      if (totaldE[0] > 0)
        h_strip0->Fill(totaldE[0]);
      if (totaldE[17] > 0)
        h_strip17->Fill(totaldE[17]);

      for (Int_t s = 1; s <= 16; s++) {
        if (leftdE[s] > 0)
          h_left[s]->Fill(leftdE[s]);
        if (rightdE[s] > 0)
          h_right[s]->Fill(rightdE[s]);
      }

      if ((j + 1) % 100000 == 0) {
        std::cout << "Processed " << j + 1 << " / " << n_entries
                  << " entries..." << std::endl;
      }
    }

    {
      std::lock_guard<std::mutex> lock(g_plot_mutex);
      hist_dir->cd();
      h_strip0->Write("", TObject::kOverwrite);
      h_strip17->Write("", TObject::kOverwrite);
      for (Int_t s = 1; s <= 16; s++) {
        h_left[s]->Write("", TObject::kOverwrite);
        h_right[s]->Write("", TObject::kOverwrite);
      }

      if (save_plots) {
        TString mode_tag = (Constants::EVENT_BUILDER_MODE ==
                            Constants::EVENT_BUILDER_NEAREST_GRID)
                               ? "nearest"
                               : "classic";
        TString plot_subdir = "strip_histograms_" + mode_tag + "/" + file_label;

        TCanvas *c0 = PlottingUtils::GetConfiguredCanvas(kFALSE);
        PlottingUtils::ConfigureAndDrawHistogram(h_strip0, kBlack);
        PlottingUtils::SaveFigure(c0, "strip_Strip0", plot_subdir,
                                  PlotSaveOptions::kLOG);
        delete c0;

        TCanvas *c17 = PlottingUtils::GetConfiguredCanvas(kFALSE);
        PlottingUtils::ConfigureAndDrawHistogram(h_strip17, kBlack);
        PlottingUtils::SaveFigure(c17, "strip_Strip17", plot_subdir,
                                  PlotSaveOptions::kLOG);
        delete c17;

        for (Int_t s = 1; s <= 16; s++) {
          TCanvas *cl = PlottingUtils::GetConfiguredCanvas(kFALSE);
          PlottingUtils::ConfigureAndDrawHistogram(h_left[s], kBlue + 1);
          PlottingUtils::SaveFigure(cl, Form("strip_L%d", s), plot_subdir,
                                    PlotSaveOptions::kLOG);
          delete cl;

          TCanvas *cr = PlottingUtils::GetConfiguredCanvas(kFALSE);
          PlottingUtils::ConfigureAndDrawHistogram(h_right[s], kRed + 1);
          PlottingUtils::SaveFigure(cr, Form("strip_R%d", s), plot_subdir,
                                    PlotSaveOptions::kLOG);
          delete cr;
        }

        std::cout << "Saved strip histograms under " << plot_subdir
                  << std::endl;
      }
    }

    delete h_strip0;
    delete h_strip17;
    for (Int_t s = 1; s <= 16; s++) {
      delete h_left[s];
      delete h_right[s];
    }

    input_output_file->Write("", TObject::kOverwrite);
    input_output_file->Close();
    std::cout << "Finished processing " << input_output_filename << std::endl;
  }
}

void StripHistograms() {
  Bool_t reprocess_initial = kTRUE;

  std::vector<TString> input_output_filenames, file_labels;

  std::vector<FileSpec> specs = BuildFileSpecs();
  for (std::size_t k = 0; k < specs.size(); k++) {
    input_output_filenames.push_back(EventsName(specs[k]));
    std::cout << "Processing file: " << std::endl;
    std::cout << input_output_filenames.back() << std::endl;
    file_labels.push_back(FileLabel(specs[k]));
  }

  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                TString(gSystem->pwd()) + "/plots",
                                TString(gSystem->pwd()) + "/root_files");
  BuildStripHistograms(input_output_filenames, file_labels, kTRUE,
                       reprocess_initial);
}
