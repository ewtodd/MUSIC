#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PipelineMutex.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <cstddef>
#include <iostream>
#include <mutex>
#include <ostream>
#include <vector>

void BuildTraces(std::vector<TString> input_output_filenames,
                 std::vector<TString> file_labels, Bool_t save_plots = kTRUE,
                 Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

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
    ULong64_t all_timestamps[36];
    UInt_t all_flags[36];
    Int_t hits[36];
    Int_t cathode, grid;

    input_output_tree->SetBranchAddress("LeftdE", leftdE);
    input_output_tree->SetBranchAddress("RightdE", rightdE);
    input_output_tree->SetBranchAddress("TotaldE", totaldE);
    input_output_tree->SetBranchAddress("AllTimestamps", all_timestamps);
    input_output_tree->SetBranchAddress("AllFlags", all_flags);
    input_output_tree->SetBranchAddress("Hits", hits);
    input_output_tree->SetBranchAddress("Cathode", &cathode);
    input_output_tree->SetBranchAddress("Grid", &grid);

    std::cout << "Loading baskets into memory..." << std::endl;
    input_output_tree->LoadBaskets();

    Long64_t n_entries = input_output_tree->GetEntries();
    std::cout << "Building traces from " << n_entries << " entries..."
              << std::endl;

    TH2D *h2_TotalE_vs_StripE =
        new TH2D("h2_TotalE_vs_StripE",
                 Form("Sum of Strips 3-4 vs Sum of Strips 1-2 (Four "
                      "Consecutive Strips with dE > %d ADC);"
                      "Sum of Strips 3-4 [ADC];"
                      "Sum of Strips 1-2 [ADC]",
                      Constants::TRIGGER_THRESHOLD),
                 250, 0, 12000, 250, 0, 12000);

    Int_t save_count = 0;
    Int_t save_count_per_strip[18] = {0};

    for (Long64_t j = 0; j < n_entries; j++) {
      input_output_tree->GetEntry(j);

      const Int_t nPoints = 18;
      Int_t strip_index[nPoints];

      for (Int_t k = 0; k < nPoints; k++) {
        strip_index[k] = k;
      }

      Bool_t save_event = kFALSE;
      Int_t triggerStrip = -1;

      for (Int_t k = 0; k < nPoints - 3; k++) {
        if (totaldE[k] > Constants::TRIGGER_THRESHOLD &&
            totaldE[k + 1] > Constants::TRIGGER_THRESHOLD &&
            totaldE[k + 2] > Constants::TRIGGER_THRESHOLD &&
            totaldE[k + 3] > Constants::TRIGGER_THRESHOLD) {

          Double_t sum_first_two = totaldE[k] + totaldE[k + 1];
          Double_t sum_second_two = totaldE[k + 2] + totaldE[k + 3];

          h2_TotalE_vs_StripE->Fill(sum_second_two, sum_first_two);

          Bool_t can_save =
              Constants::PER_STRIP_TRACE_LIMIT
                  ? (save_count_per_strip[k] < Constants::MAX_TRACE_SAVES)
                  : (save_count < Constants::MAX_TRACE_SAVES);
          if (can_save) {
            save_event = kTRUE;
            triggerStrip = k;
            if (Constants::PER_STRIP_TRACE_LIMIT)
              save_count_per_strip[k]++;
            else
              save_count++;
          }
          break;
        }
      }

      if (save_event && save_plots) {
        TGraph *TraceLeft = new TGraph(nPoints, strip_index, leftdE);
        TGraph *TraceRight = new TGraph(nPoints, strip_index, rightdE);
        TGraph *TraceTotal = new TGraph(nPoints, strip_index, totaldE);

        PlottingUtils::ConfigureGraph(
            TraceLeft, kBlue + 1,
            Form("Left dE Event %lld;Strip Index;Left dE [ADC]", j));
        TraceLeft->SetMarkerColor(kBlue + 1);
        TraceLeft->GetYaxis()->SetRangeUser(0, 3000);

        PlottingUtils::ConfigureGraph(
            TraceRight, kRed + 1,
            Form("Right dE Event %lld;Strip Index;Right dE [ADC]", j));
        TraceRight->SetMarkerColor(kRed + 1);
        TraceRight->GetYaxis()->SetRangeUser(0, 3000);

        PlottingUtils::ConfigureGraph(
            TraceTotal, kBlack,
            Form("Total dE Event %lld;Strip Index;Total dE [ADC]", j));
        TraceTotal->SetMarkerColor(kBlack);
        TraceTotal->GetYaxis()->SetRangeUser(0, 6000);

        {
          std::lock_guard<std::mutex> lock(g_plot_mutex);

          TString trace_subdir = "traces/" + file_label;
          TString trace_tag = Form("strip%d_event%lld", triggerStrip, j);

          TCanvas *c_left = PlottingUtils::GetConfiguredCanvas(kFALSE);
          TraceLeft->Draw("ALP");
          if (Constants::SAVE_LR_TRACES)
            PlottingUtils::SaveFigure(c_left, "trace_left_" + trace_tag,
                                      trace_subdir, PlotSaveOptions::kLINEAR);
          delete c_left;

          TCanvas *c_right = PlottingUtils::GetConfiguredCanvas(kFALSE);
          TraceRight->Draw("ALP");
          if (Constants::SAVE_LR_TRACES)
            PlottingUtils::SaveFigure(c_right, "trace_right_" + trace_tag,
                                      trace_subdir, PlotSaveOptions::kLINEAR);
          delete c_right;

          TCanvas *c_total = PlottingUtils::GetConfiguredCanvas(kFALSE);
          TraceTotal->Draw("ALP");
          PlottingUtils::SaveFigure(c_total, "trace_total_" + trace_tag,
                                    trace_subdir, PlotSaveOptions::kLINEAR);
          delete c_total;

          std::cout << "Saved event from strip " << triggerStrip << " (entry "
                    << j << ") under " << trace_subdir << std::endl;
        }

        delete TraceLeft;
        delete TraceRight;
        delete TraceTotal;
      }

      if ((j + 1) % 100000 == 0) {
        std::cout << "Processed " << j + 1 << " / " << n_entries
                  << " entries..." << std::endl;
      }
    }

    input_output_file->cd();
    h2_TotalE_vs_StripE->Write("", TObject::kOverwrite);
    if (save_plots) {
      std::lock_guard<std::mutex> lock(g_plot_mutex);
      TCanvas *histCanvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDraw2DHistogram(h2_TotalE_vs_StripE,
                                                 histCanvas);
      h2_TotalE_vs_StripE->GetYaxis()->SetTitleOffset(1.3);
      PlottingUtils::SaveFigure(histCanvas, "totalE_vs_stripE",
                                "traces/" + file_label,
                                PlotSaveOptions::kLINEAR);
      delete histCanvas;
    }
    delete h2_TotalE_vs_StripE;

    input_output_file->Write("", TObject::kOverwrite);
    input_output_file->Close();
    std::cout << "Finished processing " << input_output_filename << std::endl;
  }
}
void TraceCreator() {
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
  BuildTraces(input_output_filenames, file_labels, kTRUE, reprocess_initial);
}
