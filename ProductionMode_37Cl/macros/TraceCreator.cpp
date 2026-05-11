#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Normalization.hpp"
#include "Pipeline.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
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

    Baseline baseline;
    Bool_t have_baseline = LoadBaseline(input_output_file, baseline);
    if (!have_baseline) {
      std::cerr << "WARNING: no baseline tree in " << input_output_filepath
                << " -- subtracted traces will equal raw traces" << std::endl;
    }

    std::cout << "Loading baskets into memory..." << std::endl;
    input_output_tree->LoadBaskets();

    Long64_t n_entries = input_output_tree->GetEntries();
    std::cout << "Building traces from " << n_entries << " entries..."
              << std::endl;

    TH2D *h2_StripE_vs_TotalE[18];
    for (Int_t s = 0; s < 18; s++) {
      h2_StripE_vs_TotalE[s] =
          new TH2D(Form("h2_StripE_vs_TotalE_s%d", s),
                   Form("Strip %d energy vs event total energy;"
                        "Strip %d #DeltaE [ADC];"
                        "Total #DeltaE [ADC]",
                        s, s),
                   200, -100, 3500, 400, -100, 60000);
    }

    Int_t save_count = 0;

    for (Long64_t j = 0; j < n_entries; j++) {
      input_output_tree->GetEntry(j);

      Double_t normalized_total[18];
      ComputeNormalized(baseline, leftdE, rightdE, totaldE, normalized_total);

      Double_t event_total_norm = 0.0;
      for (Int_t s = 0; s < 18; s++)
        event_total_norm += normalized_total[s];

      for (Int_t s = 0; s < 18; s++)
        h2_StripE_vs_TotalE[s]->Fill(normalized_total[s], event_total_norm);

      if (save_plots && save_count < Constants::MAX_TRACE_SAVES) {
        save_count++;

        TGraph *TraceTotal = new TGraph(18);
        TGraph *TraceLeft = new TGraph(18);
        TGraph *TraceRight = new TGraph(18);
        TString total_title;
        Double_t y_lo, y_hi;
        if (have_baseline) {
          for (Int_t k = 0; k < 18; k++) {
            TraceTotal->SetPoint(k, k, normalized_total[k]);
            Double_t l = (k == 0 || k == 17)
                             ? 0.0
                             : baseline.scale[k] * Double_t(leftdE[k]);
            Double_t r = (k == 0 || k == 17)
                             ? 0.0
                             : baseline.scale[k] * Double_t(rightdE[k]);
            TraceLeft->SetPoint(k, k, l);
            TraceRight->SetPoint(k, k, r);
          }
          total_title = Form("Normalized Event %lld;Strip Index;"
                             "Normalized #DeltaE [ADC]",
                             j);
          y_lo = -200;
          y_hi = 2000;
        } else {
          for (Int_t k = 0; k < 18; k++) {
            TraceTotal->SetPoint(k, k, totaldE[k]);
            TraceLeft->SetPoint(k, k, leftdE[k]);
            TraceRight->SetPoint(k, k, rightdE[k]);
          }
          total_title =
              Form("Event %lld;Strip Index;#DeltaE [ADC]", j);
          y_lo = 0;
          y_hi = 6000;
        }
        PlottingUtils::ConfigureGraph(TraceTotal, kBlack, total_title);
        TraceTotal->SetMarkerColor(kBlack);
        TraceTotal->GetYaxis()->SetRangeUser(y_lo, y_hi);

        TraceLeft->SetLineColor(kBlue + 1);
        TraceLeft->SetMarkerColor(kBlue + 1);
        TraceLeft->SetLineWidth(2);

        TraceRight->SetLineColor(kRed + 1);
        TraceRight->SetMarkerColor(kRed + 1);
        TraceRight->SetLineWidth(2);

        {
          std::lock_guard<std::mutex> lock(g_plot_mutex);

          TString trace_subdir = "traces/" + file_label;
          TString trace_tag = Form("event%lld", j);

          TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
          TraceTotal->Draw("ALP");
          TraceLeft->Draw("LP SAME");
          TraceRight->Draw("LP SAME");
          PlottingUtils::SaveFigure(c, "trace_" + trace_tag, trace_subdir,
                                    PlotSaveOptions::kLINEAR);
          delete c;

          std::cout << "Saved event " << j << " under " << trace_subdir
                    << (have_baseline ? " (normalized)" : " (raw)")
                    << std::endl;
        }

        delete TraceTotal;
        delete TraceLeft;
        delete TraceRight;
      }

      if ((j + 1) % 100000 == 0) {
        std::cout << "Processed " << j + 1 << " / " << n_entries
                  << " entries..." << std::endl;
      }
    }

    input_output_file->cd();
    for (Int_t s = 0; s < 18; s++)
      h2_StripE_vs_TotalE[s]->Write("", TObject::kOverwrite);
    {
      std::lock_guard<std::mutex> lock(g_plot_mutex);
      TString summary_subdir = "trace_summary/" + file_label;

      for (Int_t s = 0; s < 18; s++) {
        TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
        PlottingUtils::ConfigureAndDraw2DHistogram(h2_StripE_vs_TotalE[s], c);
        h2_StripE_vs_TotalE[s]->GetYaxis()->SetTitleOffset(1.3);
        PlottingUtils::SaveFigure(c, Form("stripE_vs_totalE_s%d", s),
                                  summary_subdir, PlotSaveOptions::kLINEAR);
        delete c;
      }
    }
    for (Int_t s = 0; s < 18; s++)
      delete h2_StripE_vs_TotalE[s];

    input_output_file->Write("", TObject::kOverwrite);
    input_output_file->Close();
    std::cout << "Finished processing " << input_output_filename << std::endl;
  }
}

void TraceCreator() {
  Bool_t reprocess_initial = kTRUE;

  std::vector<TString> input_output_filenames, file_labels;

  std::vector<FileSpec> specs = BuildFileSpecs();
  for (Int_t k = 0; k < Int_t(specs.size()); k++) {
    input_output_filenames.push_back(EventsName(specs[k]));
    std::cout << "Processing file: " << input_output_filenames.back()
              << std::endl;
    file_labels.push_back(FileLabel(specs[k]));
  }

  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots",
                                project_root + "/root_files");
  BuildTraces(input_output_filenames, file_labels, kTRUE, reprocess_initial);
}
