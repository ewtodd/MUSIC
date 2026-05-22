#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Pipeline.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <mutex>
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
        static_cast<TTree *>(input_output_file->Get("events_MeV"));
    if (!input_output_tree) {
      std::cerr << "Error getting tree from: " << input_output_filepath
                << std::endl;
      input_output_file->Close();
      continue;
    }

    Float_t leftdE[18], rightdE[18], totaldE[18];
    Float_t cathode;

    input_output_tree->SetBranchAddress("LeftdE", leftdE);
    input_output_tree->SetBranchAddress("RightdE", rightdE);
    input_output_tree->SetBranchAddress("TotaldE", totaldE);
    input_output_tree->SetBranchAddress("Cathode", &cathode);

    std::cout << "Loading baskets into memory..." << std::endl;
    input_output_tree->LoadBaskets();

    Long64_t n_entries = input_output_tree->GetEntries();
    std::cout << "Building traces from " << n_entries << " entries..."
              << std::endl;

    const Int_t s_lo = Constants::FirstStrip();
    const Int_t s_hi = Constants::LastStrip();
    const Int_t n_active = Constants::NumActiveStrips();

    TH2D *h2_TotalE_vs_StripE[18] = {nullptr};
    for (Int_t s = s_lo; s <= s_hi; s++) {
      h2_TotalE_vs_StripE[s] =
          new TH2D(Form("h2_TotalE_vs_StripE_s%d", s),
                   Form(";Strip %d #DeltaE [MeV];Total #DeltaE [MeV]", s), 200,
                   Constants::STRIP_E_MIN_MEV, Constants::STRIP_E_MAX_MEV, 400,
                   Constants::TOTAL_E_MIN_MEV, Constants::TOTAL_E_MAX_MEV);
    }

    Int_t save_count = 0;

    for (Long64_t j = 0; j < n_entries; j++) {
      input_output_tree->GetEntry(j);

      Double_t event_total = 0.0;
      for (Int_t s = s_lo; s <= s_hi; s++)
        event_total += Double_t(totaldE[s]);

      for (Int_t s = s_lo; s <= s_hi; s++)
        h2_TotalE_vs_StripE[s]->Fill(Double_t(totaldE[s]), event_total);

      if (save_plots && save_count < Constants::MAX_TRACE_SAVES) {
        save_count++;

        TGraph *TraceTotal = new TGraph(n_active);
        TGraph *TraceLeft = new TGraph(n_active);
        TGraph *TraceRight = new TGraph(n_active);
        for (Int_t k = 0; k < n_active; k++) {
          Int_t s = s_lo + k;
          TraceTotal->SetPoint(k, s, Double_t(totaldE[s]));
          // Sim populates Left/Right only on indices 1..16; 0 and 17 are
          // single-ended guard strips and live in TotaldE alone.
          Double_t l = (s == 0 || s == 17) ? 0.0 : Double_t(leftdE[s]);
          Double_t r = (s == 0 || s == 17) ? 0.0 : Double_t(rightdE[s]);
          if (Constants::IGNORE_SHORT_STRIPS) {
            if (Constants::IsShortSide(s, 'L'))
              l = 0.0;
            if (Constants::IsShortSide(s, 'R'))
              r = 0.0;
          }
          TraceLeft->SetPoint(k, s, l);
          TraceRight->SetPoint(k, s, r);
        }
        TString total_title = Form("Event %lld;Strip Index;#DeltaE [MeV]", j);
        PlottingUtils::ConfigureGraph(TraceTotal, kBlack, total_title);
        TraceTotal->SetMarkerColor(kBlack);
        TraceTotal->GetYaxis()->SetRangeUser(Constants::STRIP_E_MIN_MEV,
                                             Constants::STRIP_E_MAX_MEV);

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
    for (Int_t s = s_lo; s <= s_hi; s++)
      h2_TotalE_vs_StripE[s]->Write("", TObject::kOverwrite);
    {
      std::lock_guard<std::mutex> lock(g_plot_mutex);
      TString summary_subdir = "trace_summary/" + file_label;

      for (Int_t s = s_lo; s <= s_hi; s++) {
        TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
        PlottingUtils::ConfigureAndDraw2DHistogram(h2_TotalE_vs_StripE[s], c);
        h2_TotalE_vs_StripE[s]->GetYaxis()->SetTitleOffset(1.3);
        PlottingUtils::SaveFigure(c, Form("totalE_vs_stripE_s%d", s),
                                  summary_subdir, PlotSaveOptions::kLINEAR);
        delete c;
      }
    }
    for (Int_t s = s_lo; s <= s_hi; s++)
      delete h2_TotalE_vs_StripE[s];

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
    input_output_filenames.push_back(TracesName(specs[k]));
    std::cout << "Processing file: " << input_output_filenames.back()
              << std::endl;
    file_labels.push_back(FileLabel(specs[k]));
  }

  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");
  BuildTraces(input_output_filenames, file_labels, kTRUE, reprocess_initial);
}
