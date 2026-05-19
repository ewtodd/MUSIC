#ifndef TRACE_CREATOR_HPP
#define TRACE_CREATOR_HPP

#include "Constants.hpp"
#include "IOUtils.hpp"
#include "Normalization.hpp"
#include "Pipeline.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <mutex>
#include <vector>

inline void BuildTraces(std::vector<TString> input_filenames,
                        std::vector<TString> file_labels,
                        Bool_t save_plots = kTRUE, Bool_t reprocess = kFALSE) {
  if (!reprocess || !save_plots)
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

    TTree *input_tree = static_cast<TTree *>(input_file->Get("events"));
    if (!input_tree) {
      std::cerr << "No events tree in: " << input_filepath << std::endl;
      input_file->Close();
      continue;
    }

    EnergyView ev;
    ev.Attach(input_tree);
    if (!ev.is_mev)
      std::cerr << "[" << file_label
                << "] WARNING: no MeV branches; using uncalibrated ADC values"
                << std::endl;
    const char *unit = ev.Unit();
    const Double_t strip_e_min =
        ev.is_mev ? Constants::STRIP_E_MIN_MEV : Constants::STRIP_E_MIN_ADC;
    const Double_t strip_e_max =
        ev.is_mev ? Constants::STRIP_E_MAX_MEV : Constants::STRIP_E_MAX_ADC;

    Long64_t n_entries = input_tree->GetEntries();
    Long64_t n_to_save =
        TMath::Min(Long64_t(Constants::MAX_TRACE_SAVES), n_entries);
    std::cout << "Saving " << n_to_save << " per-event traces from "
              << input_filename << "..." << std::endl;

    for (Long64_t j = 0; j < n_to_save; j++) {
      input_tree->GetEntry(j);
      ev.Decode();

      TGraph *TraceTotal = new TGraph(18);
      TGraph *TraceLeft = new TGraph(18);
      TGraph *TraceRight = new TGraph(18);
      for (Int_t k = 0; k < 18; k++) {
        TraceTotal->SetPoint(k, k, ev.total[k]);
        Double_t l = (k == 0 || k == 17) ? 0.0 : ev.left[k];
        Double_t r = (k == 0 || k == 17) ? 0.0 : ev.right[k];
        TraceLeft->SetPoint(k, k, l);
        TraceRight->SetPoint(k, k, r);
      }
      TString total_title =
          Form("Event %lld;Strip Index;#DeltaE [%s]", j, unit);
      PlottingUtils::ConfigureGraph(TraceTotal, kBlack, total_title);
      TraceTotal->SetMarkerColor(kBlack);
      TraceTotal->GetYaxis()->SetRangeUser(strip_e_min, strip_e_max);

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

    input_file->Close();
    delete input_file;
    std::cout << "Finished processing " << input_filename << std::endl;
  }
}

#endif
