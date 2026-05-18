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
#include <TH2D.h>
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
  if (!reprocess)
    return;

  Int_t n_files = input_filenames.size();
  for (Int_t i = 0; i < n_files; i++) {
    TString input_filename = input_filenames[i];
    TString input_filepath = input_filename + ".root";
    TString sidecar_filepath = input_filename + ".traces.root";
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
    const Double_t total_e_min =
        ev.is_mev ? Constants::TOTAL_E_MIN_MEV : Constants::TOTAL_E_MIN_ADC;
    const Double_t total_e_max =
        ev.is_mev ? Constants::TOTAL_E_MAX_MEV : Constants::TOTAL_E_MAX_ADC;

    input_tree->SetCacheSize(256 * 1024 * 1024);
    input_tree->AddBranchToCache("*", kTRUE);

    Long64_t n_entries = input_tree->GetEntries();
    std::cout << "Building traces from " << n_entries << " entries..."
              << std::endl;

    TH2D *h2_totalE_vs_stripE[18];
    for (Int_t s = 0; s < 18; s++) {
      h2_totalE_vs_stripE[s] = new TH2D(
          Form("h2_totalE_vs_stripE_s%d", s),
          Form(";Strip %d #DeltaE [%s];Total #DeltaE [%s]", s, unit, unit), 200,
          strip_e_min, strip_e_max, 400, total_e_min, total_e_max);
    }

    Int_t save_count = 0;

    for (Long64_t j = 0; j < n_entries; j++) {
      input_tree->GetEntry(j);
      ev.Decode();

      Double_t event_total = 0.0;
      for (Int_t s = 0; s < 18; s++)
        event_total += ev.total[s];

      for (Int_t s = 0; s < 18; s++)
        h2_totalE_vs_stripE[s]->Fill(ev.total[s], event_total);

      if (save_plots && save_count < Constants::MAX_TRACE_SAVES) {
        save_count++;

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

      if ((j + 1) % 100000 == 0) {
        std::cout << "Processed " << j + 1 << " / " << n_entries
                  << " entries..." << std::endl;
      }
    }

    input_file->Close();
    delete input_file;

    TFile *sidecar_file = IO::OpenForWriting(sidecar_filepath, "RECREATE");
    if (!sidecar_file || sidecar_file->IsZombie()) {
      std::cerr << "Error opening sidecar: " << sidecar_filepath << std::endl;
      for (Int_t s = 0; s < 18; s++)
        delete h2_totalE_vs_stripE[s];
      continue;
    }
    sidecar_file->SetCompressionAlgorithm(
        ROOT::RCompressionSetting::EAlgorithm::kZSTD);
    sidecar_file->SetCompressionLevel(5);

    {
      std::lock_guard<std::mutex> lock(g_plot_mutex);
      TString summary_subdir = "trace_summary/" + file_label;

      for (Int_t s = 0; s < 18; s++) {
        TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
        PlottingUtils::ConfigureAndDraw2DHistogram(h2_totalE_vs_stripE[s], c);
        h2_totalE_vs_stripE[s]->GetYaxis()->SetTitleOffset(1.3);
        PlottingUtils::SaveFigure(c, Form("totalE_vs_stripE_s%d", s),
                                  summary_subdir, PlotSaveOptions::kLINEAR);
        sidecar_file->cd();
        h2_totalE_vs_stripE[s]->Write(Form("h2_totalE_vs_stripE_s%d", s),
                                      TObject::kOverwrite);
        delete c;
      }
    }
    for (Int_t s = 0; s < 18; s++)
      delete h2_totalE_vs_stripE[s];
    sidecar_file->Write("", TObject::kOverwrite);
    sidecar_file->Close();
    delete sidecar_file;
    std::cout << "Finished processing " << input_filename << std::endl;
  }
}

#endif
