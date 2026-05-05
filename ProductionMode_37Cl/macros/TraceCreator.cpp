#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <ostream>
#include <vector>

void BuildTraces(std::vector<TString> input_output_filenames,
                 Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  Int_t n_files = input_output_filenames.size();
  for (Int_t i = 0; i < n_files; i++) {
    TString input_output_filename = input_output_filenames[i];
    TString input_output_filepath =
        "root_files/" + input_output_filename + ".root";
    TFile *input_output_file = new TFile(input_output_filepath, "UPDATE");

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

    TDirectory *trace_dir = input_output_file->mkdir("traces", "traces", kTRUE);

    TH2D *h2_TotalE_vs_StripE =
        new TH2D("h2_TotalE_vs_StripE",
                 Form("Sum of Strips 3-4 vs Sum of Strips 1-2 (Four "
                      "Consecutive Strips with dE > %d ADC);"
                      "Sum of Strips 3-4 [ADC];"
                      "Sum of Strips 1-2 [ADC]",
                      Constants::TRIGGER_THRESHOLD),
                 250, 0, 10000, 250, 0, 10000);

    const Int_t MAX_SAVES = 3;
    Int_t save_count = 0;

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
          if (save_count < MAX_SAVES) {
            save_event = kTRUE;
            triggerStrip = k;
            save_count++;
          }
          break;
        }
      }

      TGraph *TraceLeft = new TGraph(nPoints, strip_index, leftdE);
      TGraph *TraceRight = new TGraph(nPoints, strip_index, rightdE);
      TGraph *TraceTotal = new TGraph(nPoints, strip_index, totaldE);

      if (save_event) {

        TraceLeft->SetName(Form("TraceLeft_e%lld", j));
        TraceLeft->SetTitle(Form(
            "Left dE vs Strip Index (Event %lld);Strip Index;Left dE [ADC]",
            j));
        TraceLeft->SetLineColor(kBlue);
        TraceLeft->SetLineWidth(2);
        TraceLeft->SetMarkerColor(kBlue);
        TraceLeft->GetYaxis()->SetRangeUser(0, 2000);

        TraceRight->SetName(Form("TraceRight_e%lld", j));
        TraceRight->SetTitle(Form(
            "Right dE vs Strip Index (Event %lld);Strip Index;Right dE[ADC]",
            j));
        TraceRight->SetLineColor(kRed);
        TraceRight->SetLineWidth(2);
        TraceRight->SetMarkerColor(kRed);
        TraceRight->GetYaxis()->SetRangeUser(0, 2000);

        TraceTotal->SetName(Form("TraceTotal_e%lld", j));
        TraceTotal->SetTitle(Form(
            "Total dE vs Strip Index (Event %lld);Strip Index;Total dE [ADC]",
            j));
        TraceTotal->SetLineColor(kBlack);
        TraceTotal->SetLineWidth(2);
        TraceTotal->SetMarkerColor(kBlack);
        TraceTotal->GetYaxis()->SetRangeUser(0, 5000);

        trace_dir->cd();
        TraceLeft->Write("", TObject::kOverwrite);
        TraceRight->Write("", TObject::kOverwrite);
        TraceTotal->Write("", TObject::kOverwrite);

        TCanvas *c_traces = new TCanvas("c_traces", "Traces", 1600, 1200);
        c_traces->Divide(1, 3);

        c_traces->cd(1);
        TraceLeft->Draw("ALP");

        c_traces->cd(2);
        TraceRight->Draw("ALP");

        c_traces->cd(3);
        TraceTotal->Draw("ALP");

        TString trace_filename =
            Form("%s_Strip%d_Event%lld", input_output_filename.Data(),
                 triggerStrip, j);
        PlottingUtils::SaveFigure(c_traces, trace_filename, "", PlotSaveOptions::kLINEAR);

        delete c_traces;

        std::cout << "Saved event from strip " << triggerStrip << " (entry "
                  << j << ") as " << trace_filename << std::endl;
      }
      delete TraceLeft;
      delete TraceRight;
      delete TraceTotal;

      if ((j + 1) % 100000 == 0) {
        std::cout << "Processed " << j + 1 << " / " << n_entries
                  << " entries..." << std::endl;
      }
    }

    input_output_file->cd();
    h2_TotalE_vs_StripE->Write("", TObject::kOverwrite);
    TCanvas *histCanvas = new TCanvas("", "", 1600, 800);
    PlottingUtils::ConfigureAndDraw2DHistogram(h2_TotalE_vs_StripE, histCanvas);
    PlottingUtils::SaveFigure(histCanvas, input_output_filename + "_hist",
                              "", PlotSaveOptions::kLINEAR);
    delete histCanvas;
    delete h2_TotalE_vs_StripE;

    input_output_file->Write("", TObject::kOverwrite);
    input_output_file->Close();
    std::cout << "Finished processing " << input_output_filename << std::endl;
  }
}
void TraceCreator() {
  Bool_t reprocess_initial = kTRUE;

  std::vector<TString> input_output_filenames;

  for (Int_t i = 0; i < Constants::N_FILES; i++) {
    input_output_filenames.push_back(Form("Events_Run37_%d", i));
    std::cout << "Processing file: " << std::endl;
    std::cout << input_output_filenames[i] << std::endl;
  }

  InitUtils::SetROOTPreferences();
  BuildTraces(input_output_filenames, reprocess_initial);
}
