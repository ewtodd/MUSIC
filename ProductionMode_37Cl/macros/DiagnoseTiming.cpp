#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Pipeline.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

void ExploreOneFile(TString input_filename) {
  TString input_filepath = input_filename + ".root";
  TFile *input_file = IO::OpenForReading(input_filepath);
  if (!input_file || input_file->IsZombie()) {
    std::cerr << "Error opening file: " << input_filepath << std::endl;
    return;
  }

  TTree *input_tree = static_cast<TTree *>(input_file->Get("Data_R"));
  if (!input_tree) {
    std::cerr << "Error getting tree from: " << input_filepath << std::endl;
    input_file->Close();
    return;
  }

  UShort_t board, channel, energy;
  ULong64_t timestamp;
  input_tree->SetBranchAddress("Board", &board);
  input_tree->SetBranchAddress("Channel", &channel);
  input_tree->SetBranchAddress("Energy", &energy);
  input_tree->SetBranchAddress("Timestamp", &timestamp);

  std::cout << "Loading baskets into memory..." << std::endl;
  input_tree->LoadBaskets();

  Long64_t n_entries = input_tree->GetEntries();
  std::cout << "Total entries: " << n_entries << std::endl;

  std::map<std::pair<Int_t, Int_t>, Long64_t> counts;
  ULong64_t ts_min = 0, ts_max = 0;

  for (Long64_t j = 0; j < n_entries; j++) {
    input_tree->GetEntry(j);
    counts[{board, channel}]++;
    if (j == 0) {
      ts_min = ts_max = timestamp;
    } else {
      if (timestamp < ts_min)
        ts_min = timestamp;
      if (timestamp > ts_max)
        ts_max = timestamp;
    }
    if (j % 10000000 == 0)
      std::cout << "  Pass 1: " << j << "/" << n_entries << std::endl;
  }

  Double_t span_s = (ts_max - ts_min) / 1e12;
  std::cout << "Time span: " << span_s << " s" << std::endl;

  std::cout << "Per-channel hit counts and rates" << std::endl;
  std::cout << "Board  Chan  Name                 Count        Rate (Hz)"
            << std::endl;
  for (std::map<std::pair<Int_t, Int_t>, Long64_t>::const_iterator it =
           counts.begin();
       it != counts.end(); ++it) {
    Int_t b = it->first.first;
    Int_t c = it->first.second;
    Long64_t n = it->second;
    TString name;
    std::map<std::pair<Int_t, Int_t>, TString>::const_iterator nit =
        Constants::channelMap.find({b, c});
    name = (nit == Constants::channelMap.end()) ? "" : nit->second;
    Double_t rate = span_s > 0 ? n / span_s : 0.0;
    std::cout << "  " << b << "      " << c << "     "
              << TString::Format("%-20s", name.Data()) << " " << n << "    "
              << rate << std::endl;
  }

  std::vector<std::pair<Int_t, Int_t>> candidates = {
      {0, 0},  // Cathode
      {0, 8},  // L1   (Board 0)
      {0, 15}, // Grid (Board 0)
      {1, 0},  // L3   (Board 1)
      {2, 0},  // R4   (Board 2)
      {3, 0},  // L5   (Board 3)
  };

  std::map<std::pair<Int_t, Int_t>, TH1F *> hE;
  std::map<std::pair<Int_t, Int_t>, TH1F *> hDT;
  std::map<std::pair<Int_t, Int_t>, ULong64_t> last_ts;
  for (Int_t k = 0; k < Int_t(candidates.size()); k++) {
    std::pair<Int_t, Int_t> bc = candidates[k];
    TString tag = TString::Format("B%d_Ch%d", bc.first, bc.second);
    hE[bc] = new TH1F(PlottingUtils::GetRandomName(),
                      tag + " energy;Energy [ch];Counts", 500, 0, 16384);
    hDT[bc] = new TH1F(PlottingUtils::GetRandomName(),
                       tag + " delta-T;#Delta t [#mus];Counts", 600, 0, 3000);
    last_ts[bc] = 0;
  }

  for (Long64_t j = 0; j < n_entries; j++) {
    input_tree->GetEntry(j);
    std::pair<Int_t, Int_t> key(board, channel);
    std::map<std::pair<Int_t, Int_t>, TH1F *>::iterator eit = hE.find(key);
    if (eit == hE.end())
      continue;
    eit->second->Fill(energy);
    if (last_ts[key] != 0) {
      Double_t dt_us = (timestamp - last_ts[key]) / 1e6;
      hDT[key]->Fill(dt_us);
    }
    last_ts[key] = timestamp;
    if (j % 10000000 == 0)
      std::cout << "  Pass 2: " << j << "/" << n_entries << std::endl;
  }

  std::cout << "Per-candidate summary (energy + delta-T quantiles)"
            << std::endl;
  for (Int_t k = 0; k < Int_t(candidates.size()); k++) {
    std::pair<Int_t, Int_t> bc = candidates[k];
    TString tag = TString::Format("B%d_Ch%d", bc.first, bc.second);
    TString name;
    std::map<std::pair<Int_t, Int_t>, TString>::const_iterator cit =
        Constants::channelMap.find(bc);
    name = (cit == Constants::channelMap.end()) ? "" : cit->second;

    TH1F *he = hE[bc];
    TH1F *hd = hDT[bc];

    if (he->GetEntries() == 0) {
      std::cout << "  " << tag << "  (" << name << "): NO DATA" << std::endl;
      continue;
    }

    Double_t qE[3];
    Double_t pE[3] = {0.05, 0.5, 0.95};
    he->GetQuantiles(3, qE, pE);

    Double_t qD[5];
    Double_t pD[5] = {0.5, 0.9, 0.99, 0.999, 0.9999};
    hd->GetQuantiles(5, qD, pD);

    std::cout << "  " << tag << "  (" << name << ")" << std::endl;
    std::cout << "    E   N=" << (Long64_t)he->GetEntries()
              << "  mean=" << he->GetMean() << "  q05=" << qE[0]
              << "  q50=" << qE[1] << "  q95=" << qE[2] << std::endl;
    std::cout << "    dT[us] median=" << qD[0] << "  q90=" << qD[1]
              << "  q99=" << qD[2] << "  q99.9=" << qD[3]
              << "  q99.99=" << qD[4] << std::endl;
  }

  for (Int_t k = 0; k < Int_t(candidates.size()); k++) {
    std::pair<Int_t, Int_t> bc = candidates[k];
    TString tag = TString::Format("B%d_Ch%d", bc.first, bc.second);

    TCanvas *cE = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDrawHistogram(hE[bc], kBlue + 1,
                                             tag + ";Energy [ch];Counts");
    PlottingUtils::SaveFigure(cE, "energy_" + tag, "explore",
                              PlotSaveOptions::kLOG);
    delete cE;

    TCanvas *cD = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDrawHistogram(hDT[bc], kRed + 1,
                                             tag + ";#Delta t [#mus];Counts");
    PlottingUtils::SaveFigure(cD, "deltat_" + tag, "explore",
                              PlotSaveOptions::kLOG);
    delete cD;
  }

  input_file->Close();
  std::cout << "Plots saved under plots/explore/" << std::endl;
}

void DiagnoseTiming() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                project_root + "/plots",
                                project_root + "/root_files");

  std::vector<FileSpec> specs = BuildFileSpecs();
  for (Int_t k = 0; k < Int_t(specs.size()); k++) {
    TString name = RawRootName(specs[k]);
    if (name.EndsWith(".root"))
      name.Remove(name.Length() - 5, 5);
    std::cout << "Exploring: " << name << std::endl;
    ExploreOneFile(name);
  }
}
