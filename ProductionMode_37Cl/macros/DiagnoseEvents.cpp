#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Pipeline.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

std::pair<Int_t, Int_t> FindBoardChannel(const TString &name) {
  for (std::map<std::pair<Int_t, Int_t>, TString>::const_iterator it =
           Constants::channelMap.begin();
       it != Constants::channelMap.end(); ++it) {
    if (it->second == name)
      return it->first;
  }
  return std::pair<Int_t, Int_t>(-1, -1);
}

void DiagnoseOneFile(TString input_filename, TString file_label) {
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
  input_tree->SetBranchAddress("ShiftedTimestamp", &timestamp);

  std::cout << "Loading baskets into memory..." << std::endl;
  input_tree->LoadBaskets();

  std::vector<TString> required_names = {"Strip0", "L1", "L3", "L5", "L7",
                                         "R2",     "R4", "R6", "R8"};

  std::map<std::pair<Int_t, Int_t>, TString> required_lookup;
  for (Int_t i = 0; i < Int_t(required_names.size()); i++) {
    std::pair<Int_t, Int_t> bc = FindBoardChannel(required_names[i]);
    if (bc.first < 0) {
      std::cerr << "WARNING: " << required_names[i]
                << " not found in channelMap" << std::endl;
      continue;
    }
    required_lookup[bc] = required_names[i];
  }

  std::pair<Int_t, Int_t> grid_bc = FindBoardChannel("Grid");

  std::map<TString, std::vector<ULong64_t>> per_channel_ts;
  for (Int_t i = 0; i < Int_t(required_names.size()); i++) {
    per_channel_ts[required_names[i]] = std::vector<ULong64_t>();
    per_channel_ts[required_names[i]].reserve(3000000);
  }
  std::vector<ULong64_t> grid_ts;
  grid_ts.reserve(3000000);

  Long64_t n_entries = input_tree->GetEntries();
  std::cout << "Walking " << n_entries << " entries..." << std::endl;

  for (Long64_t j = 0; j < n_entries; j++) {
    input_tree->GetEntry(j);
    std::pair<Int_t, Int_t> bc(board, channel);

    if (bc == grid_bc) {
      grid_ts.push_back(timestamp);
    }

    std::map<std::pair<Int_t, Int_t>, TString>::const_iterator rit =
        required_lookup.find(bc);
    if (rit != required_lookup.end())
      per_channel_ts[rit->second].push_back(timestamp);

    if (j % 10000000 == 0)
      std::cout << "  Progress: " << j << "/" << n_entries << std::endl;
  }

  std::cout << "Collected:" << std::endl;
  std::cout << "  Grid: " << grid_ts.size() << " hits" << std::endl;
  for (Int_t i = 0; i < Int_t(required_names.size()); i++) {
    std::cout << "  " << required_names[i] << ": "
              << per_channel_ts[required_names[i]].size() << " hits"
              << std::endl;
  }

  if (grid_ts.empty()) {
    std::cerr << "No grid hits, abort" << std::endl;
    input_file->Close();
    return;
  }

  std::cout << "Per-channel hits-per-window summary" << std::endl;
  std::cout << "channel    P(0)     P(1)     P(2)     P(3+)" << std::endl;

  for (Int_t i = 0; i < Int_t(required_names.size()); i++) {
    TString ch_name = required_names[i];
    const std::vector<ULong64_t> &ch_ts = per_channel_ts[ch_name];

    TH1F *h_dt =
        new TH1F(PlottingUtils::GetRandomName(),
                 ch_name + " #Delta t to nearest Grid;hit_t - grid_t [#mus];"
                           "Counts",
                 400, -5, 5);

    TH1I *h_hits = new TH1I(PlottingUtils::GetRandomName(),
                            ch_name + " hits per grid window;"
                                      "Hits in (grid_{N}, grid_{N+1});"
                                      "Counts",
                            8, 0, 8);

    Int_t grid_idx = 0;
    for (Int_t k = 0; k < Int_t(ch_ts.size()); k++) {
      ULong64_t ts = ch_ts[k];
      while (grid_idx + 1 < Int_t(grid_ts.size()) &&
             grid_ts[grid_idx + 1] <= ts)
        grid_idx++;

      ULong64_t g_left = grid_ts[grid_idx];
      ULong64_t g_right = (grid_idx + 1 < Int_t(grid_ts.size()))
                              ? grid_ts[grid_idx + 1]
                              : ULong64_t(0);

      Long64_t dt_ps;
      if (g_right > 0 && (g_right - ts) < (ts - g_left))
        dt_ps = static_cast<Long64_t>(ts) - static_cast<Long64_t>(g_right);
      else
        dt_ps = static_cast<Long64_t>(ts) - static_cast<Long64_t>(g_left);

      Double_t dt_us = dt_ps / 1e6;
      h_dt->Fill(dt_us);
    }

    Int_t ch_idx = 0;
    Long64_t n_total_windows = 0;
    Long64_t n0 = 0, n1 = 0, n2 = 0, n3plus = 0;
    for (Int_t g = 0; g + 1 < Int_t(grid_ts.size()); g++) {
      ULong64_t g_ts = grid_ts[g];
      ULong64_t g_next = grid_ts[g + 1];

      while (ch_idx < Int_t(ch_ts.size()) && ch_ts[ch_idx] < g_ts)
        ch_idx++;

      Int_t hits = 0;
      Int_t kk = ch_idx;
      while (kk < Int_t(ch_ts.size()) && ch_ts[kk] < g_next) {
        hits++;
        kk++;
      }
      h_hits->Fill(TMath::Min(hits, 7));
      n_total_windows++;
      if (hits == 0)
        n0++;
      else if (hits == 1)
        n1++;
      else if (hits == 2)
        n2++;
      else
        n3plus++;
    }

    Double_t inv_total = n_total_windows > 0 ? 1.0 / n_total_windows : 0.0;
    std::cout << "  " << TString::Format("%-8s", ch_name.Data()) << " "
              << TString::Format("%.4f", n0 * inv_total) << "  "
              << TString::Format("%.4f", n1 * inv_total) << "  "
              << TString::Format("%.4f", n2 * inv_total) << "  "
              << TString::Format("%.4f", n3plus * inv_total) << std::endl;

    TCanvas *cdt = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDrawHistogram(h_dt, kBlue + 1);
    PlottingUtils::SaveFigure(cdt, "delta_t_" + ch_name,
                              "diagnostics/" + file_label,
                              PlotSaveOptions::kLOG);
    delete cdt;

    TCanvas *chw = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDrawHistogram(h_hits, kRed + 1);
    PlottingUtils::SaveFigure(chw, "hits_per_window_" + ch_name,
                              "diagnostics/" + file_label,
                              PlotSaveOptions::kLINEAR);
    delete chw;
  }

  input_file->Close();
  std::cout << "Diagnostic plots saved under plots/diagnostics/" << file_label
            << std::endl;
}

void DiagnoseEvents() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                TString(gSystem->pwd()) + "/plots",
                                TString(gSystem->pwd()) + "/root_files");

  std::vector<FileSpec> specs = BuildFileSpecs();
  for (Int_t k = 0; k < Int_t(specs.size()); k++) {
    TString name = SortedName(specs[k]);
    TString file_label = FileLabel(specs[k]);
    std::cout << "Diagnosing: " << name << std::endl;
    DiagnoseOneFile(name, file_label);
  }
}
