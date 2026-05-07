#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PipelineMutex.hpp"
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <map>
#include <ostream>
#include <utility>
#include <vector>

void ResetEvent(Int_t leftdE[18], Int_t rightdE[18], Int_t totaldE[18],
                ULong64_t timestamp_distribution[36], UInt_t all_flags[36],
                Int_t hits[36], Int_t &cathode, Int_t &grid) {

  for (Int_t k = 0; k < 18; k++) {
    leftdE[k] = rightdE[k] = totaldE[k] = 0;
  }
  for (Int_t k = 0; k < 36; k++) {
    timestamp_distribution[k] = 0;
    all_flags[k] = 0;
    hits[k] = 0;
  }
  cathode = grid = 0;
}

Int_t GetStripNumber(TString map_name) {
  TString num_str = map_name;
  num_str.Remove(0, 1);
  return num_str.Atoi();
}

Bool_t CheckEventComplete(Int_t leftdE[18], Int_t rightdE[18],
                          Int_t totaldE[18]) {
  if (totaldE[0] == 0)
    return kFALSE;

  for (Int_t strip = 1; strip <= 7; strip += 2) {
    if (leftdE[strip] == 0) {
      return kFALSE;
    }
  }

  for (Int_t strip = 2; strip <= 8; strip += 2) {
    if (rightdE[strip] == 0) {
      return kFALSE;
    }
  }

  return kTRUE;
}

void InitEvent(Int_t leftdE[18], Int_t rightdE[18], Int_t totaldE[18],
               ULong64_t all_timestamps[36], UInt_t all_flags[36],
               Int_t hits[36], Int_t &cathode, Int_t &grid,
               ULong64_t grid_ts, UShort_t grid_energy, UInt_t grid_flag) {
  ResetEvent(leftdE, rightdE, totaldE, all_timestamps, all_flags, hits, cathode,
             grid);
  grid = grid_energy;
  all_timestamps[35] = grid_ts;
  all_flags[35] = grid_flag;
  hits[35]++;
}

void FinalizeEvent(TTree *output_tree, Int_t leftdE[18], Int_t rightdE[18],
                   Int_t totaldE[18], UInt_t all_flags[36], Int_t hits[36],
                   Bool_t &is_complete, Int_t &total_events,
                   Int_t &complete_events, Int_t &incomplete_events,
                   Int_t &incomplete_with_fake,
                   Int_t &incomplete_with_saturation,
                   Int_t &incomplete_with_pileup) {
  for (Int_t s = 1; s < 17; s++) {
    totaldE[s] = leftdE[s] + rightdE[s];
  }
  is_complete = CheckEventComplete(leftdE, rightdE, totaldE);
  total_events++;

  if (is_complete) {
    output_tree->Fill();
    complete_events++;
    return;
  }

  incomplete_events++;

  Bool_t has_fake = kFALSE;
  Bool_t has_saturation = kFALSE;
  Bool_t has_pileup = kFALSE;
  for (Int_t k = 0; k < 36; k++) {
    if (hits[k] > 0) {
      UInt_t flag = all_flags[k];
      if (flag & (1 << 24))
        has_fake = kTRUE;
      if ((flag & (1 << 2)) || (flag & (1 << 3)))
        has_saturation = kTRUE;
      if (flag & (1 << 30))
        has_pileup = kTRUE;
    }
  }
  if (has_fake)
    incomplete_with_fake++;
  if (has_saturation)
    incomplete_with_saturation++;
  if (has_pileup)
    incomplete_with_pileup++;
}

void BuildEvents(std::vector<TString> input_filenames,
                 std::vector<TString> output_names, Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  Int_t n_files = input_filenames.size();
  for (Int_t i = 0; i < n_files; i++) {
    TString input_filename = input_filenames[i];
    TString input_filepath = input_filename + ".root";
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
    UInt_t flags;
    ULong64_t timestamp;
    input_tree->SetBranchAddress("Board", &board);
    input_tree->SetBranchAddress("Channel", &channel);
    input_tree->SetBranchAddress("Energy", &energy);
    input_tree->SetBranchAddress("ShiftedTimestamp", &timestamp);
    input_tree->SetBranchAddress("Flags", &flags);

    std::cout << "Loading baskets into memory..." << std::endl;
    input_tree->LoadBaskets();

    Long64_t n_entries = input_tree->GetEntries();

    std::cout << "Pass 1: collecting grid timestamps..." << std::endl;
    std::vector<ULong64_t> grid_ts;
    std::vector<UShort_t> grid_energies;
    std::vector<UInt_t> grid_flags;
    grid_ts.reserve(3000000);
    grid_energies.reserve(3000000);
    grid_flags.reserve(3000000);

    for (Long64_t j = 0; j < n_entries; j++) {
      input_tree->GetEntry(j);
      std::map<std::pair<Int_t, Int_t>, TString>::const_iterator it =
          Constants::channelMap.find(
              std::pair<Int_t, Int_t>(board, channel));
      if (it == Constants::channelMap.end())
        continue;
      if (it->second == "Grid") {
        grid_ts.push_back(timestamp);
        grid_energies.push_back(energy);
        grid_flags.push_back(flags);
      }
      if (j % 10000000 == 0)
        std::cout << "  Pass 1: " << j << "/" << n_entries << std::endl;
    }
    std::cout << "Found " << grid_ts.size() << " grid events" << std::endl;

    if (grid_ts.empty()) {
      std::cerr << "No grid events found, skipping file" << std::endl;
      input_file->Close();
      continue;
    }

    TString output_name = output_names[i];
    TString output_filepath = output_name + ".root";
    TFile *output_file = IO::OpenForWriting(output_filepath);

    Int_t leftdE[18], rightdE[18], totaldE[18];
    ULong64_t all_timestamps[36];
    UInt_t all_flags[36];
    Int_t hits[36];
    Int_t cathode, grid;
    Bool_t is_complete;

    TTree *output_tree = new TTree("event", "MUSIC events");
    output_tree->Branch("LeftdE", leftdE, "LeftdE[18]/I");
    output_tree->Branch("RightdE", rightdE, "RightdE[18]/I");
    output_tree->Branch("TotaldE", totaldE, "TotaldE[18]/I");
    output_tree->Branch("AllTimestamps", all_timestamps, "AllTimestamps[36]/l");
    output_tree->Branch("AllFlags", all_flags, "AllFlags[36]/i");
    output_tree->Branch("Hits", hits, "Hits[36]/I");
    output_tree->Branch("Cathode", &cathode, "Cathode/I");
    output_tree->Branch("Grid", &grid, "Grid/I");
    output_tree->Branch("IsComplete", &is_complete, "IsComplete/O");

    Int_t emptyChannelMapEvents = 0;
    Int_t total_events = 0;
    Int_t complete_events = 0;
    Int_t incomplete_events = 0;
    Int_t incomplete_with_fake = 0;
    Int_t incomplete_with_saturation = 0;
    Int_t incomplete_with_pileup = 0;
    Long64_t skipped_pre_first = 0;
    Long64_t skipped_post_last = 0;

    std::size_t cur_grid_idx = 0;
    InitEvent(leftdE, rightdE, totaldE, all_timestamps, all_flags, hits,
              cathode, grid, grid_ts[0], grid_energies[0], grid_flags[0]);

    ULong64_t first_grid = grid_ts.front();
    ULong64_t last_grid = grid_ts.back();

    std::cout << "Pass 2: assigning hits to nearest grid..." << std::endl;

    for (Long64_t j = 0; j < n_entries; j++) {
      input_tree->GetEntry(j);

      if (timestamp < first_grid) {
        skipped_pre_first++;
        continue;
      }
      if (timestamp > last_grid) {
        skipped_post_last++;
        continue;
      }

      std::map<std::pair<Int_t, Int_t>, TString>::const_iterator it =
          Constants::channelMap.find(
              std::pair<Int_t, Int_t>(board, channel));
      if (it == Constants::channelMap.end())
        continue;

      TString map_name = it->second;

      if (map_name == "") {
        emptyChannelMapEvents++;
        continue;
      }

      if (map_name == "Grid")
        continue;

      while (cur_grid_idx + 1 < grid_ts.size()) {
        Long64_t dist_cur =
            std::abs(static_cast<Long64_t>(timestamp) -
                     static_cast<Long64_t>(grid_ts[cur_grid_idx]));
        Long64_t dist_next =
            std::abs(static_cast<Long64_t>(timestamp) -
                     static_cast<Long64_t>(grid_ts[cur_grid_idx + 1]));
        if (dist_next >= dist_cur)
          break;

        FinalizeEvent(output_tree, leftdE, rightdE, totaldE, all_flags, hits,
                      is_complete, total_events, complete_events,
                      incomplete_events, incomplete_with_fake,
                      incomplete_with_saturation, incomplete_with_pileup);

        cur_grid_idx++;
        InitEvent(leftdE, rightdE, totaldE, all_timestamps, all_flags, hits,
                  cathode, grid, grid_ts[cur_grid_idx],
                  grid_energies[cur_grid_idx], grid_flags[cur_grid_idx]);
      }

      if (map_name == "Strip0") {
        totaldE[0] += energy;
        if (hits[0] == 0)
          all_timestamps[0] = timestamp;
        all_flags[0] |= flags;
        hits[0]++;

      } else if (map_name.BeginsWith("L")) {
        Int_t strip_num = GetStripNumber(map_name);
        leftdE[strip_num] += energy;
        if (hits[strip_num] == 0)
          all_timestamps[strip_num] = timestamp;
        all_flags[strip_num] |= flags;
        hits[strip_num]++;

      } else if (map_name.BeginsWith("R")) {
        Int_t strip_num = GetStripNumber(map_name);
        rightdE[strip_num] += energy;
        if (hits[strip_num + 16] == 0)
          all_timestamps[strip_num + 16] = timestamp;
        all_flags[strip_num + 16] |= flags;
        hits[strip_num + 16]++;

      } else if (map_name == "Strip17") {
        totaldE[17] += energy;
        if (hits[33] == 0)
          all_timestamps[33] = timestamp;
        all_flags[33] |= flags;
        hits[33]++;

      } else if (map_name == "Cathode") {
        cathode += energy;
        if (hits[34] == 0)
          all_timestamps[34] = timestamp;
        all_flags[34] |= flags;
        hits[34]++;
      }

      if (j % 10000000 == 0) {
        std::cout << "  Pass 2: " << j << "/" << n_entries << std::endl;
      }
    }

    FinalizeEvent(output_tree, leftdE, rightdE, totaldE, all_flags, hits,
                  is_complete, total_events, complete_events, incomplete_events,
                  incomplete_with_fake, incomplete_with_saturation,
                  incomplete_with_pileup);

    output_file->cd();
    output_tree->Write("event", TObject::kOverwrite);
    output_file->Close();
    input_file->Close();

    if (emptyChannelMapEvents != 0)
      std::cout << "Observed " << emptyChannelMapEvents
                << " events with empty entry in channel map." << std::endl;

    std::cout << "Skipped (before first grid): " << skipped_pre_first
              << std::endl;
    std::cout << "Skipped (after last grid): " << skipped_post_last
              << std::endl;
    std::cout << "Total events: " << total_events << std::endl;
    std::cout << "Complete events: " << complete_events << " ("
              << (100.0 * complete_events / total_events) << "%)" << std::endl;
    std::cout << "Incomplete events: " << incomplete_events << " ("
              << (100.0 * incomplete_events / total_events) << "%)"
              << std::endl;

    if (incomplete_events > 0) {
      std::cout << "\nIncomplete events with rejection-quality flags:"
                << std::endl;
      std::cout << "  Fake events: " << incomplete_with_fake << " ("
                << (100.0 * incomplete_with_fake / incomplete_events) << "%)"
                << std::endl;
      std::cout << "  Saturated: " << incomplete_with_saturation << " ("
                << (100.0 * incomplete_with_saturation / incomplete_events)
                << "%)" << std::endl;
      std::cout << "  Pileup: " << incomplete_with_pileup << " ("
                << (100.0 * incomplete_with_pileup / incomplete_events) << "%)"
                << std::endl;
    }

    std::cout << "Events saved to: " << output_filepath << std::endl;
  }
}

void EventBuilder() {
  Bool_t reprocess_initial = kTRUE;

  std::vector<TString> filenames, output_names;

  std::vector<FileSpec> specs = BuildFileSpecs();
  for (std::size_t k = 0; k < specs.size(); k++) {
    filenames.push_back(SortedName(specs[k]));
    std::cout << "Processing file: " << std::endl;
    std::cout << filenames.back() << std::endl;
    output_names.push_back(EventsName(specs[k]));
  }

  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                TString(gSystem->pwd()) + "/plots",
                                TString(gSystem->pwd()) + "/root_files");
  BuildEvents(filenames, output_names, reprocess_initial);
}
