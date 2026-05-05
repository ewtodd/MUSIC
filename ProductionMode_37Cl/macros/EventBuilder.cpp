#include "Constants.hpp"
#include "InitUtils.hpp"
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <ostream>
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

void BuildEvents(std::vector<TString> input_filenames,
                 std::vector<TString> output_names, Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  Int_t n_files = input_filenames.size();
  for (Int_t i = 0; i < n_files; i++) {
    TString input_filename = input_filenames[i];
    TString input_filepath = "root_files/" + input_filename + ".root";
    TFile *input_file = new TFile(input_filepath, "READ");

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

    TString output_name = output_names[i];
    TString output_filepath = "root_files/" + output_name + ".root";
    TFile *output_file = new TFile(output_filepath, "RECREATE");

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

    Bool_t in_event = kFALSE;

    ResetEvent(leftdE, rightdE, totaldE, all_timestamps, all_flags, hits,
               cathode, grid);

    Long64_t n_entries = input_tree->GetEntries();
    std::cout << "Building events from " << n_entries << " entries..."
              << std::endl;

    Int_t emptyChannelMapEvents = 0;
    Int_t total_events = 0;
    Int_t complete_events = 0;
    Int_t incomplete_events = 0;
    Int_t incomplete_with_fake = 0;
    Int_t incomplete_with_saturation = 0;
    Int_t incomplete_with_pileup = 0;

    for (Long64_t j = 0; j < n_entries; j++) {
      input_tree->GetEntry(j);
      TString map_name = Constants::channelMap.at({board, channel});

      if (map_name == "") {
        emptyChannelMapEvents++;
        continue;
      }

      if (map_name == "Grid") {
        if (in_event) {
          for (Int_t s = 1; s < 17; s++) {
            totaldE[s] = leftdE[s] + rightdE[s];
          }
          is_complete = CheckEventComplete(leftdE, rightdE, totaldE);
          output_tree->Fill();
          total_events++;

          if (is_complete) {
            complete_events++;
          } else {
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
        }

        ResetEvent(leftdE, rightdE, totaldE, all_timestamps, all_flags, hits,
                   cathode, grid);

        grid = energy;
        all_timestamps[35] = timestamp;
        all_flags[35] = flags;
        hits[35]++;
        in_event = kTRUE;

      } else if (in_event) {
        if (map_name == "Strip0") {
          leftdE[0] = 0;
          rightdE[0] = 0;
          totaldE[0] = energy;
          all_timestamps[0] = timestamp;
          all_flags[0] = flags;
          hits[0]++;

        } else if (map_name.BeginsWith("L")) {
          Int_t strip_num = GetStripNumber(map_name);
          leftdE[strip_num] = energy;
          all_timestamps[strip_num] = timestamp;
          all_flags[strip_num] = flags;
          hits[strip_num]++;

        } else if (map_name.BeginsWith("R")) {
          Int_t strip_num = GetStripNumber(map_name);
          rightdE[strip_num] = energy;
          all_timestamps[strip_num + 16] = timestamp;
          all_flags[strip_num + 16] = flags;
          hits[strip_num + 16]++;

        } else if (map_name == "Strip17") {
          leftdE[17] = 0;
          rightdE[17] = 0;
          totaldE[17] = energy;
          all_timestamps[33] = timestamp;
          all_flags[33] = flags;
          hits[33]++;

        } else if (map_name == "Cathode") {
          cathode = energy;
          all_timestamps[34] = timestamp;
          all_flags[34] = flags;
          hits[34]++;
        }
      }

      if (j % 10000000 == 0) {
        std::cout << "  Progress: " << j << "/" << n_entries << std::endl;
      }
    }

    output_file->cd();
    output_tree->Write("event", TObject::kOverwrite);
    output_file->Close();
    input_file->Close();

    if (emptyChannelMapEvents != 0)
      std::cout << "Observed " << emptyChannelMapEvents
                << " events with empty entry in channel map." << std::endl;

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

  for (Int_t i = 0; i < Constants::N_FILES; i++) {
    filenames.push_back(Form("Sorted_Run37_%d", i));
    std::cout << "Processing file: " << std::endl;
    std::cout << filenames[i] << std::endl;
    output_names.push_back(Form("Events_Run37_%d", i));
  }

  InitUtils::SetROOTPreferences();
  BuildEvents(filenames, output_names, reprocess_initial);
}
