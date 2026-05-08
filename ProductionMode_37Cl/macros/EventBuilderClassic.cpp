#include "BinaryUtils.hpp"
#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Pipeline.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <cstdlib>
#include <iostream>
#include <map>
#include <mutex>
#include <ostream>
#include <utility>
#include <vector>

namespace Classic {

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
               Int_t hits[36], Int_t &cathode, Int_t &grid, ULong64_t grid_ts,
               UShort_t grid_energy, UInt_t grid_flag) {
  ResetEvent(leftdE, rightdE, totaldE, all_timestamps, all_flags, hits, cathode,
             grid);
  grid = grid_energy;
  all_timestamps[35] = grid_ts;
  all_flags[35] = grid_flag;
  hits[35]++;
}

void FillEventHistograms(Int_t totaldE[18], Int_t hits[36], Bool_t has_any_flag,
                         TH2F *h_music, TH2F *h_music_clean,
                         TH2F *h_music_flagged, TH1F *h_mult) {
  for (Int_t s = 0; s < 18; s++) {
    h_music->Fill(Double_t(s), Double_t(totaldE[s]));
    if (has_any_flag)
      h_music_flagged->Fill(Double_t(s), Double_t(totaldE[s]));
    else
      h_music_clean->Fill(Double_t(s), Double_t(totaldE[s]));
  }
  Int_t mult = 0;
  for (Int_t k = 0; k < 36; k++) {
    mult += hits[k];
  }
  h_mult->Fill(Double_t(mult));
}

void GetFlagSummary(UInt_t all_flags[36], Int_t hits[36], Bool_t &has_fake,
                    Bool_t &has_saturation, Bool_t &has_pileup) {
  has_fake = kFALSE;
  has_saturation = kFALSE;
  has_pileup = kFALSE;
  for (Int_t k = 0; k < 36; k++) {
    if (hits[k] > 0) {
      UInt_t flag = all_flags[k];
      if (flag & CoMPASSData::FAKE_EVENT)
        has_fake = kTRUE;
      if ((flag & CoMPASSData::SATURATION_IN_GATE) ||
          (flag & CoMPASSData::INPUT_SATURATING))
        has_saturation = kTRUE;
      if (flag & CoMPASSData::PILEUP)
        has_pileup = kTRUE;
    }
  }
}

void FinalizeEvent(TTree *output_tree, Int_t leftdE[18], Int_t rightdE[18],
                   Int_t totaldE[18], UInt_t all_flags[36], Int_t hits[36],
                   Bool_t &is_complete, Int_t &total_events,
                   Int_t &complete_events, Int_t &complete_with_fake,
                   Int_t &complete_with_saturation, Int_t &complete_with_pileup,
                   Int_t &complete_rejected, Int_t &incomplete_events,
                   Int_t &incomplete_with_fake,
                   Int_t &incomplete_with_saturation,
                   Int_t &incomplete_with_pileup, TH2F *h_music,
                   TH2F *h_music_clean, TH2F *h_music_flagged, TH1F *h_mult) {
  for (Int_t s = 1; s < 17; s++) {
    totaldE[s] = leftdE[s] + rightdE[s];
  }
  is_complete = CheckEventComplete(leftdE, rightdE, totaldE);
  total_events++;

  Bool_t has_fake = kFALSE;
  Bool_t has_saturation = kFALSE;
  Bool_t has_pileup = kFALSE;
  GetFlagSummary(all_flags, hits, has_fake, has_saturation, has_pileup);
  Bool_t has_any_flag = has_fake || has_saturation || has_pileup;

  if (is_complete) {
    Bool_t reject = Constants::REJECT_FLAGGED_EVENTS && has_any_flag;
    if (!reject) {
      output_tree->Fill();
      FillEventHistograms(totaldE, hits, has_any_flag, h_music, h_music_clean,
                          h_music_flagged, h_mult);
    } else {
      complete_rejected++;
    }
    complete_events++;
    if (has_fake)
      complete_with_fake++;
    if (has_saturation)
      complete_with_saturation++;
    if (has_pileup)
      complete_with_pileup++;
    return;
  }

  incomplete_events++;
  if (has_fake)
    incomplete_with_fake++;
  if (has_saturation)
    incomplete_with_saturation++;
  if (has_pileup)
    incomplete_with_pileup++;
}

} // namespace Classic

void BuildEvents(std::vector<TString> input_filenames,
                 std::vector<TString> output_names,
                 std::vector<TString> file_labels, Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  Int_t n_files = input_filenames.size();
  for (Int_t i = 0; i < n_files; i++) {
    TString input_filename = input_filenames[i];
    TString file_label = file_labels[i];
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

    TH2F *h_music = new TH2F("hMUSIC",
                             "MUSIC strip energies (complete events);"
                             "Strip;Energy [ADC]",
                             18, -0.5, 17.5, 2000, 0, 16384);
    TH2F *h_music_clean =
        new TH2F("hMUSICClean",
                 "MUSIC strip energies (complete, no rejection flags);"
                 "Strip;Energy [ADC]",
                 18, -0.5, 17.5, 2000, 0, 16384);
    TH2F *h_music_flagged =
        new TH2F("hMUSICFlagged",
                 "MUSIC strip energies (complete, with rejection flags);"
                 "Strip;Energy [ADC]",
                 18, -0.5, 17.5, 2000, 0, 16384);
    TH1F *h_mult = new TH1F("hMult",
                            "Event multiplicity (complete events);"
                            "Multiplicity;Counts",
                            36, -0.5, 35.5);

    Int_t emptyChannelMapEvents = 0;
    Int_t total_events = 0;
    Int_t complete_events = 0;
    Int_t complete_with_fake = 0;
    Int_t complete_with_saturation = 0;
    Int_t complete_with_pileup = 0;
    Int_t complete_rejected = 0;
    Int_t incomplete_events = 0;
    Int_t incomplete_with_fake = 0;
    Int_t incomplete_with_saturation = 0;
    Int_t incomplete_with_pileup = 0;
    Long64_t grid_hits_total = 0;
    Long64_t cathode_hits_total = 0;
    Long64_t events_with_cathode = 0;
    Bool_t cur_event_has_cathode = kFALSE;

    Bool_t in_event = kFALSE;
    Classic::ResetEvent(leftdE, rightdE, totaldE, all_timestamps, all_flags,
                        hits, cathode, grid);

    std::cout << "Building events from " << n_entries << " entries..."
              << std::endl;

    input_file->cd();

    for (Long64_t j = 0; j < n_entries; j++) {
      input_tree->GetEntry(j);

      std::map<std::pair<Int_t, Int_t>, TString>::const_iterator it =
          Constants::channelMap.find(std::pair<Int_t, Int_t>(board, channel));
      if (it == Constants::channelMap.end())
        continue;

      TString map_name = it->second;

      if (map_name == "") {
        emptyChannelMapEvents++;
        continue;
      }

      if (map_name == "Cathode")
        cathode_hits_total++;

      if (map_name == "Grid") {
        grid_hits_total++;
        if (in_event) {
          if (cur_event_has_cathode)
            events_with_cathode++;
          Classic::FinalizeEvent(
              output_tree, leftdE, rightdE, totaldE, all_flags, hits,
              is_complete, total_events, complete_events, complete_with_fake,
              complete_with_saturation, complete_with_pileup, complete_rejected,
              incomplete_events, incomplete_with_fake,
              incomplete_with_saturation, incomplete_with_pileup, h_music,
              h_music_clean, h_music_flagged, h_mult);
        }
        Classic::InitEvent(leftdE, rightdE, totaldE, all_timestamps, all_flags,
                           hits, cathode, grid, timestamp, energy, flags);
        in_event = kTRUE;
        cur_event_has_cathode = kFALSE;

      } else if (in_event) {
        if (map_name == "Strip0") {
          totaldE[0] += energy;
          if (hits[0] == 0)
            all_timestamps[0] = timestamp;
          all_flags[0] |= flags;
          hits[0]++;

        } else if (map_name.BeginsWith("L")) {
          Int_t strip_num = Classic::GetStripNumber(map_name);
          leftdE[strip_num] += energy;
          if (hits[strip_num] == 0)
            all_timestamps[strip_num] = timestamp;
          all_flags[strip_num] |= flags;
          hits[strip_num]++;

        } else if (map_name.BeginsWith("R")) {
          Int_t strip_num = Classic::GetStripNumber(map_name);
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
          cur_event_has_cathode = kTRUE;
        }
      }

      if (j % 10000000 == 0) {
        std::cout << "  Progress: " << j << "/" << n_entries << std::endl;
      }
    }

    if (in_event) {
      if (cur_event_has_cathode)
        events_with_cathode++;
      Classic::FinalizeEvent(output_tree, leftdE, rightdE, totaldE, all_flags,
                             hits, is_complete, total_events, complete_events,
                             complete_with_fake, complete_with_saturation,
                             complete_with_pileup, complete_rejected,
                             incomplete_events, incomplete_with_fake,
                             incomplete_with_saturation, incomplete_with_pileup,
                             h_music, h_music_clean, h_music_flagged, h_mult);
    }

    output_file->cd();
    output_tree->Write("event", TObject::kOverwrite);
    h_music->Write("", TObject::kOverwrite);
    h_music_clean->Write("", TObject::kOverwrite);
    h_music_flagged->Write("", TObject::kOverwrite);
    h_mult->Write("", TObject::kOverwrite);

    {
      std::lock_guard<std::mutex> lock(g_plot_mutex);
      TString subdir = "events_classic/" + file_label;

      TCanvas *c_music = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDraw2DHistogram(h_music, c_music);
      h_music->GetYaxis()->SetTitleOffset(1.4);
      c_music->SetLeftMargin(0.18);
      PlottingUtils::SaveFigure(c_music, "music_strip_energies", subdir,
                                PlotSaveOptions::kLINEAR);
      delete c_music;

      TCanvas *c_music_clean = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDraw2DHistogram(h_music_clean, c_music_clean);
      h_music_clean->GetYaxis()->SetTitleOffset(1.4);
      c_music_clean->SetLeftMargin(0.18);
      PlottingUtils::SaveFigure(c_music_clean, "music_strip_energies_clean",
                                subdir, PlotSaveOptions::kLINEAR);
      delete c_music_clean;

      TCanvas *c_music_flagged = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDraw2DHistogram(h_music_flagged,
                                                 c_music_flagged);
      h_music_flagged->GetYaxis()->SetTitleOffset(1.4);
      c_music_flagged->SetLeftMargin(0.18);
      PlottingUtils::SaveFigure(c_music_flagged, "music_strip_energies_flagged",
                                subdir, PlotSaveOptions::kLINEAR);
      delete c_music_flagged;

      TCanvas *c_mult = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDrawHistogram(h_mult, kBlue + 1);
      PlottingUtils::SaveFigure(c_mult, "multiplicity", subdir,
                                PlotSaveOptions::kLOG);
      delete c_mult;
    }

    output_file->Close();
    input_file->Close();

    if (emptyChannelMapEvents != 0)
      std::cout << "Observed " << emptyChannelMapEvents
                << " events with empty entry in channel map." << std::endl;

    std::cout << "Grid hits total: " << grid_hits_total << std::endl;
    std::cout << "Cathode hits total: " << cathode_hits_total
              << " (cathode/grid = "
              << (grid_hits_total > 0
                      ? Double_t(cathode_hits_total) / grid_hits_total
                      : 0.0)
              << ")" << std::endl;
    std::cout << "Total events: " << total_events << std::endl;
    if (total_events > 0) {
      std::cout << "Events with cathode hit: " << events_with_cathode << " ("
                << (100.0 * events_with_cathode / total_events) << "%)"
                << std::endl;
    }
    std::cout << "Complete events: " << complete_events << " ("
              << (100.0 * complete_events / total_events) << "%)" << std::endl;
    std::cout << "Incomplete events: " << incomplete_events << " ("
              << (100.0 * incomplete_events / total_events) << "%)"
              << std::endl;
    if (Constants::REJECT_FLAGGED_EVENTS) {
      Int_t stored = complete_events - complete_rejected;
      std::cout << "Stored events (REJECT_FLAGGED_EVENTS=true): " << stored
                << " (" << (complete_events > 0
                                ? 100.0 * stored / complete_events
                                : 0.0)
                << "% of complete; " << complete_rejected << " rejected)"
                << std::endl;
    }

    if (complete_events > 0) {
      std::cout << "Complete events with rejection-quality flags:"
                << std::endl;
      std::cout << "  Fake events: " << complete_with_fake << " ("
                << (100.0 * complete_with_fake / complete_events) << "%)"
                << std::endl;
      std::cout << "  Saturated: " << complete_with_saturation << " ("
                << (100.0 * complete_with_saturation / complete_events) << "%)"
                << std::endl;
      std::cout << "  Pileup: " << complete_with_pileup << " ("
                << (100.0 * complete_with_pileup / complete_events) << "%)"
                << std::endl;
    }

    if (incomplete_events > 0) {
      std::cout << "Incomplete events with rejection-quality flags:"
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

void EventBuilderClassic() {
  Bool_t reprocess_initial = kTRUE;

  std::vector<TString> filenames, output_names, file_labels;

  std::vector<FileSpec> specs = BuildFileSpecs();
  for (Int_t k = 0; k < Int_t(specs.size()); k++) {
    filenames.push_back(SortedName(specs[k]));
    std::cout << "Processing file: " << filenames.back() << std::endl;
    output_names.push_back(EventsClassicName(specs[k]));
    file_labels.push_back(FileLabel(specs[k]));
  }

  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                TString(gSystem->pwd()) + "/plots",
                                TString(gSystem->pwd()) + "/root_files");
  BuildEvents(filenames, output_names, file_labels, reprocess_initial);
}
