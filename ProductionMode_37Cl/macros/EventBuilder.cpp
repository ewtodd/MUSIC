#include "EventBuilder.hpp"
#include "BinaryUtils.hpp"
#include "Constants.hpp"
#include "InitUtils.hpp"
#include "Pipeline.hpp"
#include "PlottingUtils.hpp"
#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TObject.h>
#include <TParameter.h>
#include <TString.h>
#include <TTree.h>
#include <iostream>
#include <mutex>
#include <vector>

struct EventCounters {
  Int_t total_events;
  Int_t complete_events;
  Int_t complete_with_fake;
  Int_t complete_with_saturation;
  Int_t complete_with_pileup;
  Int_t complete_rejected;
  Int_t incomplete_events;
  Int_t incomplete_with_fake;
  Int_t incomplete_with_saturation;
  Int_t incomplete_with_pileup;
  Long64_t events_with_cathode;
  Long64_t events_with_multi_cathode;
  Long64_t events_with_multi_anode_hit;
  Long64_t dropped_anode_hits_total;
  Long64_t dropped_cathode_hits_total;
};

static void FinalizeEvent(
    NearestGrid::EventState &e, NearestGrid::PerChannelData *pc,
    Bool_t write_per_channel, TTree *output_tree, Int_t *leftdE_branch,
    Int_t *rightdE_branch, Int_t *totaldE_branch,
    ULong64_t *all_timestamps_branch, UInt_t *all_flags_branch,
    Int_t *hits_branch, Int_t &cathode_branch, Int_t &grid_branch,
    UInt_t &flags_or_branch, Bool_t &is_complete_branch, TH2F *h_music,
    TH2F *h_music_clean, TH2F *h_music_flagged, TH1F *h_mult,
    TH2F *h2_totalE_vs_stripE[18], TH2F *h2_totalE_vs_short[18],
    TH2F *h2_totalE_vs_cathode, EventCounters &c) {
  for (Int_t s = 1; s < 17; s++)
    e.totaldE[s] = e.leftdE[s] + e.rightdE[s];

  if (e.had_cathode)
    c.events_with_cathode++;
  if (e.hits[34] > 1) {
    c.events_with_multi_cathode++;
    c.dropped_cathode_hits_total += (e.hits[34] - 1);
  }
  Bool_t any_anode_multi = kFALSE;
  for (Int_t k = 0; k < 34; k++) {
    if (e.hits[k] > 1) {
      any_anode_multi = kTRUE;
      c.dropped_anode_hits_total += (e.hits[k] - 1);
    }
  }
  if (any_anode_multi)
    c.events_with_multi_anode_hit++;

  is_complete_branch = NearestGrid::CheckEventComplete(e);
  c.total_events++;

  Bool_t has_fake = kFALSE, has_saturation = kFALSE, has_pileup = kFALSE;
  NearestGrid::GetFlagSummary(e, has_fake, has_saturation, has_pileup);
  Bool_t has_any_flag = has_fake || has_saturation || has_pileup;

  if (is_complete_branch) {
    Bool_t reject = Constants::REJECT_FLAGGED_EVENTS && has_any_flag;
    if (!reject) {
      for (Int_t k = 0; k < 18; k++) {
        leftdE_branch[k] = e.leftdE[k];
        rightdE_branch[k] = e.rightdE[k];
        totaldE_branch[k] = e.totaldE[k];
      }
      for (Int_t k = 0; k < 36; k++)
        hits_branch[k] = e.hits[k];
      if (write_per_channel && pc) {
        for (Int_t k = 0; k < 36; k++) {
          all_timestamps_branch[k] = pc->timestamps[k];
          all_flags_branch[k] = pc->flags[k];
        }
      }
      cathode_branch = e.cathode;
      grid_branch = e.grid;
      flags_or_branch = e.flags_or;
      output_tree->Fill();

      Double_t event_total = 0.0;
      for (Int_t s = 0; s < 18; s++)
        event_total += Double_t(e.totaldE[s]);
      for (Int_t s = 0; s < 18; s++) {
        h_music->Fill(Double_t(s), Double_t(e.totaldE[s]));
        if (has_any_flag)
          h_music_flagged->Fill(Double_t(s), Double_t(e.totaldE[s]));
        else
          h_music_clean->Fill(Double_t(s), Double_t(e.totaldE[s]));
        h2_totalE_vs_stripE[s]->Fill(Double_t(e.totaldE[s]), event_total);
      }
      for (Int_t s = 1; s <= 16; s++) {
        // Short side: opposite of LongAnodeSide. Odd strips -> R short,
        // even strips -> L short.
        Int_t short_val = (s % 2 == 1) ? e.rightdE[s] : e.leftdE[s];
        if (short_val > 0)
          h2_totalE_vs_short[s]->Fill(Double_t(short_val), event_total);
      }
      if (e.had_cathode)
        h2_totalE_vs_cathode->Fill(Double_t(e.cathode), event_total);
      Int_t mult = 0;
      for (Int_t k = 0; k < 36; k++)
        mult += e.hits[k];
      h_mult->Fill(Double_t(mult));
    } else {
      c.complete_rejected++;
    }
    c.complete_events++;
    if (has_fake)
      c.complete_with_fake++;
    if (has_saturation)
      c.complete_with_saturation++;
    if (has_pileup)
      c.complete_with_pileup++;
  } else {
    c.incomplete_events++;
    if (has_fake)
      c.incomplete_with_fake++;
    if (has_saturation)
      c.incomplete_with_saturation++;
    if (has_pileup)
      c.incomplete_with_pileup++;
  }
}

Bool_t BuildEventsFromSortedHits(const std::vector<RawHit> &hits,
                                 const NearestGrid::SlotMap &slot_map,
                                 const TString &output_name,
                                 const TString &file_label) {
  TString output_filepath = output_name + ".root";
  TFile *output_file = IO::OpenForWriting(output_filepath);
  if (!output_file || output_file->IsZombie()) {
    std::cerr << "Error opening output: " << output_filepath << std::endl;
    return kFALSE;
  }

  Int_t leftdE[18], rightdE[18], totaldE[18];
  ULong64_t all_timestamps[36];
  UInt_t all_flags[36];
  Int_t hits_arr[36];
  Int_t cathode, grid;
  UInt_t flags_or;
  Bool_t is_complete;

  TTree *output_tree = new TTree("events", "MUSIC events");
  output_tree->Branch("LeftdE", leftdE, "LeftdE[18]/I");
  output_tree->Branch("RightdE", rightdE, "RightdE[18]/I");
  output_tree->Branch("TotaldE", totaldE, "TotaldE[18]/I");
  Bool_t write_per_channel = Constants::SAVE_PER_CHANNEL_TIMESTAMPS_FLAGS;
  Bool_t use_per_channel = Constants::USE_NEAREST_TO_GRID || write_per_channel;
  if (write_per_channel) {
    output_tree->Branch("AllTimestamps", all_timestamps, "AllTimestamps[36]/l");
    output_tree->Branch("AllFlags", all_flags, "AllFlags[36]/i");
  }
  output_tree->Branch("Hits", hits_arr, "Hits[36]/I");
  output_tree->Branch("Cathode", &cathode, "Cathode/I");
  output_tree->Branch("Grid", &grid, "Grid/I");
  output_tree->Branch("FlagsOR", &flags_or, "FlagsOR/i");
  output_tree->Branch("IsComplete", &is_complete, "IsComplete/O");

  // Large baskets + disable auto-flush so ZSTD compresses in big chunks
  // instead of many small basket flushes during Fill().
  output_tree->SetAutoFlush(0);
  for (Int_t b = 0; b < output_tree->GetListOfBranches()->GetEntries(); b++) {
    TBranch *branch =
        static_cast<TBranch *>(output_tree->GetListOfBranches()->At(b));
    branch->SetBasketSize(128 * 1024 * 1024);
  }

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

  TH2F *h2_totalE_vs_stripE[18];
  for (Int_t s = 0; s < 18; s++) {
    // Strip 0 dynamic range exceeds the L/R-tuned default; widen so its
    // beam + pileup peaks fit on the calibration h2.
    Double_t strip_e_max = (s == 0) ? 10000.0 : Constants::STRIP_E_MAX_ADC;
    h2_totalE_vs_stripE[s] =
        new TH2F(Form("h2_totalE_vs_stripE_s%d", s),
                 Form(";Strip %d #DeltaE [ADC];Total #DeltaE [ADC]", s), 200,
                 Constants::STRIP_E_MIN_ADC, strip_e_max, 400,
                 Constants::TOTAL_E_MIN_ADC, Constants::TOTAL_E_MAX_ADC);
  }
  TH2F *h2_totalE_vs_cathode =
      new TH2F("h2_totalE_vs_cathode",
               ";Cathode #DeltaE [ADC];Total #DeltaE [ADC]", 200, 0.0, 16384.0,
               400, Constants::TOTAL_E_MIN_ADC, Constants::TOTAL_E_MAX_ADC);
  // Short-side per-strip h2s for split strips 1..16. Indices 0 and 17 unused.
  TH2F *h2_totalE_vs_short[18] = {nullptr};
  for (Int_t s = 1; s <= 16; s++) {
    h2_totalE_vs_short[s] = new TH2F(
        Form("h2_totalE_vs_short_s%d", s),
        Form(";Strip %d short-side #DeltaE [ADC];Total #DeltaE [ADC]", s), 200,
        Constants::STRIP_E_MIN_ADC, Constants::STRIP_E_MAX_ADC, 400,
        Constants::TOTAL_E_MIN_ADC, Constants::TOTAL_E_MAX_ADC);
  }

  Long64_t n_entries = Long64_t(hits.size());

  struct PendingHit {
    Int_t slot;
    UShort_t energy;
    ULong64_t timestamp;
    UInt_t flags;
  };

  NearestGrid::EventState cur_event;
  NearestGrid::PerChannelData cur_per_channel;
  ULong64_t cur_grid_ts = 0;
  Bool_t have_cur = kFALSE;

  std::vector<PendingHit> pending;
  pending.reserve(4096);

  EventCounters cnt;
  cnt.total_events = 0;
  cnt.complete_events = 0;
  cnt.complete_with_fake = 0;
  cnt.complete_with_saturation = 0;
  cnt.complete_with_pileup = 0;
  cnt.complete_rejected = 0;
  cnt.incomplete_events = 0;
  cnt.incomplete_with_fake = 0;
  cnt.incomplete_with_saturation = 0;
  cnt.incomplete_with_pileup = 0;
  cnt.events_with_cathode = 0;
  cnt.events_with_multi_cathode = 0;
  cnt.events_with_multi_anode_hit = 0;
  cnt.dropped_anode_hits_total = 0;
  cnt.dropped_cathode_hits_total = 0;
  Int_t &total_events = cnt.total_events;
  Int_t &complete_events = cnt.complete_events;
  Int_t &complete_with_fake = cnt.complete_with_fake;
  Int_t &complete_with_saturation = cnt.complete_with_saturation;
  Int_t &complete_with_pileup = cnt.complete_with_pileup;
  Int_t &complete_rejected = cnt.complete_rejected;
  Int_t &incomplete_events = cnt.incomplete_events;
  Int_t &incomplete_with_fake = cnt.incomplete_with_fake;
  Int_t &incomplete_with_saturation = cnt.incomplete_with_saturation;
  Int_t &incomplete_with_pileup = cnt.incomplete_with_pileup;
  Long64_t &events_with_cathode = cnt.events_with_cathode;
  Long64_t &events_with_multi_cathode = cnt.events_with_multi_cathode;
  Long64_t &events_with_multi_anode_hit = cnt.events_with_multi_anode_hit;
  Long64_t &dropped_anode_hits_total = cnt.dropped_anode_hits_total;
  Long64_t &dropped_cathode_hits_total = cnt.dropped_cathode_hits_total;

  Int_t n_grids = 0;
  ULong64_t first_grid_ts = 0;
  ULong64_t last_grid_ts = 0;
  Int_t emptyChannelMapEvents = 0;
  Long64_t cathode_hits_total = 0;

  std::cout << "[" << file_label << "] Streaming pass over " << n_entries
            << " sorted hits..." << std::endl;
  for (Int_t i = 0; i < Int_t(n_entries); i++) {
    const RawHit &h = hits[i];

    if (h.board >= Constants::N_BOARDS || h.channel >= Constants::N_CHANNELS) {
      emptyChannelMapEvents++;
      continue;
    }
    Int_t slot = slot_map[h.board * Constants::N_CHANNELS + h.channel];
    if (slot < 0) {
      emptyChannelMapEvents++;
      continue;
    }

    if (slot == NearestGrid::kSlotGrid) {
      if (have_cur) {
        // Snapshot the current grid as prev_event so we can split-flush
        // pending and finalize it without touching the new grid's state.
        NearestGrid::EventState prev_event = cur_event;
        NearestGrid::PerChannelData prev_per_channel = cur_per_channel;
        ULong64_t prev_grid_ts = cur_grid_ts;

        // Initialise cur from the new grid hit.
        NearestGrid::InitEventState(cur_event, h.energy, h.flags);
        if (use_per_channel)
          NearestGrid::InitPerChannelData(cur_per_channel, h.timestamp,
                                          h.flags);
        cur_grid_ts = h.timestamp;

        // Split pending at the midpoint between prev and the new cur.
        ULong64_t midpoint = (prev_grid_ts + cur_grid_ts) / 2;
        NearestGrid::PerChannelData *pc_prev =
            use_per_channel ? &prev_per_channel : nullptr;
        NearestGrid::PerChannelData *pc_cur =
            use_per_channel ? &cur_per_channel : nullptr;
        Int_t split = 0;
        while (split < Int_t(pending.size()) &&
               pending[split].timestamp < midpoint)
          split++;
        for (Int_t p = 0; p < split; p++) {
          NearestGrid::AssignHit(prev_event, pc_prev, prev_grid_ts,
                                 pending[p].slot, pending[p].energy,
                                 pending[p].timestamp, pending[p].flags);
        }
        for (Int_t p = split; p < Int_t(pending.size()); p++) {
          NearestGrid::AssignHit(cur_event, pc_cur, cur_grid_ts,
                                 pending[p].slot, pending[p].energy,
                                 pending[p].timestamp, pending[p].flags);
        }
        pending.clear();

        // prev_event is now complete: finalize, Fill(), and discard.
        FinalizeEvent(prev_event, pc_prev, write_per_channel, output_tree,
                      leftdE, rightdE, totaldE, all_timestamps, all_flags,
                      hits_arr, cathode, grid, flags_or, is_complete, h_music,
                      h_music_clean, h_music_flagged, h_mult,
                      h2_totalE_vs_stripE, h2_totalE_vs_short,
                      h2_totalE_vs_cathode, cnt);
      } else {
        // First grid: initialise cur and flush all pending into it. Pending
        // here only contains hits with timestamps before the first grid;
        // they all belong to the first grid.
        NearestGrid::InitEventState(cur_event, h.energy, h.flags);
        if (use_per_channel)
          NearestGrid::InitPerChannelData(cur_per_channel, h.timestamp,
                                          h.flags);
        cur_grid_ts = h.timestamp;
        have_cur = kTRUE;

        NearestGrid::PerChannelData *pc_cur =
            use_per_channel ? &cur_per_channel : nullptr;
        for (Int_t p = 0; p < Int_t(pending.size()); p++) {
          NearestGrid::AssignHit(cur_event, pc_cur, cur_grid_ts,
                                 pending[p].slot, pending[p].energy,
                                 pending[p].timestamp, pending[p].flags);
        }
        pending.clear();
      }
      if (n_grids == 0)
        first_grid_ts = h.timestamp;
      last_grid_ts = h.timestamp;
      n_grids++;
    } else {
      if (slot == NearestGrid::kSlotCathode)
        cathode_hits_total++;
      PendingHit ph;
      ph.slot = slot;
      ph.energy = h.energy;
      ph.timestamp = h.timestamp;
      ph.flags = h.flags;
      pending.push_back(ph);
    }

    if (i % 10000000 == 0)
      std::cout << "  Stream progress: " << i << "/" << n_entries << std::endl;
  }

  if (have_cur) {
    NearestGrid::PerChannelData *pc_cur =
        use_per_channel ? &cur_per_channel : nullptr;
    for (Int_t p = 0; p < Int_t(pending.size()); p++) {
      NearestGrid::AssignHit(cur_event, pc_cur, cur_grid_ts, pending[p].slot,
                             pending[p].energy, pending[p].timestamp,
                             pending[p].flags);
    }
    pending.clear();
    FinalizeEvent(cur_event, pc_cur, write_per_channel, output_tree, leftdE,
                  rightdE, totaldE, all_timestamps, all_flags, hits_arr,
                  cathode, grid, flags_or, is_complete, h_music, h_music_clean,
                  h_music_flagged, h_mult, h2_totalE_vs_stripE,
                  h2_totalE_vs_short, h2_totalE_vs_cathode, cnt);
  }

  std::cout << "Found " << n_grids << " grid hits." << std::endl;

  if (n_grids == 0) {
    std::cerr << "No grid hits in file, skipping." << std::endl;
    output_file->Close();
    delete output_file;
    delete h_music;
    delete h_music_clean;
    delete h_music_flagged;
    delete h_mult;
    for (Int_t s = 0; s < 18; s++)
      delete h2_totalE_vs_stripE[s];
    for (Int_t s = 1; s <= 16; s++)
      delete h2_totalE_vs_short[s];
    delete h2_totalE_vs_cathode;
    return kFALSE;
  }

  Double_t span_s = (last_grid_ts > first_grid_ts)
                        ? Double_t(last_grid_ts - first_grid_ts) / 1e12
                        : 0.0;
  Double_t grid_rate_hz = (span_s > 0.0) ? Double_t(n_grids) / span_s : 0.0;

  output_file->cd();
  TParameter<Double_t>("grid_rate_hz", grid_rate_hz).Write();
  output_tree->Write("events", TObject::kOverwrite);
  h_mult->Write("", TObject::kOverwrite);

  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    TString subdir = "events_nearest/" + file_label;

    TCanvas *c_music = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(h_music, c_music);
    h_music->GetYaxis()->SetTitleOffset(1.4);
    c_music->SetLeftMargin(0.18);
    if (Constants::SAVE_PLOTS)
      PlottingUtils::SaveFigure(c_music, "music_strip_energies", subdir,
                                PlotSaveOptions::kLINEAR);
    output_file->cd();
    c_music->Write(h_music->GetName(), TObject::kOverwrite);
    delete c_music;

    TCanvas *c_music_clean = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(h_music_clean, c_music_clean);
    h_music_clean->GetYaxis()->SetTitleOffset(1.4);
    c_music_clean->SetLeftMargin(0.18);
    if (Constants::SAVE_PLOTS)
      PlottingUtils::SaveFigure(c_music_clean, "music_strip_energies_clean",
                                subdir, PlotSaveOptions::kLINEAR);
    output_file->cd();
    c_music_clean->Write(h_music_clean->GetName(), TObject::kOverwrite);
    delete c_music_clean;

    TCanvas *c_music_flagged = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(h_music_flagged,
                                               c_music_flagged);
    h_music_flagged->GetYaxis()->SetTitleOffset(1.4);
    c_music_flagged->SetLeftMargin(0.18);
    if (Constants::SAVE_PLOTS)
      PlottingUtils::SaveFigure(c_music_flagged, "music_strip_energies_flagged",
                                subdir, PlotSaveOptions::kLINEAR);
    output_file->cd();
    c_music_flagged->Write(h_music_flagged->GetName(), TObject::kOverwrite);
    delete c_music_flagged;

    TCanvas *c_mult = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDrawHistogram(h_mult, kBlue + 1);
    if (Constants::SAVE_PLOTS)
      PlottingUtils::SaveFigure(c_mult, "multiplicity", subdir,
                                PlotSaveOptions::kLOG);
    delete c_mult;

    TString trace_summary_subdir = "trace_summary/" + file_label;
    for (Int_t s = 0; s < 18; s++) {
      TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDraw2DHistogram(h2_totalE_vs_stripE[s], c);
      h2_totalE_vs_stripE[s]->GetYaxis()->SetTitleOffset(1.3);
      if (Constants::SAVE_PLOTS)
        PlottingUtils::SaveFigure(c, Form("totalE_vs_stripE_s%d", s),
                                  trace_summary_subdir,
                                  PlotSaveOptions::kLINEAR);
      output_file->cd();
      c->Write(h2_totalE_vs_stripE[s]->GetName(), TObject::kOverwrite);
      delete c;
    }

    TCanvas *c_cath = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(h2_totalE_vs_cathode, c_cath);
    h2_totalE_vs_cathode->GetYaxis()->SetTitleOffset(1.3);
    if (Constants::SAVE_PLOTS)
      PlottingUtils::SaveFigure(c_cath, "totalE_vs_cathode",
                                trace_summary_subdir, PlotSaveOptions::kLINEAR);
    output_file->cd();
    c_cath->Write(h2_totalE_vs_cathode->GetName(), TObject::kOverwrite);
    delete c_cath;

    for (Int_t s = 1; s <= 16; s++) {
      TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDraw2DHistogram(h2_totalE_vs_short[s], c);
      h2_totalE_vs_short[s]->GetYaxis()->SetTitleOffset(1.3);
      if (Constants::SAVE_PLOTS)
        PlottingUtils::SaveFigure(c, Form("totalE_vs_short_s%d", s),
                                  trace_summary_subdir,
                                  PlotSaveOptions::kLINEAR);
      output_file->cd();
      c->Write(h2_totalE_vs_short[s]->GetName(), TObject::kOverwrite);
      delete c;
    }
  }

  for (Int_t s = 0; s < 18; s++)
    delete h2_totalE_vs_stripE[s];
  for (Int_t s = 1; s <= 16; s++)
    delete h2_totalE_vs_short[s];
  delete h2_totalE_vs_cathode;

  output_file->Close();
  delete output_file;

  if (emptyChannelMapEvents != 0)
    std::cout << "Observed " << emptyChannelMapEvents
              << " hits with empty entry in channel map." << std::endl;

  std::cout << "Grid hits total: " << n_grids << std::endl;
  std::cout << "Grid trigger rate: " << grid_rate_hz << " Hz over " << span_s
            << " s" << std::endl;
  std::cout << "Cathode hits total: " << cathode_hits_total
            << " (cathode/grid = "
            << (n_grids > 0 ? Double_t(cathode_hits_total) / n_grids : 0.0)
            << ")" << std::endl;
  std::cout << "Total events: " << total_events << std::endl;
  if (total_events > 0) {
    std::cout << "Events with cathode hit: " << events_with_cathode << " ("
              << (100.0 * events_with_cathode / total_events) << "%)"
              << std::endl;
    std::cout << "Events with multiple cathode hits: "
              << events_with_multi_cathode << " ("
              << (100.0 * events_with_multi_cathode / total_events)
              << "% of all, "
              << (events_with_cathode > 0
                      ? 100.0 * events_with_multi_cathode / events_with_cathode
                      : 0.0)
              << "% of cathode events)" << std::endl;
    std::cout << "Events with multi-hit on any anode: "
              << events_with_multi_anode_hit << " ("
              << (100.0 * events_with_multi_anode_hit / total_events) << "%)"
              << std::endl;
    std::cout << "Dropped hits (kept closest-to-grid): anode="
              << dropped_anode_hits_total
              << ", cathode=" << dropped_cathode_hits_total << std::endl;
  }
  std::cout << "Complete events: " << complete_events << " ("
            << (100.0 * complete_events / total_events) << "%)" << std::endl;
  std::cout << "Incomplete events: " << incomplete_events << " ("
            << (100.0 * incomplete_events / total_events) << "%)" << std::endl;
  if (Constants::REJECT_FLAGGED_EVENTS) {
    Int_t stored = complete_events - complete_rejected;
    std::cout << "Stored events (REJECT_FLAGGED_EVENTS=true): " << stored
              << " ("
              << (complete_events > 0 ? 100.0 * stored / complete_events : 0.0)
              << "% of complete; " << complete_rejected << " rejected)"
              << std::endl;
  }
  if (complete_events > 0) {
    std::cout << "Complete events with rejection-quality flags:" << std::endl;
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
    std::cout << "Incomplete events with rejection-quality flags:" << std::endl;
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
  return kTRUE;
}
