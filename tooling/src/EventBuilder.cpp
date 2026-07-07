#include "EventBuilder.hpp"

void EventBuilder::ResetEventState(EventState &e) {
  for (Int_t k = 0; k < 18; k++) {
    e.leftdE[k] = 0;
    e.rightdE[k] = 0;
    e.totaldE[k] = 0;
  }
  for (Int_t k = 0; k < Constants::N_ARR_SLOTS; k++)
    e.hits[k] = 0;
  e.cathode = -1;
  e.grid = 0;
  e.flags_or = 0;
  e.had_cathode = kFALSE;
}

void EventBuilder::ResetPerChannelData(PerChannelData &p) {
  for (Int_t k = 0; k < Constants::N_ARR_SLOTS; k++) {
    p.timestamps[k] = 0;
    p.energies[k] = 0;
    p.flags[k] = 0;
  }
}

Bool_t EventBuilder::ShouldKeepHit(ULong64_t cand_ts, ULong64_t prev_ts,
                                   UShort_t cand_energy, UShort_t prev_energy,
                                   ULong64_t ref_ts, DedupStrategy strategy) {
  if (strategy == kSMALLEST_ENERGY)
    return cand_energy < prev_energy;
  if (strategy == kLARGEST_ENERGY)
    return cand_energy > prev_energy;
  if (strategy == kEARLIEST_TIMESTAMP)
    return cand_ts < prev_ts;
  if (strategy == kLATEST_TIMESTAMP)
    return cand_ts > prev_ts;

  ULong64_t cand_dt = (cand_ts > ref_ts) ? cand_ts - ref_ts : ref_ts - cand_ts;
  ULong64_t prev_dt = (prev_ts > ref_ts) ? prev_ts - ref_ts : ref_ts - prev_ts;
  return cand_dt < prev_dt;
}

EventBuilder::SlotMap EventBuilder::BuildSlotMap() {
  SlotMap sm;
  sm.resize(Constants::cfg.N_BOARDS * Constants::cfg.N_CHANNELS, -1);
  for (std::map<std::pair<Int_t, Int_t>, TString>::const_iterator it =
           Constants::ActiveChannelMap().begin();
       it != Constants::ActiveChannelMap().end(); ++it) {
    Int_t board = it->first.first;
    Int_t channel = it->first.second;
    if (board < 0 || board >= Constants::cfg.N_BOARDS)
      continue;
    if (channel < 0 || channel >= Constants::cfg.N_CHANNELS)
      continue;
    Int_t idx = board * Constants::cfg.N_CHANNELS + channel;
    const TString &name = it->second;
    if (name == "Strip0") {
      sm[idx] = Constants::ARR_SLOT_STRIP_0;
    } else if (name == "Strip17") {
      sm[idx] = Constants::ARR_SLOT_STRIP_17;
    } else if (name == "Cathode") {
      sm[idx] = Constants::ARR_SLOT_CATHODE;
    } else if (name == "Grid") {
      sm[idx] = Constants::ARR_SLOT_GRID;
    } else if (name.Length() >= 2 && (name[0] == 'L' || name[0] == 'R')) {
      TString numstr = name(1, name.Length() - 1);
      if (numstr.IsDigit()) {
        Int_t s = numstr.Atoi();
        if (s >= 1 && s <= 16)
          sm[idx] = (name[0] == 'L') ? s : (s + 16);
      }
    }
  }
  return sm;
}

void EventBuilder::AssignHit(EventState &e, PerChannelData *pc,
                             ULong64_t ref_ts, Int_t slot, UShort_t energy,
                             ULong64_t timestamp, UInt_t flags,
                             DedupStrategy strategy) {
  // Slots are 0..ARR_SLOT_GRID by construction (BuildSlotMap); guard
  // defensively so an out-of-range slot can never index the arrays below.
  if (slot < 0 || slot >= Constants::N_ARR_SLOTS)
    return;
  e.flags_or |= flags;

  Bool_t wins;
  if (e.hits[slot] == 0) {
    wins = kTRUE;
  } else if (pc) {
    wins = ShouldKeepHit(timestamp, pc->timestamps[slot], energy,
                         pc->energies[slot], ref_ts, strategy);
  } else {
    wins = kFALSE;
  }

  if (wins) {
    if (slot == Constants::ARR_SLOT_STRIP_0)
      e.totaldE[0] = energy;
    else if (slot <= 16)
      e.leftdE[slot] = energy;
    else if (slot <= 32)
      e.rightdE[slot - 16] = energy;
    else if (slot == Constants::ARR_SLOT_STRIP_17)
      e.totaldE[17] = energy;
    else if (slot == Constants::ARR_SLOT_CATHODE)
      e.cathode = energy;
    else if (slot == Constants::ARR_SLOT_GRID)
      e.grid = energy;
    if (pc) {
      pc->timestamps[slot] = timestamp;
      pc->energies[slot] = energy;
      pc->flags[slot] = flags;
    }
  }
  if (slot == Constants::ARR_SLOT_CATHODE)
    e.had_cathode = kTRUE;
  e.hits[slot]++;
}

Bool_t EventBuilder::CheckEventComplete(const EventState &e) {
  if (e.totaldE[0] == 0)
    return kFALSE;
  for (Int_t strip = 1; strip <= 14; strip += 2) {
    if (e.leftdE[strip] == 0)
      return kFALSE;
  }
  for (Int_t strip = 2; strip <= 15; strip += 2) {
    if (e.rightdE[strip] == 0)
      return kFALSE;
  }
  return kTRUE;
}

void EventBuilder::GetFlagSummary(const EventState &e, Bool_t &has_fake,
                                  Bool_t &has_saturation, Bool_t &has_pileup) {
  has_fake = (e.flags_or & CoMPASSData::FAKE_EVENT) != 0;
  has_saturation = ((e.flags_or & CoMPASSData::SATURATION_IN_GATE) ||
                    (e.flags_or & CoMPASSData::INPUT_SATURATING)) != 0;
  has_pileup = (e.flags_or & CoMPASSData::PILEUP) != 0;
}

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
  // Per-required-channel miss counts over all events: how often each
  // completeness-required long strip / guard strip was zero. miss_long is
  // indexed by strip 1..16 (counts L==0 for odd strips, R==0 for even).
  Long64_t miss_long[17];
  Long64_t miss_strip0;
  Long64_t miss_strip17;
};

void FinalizeEvent(EventState &e, PerChannelData *pc, TTree *output_tree,
                   UShort_t *left_0_17_branch, UShort_t *rightdE_branch,
                   UShort_t *hits_branch, Short_t &cathode_branch,
                   Short_t &grid_branch, UInt_t &flags_or_branch, TH2F *h_music,
                   TH2F *h_music_clean, TH2F *h_music_flagged, TH1F *h_mult,
                   TH2F *h2_totalE_vs_stripE[18], TH2F *h2_totalE_vs_L[18],
                   TH2F *h2_totalE_vs_R[18], TH2F *h2_totalE_vs_cathode,
                   TH1F *h1_stripE[18], TH1F *h1_L[18], TH1F *h1_R[18],
                   TH1F *h1_cathode, TH2F *h2_long_back_vs_front,
                   EventCounters &c) {
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

  Bool_t is_complete = EventBuilder::CheckEventComplete(e);
  c.total_events++;

  // Per-required-channel miss tally (over all events): which completeness
  // channel was zero. Long side = L for odd strips, R for even strips.
  if (e.totaldE[0] == 0)
    c.miss_strip0++;
  if (e.totaldE[17] == 0)
    c.miss_strip17++;
  for (Int_t s = 1; s <= 16; s++) {
    Int_t v = (s % 2 == 1) ? e.leftdE[s] : e.rightdE[s];
    if (v == 0)
      c.miss_long[s]++;
  }

  Bool_t has_fake = kFALSE, has_saturation = kFALSE, has_pileup = kFALSE;
  EventBuilder::GetFlagSummary(e, has_fake, has_saturation, has_pileup);
  Bool_t has_any_flag = has_fake || has_saturation || has_pileup;

  if (is_complete) {
    Bool_t reject = Constants::cfg.REJECT_FLAGGED_EVENTS && has_any_flag;
    if (!reject) {
      for (Int_t k = 0; k < 18; k++) {
        left_0_17_branch[k] = UShort_t(e.leftdE[k]);
        rightdE_branch[k] = UShort_t(e.rightdE[k]);
      }
      // Guard strips are single-ended: their energy lives in totaldE[0]/[17]
      // (leftdE/rightdE are 0 there). Park them in the left array at 0/17 so
      // the recompute total[s] = left_0_17[s] + rightdE[s] reproduces them too.
      left_0_17_branch[0] = UShort_t(e.totaldE[0]);
      left_0_17_branch[17] = UShort_t(e.totaldE[17]);
      for (Int_t k = 0; k < Constants::N_ARR_SLOTS; k++)
        hits_branch[k] = UShort_t(e.hits[k]);
      cathode_branch = Short_t(e.cathode);
      grid_branch = Short_t(e.grid);
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
        h1_stripE[s]->Fill(Double_t(e.totaldE[s]));
      }
      for (Int_t s = 1; s <= 16; s++) {
        if (e.leftdE[s] > 0) {
          h2_totalE_vs_L[s]->Fill(Double_t(e.leftdE[s]), event_total);
          h1_L[s]->Fill(Double_t(e.leftdE[s]));
        }
        if (e.rightdE[s] > 0) {
          h2_totalE_vs_R[s]->Fill(Double_t(e.rightdE[s]), event_total);
          h1_R[s]->Fill(Double_t(e.rightdE[s]));
        }
      }
      if (e.had_cathode) {
        h2_totalE_vs_cathode->Fill(Double_t(e.cathode), event_total);
        h1_cathode->Fill(Double_t(e.cathode));
      }
      // Density diagnostic: front-long-L sum vs back-long-L sum, gated by
      // (count of LeftdE==0 across strips 1..16) < 8. Beam blob plus its
      // first pileup separate cleanly along the diagonal.
      Int_t n_left_zero = 0;
      for (Int_t s = 1; s <= 16; s++) {
        if (e.leftdE[s] == 0)
          n_left_zero++;
      }
      if (n_left_zero < 8) {
        Double_t front_long_l =
            Double_t(e.leftdE[1] + e.leftdE[3] + e.leftdE[5] + e.leftdE[7]);
        Double_t back_long_l =
            Double_t(e.leftdE[9] + e.leftdE[11] + e.leftdE[13] + e.leftdE[15]);
        h2_long_back_vs_front->Fill(back_long_l, front_long_l);
      }
      Int_t mult = 0;
      for (Int_t k = 0; k < Constants::N_ARR_SLOTS; k++)
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

Bool_t EventBuilder::BuildEventsFromSortedHits(const std::vector<RawHit> &hits,
                                               const SlotMap &slot_map,
                                               const TString &output_name,
                                               const TString &file_label) {
  TString output_filepath = output_name + ".root";
  TFile *output_file = IO::OpenForWriting(output_filepath);
  if (!output_file || output_file->IsZombie()) {
    std::cerr << "Error opening output: " << output_filepath << std::endl;
    return kFALSE;
  }

  // ADC energies are 14-bit unsigned at the source; store them as UShort_t, not
  // Int_t. Index 0/17 of left_0_17_dE hold the single-ended guard strips
  // (Strip0/Strip17), indices 1..16 the left ends of the split anodes; the
  // right ends live in rightdE (which is 0 at 0/17). The full per-strip deposit
  // is recomputed as left_0_17_dE[s] + rightdE[s] on read -- it holds for all
  // 18 indices since rightdE is 0 at the guards -- so no TotaldE branch is
  // stored.
  UShort_t left_0_17_dE[18], rightdE[18];
  UShort_t hits_arr[36];
  // 14-bit ADC (<=16383), so Short_t holds every value with room to spare while
  // preserving Cathode's -1 "no cathode hit" sentinel (Grid is non-negative but
  // shares the type for symmetry). Unsplit anode/guard values are non-negative,
  // hence the UShort_t arrays above.
  Short_t cathode, grid;
  UInt_t flags_or;

  TTree *output_tree = new TTree("events", "MUSIC events");
  output_tree->Branch("Left_0_17_dE", left_0_17_dE, "Left_0_17_dE[18]/s");
  output_tree->Branch("RightdE", rightdE, "RightdE[18]/s");
  output_tree->Branch("Hits", hits_arr, "Hits[36]/s");
  output_tree->Branch("Cathode", &cathode, "Cathode/S");
  output_tree->Branch("Grid", &grid, "Grid/S");
  output_tree->Branch("FlagsOR", &flags_or, "FlagsOR/i");

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
  TH1F *h1_stripE[18];
  for (Int_t s = 0; s < 18; s++) {
    // Strip 0 dynamic range exceeds the L/R-tuned default; widen so its
    // beam + pileup peaks fit on the calibration h2.
    Double_t strip_e_max = (s == 0) ? 10000.0 : Constants::cfg.STRIP_E_MAX_ADC;
    h2_totalE_vs_stripE[s] = new TH2F(
        Form("h2_totalE_vs_stripE_s%d", s),
        Form(";Strip %d #DeltaE [ADC];Total #DeltaE [ADC]", s), 200,
        Constants::cfg.STRIP_E_MIN_ADC, strip_e_max, 400,
        Constants::cfg.TOTAL_E_MIN_ADC, Constants::cfg.TOTAL_E_MAX_ADC);
    h1_stripE[s] = new TH1F(Form("h1_stripE_s%d", s),
                            Form(";Strip %d #DeltaE [ADC];Counts", s), 400,
                            Constants::cfg.STRIP_E_MIN_ADC, strip_e_max);
  }
  TH2F *h2_totalE_vs_cathode = new TH2F(
      "h2_totalE_vs_cathode", ";Cathode #DeltaE [ADC];Total #DeltaE [ADC]", 200,
      0.0, 16384.0, 400, Constants::cfg.TOTAL_E_MIN_ADC,
      Constants::cfg.TOTAL_E_MAX_ADC);
  TH1F *h1_cathode = new TH1F("h1_cathode", ";Cathode #DeltaE [ADC];Counts",
                              400, 0.0, 16384.0);
  // Per-side per-strip h2s/h1s for split strips 1..16. Indices 0 and 17
  // unused (guard strips are unsplit).
  TH2F *h2_totalE_vs_L[18] = {nullptr};
  TH2F *h2_totalE_vs_R[18] = {nullptr};
  TH1F *h1_L[18] = {nullptr};
  TH1F *h1_R[18] = {nullptr};
  for (Int_t s = 1; s <= 16; s++) {
    h2_totalE_vs_L[s] = new TH2F(
        Form("h2_totalE_vs_L_s%d", s),
        Form(";Strip %d L #DeltaE [ADC];Total #DeltaE [ADC]", s), 200,
        Constants::cfg.STRIP_E_MIN_ADC, Constants::cfg.STRIP_E_MAX_ADC, 400,
        Constants::cfg.TOTAL_E_MIN_ADC, Constants::cfg.TOTAL_E_MAX_ADC);
    h2_totalE_vs_R[s] = new TH2F(
        Form("h2_totalE_vs_R_s%d", s),
        Form(";Strip %d R #DeltaE [ADC];Total #DeltaE [ADC]", s), 200,
        Constants::cfg.STRIP_E_MIN_ADC, Constants::cfg.STRIP_E_MAX_ADC, 400,
        Constants::cfg.TOTAL_E_MIN_ADC, Constants::cfg.TOTAL_E_MAX_ADC);
    h1_L[s] = new TH1F(
        Form("h1_L_s%d", s), Form(";Strip %d L #DeltaE [ADC];Counts", s), 400,
        Constants::cfg.STRIP_E_MIN_ADC, Constants::cfg.STRIP_E_MAX_ADC);
    h1_R[s] = new TH1F(
        Form("h1_R_s%d", s), Form(";Strip %d R #DeltaE [ADC];Counts", s), 400,
        Constants::cfg.STRIP_E_MIN_ADC, Constants::cfg.STRIP_E_MAX_ADC);
  }
  // Density diagnostic for beam-event isolation. Axes are the long-L sums
  // across the back half (strips 9,11,13,15) and front half (strips
  // 1,3,5,7). With the per-event gate applied at fill time, the beam and
  // first-pileup blobs sit along the diagonal.
  TH2F *h2_long_back_vs_front =
      new TH2F("h2_long_back_vs_front",
               ";Long-L back sum (s9+s11+s13+s15) [ADC];"
               "Long-L front sum (s1+s3+s5+s7) [ADC]",
               400, 0.0, 10000.0, 400, 0.0, 10000.0);

  // Determine operating mode: reference-channel or blind time-window.
  Bool_t ref_mode = Constants::cfg.REFERENCE_CHANNEL != "NONE";

  // Resolve reference channel name to slot ID by scanning active channel map.
  Int_t ref_slot = -1;
  if (ref_mode) {
    for (std::map<std::pair<Int_t, Int_t>, TString>::const_iterator it =
             Constants::ActiveChannelMap().begin();
         it != Constants::ActiveChannelMap().end(); ++it) {
      if (it->second == Constants::cfg.REFERENCE_CHANNEL) {
        Int_t board = it->first.first;
        Int_t channel = it->first.second;
        if (board >= 0 && board < Constants::cfg.N_BOARDS && channel >= 0 &&
            channel < Constants::cfg.N_CHANNELS) {
          ref_slot = slot_map[board * Constants::cfg.N_CHANNELS + channel];
          break;
        }
      }
    }
    if (ref_slot < 0) {
      std::cerr << "FATAL: reference channel '"
                << Constants::cfg.REFERENCE_CHANNEL
                << "' not found in channel map. Cannot build events."
                << std::endl;
      output_file->Close();
      delete output_file;
      delete h_music;
      delete h_music_clean;
      delete h_music_flagged;
      delete h_mult;
      for (Int_t s = 0; s < 18; s++) {
        delete h2_totalE_vs_stripE[s];
        delete h1_stripE[s];
      }
      for (Int_t s = 1; s <= 16; s++) {
        delete h2_totalE_vs_L[s];
        delete h2_totalE_vs_R[s];
        delete h1_L[s];
        delete h1_R[s];
      }
      delete h2_totalE_vs_cathode;
      delete h1_cathode;
      delete h2_long_back_vs_front;
      return kFALSE;
    }
  }

  Long64_t n_entries = Long64_t(hits.size());

  struct PendingHit {
    Int_t slot;
    UShort_t energy;
    ULong64_t timestamp;
    UInt_t flags;
  };

  EventState cur_event;
  PerChannelData cur_per_channel;
  ULong64_t cur_ref_ts = 0;
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
  for (Int_t s = 0; s < 17; s++)
    cnt.miss_long[s] = 0;
  cnt.miss_strip0 = 0;
  cnt.miss_strip17 = 0;
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

  Int_t n_ref = 0;
  ULong64_t first_ref_ts = 0;
  ULong64_t last_ref_ts = 0;
  Int_t emptyChannelMapEvents = 0;
  Long64_t cathode_hits_total = 0;
  Long64_t dropped_outside_window = 0;

  ULong64_t window_ps = ULong64_t(Constants::cfg.EVENT_TIME_WINDOW_US * 1.0e6);
  DedupStrategy dedup_strat = Constants::cfg.DEDUP_STRATEGY;
  PerChannelData *pc_cur = &cur_per_channel;

  std::cout << "[" << file_label << "] Streaming pass over " << n_entries
            << " sorted hits (reference: " << Constants::cfg.REFERENCE_CHANNEL
            << ", window: " << Constants::cfg.EVENT_TIME_WINDOW_US << " us)"
            << std::endl;

  for (Int_t i = 0; i < Int_t(n_entries); i++) {
    const RawHit &h = hits[i];

    if (h.board >= Constants::cfg.N_BOARDS ||
        h.channel >= Constants::cfg.N_CHANNELS) {
      emptyChannelMapEvents++;
      continue;
    }
    Int_t slot = slot_map[h.board * Constants::cfg.N_CHANNELS + h.channel];
    if (slot < 0) {
      emptyChannelMapEvents++;
      continue;
    }

    if (slot == Constants::ARR_SLOT_CATHODE)
      cathode_hits_total++;

    if (ref_mode) {
      // --- Reference-channel mode ---
      // Reference hit seeds a new event; non-reference hits queue until
      // next reference hit (or end of stream).
      if (slot == ref_slot) {
        if (n_ref == 0)
          first_ref_ts = h.timestamp;
        last_ref_ts = h.timestamp;
        n_ref++;

        if (have_cur) {
          // Flush pending hits that fall within the current window.
          for (Int_t p = 0; p < Int_t(pending.size()); p++) {
            if (pending[p].timestamp - cur_ref_ts <= window_ps) {
              AssignHit(cur_event, pc_cur, cur_ref_ts, pending[p].slot,
                        pending[p].energy, pending[p].timestamp,
                        pending[p].flags, dedup_strat);
            } else {
              dropped_outside_window++;
            }
          }
          pending.clear();

          // Assign the reference hit itself.
          AssignHit(cur_event, pc_cur, cur_ref_ts, ref_slot, h.energy,
                    h.timestamp, h.flags, dedup_strat);

          // Finalize the completed event.
          FinalizeEvent(cur_event, pc_cur, output_tree, left_0_17_dE, rightdE,
                        hits_arr, cathode, grid, flags_or, h_music,
                        h_music_clean, h_music_flagged, h_mult,
                        h2_totalE_vs_stripE, h2_totalE_vs_L, h2_totalE_vs_R,
                        h2_totalE_vs_cathode, h1_stripE, h1_L, h1_R, h1_cathode,
                        h2_long_back_vs_front, cnt);
        }

        // Seed new event from this reference hit.
        ResetEventState(cur_event);
        ResetPerChannelData(cur_per_channel);
        cur_ref_ts = h.timestamp;
        have_cur = kTRUE;
      } else {
        if (have_cur && h.timestamp - cur_ref_ts <= window_ps) {
          PendingHit ph;
          ph.slot = slot;
          ph.energy = h.energy;
          ph.timestamp = h.timestamp;
          ph.flags = h.flags;
          pending.push_back(ph);
        } else if (have_cur) {
          dropped_outside_window++;
        }
      }
    } else {
      // --- Time-window mode (REFERENCE_CHANNEL == "NONE") ---
      // First unassigned hit opens a window; all hits within window_ps
      // belong to that event. Dedup still applies, ref_ts = window anchor.
      if (!have_cur) {
        ResetEventState(cur_event);
        ResetPerChannelData(cur_per_channel);
        cur_ref_ts = h.timestamp;
        have_cur = kTRUE;
        AssignHit(cur_event, pc_cur, cur_ref_ts, slot, h.energy, h.timestamp,
                  h.flags, dedup_strat);
      } else if (h.timestamp - cur_ref_ts <= window_ps) {
        AssignHit(cur_event, pc_cur, cur_ref_ts, slot, h.energy, h.timestamp,
                  h.flags, dedup_strat);
      } else {
        // Window exceeded — finalize and start new window.
        FinalizeEvent(cur_event, pc_cur, output_tree, left_0_17_dE, rightdE,
                      hits_arr, cathode, grid, flags_or, h_music, h_music_clean,
                      h_music_flagged, h_mult, h2_totalE_vs_stripE,
                      h2_totalE_vs_L, h2_totalE_vs_R, h2_totalE_vs_cathode,
                      h1_stripE, h1_L, h1_R, h1_cathode, h2_long_back_vs_front,
                      cnt);
        ResetEventState(cur_event);
        ResetPerChannelData(cur_per_channel);
        cur_ref_ts = h.timestamp;
        AssignHit(cur_event, pc_cur, cur_ref_ts, slot, h.energy, h.timestamp,
                  h.flags, dedup_strat);
      }
    }

    if (i % 10000000 == 0)
      std::cout << "  Stream progress: " << i << "/" << n_entries << std::endl;
  }

  // Finalize the last event.
  if (have_cur) {
    if (ref_mode) {
      for (Int_t p = 0; p < Int_t(pending.size()); p++) {
        if (pending[p].timestamp - cur_ref_ts <= window_ps) {
          AssignHit(cur_event, pc_cur, cur_ref_ts, pending[p].slot,
                    pending[p].energy, pending[p].timestamp, pending[p].flags,
                    dedup_strat);
        } else {
          dropped_outside_window++;
        }
      }
      pending.clear();
    }
    FinalizeEvent(cur_event, pc_cur, output_tree, left_0_17_dE, rightdE,
                  hits_arr, cathode, grid, flags_or, h_music, h_music_clean,
                  h_music_flagged, h_mult, h2_totalE_vs_stripE, h2_totalE_vs_L,
                  h2_totalE_vs_R, h2_totalE_vs_cathode, h1_stripE, h1_L, h1_R,
                  h1_cathode, h2_long_back_vs_front, cnt);
  }

  std::cout << "Found " << n_ref << " " << Constants::cfg.REFERENCE_CHANNEL
            << " hits." << std::endl;

  if (ref_mode && n_ref == 0) {
    std::cerr << "No " << Constants::cfg.REFERENCE_CHANNEL
              << " hits in file, skipping." << std::endl;
    output_file->Close();
    delete output_file;
    delete h_music;
    delete h_music_clean;
    delete h_music_flagged;
    delete h_mult;
    for (Int_t s = 0; s < 18; s++) {
      delete h2_totalE_vs_stripE[s];
      delete h1_stripE[s];
    }
    for (Int_t s = 1; s <= 16; s++) {
      delete h2_totalE_vs_L[s];
      delete h2_totalE_vs_R[s];
      delete h1_L[s];
      delete h1_R[s];
    }
    delete h2_totalE_vs_cathode;
    delete h1_cathode;
    delete h2_long_back_vs_front;
    return kFALSE;
  }

  Double_t span_s = (last_ref_ts > first_ref_ts)
                        ? Double_t(last_ref_ts - first_ref_ts) / 1e12
                        : 0.0;
  Double_t ref_rate_hz = (span_s > 0.0) ? Double_t(n_ref) / span_s : 0.0;

  output_file->cd();
  TParameter<Double_t>("grid_rate_hz", ref_rate_hz).Write();
  output_tree->Write("events", TObject::kOverwrite);
  h_mult->Write("", TObject::kOverwrite);

  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    TString subdir = "events_nearest/" + file_label;

    TCanvas *c_music = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(h_music, c_music);
    h_music->GetYaxis()->SetTitleOffset(1.4);
    c_music->SetLeftMargin(0.18);
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(c_music, "music_strip_energies", subdir,
                                PlotSaveOptions::kLINEAR);
    output_file->cd();
    c_music->Write(h_music->GetName(), TObject::kOverwrite);
    delete c_music;

    TCanvas *c_music_clean = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(h_music_clean, c_music_clean);
    h_music_clean->GetYaxis()->SetTitleOffset(1.4);
    c_music_clean->SetLeftMargin(0.18);
    if (Constants::cfg.SAVE_PLOTS)
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
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(c_music_flagged, "music_strip_energies_flagged",
                                subdir, PlotSaveOptions::kLINEAR);
    output_file->cd();
    c_music_flagged->Write(h_music_flagged->GetName(), TObject::kOverwrite);
    delete c_music_flagged;

    TCanvas *c_mult = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDrawHistogram(h_mult, kBlue + 1);
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(c_mult, "multiplicity", subdir,
                                PlotSaveOptions::kLOG);
    delete c_mult;

    TString trace_summary_subdir = "trace_summary/" + file_label;
    for (Int_t s = 0; s < 18; s++) {
      TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDraw2DHistogram(h2_totalE_vs_stripE[s], c);
      h2_totalE_vs_stripE[s]->GetYaxis()->SetTitleOffset(1.3);
      if (Constants::cfg.SAVE_PLOTS)
        PlottingUtils::SaveFigure(c, Form("totalE_vs_stripE_s%d", s),
                                  trace_summary_subdir,
                                  PlotSaveOptions::kLINEAR);
      output_file->cd();
      c->Write(h2_totalE_vs_stripE[s]->GetName(), TObject::kOverwrite);
      delete c;

      TCanvas *c1 = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDrawHistogram(h1_stripE[s], kBlue + 1);
      if (Constants::cfg.SAVE_PLOTS)
        PlottingUtils::SaveFigure(c1, Form("stripE_s%d", s),
                                  trace_summary_subdir, PlotSaveOptions::kLOG);
      output_file->cd();
      h1_stripE[s]->Write(h1_stripE[s]->GetName(), TObject::kOverwrite);
      delete c1;
    }

    TCanvas *c_cath = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(h2_totalE_vs_cathode, c_cath);
    h2_totalE_vs_cathode->GetYaxis()->SetTitleOffset(1.3);
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(c_cath, "totalE_vs_cathode",
                                trace_summary_subdir, PlotSaveOptions::kLINEAR);
    output_file->cd();
    c_cath->Write(h2_totalE_vs_cathode->GetName(), TObject::kOverwrite);
    delete c_cath;

    TCanvas *c_cath1 = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDrawHistogram(h1_cathode, kBlue + 1);
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(c_cath1, "cathode", trace_summary_subdir,
                                PlotSaveOptions::kLOG);
    output_file->cd();
    h1_cathode->Write(h1_cathode->GetName(), TObject::kOverwrite);
    delete c_cath1;

    for (Int_t s = 1; s <= 16; s++) {
      TCanvas *cL = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDraw2DHistogram(h2_totalE_vs_L[s], cL);
      h2_totalE_vs_L[s]->GetYaxis()->SetTitleOffset(1.3);
      if (Constants::cfg.SAVE_PLOTS)
        PlottingUtils::SaveFigure(cL, Form("totalE_vs_L_s%d", s),
                                  trace_summary_subdir,
                                  PlotSaveOptions::kLINEAR);
      output_file->cd();
      cL->Write(h2_totalE_vs_L[s]->GetName(), TObject::kOverwrite);
      delete cL;

      TCanvas *cR = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDraw2DHistogram(h2_totalE_vs_R[s], cR);
      h2_totalE_vs_R[s]->GetYaxis()->SetTitleOffset(1.3);
      if (Constants::cfg.SAVE_PLOTS)
        PlottingUtils::SaveFigure(cR, Form("totalE_vs_R_s%d", s),
                                  trace_summary_subdir,
                                  PlotSaveOptions::kLINEAR);
      output_file->cd();
      cR->Write(h2_totalE_vs_R[s]->GetName(), TObject::kOverwrite);
      delete cR;

      TCanvas *c1L = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDrawHistogram(h1_L[s], kBlue + 1);
      if (Constants::cfg.SAVE_PLOTS)
        PlottingUtils::SaveFigure(c1L, Form("L_s%d", s), trace_summary_subdir,
                                  PlotSaveOptions::kLOG);
      output_file->cd();
      h1_L[s]->Write(h1_L[s]->GetName(), TObject::kOverwrite);
      delete c1L;

      TCanvas *c1R = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDrawHistogram(h1_R[s], kBlue + 1);
      if (Constants::cfg.SAVE_PLOTS)
        PlottingUtils::SaveFigure(c1R, Form("R_s%d", s), trace_summary_subdir,
                                  PlotSaveOptions::kLOG);
      output_file->cd();
      h1_R[s]->Write(h1_R[s]->GetName(), TObject::kOverwrite);
      delete c1R;
    }

    TCanvas *c_dens = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(h2_long_back_vs_front, c_dens);
    h2_long_back_vs_front->GetYaxis()->SetTitleOffset(1.3);
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(c_dens, "long_back_vs_front",
                                trace_summary_subdir, PlotSaveOptions::kLINEAR);
    output_file->cd();
    c_dens->Write(h2_long_back_vs_front->GetName(), TObject::kOverwrite);
    delete c_dens;
  }

  for (Int_t s = 0; s < 18; s++) {
    delete h2_totalE_vs_stripE[s];
    delete h1_stripE[s];
  }
  for (Int_t s = 1; s <= 16; s++) {
    delete h2_totalE_vs_L[s];
    delete h2_totalE_vs_R[s];
    delete h1_L[s];
    delete h1_R[s];
  }
  delete h2_totalE_vs_cathode;
  delete h1_cathode;
  delete h2_long_back_vs_front;

  output_file->Close();
  delete output_file;

  if (emptyChannelMapEvents != 0)
    std::cout << "Observed " << emptyChannelMapEvents
              << " hits with empty entry in channel map." << std::endl;

  if (dropped_outside_window > 0)
    std::cout << "Dropped " << dropped_outside_window
              << " hits outside coincidence window." << std::endl;

  std::cout << Constants::cfg.REFERENCE_CHANNEL << " hits total: " << n_ref
            << std::endl;
  std::cout << Constants::cfg.REFERENCE_CHANNEL
            << " trigger rate: " << ref_rate_hz << " Hz over " << span_s << " s"
            << std::endl;
  std::cout << "Cathode hits total: " << cathode_hits_total << " (cathode/"
            << Constants::cfg.REFERENCE_CHANNEL << " = "
            << (n_ref > 0 ? Double_t(cathode_hits_total) / n_ref : 0.0) << ")"
            << std::endl;
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
    std::cout << "Dropped hits (dedup strategy): anode="
              << dropped_anode_hits_total
              << ", cathode=" << dropped_cathode_hits_total << std::endl;
  }
  std::cout << "Complete events: " << complete_events << " ("
            << (100.0 * complete_events / total_events) << "%)" << std::endl;
  std::cout << "Incomplete events: " << incomplete_events << " ("
            << (100.0 * incomplete_events / total_events) << "%)" << std::endl;
  if (total_events > 0) {
    std::cout << "Per-channel miss rate (% of all events where the "
                 "completeness-required channel was zero):"
              << std::endl;
    std::cout << "  Strip0=" << (100.0 * cnt.miss_strip0 / total_events)
              << "%   Strip17=" << (100.0 * cnt.miss_strip17 / total_events)
              << "%" << std::endl;
    for (Int_t s = 1; s <= 16; s++) {
      const Char_t *side = (s % 2 == 1) ? "L" : "R";
      std::cout << "  " << side << s << "="
                << (100.0 * cnt.miss_long[s] / total_events) << "%";
      if (s % 4 == 0)
        std::cout << std::endl;
    }
  }
  if (Constants::cfg.REJECT_FLAGGED_EVENTS) {
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
