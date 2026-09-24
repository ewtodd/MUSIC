#include "EventBuilder.hpp"
#include "EventsSummary.hpp"

void EventBuilder::ResetEventState(EventState &e) {
  for (Int_t k = 0; k < 16; k++) {
    e.leftdE[k] = 0;
    e.rightdE[k] = 0;
  }
  e.strip0dE = 0;
  e.strip17dE = 0;
  for (Int_t k = 0; k < Constants::N_ARR_SLOTS; k++)
    e.hits[k] = 0;
  e.cathode = -1;
  e.grid = 0;
  e.flags_or = 0;
  e.had_cathode = kFALSE;
  e.ref_ts = 0;
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
  if (strategy == kLATEST_TIMESTAMP)
    return cand_ts > prev_ts;

  ULong64_t cand_dt = (cand_ts > ref_ts) ? cand_ts - ref_ts : ref_ts - cand_ts;
  ULong64_t prev_dt = (prev_ts > ref_ts) ? prev_ts - ref_ts : ref_ts - prev_ts;
  return cand_dt < prev_dt;
}

EventBuilder::SlotMap EventBuilder::BuildSlotMap() {
  SlotMap sm;
  sm.resize(Constants::ActiveNBoards() * Constants::ActiveNChannels(), -1);
  for (std::map<std::pair<Int_t, Int_t>, TString>::const_iterator it =
           Constants::ActiveChannelMap().begin();
       it != Constants::ActiveChannelMap().end(); ++it) {
    Int_t board = it->first.first;
    Int_t channel = it->first.second;
    if (board < 0 || board >= Constants::ActiveNBoards())
      continue;
    if (channel < 0 || channel >= Constants::ActiveNChannels())
      continue;
    Int_t idx = board * Constants::ActiveNChannels() + channel;
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
    // Slots 1-16 are the left ends of strips 1-16, 17-32 the right ends.
    if (slot == Constants::ARR_SLOT_STRIP_0)
      e.strip0dE = energy;
    else if (slot <= 16)
      e.leftdE[slot - 1] = energy;
    else if (slot <= 32)
      e.rightdE[slot - 17] = energy;
    else if (slot == Constants::ARR_SLOT_STRIP_17)
      e.strip17dE = energy;
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

// The same first event-level cut as the analysis (AllStripsFired), so no
// half-read event reaches the beam-gate fits.
Bool_t EventBuilder::CheckEventComplete(const EventState &e) {
  if (!Constants::cfg.IGNORE_STRIP_0 && Constants::cfg.REQUIRE_STRIP_0 &&
      e.strip0dE == 0)
    return kFALSE;
  if (!Constants::cfg.IGNORE_STRIP_17 && e.strip17dE == 0)
    return kFALSE;
  for (Int_t strip = 1; strip <= 16; strip++) {
    const Int_t v =
        (strip % 2) != 0 ? e.leftdE[strip - 1] : e.rightdE[strip - 1];
    if (v == 0)
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
  /// Per-required-channel miss counts over all events: how often each
  /// completeness-required long end or unsegmented strip was zero. miss_long is
  /// indexed by strip 1..16 (counts L==0 for odd strips, R==0 for even).
  Long64_t miss_long[17];
  Long64_t miss_strip0;
  Long64_t miss_strip17;
};

void FinalizeEvent(EventState &e, TTree *output_tree, UShort_t *leftdE_branch,
                   UShort_t *rightdE_branch, UShort_t &strip0_branch,
                   UShort_t &strip17_branch, UShort_t *hits_branch,
                   Short_t &cathode_branch, Short_t &grid_branch,
                   UInt_t &flags_or_branch, ULong64_t &seed_ts_branch,
                   SummaryHistograms &hSum, EventCounters &c,
                   Long64_t event_idx = -1,
                   std::vector<TGraph *> *sample_traces = nullptr,
                   Long64_t sample_stride = 0, Int_t *n_sampled = nullptr) {
  // Collect sample traces for overlay plot
  if (sample_traces && event_idx >= 0 && sample_stride > 0 &&
      event_idx % sample_stride == 0 &&
      Int_t(sample_traces->size()) < Constants::cfg.SAVE_SAMPLE_TRACES) {
    Double_t total[18];
    for (Int_t s = 0; s < 18; s++)
      total[s] = Double_t(e.Total(s));
    sample_traces->push_back(EventsSummary::BuildTraceFromTotals(total));
    if (n_sampled)
      (*n_sampled)++;
  }

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
  if (e.strip0dE == 0)
    c.miss_strip0++;
  if (e.strip17dE == 0)
    c.miss_strip17++;
  for (Int_t s = 1; s <= 16; s++) {
    Int_t v = (s % 2 == 1) ? e.leftdE[s - 1] : e.rightdE[s - 1];
    if (v == 0)
      c.miss_long[s]++;
  }

  Bool_t has_fake = kFALSE, has_saturation = kFALSE, has_pileup = kFALSE;
  EventBuilder::GetFlagSummary(e, has_fake, has_saturation, has_pileup);
  Bool_t has_any_flag = has_fake || has_saturation || has_pileup;

  if (is_complete) {
    Bool_t reject = Constants::cfg.REJECT_FLAGGED_EVENTS && has_any_flag;
    if (Constants::ActiveDedupStrategy() == kDISCARD && any_anode_multi)
      reject = kTRUE;
    if (!reject) {
      for (Int_t k = 0; k < 16; k++) {
        leftdE_branch[k] = UShort_t(e.leftdE[k]);
        rightdE_branch[k] = UShort_t(e.rightdE[k]);
      }
      strip0_branch = UShort_t(e.strip0dE);
      strip17_branch = UShort_t(e.strip17dE);
      for (Int_t k = 0; k < Constants::N_ARR_SLOTS; k++)
        hits_branch[k] = UShort_t(e.hits[k]);
      cathode_branch = Short_t(e.cathode);
      grid_branch = Short_t(e.grid);
      flags_or_branch = e.flags_or;
      seed_ts_branch = e.ref_ts;
      output_tree->Fill();

      for (Int_t s = 0; s < 18; s++)
        hSum.h_music->Fill(Double_t(s), Double_t(e.Total(s)));

      for (Int_t s = 1; s <= 16; s++) {
        const Bool_t lIsLong = (s % 2) != 0;
        hSum.h2_long_vs_short[s - 1]->Fill(
            Double_t(lIsLong ? e.rightdE[s - 1] : e.leftdE[s - 1]),
            Double_t(lIsLong ? e.leftdE[s - 1] : e.rightdE[s - 1]));
      }

      if (hSum.h1_cathode && e.had_cathode)
        hSum.h1_cathode->Fill(Double_t(e.cathode));

      if (hSum.h1_strip17)
        hSum.h1_strip17->Fill(Double_t(e.strip17dE));

      if (hSum.h2_strip0_vs_grid)
        hSum.h2_strip0_vs_grid->Fill(
            Double_t(e.grid),
            Double_t(SummaryGridPartnerIsStrip1() ? e.Total(1) : e.strip0dE));

      if (hSum.h1_strip0)
        hSum.h1_strip0->Fill(Double_t(e.strip0dE));

      if (hSum.h1_grid)
        hSum.h1_grid->Fill(Double_t(e.grid));

      Int_t mult = 0;
      for (Int_t k = 0; k < Constants::N_ARR_SLOTS; k++)
        mult += e.hits[k];
      hSum.h_mult->Fill(Double_t(mult));
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

// A hit a reference-channel event has not yet assigned, queued until the
// next reference hit (or the stream end) decides its window.
struct PendingHit {
  Int_t slot;
  UShort_t energy;
  ULong64_t timestamp;
  UInt_t flags;
};

// The events tree's branch variables. The tree holds their addresses, so
// the instance must outlive the tree's Fills.
struct EventBranches {
  /// ADC energies are 14-bit unsigned at the source; store them as UShort_t,
  /// not Int_t. LeftdE and RightdE hold the two ends of the segmented strips
  /// 1-16, indexed by strip - 1; Strip0dE and Strip17dE the two unsegmented
  /// strips. A strip's deposit is left + right, computed on read, and never
  /// stored.
  UShort_t leftdE[16];
  UShort_t rightdE[16];
  UShort_t strip0dE;
  UShort_t strip17dE;
  UShort_t hits_arr[36];
  /// 14-bit ADC (<=16383), so Short_t holds every value with room to spare
  /// while preserving Cathode's -1 "no cathode hit" sentinel (Grid is
  /// non-negative but shares the type for symmetry). Anode values are
  /// non-negative, hence the UShort_t above.
  Short_t cathode;
  Short_t grid;
  UInt_t flags_or;
  ULong64_t seed_ts;
};

// The build's configuration, resolved once from the dataset.
struct BuildConfig {
  Bool_t ref_mode;
  Int_t ref_slot;
  ULong64_t window_ps;
  DedupStrategy dedup_strat;
  Bool_t holdoff_on;
  Double_t holdoff_ratio;
  ULong64_t holdoff_ps;
};

// The event currently being built and the hits queued for it.
struct OpenEvent {
  EventState event;
  PerChannelData per_channel;
  PerChannelData *pc;
  ULong64_t ref_ts;
  ULong64_t window_end;
  UShort_t seed_energy;
  Bool_t have;
  std::vector<PendingHit> pending;
};

// What the build writes into: tree, branch variables, histograms, counters.
struct EventSink {
  TTree *tree;
  EventBranches *br;
  SummaryHistograms *hSum;
  EventCounters *cnt;
  std::vector<TGraph *> *sample_traces;
  Long64_t sample_stride;
  Int_t n_sampled;
  Long64_t event_idx;
};

// The run's summary counters: the reference hits and their first and last
// stamps, the rejections, the cathode total, and the seed holdoff tallies.
struct BuildTallies {
  Int_t n_ref;
  ULong64_t first_ref_ts;
  ULong64_t last_ref_ts;
  Int_t empty_channel_map_events;
  Long64_t cathode_hits_total;
  Long64_t dropped_outside_window;
  Long64_t seeds_pretrigger;
  Long64_t seeds_merged;
  Long64_t seeds_close_pair;
  Long64_t seeds_close_empty;
  Long64_t seeds_close_anodes;
};

// A pre-trigger has no anodes pending: the strips fire after the real grid
// trigger. Size alone cannot tell it from a first particle read low.
static const Int_t kPreTriggerMaxPending = 2;

// The events tree; the large baskets have auto-flush disabled so ZSTD
// compresses in big chunks instead of many small basket flushes.
static TTree *BookEventBranches(EventBranches &br) {
  TTree *output_tree = new TTree("events", "MUSIC events");
  output_tree->Branch("LeftdE", br.leftdE, "LeftdE[16]/s");
  output_tree->Branch("RightdE", br.rightdE, "RightdE[16]/s");
  output_tree->Branch("Strip0dE", &br.strip0dE, "Strip0dE/s");
  output_tree->Branch("Strip17dE", &br.strip17dE, "Strip17dE/s");
  output_tree->Branch("Hits", br.hits_arr, "Hits[36]/s");
  output_tree->Branch("Cathode", &br.cathode, "Cathode/S");
  output_tree->Branch("Grid", &br.grid, "Grid/S");
  output_tree->Branch("FlagsOR", &br.flags_or, "FlagsOR/i");
  output_tree->Branch("SeedTs", &br.seed_ts, "SeedTs/l");

  output_tree->SetAutoFlush(0);
  for (Int_t b = 0; b < output_tree->GetListOfBranches()->GetEntries(); b++) {
    TBranch *branch =
        static_cast<TBranch *>(output_tree->GetListOfBranches()->At(b));
    branch->SetBasketSize(128 * 1024 * 1024);
  }
  return output_tree;
}

// The summary histograms, ranged for the active dataset.
static void ConfigureSummaryHistograms(SummaryHistograms &hSum) {
  SummaryHistConfig cfg;
  cfg.unit_label = "ADC";
  cfg.strip_e_min = Constants::ActiveStripEMinAdc();
  cfg.strip_e_max = Constants::ActiveStripEMaxAdc();
  cfg.odd_even_split = kTRUE;
  cfg.left_odd_max = Constants::ActiveLeftOddMaxAdc();
  cfg.left_even_max = Constants::ActiveLeftEvenMaxAdc();
  cfg.right_odd_max = Constants::ActiveRightOddMaxAdc();
  cfg.right_even_max = Constants::ActiveRightEvenMaxAdc();
  cfg.cathode_max = Constants::ActiveCathodeMaxAdc();
  cfg.grid_max = Constants::ActiveGridMaxAdc();
  cfg.strip0_max = Constants::ActiveStrip0MaxAdc();
  cfg.music_energy_bins = 2000;
  CreateSummaryHistograms(hSum, cfg);
}

// The summary histograms' deletion, in the order the allocation made.
static void DeleteSummaryHistograms(SummaryHistograms &hSum) {
  for (Int_t s = 1; s <= 16; s++)
    delete hSum.h2_long_vs_short[s - 1];
  delete hSum.h_music;
  delete hSum.h_mult;
  delete hSum.h1_cathode;
  delete hSum.h1_strip17;
  delete hSum.h2_strip0_vs_grid;
  delete hSum.h1_strip0;
  delete hSum.h1_grid;
}

// The sample-trace stride: a rough event estimate, each complete event
// having ~20 hits (18 strips + cathode + grid).
static Long64_t SampleTraceStride(Long64_t n_hits) {
  Long64_t sample_stride = 0;
  if (Constants::cfg.SAVE_SAMPLE_TRACES > 0) {
    Long64_t est_events = n_hits / 20;
    sample_stride = est_events / Long64_t(Constants::cfg.SAVE_SAMPLE_TRACES);
    if (sample_stride < 1)
      sample_stride = 1;
  }
  return sample_stride;
}

// The reference channel's slot, the active channel map scanned for its
// name; -1 when not found.
static Int_t ResolveReferenceSlot(const EventBuilder::SlotMap &slot_map) {
  Int_t ref_slot = -1;
  for (std::map<std::pair<Int_t, Int_t>, TString>::const_iterator it =
           Constants::ActiveChannelMap().begin();
       it != Constants::ActiveChannelMap().end(); ++it) {
    if (it->second == Constants::ActiveReferenceChannel()) {
      Int_t board = it->first.first;
      Int_t channel = it->first.second;
      if (board >= 0 && board < Constants::ActiveNBoards() && channel >= 0 &&
          channel < Constants::ActiveNChannels()) {
        ref_slot = slot_map[board * Constants::ActiveNChannels() + channel];
        break;
      }
    }
  }
  return ref_slot;
}

// The event counters, zeroed.
static void ZeroEventCounters(EventCounters &cnt) {
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
}

// A new event, seeded by `h` on `slot`: the state reset, the reference
// stamp, and in reference mode the window end and the seed's energy.
static inline void OpenNewEvent(OpenEvent &oe, const RawHit &h, Int_t slot,
                                const BuildConfig &cfg, Bool_t is_ref) {
  EventBuilder::ResetEventState(oe.event);
  EventBuilder::ResetPerChannelData(oe.per_channel);
  oe.ref_ts = h.timestamp;
  oe.event.ref_ts = oe.ref_ts;
  if (is_ref) {
    oe.window_end = oe.ref_ts + cfg.window_ps;
    oe.seed_energy = h.energy;
  }
  oe.have = kTRUE;
  EventBuilder::AssignHit(oe.event, oe.pc, oe.ref_ts, slot, h.energy,
                          h.timestamp, h.flags, cfg.dedup_strat);
}

// The queued hits assigned to the open event, the out-of-window ones
// counted.
static inline void FlushPending(OpenEvent &oe, const BuildConfig &cfg,
                                Long64_t &dropped_outside_window) {
  for (Int_t p = 0; p < Int_t(oe.pending.size()); p++) {
    if (oe.pending[p].timestamp <= oe.window_end) {
      EventBuilder::AssignHit(oe.event, oe.pc, oe.ref_ts, oe.pending[p].slot,
                              oe.pending[p].energy, oe.pending[p].timestamp,
                              oe.pending[p].flags, cfg.dedup_strat);
    } else {
      dropped_outside_window++;
    }
  }
  oe.pending.clear();
}

// The open event finalized into the sink; the event index advances.
static inline void CloseEvent(const OpenEvent &oe, EventSink &sink) {
  FinalizeEvent(oe.event, sink.tree, sink.br->leftdE, sink.br->rightdE,
                sink.br->strip0dE, sink.br->strip17dE, sink.br->hits_arr,
                sink.br->cathode, sink.br->grid, sink.br->flags_or,
                sink.br->seed_ts, *sink.hSum, *sink.cnt, sink.event_idx,
                sink.sample_traces, sink.sample_stride, &sink.n_sampled);
  sink.event_idx++;
}

// The seed holdoff: the pulse behind a pre-trigger re-seeds the open event
// rather than splitting the anodes. kTRUE when the hit was consumed.
static inline Bool_t MergePreTriggerSeed(const RawHit &h, Int_t ref_slot,
                                         OpenEvent &oe, const BuildConfig &cfg,
                                         BuildTallies &t) {
  if (oe.have && h.timestamp - oe.ref_ts <= cfg.holdoff_ps) {
    const Bool_t empty = Int_t(oe.pending.size()) <= kPreTriggerMaxPending;
    if (empty)
      t.seeds_close_empty++;
    else
      t.seeds_close_anodes++;
    if (empty &&
        Double_t(oe.seed_energy) < cfg.holdoff_ratio * Double_t(h.energy)) {
      t.seeds_pretrigger++;
      if (cfg.holdoff_on) {
        t.seeds_merged++;
        oe.ref_ts = h.timestamp;
        oe.event.ref_ts = oe.ref_ts;
        oe.window_end = oe.ref_ts + cfg.window_ps;
        // The larger energy is the real pulse; the stamp is the later hit.
        const UShort_t keep_e = std::max(oe.seed_energy, h.energy);
        oe.seed_energy = keep_e;
        oe.event.hits[ref_slot] = 0;
        EventBuilder::AssignHit(oe.event, oe.pc, oe.ref_ts, ref_slot, keep_e,
                                h.timestamp, h.flags, cfg.dedup_strat);
        return kTRUE;
      }
    } else {
      t.seeds_close_pair++;
    }
  }
  return kFALSE;
}

// Reference-channel mode: a reference hit seeds a new event, the others
// queue until the next one. kTRUE when the hit was skipped.
static inline Bool_t ReferenceModeHit(const RawHit &h, Int_t slot,
                                      const BuildConfig &cfg, OpenEvent &oe,
                                      EventSink &sink, BuildTallies &t) {
  if (slot == cfg.ref_slot) {
    // Grid ADC window filter: skip reference hits outside the accepted
    // range so they don't seed an event.
    if (Double_t(h.energy) < Constants::ActiveReferenceChannelMinAdc() ||
        Double_t(h.energy) > Constants::ActiveReferenceChannelMaxAdc()) {
      return kTRUE;
    }
    if (t.n_ref == 0)
      t.first_ref_ts = h.timestamp;
    t.last_ref_ts = h.timestamp;
    t.n_ref++;

    if (MergePreTriggerSeed(h, cfg.ref_slot, oe, cfg, t))
      return kTRUE;

    if (oe.have) {
      // Flush pending hits that fall within the current window.
      FlushPending(oe, cfg, t.dropped_outside_window);

      // Finalize the completed event.
      CloseEvent(oe, sink);
    }

    // This hit seeds the new event, keeping its energy; the outgoing event
    // gets a different ion's amplitude, uncorrelated with its anodes.
    OpenNewEvent(oe, h, cfg.ref_slot, cfg, kTRUE);
  } else {
    if (oe.have && h.timestamp <= oe.window_end) {
      PendingHit ph;
      ph.slot = slot;
      ph.energy = h.energy;
      ph.timestamp = h.timestamp;
      ph.flags = h.flags;
      oe.pending.push_back(ph);
    } else if (oe.have) {
      t.dropped_outside_window++;
    }
  }
  return kFALSE;
}

// Time-window mode (REFERENCE_CHANNEL == "NONE"): the first hit opens a
// window; its passing closes the event and opens the next.
static inline void WindowModeHit(const RawHit &h, Int_t slot,
                                 const BuildConfig &cfg, OpenEvent &oe,
                                 EventSink &sink) {
  if (!oe.have) {
    OpenNewEvent(oe, h, slot, cfg, kFALSE);
  } else if (h.timestamp - oe.ref_ts <= cfg.window_ps) {
    EventBuilder::AssignHit(oe.event, oe.pc, oe.ref_ts, slot, h.energy,
                            h.timestamp, h.flags, cfg.dedup_strat);
  } else {
    // Window exceeded — finalize and start new window.
    CloseEvent(oe, sink);
    OpenNewEvent(oe, h, slot, cfg, kFALSE);
  }
}

// The one-line summary outside the per-file sample: what the build did.
static void PrintBuildOneLiner(const TString &file_label,
                               const EventCounters &cnt, Int_t n_ref,
                               Double_t ref_rate_hz,
                               Long64_t dropped_outside_window) {
  // One string, one write: workers share the stream and a line built
  // piecewise interleaves with another worker's.
  const TString line = Form(
      "[events] %s: %lld events, %lld complete (%.1f%%), %lld %s hits at "
      "%.0f Hz, %lld outside window, dedup dropped %lld anode",
      file_label.Data(), Long64_t(cnt.total_events),
      Long64_t(cnt.complete_events),
      cnt.total_events > 0 ? 100.0 * cnt.complete_events / cnt.total_events
                           : 0.0,
      Long64_t(n_ref), Constants::ActiveReferenceChannel().Data(), ref_rate_hz,
      Long64_t(dropped_outside_window), Long64_t(cnt.dropped_anode_hits_total));
  std::cout << line << std::endl;
}

// The full summary block, for the sample files.
static void PrintBuildSummary(const BuildConfig &cfg, const BuildTallies &t,
                              const EventCounters &cnt, Double_t span_s,
                              Double_t ref_rate_hz, Double_t holdoff_us,
                              const TString &output_filepath) {
  if (t.empty_channel_map_events != 0)
    std::cout << "Observed " << t.empty_channel_map_events
              << " hits with empty entry in channel map." << std::endl;

  if (t.dropped_outside_window > 0)
    std::cout << "Dropped " << t.dropped_outside_window
              << " hits outside coincidence window." << std::endl;
  if (cfg.ref_mode && t.n_ref > 0) {
    // What the seed holdoff did, or would do.
    const Double_t pct = 100.0 / Double_t(t.n_ref);
    std::cout << "Seed holdoff " << (cfg.holdoff_on ? "ON" : "off") << " ("
              << (cfg.holdoff_on ? holdoff_us : 10.0) << " us"
              << (cfg.holdoff_on ? "" : " diagnostic") << ", guard ratio "
              << cfg.holdoff_ratio << "x): close seed pairs "
              << t.seeds_close_empty + t.seeds_close_anodes << " ("
              << Form("%.2f", pct * Double_t(t.seeds_close_empty +
                                             t.seeds_close_anodes))
              << "% of seeds): " << t.seeds_close_empty
              << " with no anodes pending ("
              << Form("%.2f", pct * Double_t(t.seeds_close_empty)) << "%), "
              << t.seeds_close_anodes << " with anodes pending ("
              << Form("%.2f", pct * Double_t(t.seeds_close_anodes))
              << "%). Pre-triggers (no anodes, current seed under the guard): "
              << t.seeds_pretrigger << " ("
              << Form("%.2f", pct * Double_t(t.seeds_pretrigger)) << "%), "
              << t.seeds_merged << " merged; " << t.seeds_close_pair
              << " kept apart." << std::endl;
  }

  std::cout << Constants::ActiveReferenceChannel() << " hits total: " << t.n_ref
            << std::endl;
  std::cout << Constants::ActiveReferenceChannel()
            << " trigger rate: " << ref_rate_hz << " Hz over " << span_s << " s"
            << std::endl;
  std::cout << "Cathode hits total: " << t.cathode_hits_total << " (cathode/"
            << Constants::ActiveReferenceChannel() << " = "
            << (t.n_ref > 0 ? Double_t(t.cathode_hits_total) / t.n_ref : 0.0)
            << ")" << std::endl;
  std::cout << "Total events: " << cnt.total_events << std::endl;
  if (cnt.total_events > 0) {
    std::cout << "Events with cathode hit: " << cnt.events_with_cathode << " ("
              << (100.0 * cnt.events_with_cathode / cnt.total_events) << "%)"
              << std::endl;
    std::cout << "Events with multiple cathode hits: "
              << cnt.events_with_multi_cathode << " ("
              << (100.0 * cnt.events_with_multi_cathode / cnt.total_events)
              << "% of all, "
              << (cnt.events_with_cathode > 0
                      ? 100.0 * cnt.events_with_multi_cathode /
                            cnt.events_with_cathode
                      : 0.0)
              << "% of cathode events)" << std::endl;
    std::cout << "Events with multi-hit on any anode: "
              << cnt.events_with_multi_anode_hit << " ("
              << (100.0 * cnt.events_with_multi_anode_hit / cnt.total_events)
              << "%)" << std::endl;
    std::cout << "Dropped hits (dedup strategy): anode="
              << cnt.dropped_anode_hits_total
              << ", cathode=" << cnt.dropped_cathode_hits_total << std::endl;
  }
  std::cout << "Complete events: " << cnt.complete_events << " ("
            << (100.0 * cnt.complete_events / cnt.total_events) << "%)"
            << std::endl;
  std::cout << "Incomplete events: " << cnt.incomplete_events << " ("
            << (100.0 * cnt.incomplete_events / cnt.total_events) << "%)"
            << std::endl;
  if (cnt.total_events > 0) {
    std::cout << "Per-channel miss rate (% of all events where the "
                 "completeness-required channel was zero):"
              << std::endl;
    std::cout << "  Strip0=" << (100.0 * cnt.miss_strip0 / cnt.total_events)
              << "%   Strip17=" << (100.0 * cnt.miss_strip17 / cnt.total_events)
              << "%" << std::endl;
    for (Int_t s = 1; s <= 16; s++) {
      const Char_t *side = (s % 2 == 1) ? "L" : "R";
      std::cout << "  " << side << s << "="
                << (100.0 * cnt.miss_long[s] / cnt.total_events) << "%";
      if (s % 4 == 0)
        std::cout << std::endl;
    }
  }
  if (Constants::cfg.REJECT_FLAGGED_EVENTS) {
    Int_t stored = cnt.complete_events - cnt.complete_rejected;
    std::cout << "Stored events (REJECT_FLAGGED_EVENTS=true): " << stored
              << " ("
              << (cnt.complete_events > 0 ? 100.0 * stored / cnt.complete_events
                                          : 0.0)
              << "% of complete; " << cnt.complete_rejected << " rejected)"
              << std::endl;
  }
  if (cnt.complete_events > 0) {
    std::cout << "Complete events with rejection-quality flags:" << std::endl;
    std::cout << "  Fake events: " << cnt.complete_with_fake << " ("
              << (100.0 * cnt.complete_with_fake / cnt.complete_events) << "%)"
              << std::endl;
    std::cout << "  Saturated: " << cnt.complete_with_saturation << " ("
              << (100.0 * cnt.complete_with_saturation / cnt.complete_events)
              << "%)" << std::endl;
    std::cout << "  Pileup: " << cnt.complete_with_pileup << " ("
              << (100.0 * cnt.complete_with_pileup / cnt.complete_events)
              << "%)" << std::endl;
  }
  if (cnt.incomplete_events > 0) {
    std::cout << "Incomplete events with rejection-quality flags:" << std::endl;
    std::cout << "  Fake events: " << cnt.incomplete_with_fake << " ("
              << (100.0 * cnt.incomplete_with_fake / cnt.incomplete_events)
              << "%)" << std::endl;
    std::cout << "  Saturated: " << cnt.incomplete_with_saturation << " ("
              << (100.0 * cnt.incomplete_with_saturation /
                  cnt.incomplete_events)
              << "%)" << std::endl;
    std::cout << "  Pileup: " << cnt.incomplete_with_pileup << " ("
              << (100.0 * cnt.incomplete_with_pileup / cnt.incomplete_events)
              << "%)" << std::endl;
  }

  std::cout << "Events saved to: " << output_filepath << std::endl;
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

  EventBranches br;
  TTree *output_tree = BookEventBranches(br);

  SummaryHistograms hSum;
  ConfigureSummaryHistograms(hSum);

  // Sample traces for overlay plot
  std::vector<TGraph *> sample_traces;
  Long64_t sample_stride = SampleTraceStride(Long64_t(hits.size()));
  Int_t n_sampled = 0;

  BuildConfig cfg;
  // Determine operating mode: reference-channel or blind time-window.
  cfg.ref_mode = Constants::ActiveReferenceChannel() != "NONE";
  // Resolve reference channel name to slot ID by scanning active channel map.
  cfg.ref_slot = cfg.ref_mode ? ResolveReferenceSlot(slot_map) : -1;
  if (cfg.ref_mode && cfg.ref_slot < 0) {
    std::cerr << "FATAL: reference channel '"
              << Constants::ActiveReferenceChannel()
              << "' not found in channel map. Cannot build events."
              << std::endl;
    output_file->Close();
    delete output_file;
    DeleteSummaryHistograms(hSum);
    for (Int_t i = 0; i < Int_t(sample_traces.size()); i++)
      delete sample_traces[i];
    return kFALSE;
  }
  cfg.window_ps = ULong64_t(Constants::ActiveEventTimeWindowUs() * 1.0e6);
  cfg.dedup_strat = Constants::ActiveDedupStrategy();
  // Seed holdoff (SEED_HOLDOFF_US); off, the candidates are still counted in
  // a 10 us window so the summary says what turning it on would do.
  const Double_t holdoff_us = Constants::ActiveSeedHoldoffUs();
  cfg.holdoff_ratio = Constants::ActiveSeedHoldoffMaxRatio();
  cfg.holdoff_on = holdoff_us > 0.0;
  cfg.holdoff_ps = ULong64_t((cfg.holdoff_on ? holdoff_us : 10.0) * 1.0e6);

  Long64_t n_entries = Long64_t(hits.size());

  OpenEvent oe;
  oe.pc = &oe.per_channel;
  oe.ref_ts = 0;
  oe.window_end = 0;
  oe.seed_energy = 0;
  oe.have = kFALSE;
  oe.pending.reserve(4096);

  EventCounters cnt;
  ZeroEventCounters(cnt);

  BuildTallies t;
  t.n_ref = 0;
  t.first_ref_ts = 0;
  t.last_ref_ts = 0;
  t.empty_channel_map_events = 0;
  t.cathode_hits_total = 0;
  t.dropped_outside_window = 0;
  t.seeds_pretrigger = 0;
  t.seeds_merged = 0;
  t.seeds_close_pair = 0;
  t.seeds_close_empty = 0;
  t.seeds_close_anodes = 0;

  EventSink sink;
  sink.tree = output_tree;
  sink.br = &br;
  sink.hSum = &hSum;
  sink.cnt = &cnt;
  sink.sample_traces = &sample_traces;
  sink.sample_stride = sample_stride;
  sink.n_sampled = n_sampled;
  sink.event_idx = 0;

  if (Constants::FileInSample())
    std::cout << "[" << file_label << "] Streaming pass over " << n_entries
              << " sorted hits (reference: "
              << Constants::ActiveReferenceChannel()
              << ", window: " << Constants::ActiveEventTimeWindowUs() << " us)"
              << std::endl;

  for (Long64_t i = 0; i < n_entries; i++) {
    const RawHit &h = hits[i];

    if (h.board >= Constants::ActiveNBoards() ||
        h.channel >= Constants::ActiveNChannels()) {
      t.empty_channel_map_events++;
      continue;
    }
    Int_t slot = slot_map[h.board * Constants::ActiveNChannels() + h.channel];
    if (slot < 0) {
      t.empty_channel_map_events++;
      continue;
    }

    if (slot == Constants::ARR_SLOT_CATHODE)
      t.cathode_hits_total++;

    if (cfg.ref_mode) {
      if (ReferenceModeHit(h, slot, cfg, oe, sink, t))
        continue;
    } else {
      WindowModeHit(h, slot, cfg, oe, sink);
    }

#if MUSIC_HOT_PATH_LOGGING
    if (i % 10000000 == 0)
      std::cout << "  Stream progress: " << i << "/" << n_entries << std::endl;
#endif
  }

  // Finalize the last event.
  if (oe.have) {
    if (cfg.ref_mode)
      FlushPending(oe, cfg, t.dropped_outside_window);
    CloseEvent(oe, sink);
  }

  if (Constants::FileInSample())
    std::cout << "Found " << t.n_ref << " "
              << Constants::ActiveReferenceChannel() << " hits." << std::endl;

  if (cfg.ref_mode && t.n_ref == 0) {
    std::cerr << "No " << Constants::ActiveReferenceChannel()
              << " hits in file, skipping." << std::endl;
    output_file->Close();
    delete output_file;
    DeleteSummaryHistograms(hSum);
    return kFALSE;
  }

  Double_t span_s = (t.last_ref_ts > t.first_ref_ts)
                        ? Double_t(t.last_ref_ts - t.first_ref_ts) / 1e12
                        : 0.0;
  Double_t ref_rate_hz = (span_s > 0.0) ? Double_t(t.n_ref) / span_s : 0.0;

  output_file->cd();
  TParameter<Double_t>("grid_rate_hz", ref_rate_hz).Write();
  output_tree->Write("events", TObject::kOverwrite);
  hSum.h_mult->Write("", TObject::kOverwrite);

  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    TString subdir = "events_summary/" + file_label;
    SaveAndDeleteSummaryHistograms(hSum, output_file, subdir, "");

    // Sample traces overlay
    if (!sample_traces.empty()) {
      EventsSummary::SaveSampleTraces(sample_traces, "sample_traces", subdir,
                                      0.0, Constants::ActiveStripEMaxAdc(),
                                      "Energy [ADC]");
    }
  }

  output_file->Close();
  delete output_file;

  // Outside the per-file sample one line says what the build did; the full
  // block below is for the sample files.
  if (!Constants::FileInSample()) {
    PrintBuildOneLiner(file_label, cnt, t.n_ref, ref_rate_hz,
                       t.dropped_outside_window);
    return kTRUE;
  }
  PrintBuildSummary(cfg, t, cnt, span_s, ref_rate_hz, holdoff_us,
                    output_filepath);
  return kTRUE;
}
