#ifndef EVENT_BUILDER_HPP
#define EVENT_BUILDER_HPP

#include "BinaryUtils.hpp"
#include "DedupStrategy.hpp"
#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "PlottingUtils.hpp"
#include <Rtypes.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TObject.h>
#include <TParameter.h>
#include <TString.h>
#include <TTree.h>
#include <array>
#include <iostream>
#include <map>
#include <mutex>
#include <utility>
#include <vector>

// Slot layout of the per-event arrays, fixed by the MUSIC detector geometry:
// 0 = Strip0, 1..16 = L1..L16, 17..32 = R1..R16, 33 = Strip17, 34 = Cathode,
// 35 = Grid.
struct EventState {
  Int_t leftdE[18];
  Int_t rightdE[18];
  Int_t totaldE[18];
  Int_t hits[36];
  Int_t cathode;
  Int_t grid;
  UInt_t flags_or;
  Bool_t had_cathode;
};

struct EventCounters {
  Int_t total_events;
  Int_t complete_events;
  Int_t complete_with_fake;
  Int_t complete_with_saturation;
  Int_t complete_with_pileup;
  Int_t complete_rejected;
  Int_t complete_rejected_multi;
  Int_t incomplete_events;
  Int_t incomplete_with_fake;
  Int_t incomplete_with_saturation;
  Int_t incomplete_with_pileup;
  Int_t events_with_cathode;
  Int_t events_with_multi_cathode;
  Int_t events_with_multi_anode_hit;
  Int_t dropped_anode_hits_total;
  Int_t dropped_cathode_hits_total;
  Int_t missing_long_0_17[18];
};

struct PerChannelData {
  ULong64_t timestamps[36];
  UShort_t energies[36];
  UInt_t flags[36];
};

class EventBuilder {
public:
  typedef std::vector<Int_t> SlotMap;

  static void ResetEventState(EventState &e);
  static void ResetPerChannelData(PerChannelData &p);
  static Bool_t ShouldKeepHit(ULong64_t cand_ts, ULong64_t prev_ts,
                              UShort_t cand_energy, UShort_t prev_energy,
                              ULong64_t ref_ts, DedupStrategy strategy);
  static SlotMap BuildSlotMap();
  static void AssignHit(EventState &e, PerChannelData *pc, ULong64_t ref_ts,
                        Int_t slot, UShort_t energy, ULong64_t timestamp,
                        UInt_t flags, DedupStrategy strategy);
  static Bool_t CheckEventComplete(const EventState &e);
  static void GetFlagSummary(const EventState &e, Bool_t &has_fake,
                             Bool_t &has_saturation, Bool_t &has_pileup);
  static Bool_t BuildEventsFromSortedHits(const std::vector<RawHit> &hits,
                                          const SlotMap &slot_map,
                                          const TString &output_name,
                                          const TString &file_label);
};

#endif
