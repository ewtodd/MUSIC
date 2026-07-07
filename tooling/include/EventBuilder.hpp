#ifndef EVENT_BUILDER_HPP
#define EVENT_BUILDER_HPP

#include "BinaryUtils.hpp"
#include "DedupStrategy.hpp"
#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "PlottingUtils.hpp"
#include "SlotLayout.hpp"
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

struct EventState {
  Int_t leftdE[18];
  Int_t rightdE[18];
  Int_t totaldE[18];
  Int_t hits[Constants::N_ARR_SLOTS];
  Int_t cathode;
  Int_t grid;
  UInt_t flags_or;
  Bool_t had_cathode;
};

struct PerChannelData {
  ULong64_t timestamps[Constants::N_ARR_SLOTS];
  UShort_t energies[Constants::N_ARR_SLOTS];
  UInt_t flags[Constants::N_ARR_SLOTS];
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
