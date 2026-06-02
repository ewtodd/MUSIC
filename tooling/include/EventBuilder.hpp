#ifndef EVENT_BUILDER_HPP
#define EVENT_BUILDER_HPP

#include "BinaryUtils.hpp"
#include "Constants.hpp"
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

struct PerChannelData {
  ULong64_t timestamps[36];
  UInt_t flags[36];
};

class EventBuilder {
public:
  typedef std::array<Int_t, Constants::N_BOARDS * Constants::N_CHANNELS>
      SlotMap;

  static void InitEventState(EventState &e, UShort_t grid_energy,
                             UInt_t grid_flag);
  static void InitPerChannelData(PerChannelData &p, ULong64_t grid_ts,
                                 UInt_t grid_flag);
  static void ResetEventState(EventState &e);
  static void ResetPerChannelData(PerChannelData &p);
  static Bool_t CloserToGrid(ULong64_t cand_ts, ULong64_t prev_ts,
                             ULong64_t grid_ts);
  static SlotMap BuildSlotMap();
  static void AssignHit(EventState &e, PerChannelData *pc, ULong64_t grid_ts,
                        Int_t slot, UShort_t energy, ULong64_t timestamp,
                        UInt_t flags);
  static Bool_t CheckEventComplete(const EventState &e);
  static void GetFlagSummary(const EventState &e, Bool_t &has_fake,
                             Bool_t &has_saturation, Bool_t &has_pileup);
  static Bool_t BuildEventsFromSortedHits(const std::vector<RawHit> &hits,
                                          const SlotMap &slot_map,
                                          const TString &output_name,
                                          const TString &file_label);
};

#endif
