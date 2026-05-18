#ifndef EVENT_BUILDER_HPP
#define EVENT_BUILDER_HPP

#include "BinaryUtils.hpp"
#include "Constants.hpp"
#include <TString.h>
#include <array>
#include <map>
#include <utility>
#include <vector>

namespace NearestGrid {

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

inline void InitEventState(EventState &e, UShort_t grid_energy,
                           UInt_t grid_flag) {
  for (Int_t k = 0; k < 18; k++) {
    e.leftdE[k] = 0;
    e.rightdE[k] = 0;
    e.totaldE[k] = 0;
  }
  for (Int_t k = 0; k < 36; k++)
    e.hits[k] = 0;
  e.cathode = -1;
  e.grid = grid_energy;
  e.flags_or = grid_flag;
  e.had_cathode = kFALSE;
  e.hits[35] = 1;
}

inline void InitPerChannelData(PerChannelData &p, ULong64_t grid_ts,
                               UInt_t grid_flag) {
  for (Int_t k = 0; k < 36; k++) {
    p.timestamps[k] = 0;
    p.flags[k] = 0;
  }
  p.timestamps[35] = grid_ts;
  p.flags[35] = grid_flag;
}

inline Bool_t CloserToGrid(ULong64_t cand_ts, ULong64_t prev_ts,
                           ULong64_t grid_ts) {
  ULong64_t cand_dt =
      (cand_ts > grid_ts) ? cand_ts - grid_ts : grid_ts - cand_ts;
  ULong64_t prev_dt =
      (prev_ts > grid_ts) ? prev_ts - grid_ts : grid_ts - prev_ts;
  return cand_dt < prev_dt;
}

constexpr Int_t kSlotStrip0 = 0;
constexpr Int_t kSlotStrip17 = 33;
constexpr Int_t kSlotCathode = 34;
constexpr Int_t kSlotGrid = 35;

typedef std::array<Int_t, Constants::N_BOARDS * Constants::N_CHANNELS> SlotMap;

inline SlotMap BuildSlotMap() {
  SlotMap sm;
  sm.fill(-1);
  for (std::map<std::pair<Int_t, Int_t>, TString>::const_iterator it =
           Constants::channelMap.begin();
       it != Constants::channelMap.end(); ++it) {
    Int_t board = it->first.first;
    Int_t channel = it->first.second;
    if (board < 0 || board >= Constants::N_BOARDS)
      continue;
    if (channel < 0 || channel >= Constants::N_CHANNELS)
      continue;
    Int_t idx = board * Constants::N_CHANNELS + channel;
    const TString &name = it->second;
    if (name == "Strip0") {
      sm[idx] = kSlotStrip0;
    } else if (name == "Strip17") {
      sm[idx] = kSlotStrip17;
    } else if (name == "Cathode") {
      sm[idx] = kSlotCathode;
    } else if (name == "Grid") {
      sm[idx] = kSlotGrid;
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

inline void AssignHit(EventState &e, PerChannelData *pc, ULong64_t grid_ts,
                      Int_t slot, UShort_t energy, ULong64_t timestamp,
                      UInt_t flags) {
  e.flags_or |= flags;

  Bool_t wins;
  if (e.hits[slot] == 0) {
    wins = kTRUE;
  } else if (pc) {
    wins = CloserToGrid(timestamp, pc->timestamps[slot], grid_ts);
  } else {
    wins = kFALSE;
  }

  if (wins) {
    if (slot == kSlotStrip0)
      e.totaldE[0] = energy;
    else if (slot <= 16)
      e.leftdE[slot] = energy;
    else if (slot <= 32)
      e.rightdE[slot - 16] = energy;
    else if (slot == kSlotStrip17)
      e.totaldE[17] = energy;
    else if (slot == kSlotCathode)
      e.cathode = energy;
    if (pc) {
      pc->timestamps[slot] = timestamp;
      pc->flags[slot] = flags;
    }
  }
  if (slot == kSlotCathode)
    e.had_cathode = kTRUE;
  e.hits[slot]++;
}

inline Bool_t CheckEventComplete(const EventState &e) {
  if (e.totaldE[0] == 0)
    return kFALSE;
  for (Int_t strip = 1; strip <= 7; strip += 2) {
    if (e.leftdE[strip] == 0)
      return kFALSE;
  }
  for (Int_t strip = 2; strip <= 8; strip += 2) {
    if (e.rightdE[strip] == 0)
      return kFALSE;
  }
  return kTRUE;
}

inline void GetFlagSummary(const EventState &e, Bool_t &has_fake,
                           Bool_t &has_saturation, Bool_t &has_pileup) {
  has_fake = (e.flags_or & CoMPASSData::FAKE_EVENT) != 0;
  has_saturation = ((e.flags_or & CoMPASSData::SATURATION_IN_GATE) ||
                    (e.flags_or & CoMPASSData::INPUT_SATURATING)) != 0;
  has_pileup = (e.flags_or & CoMPASSData::PILEUP) != 0;
}

} // namespace NearestGrid

Bool_t BuildEventsFromSortedHits(const std::vector<RawHit> &hits,
                                 const NearestGrid::SlotMap &slot_map,
                                 const TString &output_name,
                                 const TString &file_label);

#endif
