#ifndef HIT_ADAPTER_HPP
#define HIT_ADAPTER_HPP

#include "BinaryUtils.hpp"
#include <Rtypes.h>

// SOLARIS flag bit definitions (from SOLARIS DAQ source)
namespace SOLFlags {
// High priority (8 bits)
const UInt_t PILEUP = 0x01;           // bit 0
const UInt_t EVENT_SATURATION = 0x04; // bit 2
const UInt_t POST_SATURATION = 0x08;  // bit 3
const UInt_t CHARGE_OVERFLOW = 0x10;  // bit 4
const UInt_t SCA_SELECTED = 0x20;     // bit 5
const UInt_t FINE_TS_VALID = 0x40;    // bit 6

// Low priority (12 bits)
const UInt_t EXT_INHIBIT = 0x0001;    // bit 0
const UInt_t UNDER_SAT = 0x0002;      // bit 1
const UInt_t OVER_SAT = 0x0004;       // bit 2
const UInt_t EXT_TRIGGER = 0x0008;    // bit 3
const UInt_t GLOBAL_TRIGGER = 0x0010; // bit 4
const UInt_t SW_TRIGGER = 0x0020;     // bit 5
const UInt_t SELF_TRIGGER = 0x0040;   // bit 6
const UInt_t LVDS_TRIGGER = 0x0080;   // bit 7
const UInt_t TRIGGER_64CH = 0x0100;   // bit 8
const UInt_t ITLA_TRIGGER = 0x0200;   // bit 9
const UInt_t ITLB_TRIGGER = 0x0400;   // bit 10
} // namespace SOLFlags

// Map SOLARIS flags to CoMPASS-compatible flag bits so downstream code
// (EventBuilder, timing filters, etc.) can use the same flag checks.
inline UInt_t MapSOLFlagsToCoMPASS(UShort_t sol_flags_high,
                                   UShort_t sol_flags_low) {
  UInt_t mapped = 0;

  // SOLARIS high-priority pileup -> CoMPASS PILEUP
  if (sol_flags_high & SOLFlags::PILEUP) {
    mapped |= CoMPASSData::PILEUP;
  }

  // SOLARIS event saturation -> CoMPASS SATURATION_IN_GATE
  if (sol_flags_high & SOLFlags::EVENT_SATURATION) {
    mapped |= CoMPASSData::SATURATION_IN_GATE;
  }

  // SOLARIS post-saturation -> CoMPASS INPUT_SATURATING
  if (sol_flags_high & SOLFlags::POST_SATURATION) {
    mapped |= CoMPASSData::INPUT_SATURATING;
  }

  // Store remaining SOLARIS flags in upper bits for downstream inspection
  // (bits 16-31 are unused by CoMPASS)
  mapped |= (static_cast<UInt_t>(sol_flags_high) << 16);
  mapped |= (static_cast<UInt_t>(sol_flags_low) << 20);

  return mapped;
}

// Convert a SOLHit to a RawHit for the pipeline: board = 0 (flat-channel
// SOLARIS digitizer); the all-unsigned SOLHit needs no casts.
inline RawHit SOLHitToRawHit(const SOLHit &sol) {
  RawHit raw;
  raw.board = 0;
  raw.channel = sol.channel;
  raw.energy = sol.energy;
  raw.timestamp = sol.timestamp * 1000;
  raw.flags = MapSOLFlagsToCoMPASS(sol.flags_high, sol.flags_low);
  return raw;
}

// Convert a vector of SOLHit to RawHit.
inline std::vector<RawHit>
SOLHitsToRawHits(const std::vector<SOLHit> &sol_hits) {
  std::vector<RawHit> raw_hits;
  raw_hits.reserve(sol_hits.size());
  for (Int_t i = 0; i < Int_t(sol_hits.size()); i++) {
    raw_hits.push_back(SOLHitToRawHit(sol_hits[i]));
  }
  return raw_hits;
}

#endif
