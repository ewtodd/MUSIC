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

// Layout of the word returned by MapSOLFlagsToCoMPASS.
//
// The two SOLARIS words carry 18 defined bits between them -- flags_high 0-6
// and flags_low 0-10 -- and only 16 bits are free above the CoMPASS-compatible
// half, so they cannot both be kept whole. What makes it fit is that three of
// the flags_high bits (PILEUP, EVENT_SATURATION, POST_SATURATION) are already
// represented as their CoMPASS equivalents below, so carrying them again would
// be redundant: bits 2-6 of flags_high plus all 11 of flags_low come to exactly
// 16. Only flags_high bit 0 is dropped, and CoMPASSData::PILEUP is that bit.
//
// clang-format off
//   bits  0-15   CoMPASS-compatible flags (shared downstream code unchanged)
//   bits 16-26   SOLARIS flags_low, all 11 defined bits
//   bits 27-31   SOLARIS flags_high bits 2-6
// clang-format on
//
// Do not widen these without checking CoMPASSData: it defines PLL_LOCK_LOSS,
// OVER_TEMPERATURE and ADC_SHUTDOWN at bits 19-21, which is why the upper half
// is only safe to reuse for hits that came from a SOLARIS digitizer.
namespace SOLPack {
const Int_t LOW_SHIFT = 16;
const UInt_t LOW_MASK = 0x7FF; // flags_low bits 0-10
const Int_t HIGH_SHIFT = 27;
const Int_t HIGH_SKIP = 2;     // flags_high bits 0-1 are not carried here
const UInt_t HIGH_MASK = 0x1F; // flags_high bits 2-6
} // namespace SOLPack

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

  // The rest, packed per the layout above. Each word is masked to the width it
  // actually defines; the previous shifts (16 and 20) overlapped in bits 20-23,
  // so flags_high bits 4-7 and flags_low bits 0-3 were ORed on top of each
  // other and neither could be read back.
  mapped |= (static_cast<UInt_t>(sol_flags_low) & SOLPack::LOW_MASK)
            << SOLPack::LOW_SHIFT;
  mapped |= ((static_cast<UInt_t>(sol_flags_high) >> SOLPack::HIGH_SKIP) &
             SOLPack::HIGH_MASK)
            << SOLPack::HIGH_SHIFT;

  return mapped;
}

// Read the packed SOLARIS words back out. Use these rather than open-coding the
// shifts, so the layout stays defined in one place.
inline UShort_t SOLFlagsLowFrom(UInt_t mapped) {
  return UShort_t((mapped >> SOLPack::LOW_SHIFT) & SOLPack::LOW_MASK);
}

// flags_high bits 2-6. Bit 0 (pile-up) is not carried here; read
// CoMPASSData::PILEUP for it. Bits 1 and 7 are undefined by the format.
inline UShort_t SOLFlagsHighFrom(UInt_t mapped) {
  return UShort_t(((mapped >> SOLPack::HIGH_SHIFT) & SOLPack::HIGH_MASK)
                  << SOLPack::HIGH_SKIP);
}

// Convert a SOLHit to a RawHit for compatibility with the existing pipeline.
// Board is set to 0 (flat-channel SOLARIS digitizer). With the updated
// SOLHit (all unsigned fields), no casts are needed.
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
