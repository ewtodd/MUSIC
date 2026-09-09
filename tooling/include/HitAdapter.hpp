#ifndef HIT_ADAPTER_HPP
#define HIT_ADAPTER_HPP

#include "BinaryUtils.hpp"
#include <Rtypes.h>

/**
 * @file HitAdapter.hpp
 * @brief Converting SOLARIS hits into the CoMPASS-shaped form the pipeline
 * uses.
 *
 * The pipeline was written against CoMPASS. Rather than fork every downstream
 * stage for the SOLARIS era, SOLARIS hits are adapted into `RawHit` with
 * CoMPASS-compatible flag bits, so EventBuilder, the timing filters and
 * everything after them work unchanged.
 */

/// @brief SOLARIS status flag bits, from the SOLARIS DAQ source.
namespace SOLFlags {
/// @name High-priority word, 8 bits
/// @{
const UInt_t PILEUP = 0x01;           ///< bit 0
const UInt_t EVENT_SATURATION = 0x04; ///< bit 2
const UInt_t POST_SATURATION = 0x08;  ///< bit 3
const UInt_t CHARGE_OVERFLOW = 0x10;  ///< bit 4
const UInt_t SCA_SELECTED = 0x20;     ///< bit 5
const UInt_t FINE_TS_VALID = 0x40;    ///< bit 6
/// @}

/// @name Low-priority word, 12 bits
/// @{
const UInt_t EXT_INHIBIT = 0x0001;    ///< bit 0
const UInt_t UNDER_SAT = 0x0002;      ///< bit 1
const UInt_t OVER_SAT = 0x0004;       ///< bit 2
const UInt_t EXT_TRIGGER = 0x0008;    ///< bit 3
const UInt_t GLOBAL_TRIGGER = 0x0010; ///< bit 4
const UInt_t SW_TRIGGER = 0x0020;     ///< bit 5
const UInt_t SELF_TRIGGER = 0x0040;   ///< bit 6
const UInt_t LVDS_TRIGGER = 0x0080;   ///< bit 7
const UInt_t TRIGGER_64CH = 0x0100;   ///< bit 8
const UInt_t ITLA_TRIGGER = 0x0200;   ///< bit 9
const UInt_t ITLB_TRIGGER = 0x0400;   ///< bit 10
/// @}
} // namespace SOLFlags

/**
 * @brief Layout of the packed word MapSOLFlagsToCoMPASS() returns.
 *
 * The two SOLARIS words carry 18 defined bits between them — `flags_high` 0-6
 * and `flags_low` 0-10 — while only 16 bits are free above the
 * CoMPASS-compatible half. They cannot both be kept whole.
 *
 * What makes it fit: three of the `flags_high` bits (pileup, event saturation,
 * post saturation) are already represented as their CoMPASS equivalents in the
 * lower half, so carrying them again would be redundant. Bits 2-6 of
 * `flags_high` plus all 11 of `flags_low` come to exactly 16. Only
 * `flags_high` bit 0 is dropped, and `CoMPASSData::PILEUP` is that bit.
 *
 * ```
 *   bits  0-15   CoMPASS-compatible flags (downstream code unchanged)
 *   bits 16-26   SOLARIS flags_low, all 11 defined bits
 *   bits 27-31   SOLARIS flags_high bits 2-6
 * ```
 *
 * @warning Do not widen these without checking `CoMPASSData`. It defines
 *          `PLL_LOCK_LOSS`, `OVER_TEMPERATURE` and `ADC_SHUTDOWN` at bits
 *          19-21, which is why the upper half is only safe to reuse for hits
 *          that came from a SOLARIS digitiser in the first place.
 */
namespace SOLPack {
const Int_t LOW_SHIFT = 16;    ///< Where flags_low starts.
const UInt_t LOW_MASK = 0x7FF; ///< flags_low bits 0-10.
const Int_t HIGH_SHIFT = 27;   ///< Where the flags_high remnant starts.
const Int_t HIGH_SKIP = 2;     ///< flags_high bits 0-1 are not carried here.
const UInt_t HIGH_MASK = 0x1F; ///< flags_high bits 2-6.
} // namespace SOLPack

/**
 * @brief Pack SOLARIS flags into a CoMPASS-compatible word.
 *
 * The three conditions that have CoMPASS equivalents are mapped onto them, so
 * downstream flag checks work unchanged. The remainder is packed above per
 * SOLPack.
 *
 * @param sol_flags_high SOLARIS high-priority flag word.
 * @param sol_flags_low  SOLARIS low-priority flag word.
 * @return The packed word, for `RawHit::flags`.
 *
 * @note An earlier layout used shifts of 16 and 20, which overlapped in bits
 *       20-23: `flags_high` bits 4-7 and `flags_low` bits 0-3 were ORed on top
 *       of each other and neither could be read back. Each word is now masked
 *       to the width it actually defines.
 */
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

/**
 * @brief Recover the SOLARIS low-priority word from a packed value.
 * @param mapped Word from MapSOLFlagsToCoMPASS().
 * @return `flags_low`, bits 0-10.
 * @note Use this rather than open-coding the shift, so the layout stays defined
 *       in one place.
 */
inline UShort_t SOLFlagsLowFrom(UInt_t mapped) {
  return UShort_t((mapped >> SOLPack::LOW_SHIFT) & SOLPack::LOW_MASK);
}

/**
 * @brief Recover the SOLARIS high-priority word from a packed value.
 * @param mapped Word from MapSOLFlagsToCoMPASS().
 * @return `flags_high` bits 2-6, in their original positions.
 * @note Bit 0, pileup, is not carried here — read `CoMPASSData::PILEUP`
 *       instead. Bits 1 and 7 are undefined by the format.
 */
inline UShort_t SOLFlagsHighFrom(UInt_t mapped) {
  return UShort_t(((mapped >> SOLPack::HIGH_SHIFT) & SOLPack::HIGH_MASK)
                  << SOLPack::HIGH_SKIP);
}

/**
 * @brief Adapt one SOLARIS hit into the pipeline's `RawHit`.
 *
 * @param sol Hit to convert.
 * @return The adapted hit.
 *
 * @note Board is always `0`: the SOLARIS digitiser presents a flat channel
 *       space with no board dimension. Timestamps are scaled from nanoseconds
 *       to the picoseconds the pipeline works in.
 */
inline RawHit SOLHitToRawHit(const SOLHit &sol) {
  RawHit raw;
  raw.board = 0;
  raw.channel = sol.channel;
  raw.energy = sol.energy;
  raw.timestamp = sol.timestamp * 1000;
  raw.flags = MapSOLFlagsToCoMPASS(sol.flags_high, sol.flags_low);
  return raw;
}

/// @brief Adapt a whole run's SOLARIS hits.
/// @param sol_hits Hits to convert.
/// @return The adapted hits, in the same order.
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
