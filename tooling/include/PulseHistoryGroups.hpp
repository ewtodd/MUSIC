#ifndef PULSE_HISTORY_GROUPS_HPP
#define PULSE_HISTORY_GROUPS_HPP

#include <Rtypes.h>

/**
 * @file PulseHistoryGroups.hpp
 * @brief The per-group settings of the pulse-history correction.
 *
 * Its own header because both the flat `DatasetConfig` and `RunEpoch` carry
 * a set: MakeEpoch() copies the flat one into the epoch, and an epoch
 * overrides what its electronics need (a known preamp decay for one era, a
 * chain that needs no correction in another).
 */

/**
 * @brief Shape of the pulse-history kernel applied to a channel group.
 *
 * Both are fitted whenever any group asks for the form, so the report can
 * compare them; this chooses which one corrects the group's energies.
 */
enum PulseHistoryKernel {
  /// One free coefficient per logarithmic dt bin: tracks whatever the
  /// electronics did, at the price of bin-to-bin noise on sparse groups.
  kPulseHistoryBinned,
  /// The pole-zero form (PoleZeroMUSIC, doc/polezero.tex Eq. kernel): the
  /// direct tail of the previous pulse under the read point plus that pulse's
  /// share of the trapezoid baseline,
  ///   k(dt) = -c exp(-(dt + t_m)/tau) + c_B (exp(-t_f/tau) - exp(-dt/tau)),
  /// the second term only once the baseline unfreezes (dt > t_f). tau is the
  /// preamp decay: given per group (`PulseHistoryGroupOption::tau_us`) or,
  /// when not, profiled on a grid; the profile runs either way so a given
  /// value can be checked against what the data prefer. c and c_B are fitted
  /// per amplitude band, c_B free because the restorer's effective window is
  /// not known from the settings. t_m and t_f come from the trapezoid
  /// settings below.
  /// Gaps shorter than the trapezoid, where the previous pulse's own top and
  /// fall sit under the read point, keep one free coefficient per bin.
  kPulseHistoryForm
};

/// @brief Per-group pulse-history setting: whether the group is corrected and
/// with which kernel shape.
struct PulseHistoryGroupOption {
  Bool_t enabled = kFALSE;
  PulseHistoryKernel kernel = kPulseHistoryBinned;
  /// Known preamp decay time in microseconds, for the form. Positive: the
  /// applied form uses this tau and the report says where the free profile
  /// lands relative to it. Zero: the profiled tau is applied.
  Double_t tau_us = 0.0;
};

/**
 * @brief Which channel groups the pulse-history correction fits and applies
 * to, one setting per kernel group of PulseHistory (see PulseHistory::Group).
 *
 * Every group defaults to off. A group that is off keeps its channels in the
 * beam-like selection the kernels are fitted on (that selection needs every
 * long end and is not changed by this mask) but gets no kernel and no
 * correction, so its energies are left as read out. The set in force is
 * Constants::ActivePulseHistoryGroups(): the active epoch's copy when one is
 * set, the flat `PULSE_HISTORY_GROUPS` otherwise.
 */
struct PulseHistoryGroups {
  PulseHistoryGroupOption long_left;   ///< Long end of the odd strips (L).
  PulseHistoryGroupOption long_right;  ///< Long end of the even strips (R).
  PulseHistoryGroupOption short_left;  ///< Short end of the even strips (L).
  PulseHistoryGroupOption short_right; ///< Short end of the odd strips (R).
  PulseHistoryGroupOption strip0;      ///< Strip 0, unsegmented.
  PulseHistoryGroupOption strip17;     ///< Strip 17, unsegmented.
  /// The Frisch grid, the channel whose hits seed the events. Its height per
  /// event is the seed hit's own; the correction is applied before event
  /// building, so the seed's ADC gate reads the corrected value.
  PulseHistoryGroupOption grid;

  /// @brief Whether any group is on: an epoch with none skips the
  ///        correction altogether.
  Bool_t AnyEnabled() const {
    return long_left.enabled || long_right.enabled || short_left.enabled ||
           short_right.enabled || strip0.enabled || strip17.enabled ||
           grid.enabled;
  }
};

#endif
