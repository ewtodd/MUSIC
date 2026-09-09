#ifndef PULSE_HISTORY_HPP
#define PULSE_HISTORY_HPP

/**
 * @file PulseHistory.hpp
 * @brief Pole-zero pulse-history correction on the raw hit stream.
 *
 * A preamplifier's response to a pulse does not settle instantly, so a hit
 * arriving soon after a previous one sits on the tail of that earlier pulse and
 * reads high. This measures that shift per channel group as a kernel over the
 * time since the previous pulse, then subtracts it.
 *
 * The kernel is measured on each subfile's own beam-like events rather than
 * assumed, so it tracks whatever the electronics were doing during that run.
 */
#include "BinaryUtils.hpp"
#include <Rtypes.h>
#include <TString.h>
#include <vector>

class TH1D;
class TH2D;
class TProfile;

namespace PulseHistory {

/// Flag set on a hit whose correction was clamped rather than applied in full.
const UInt_t kFlagClamped = 0x2000;
/// Bins in the kernel, spanning #kLogLo to #kLogHi in log10 of dt.
const Int_t kNBins = 12;
const Double_t kLogLo = -6.0; ///< Lowest dt bin edge: log10 of 1 us in seconds.
const Double_t kLogHi =
    -4.0; ///< Highest dt bin edge: log10 of 100 us in seconds.
/// Maximum amplitude bands a kernel may split the previous pulse into.
const Int_t kMaxAmpBins = 6;

/**
 * @brief Channel groups, one kernel each.
 *
 * The long and short end of the split strips on either chain, plus the two
 * single-pad guards.
 *
 * @note A short end fires in only a fraction of beam events and its height
 *       follows the track position, so its kernel is fitted on the fired hits
 *       alone and describes a smaller share of its variance than a long end's.
 */
enum Group {
  kNone = 0,       ///< Not part of any corrected group.
  kLongLeft = 1,   ///< Long end, left chain.
  kLongRight = 2,  ///< Long end, right chain.
  kShortLeft = 3,  ///< Short end, left chain.
  kShortRight = 4, ///< Short end, right chain.
  kGuard0 = 5,     ///< Single-pad guard before strip 0.
  kGuard17 = 6,    ///< Single-pad guard after strip 17.
  kNGroups = 7     ///< Count of groups; not a group itself.
};
/// @brief Human-readable name of a Group.
const char *GroupName(Int_t g);
/// @brief Short tag for filenames: `L`, `R`, `Ls`, `Rs`, `S0`, `S17`.
const char *GroupTag(Int_t g);
/// @brief Which readout chain a group belongs to.
/// @return `0` for the left chain, `1` for the right, `-1` for the guards.
Int_t ChainOf(Int_t g);
/// @brief Whether a group is a long end rather than a short end or guard.
Bool_t IsLongGroup(Int_t g);

/**
 * @brief The fitted correction for one channel group.
 *
 * Check #ok before use; a group that could not be fitted leaves the
 * coefficients at zero.
 */
struct Kernel {
  Bool_t ok = kFALSE; ///< Whether this group was successfully fitted.
  Int_t n_amp = 1;    ///< Amplitude bands actually used, at most #kMaxAmpBins.
  Double_t k[kMaxAmpBins]
            [kNBins];        ///< Coefficients, `[amplitude band][dt bin]`.
  Double_t intercept = 0.0;  ///< Fit intercept, in ADC.
  Double_t r2 = 0.0;         ///< Coefficient of determination.
  Double_t rms_before = 0.0; ///< Channel-deviation RMS before correction.
  Double_t rms_after = 0.0;  ///< And after, so the gain is visible.
  Long64_t n = 0;            ///< (event, channel) pairs entering the fit.
  Long64_t n_beam = 0; ///< Beam-like events passing this group's selection.

  /// @brief Construct unfitted, with zeroed coefficients.
  Kernel();
};

/**
 * @brief Everything one subfile's pulse-history pass produced.
 *
 * @note Owns the diagnostic histograms. SavePlots() draws and frees them; the
 *       destructor frees any that survive.
 */
struct Result {
  Kernel kernel[kNGroups];    ///< One kernel per Group.
  Long64_t n_hits = 0;        ///< Hits examined.
  Long64_t n_seeds = 0;       ///< Hits usable as a previous pulse.
  Long64_t n_beam_events = 0; ///< Beam-like events the fit drew on.
  Long64_t n_corrected = 0;   ///< Hits that received a correction.
  Long64_t n_clamped = 0; ///< Corrections clamped rather than applied in full.
  Long64_t n_clamped_group[kNGroups]; ///< Clamps per group.
  Double_t mean_shift[kNGroups];      ///< Mean subtracted term per group, ADC.
  Bool_t sorted_input = kTRUE; ///< Whether the input arrived time-ordered.
  /// Beam peak per (board, channel) index, in ADC. Zero where the channel is
  /// not a long end, which is what AmpBinOf() keys off.
  std::vector<Double_t> mode;
  /// @name Diagnostics
  /// Filled when `SAVE_PLOTS` is on, drawn and freed by SavePlots(). Null
  /// otherwise.
  /// @{
  TH2D *dev_vs_pred[kNGroups]; ///< Channel deviation against predicted shift.
  TH1D *dev_before[kNGroups];  ///< Deviation distribution before correction.
  TH1D *dev_after[kNGroups];   ///< And after.
  TH1D *shift[kNGroups];       ///< Shift applied per hit.
  TProfile
      *dtprev_before[kNGroups]; ///< Deviation against time to previous pulse.
  TProfile *dtprev_after[kNGroups]; ///< And after.
  /// @}

  Result();  ///< @brief Construct with zeroed counters and null diagnostics.
  ~Result(); ///< @brief Frees any diagnostic histograms still held.
};

/// @brief Kernel bin for a time since the previous pulse.
/// @param dt_s Time difference, in seconds.
/// @return The bin index, or `-1` outside the kernel's range.
Int_t BinOf(Double_t dt_s);
/// @brief Centre of a kernel bin, in microseconds, for reports.
/// @param b Bin index.
Double_t BinCentreUs(Int_t b);
/// @brief Amplitude band of a previous pulse.
/// @param e_prev Previous pulse energy, in ADC.
/// @param mode   That channel's beam peak, in ADC.
/// @param n_amp  Bands in use for this kernel.
/// @return The band index.
Int_t AmpBinOf(Double_t e_prev, Double_t mode, Int_t n_amp);

/// @brief Group lookup for every (board, channel) under the active map.
/// @return A vector indexed as the channel map is, holding a Group per channel.
std::vector<Int_t> BuildGroupMap();

/**
 * @brief Measure the kernels on this subfile's beam-like events.
 *
 * @param[in,out] hits Raw hits for the subfile. **Sorted in place** if not
 *                     already time-ordered, which the event builder needs
 *                     anyway.
 * @param group_of     Group per channel, from BuildGroupMap().
 * @param[out] res     Fitted kernels and diagnostics.
 * @param file_label   Label from FileSet::FileLabel(), used in plots and logs.
 *
 * @return `kFALSE` when no group could be fitted.
 *
 * @todo Always pass through pulse history, and sort here before event building
 *       by default, so the CUDA sort can be used.
 */
Bool_t Measure(std::vector<RawHit> &hits, const std::vector<Int_t> &group_of,
               Result &res, const TString &file_label);

/**
 * @brief Apply the measured kernels to the hit stream, in place.
 *
 * @param[in,out] hits Hits to correct.
 * @param group_of     Group per channel, from BuildGroupMap().
 * @param[in,out] res  Kernels to apply; its clamp and shift counters are
 *                     updated as it goes.
 *
 * @note Only bins up to `Constants::cfg.PULSE_HISTORY_APPLY_MAX_US` are used.
 *       Beyond that the fitted coefficients are degenerate with the intercept
 *       and physically zero, so applying them would add noise rather than
 *       remove it.
 */
void Apply(std::vector<RawHit> &hits, const std::vector<Int_t> &group_of,
           Result &res);

/// @brief Format the pass as a human-readable report.
/// @param res        Result to summarise.
/// @param file_label Subfile label.
TString Report(const Result &res, const TString &file_label);
/// @brief Draw and save the diagnostics, then free them.
/// @param[in,out] res Result whose histograms are drawn and then deleted.
/// @param file_label  Subfile label, used in the plot names.
void SavePlots(Result &res, const TString &file_label);
/// @brief Record the kernels and counters alongside a subfile's events.
/// @param events_subpath Events file to write into.
/// @param res            Result to record.
void WriteToEventsFile(const TString &events_subpath, const Result &res);

} // namespace PulseHistory

#endif
