#ifndef PULSE_HISTORY_HPP
#define PULSE_HISTORY_HPP

/**
 * @file PulseHistory.hpp
 * @brief Pole-zero pulse-history correction on the raw hit stream.
 *
 * A preamplifier's response to a pulse does not settle instantly, so a hit
 * arriving soon after a previous one sits on the tail of that earlier pulse and
 * reads high or low. This measures that shift per channel as a kernel over
 * the time since the previous pulses, then subtracts it.
 *
 * The kernels are measured on each subfile's own beam-like events rather
 * than assumed, so they track whatever the electronics were doing during that
 * run. Channels are configured in groups (the long and short ends of either
 * chain, the two unsegmented strips): a group decides whether its channels
 * are corrected, with which kernel shape and with which given decay time,
 * and its own fit over all its channels is the fallback for a channel with
 * too few pairs for a kernel of its own.
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
/// Bins in the kernel, spanning #kLogLo to #kLogHi in log10 of dt: twelve
/// per decade, 0.083 dex each.
const Int_t kNBins = 30;
const Double_t kLogLo = -6.0; ///< Lowest dt bin edge: log10 of 1 us in seconds.
const Double_t kLogHi =
    -3.5; ///< Highest dt bin edge: log10 of 316 us in seconds.
/// Maximum amplitude bands a kernel may split the previous pulse into.
const Int_t kMaxAmpBins = 6;
/// Fewest (event, channel) pairs a channel needs for a kernel of its own;
/// below it the channel is corrected with its group's kernel.
const Long64_t kMinPairsPerChannel = 2000;
/// @name Pole-zero kernel form
/// The preamp decay time of the form (see `kPulseHistoryForm`) is not linear
/// in the fit, so it is profiled: the linear coefficients are solved for every
/// tau of a logarithmic grid from #kTauGridLoUs to #kTauGridHiUs and the tau
/// with the smallest residual is kept. The grid step is the resolution of the
/// profiled tau, about 15 percent at #kNTauGrid = 21. A tau given in the
/// configuration is solved for as well, at its exact value. The grid is the
/// main cost of the fit scan: two features per tau, band and past pulse, so
/// halving it roughly halves the pulse-history time when the form is on.
/// @{
const Int_t kNTauGrid = 21;
const Double_t kTauGridLoUs = 4.0;
const Double_t kTauGridHiUs = 64.0;
/// @}
/// @name Decay-time fit
/// Every fitted kernel's group additionally summarises the pre-correction
/// deviation-against-dt profile with `p0 + p1 exp(-dt / p2) + p3 dt` over
/// this window in microseconds. p2 is the preamp decay time the pole-zero
/// setting failed to cancel. The window starts after the trapezoid freeze
/// (holdoff plus margin, where the deviation follows the residual itself)
/// and ends inside the kernel reach, where the profile still carries the
/// tail and the linear term has not drowned it.
/// @{
const Double_t kTauFitLoUs = 8.0;
const Double_t kTauFitHiUs = 130.0;
/// Below this many (event, channel) pairs in the window the fit is skipped.
const Int_t kTauFitMinEntries = 100;
/// @}

/**
 * @brief Channel groups, one kernel each.
 *
 * The long and short end of the segmented strips on either chain, plus the
 * two unsegmented strips 0 and 17.
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
  kStrip0 = 5,     ///< Strip 0, unsegmented.
  kStrip17 = 6,    ///< Strip 17, unsegmented.
  kNGroups = 7     ///< Count of groups; not a group itself.
};
/// @brief Human-readable name of a Group.
const char *GroupName(Int_t g);
/// @brief Short tag for filenames: `L`, `R`, `Ls`, `Rs`, `S0`, `S17`.
const char *GroupTag(Int_t g);
/// @brief Which readout chain a group belongs to.
/// @return `0` for the left chain, `1` for the right, `-1` for strips 0/17.
Int_t ChainOf(Int_t g);
/// @brief Whether a group is a long end rather than a short end or an
///        unsegmented strip.
Bool_t IsLongGroup(Int_t g);
/// @brief Whether the configuration (`PULSE_HISTORY_GROUPS`) has this group
///        corrected. A disabled group gets no kernel and no correction but
///        still takes part in the beam-like selection.
Bool_t GroupEnabled(Int_t g);
/// @brief Whether the configuration asks for the pole-zero form on this
///        group rather than the binned kernel.
Bool_t GroupWantsForm(Int_t g);
/// @brief The preamp decay time the configuration gives for this group, in
///        microseconds; 0 when it is to be profiled.
Double_t GroupTauUs(Int_t g);

/**
 * @brief One fitted correction: a channel's own, or a group's.
 *
 * Every channel with enough pairs gets its own kernel, since the preamps
 * differ within a chain; the group's kernel, fitted on the sum of its
 * channels, is the fallback for the rest (#from_group) and the reference the
 * channels are compared against. Check #ok before use; a kernel that could
 * not be fitted leaves the coefficients at zero.
 */
struct Kernel {
  Bool_t ok = kFALSE; ///< Whether this kernel was successfully fitted.
  /// Channel kernels only: this is a copy of the group's kernel because the
  /// channel's own could not be fitted; #why_group says why.
  Bool_t from_group = kFALSE;
  TString why_group;
  Int_t n_amp = 1; ///< Amplitude bands actually used, at most #kMaxAmpBins.
  /// The applied kernel at the bin centres, `[amplitude band][dt bin]`: the
  /// binned coefficients themselves, or the form evaluated there when #form.
  /// What the report prints and the plots draw.
  Double_t k[kMaxAmpBins][kNBins];
  /// The binned fit's coefficients, always kept so the form can be judged
  /// against them. Equal to #k unless #form.
  Double_t k_binned[kMaxAmpBins][kNBins];

  /// @name Pole-zero form
  /// Fitted alongside the bins whenever any group asks for the form
  /// (`kPulseHistoryForm`):
  ///   k(dt) = -c exp(-(dt + t_m)/tau) + c_B (exp(-t_f/tau) - exp(-dt/tau))
  /// for dt >= t_lo, the baseline term only for dt > t_f; below t_lo one free
  /// coefficient per bin (#k_low). #form says whether it is the kernel that
  /// is applied; #form_ok whether it could be fitted at all.
  /// @{
  Bool_t form = kFALSE;
  Bool_t form_ok = kFALSE;
  Int_t j_tau = -1;      ///< Grid index of the applied tau; internal.
  Double_t tau_us = 0.0; ///< The applied form's decay time: given or profiled.
  Bool_t tau_given = kFALSE; ///< Whether #tau_us came from the configuration.
  /// The profile's own preference, on the grid, and its fit quality: the check
  /// on a given #tau_us. Equal to #tau_us and #r2_form when nothing is given.
  Double_t tau_free_us = 0.0;
  Double_t r2_form_free = 0.0;
  Double_t c[kMaxAmpBins];  ///< Direct-tail amplitude per band; > 0 undershoot.
  Double_t cb[kMaxAmpBins]; ///< Baseline-term amplitude per band.
  Double_t k_low[kMaxAmpBins][kNBins]; ///< Free coefficients, bins < #n_free.
  Int_t n_free = 0;       ///< Bins below t_lo that stay free under the form.
  Double_t t_m_us = 0.0;  ///< Read point after a pulse: rise + peaking.
  Double_t t_f_us = 0.0;  ///< Baseline freeze after a trigger: 2 rise + flat.
  Double_t t_lo_us = 0.0; ///< From where the form holds: (2 rise + flat) - t_m.
  Double_t intercept_form = 0.0; ///< The form fit's intercept, ADC.
  Double_t r2_form = 0.0;        ///< Coefficient of determination of the form.
  Double_t rms_after_form = 0.0; ///< Deviation RMS after the form.
  Double_t r2_binned = 0.0;      ///< And the binned fit's, for comparison.
  /// Channel kernels only: the R^2 the group's kernel reaches on this
  /// channel's pairs, so the gain from a kernel of its own is `r2 - r2_group`.
  Double_t r2_group = 0.0;
  Double_t rms_after_binned = 0.0;
  /// @brief The form at one gap for one band, per unit previous height.
  /// @param a     Amplitude band.
  /// @param dt_us Gap to the previous pulse, microseconds.
  Double_t FormAt(Int_t a, Double_t dt_us) const;
  /// @brief Effective restorer window implied by the fitted amplitudes,
  ///        `c tau / c_B`, in microseconds; 0 when c_B is 0.
  Double_t EffectiveWindowUs(Int_t a) const;
  /// @}
  /// @brief Single-exponential summary of the dt profile for this group.
  ///
  /// Fitted to the pre-correction mean deviation against dt over
  /// `p0 + p1 exp(-dt / p2) + p3 dt`, in microseconds and ADC. The linear
  /// term is not cosmetic: the trapezoid baseline is an average over a few
  /// hundred microseconds, and a hit with a long gap behind it finds that
  /// average less depressed by pile-up, lifting the profile roughly in
  /// proportion to the gap. Without p3 the fitted time constant comes out
  /// about 15 percent short on every group with a pole-zero mismatch.
  /// p2 is the preamp decay time the pole-zero failed to cancel; #flat flags
  /// the profiles with no decay to measure (see the flat rule in
  /// PulseHistory.cpp).
  struct TauFit {
    Bool_t ok = kFALSE;   ///< Whether the fit converged.
    Bool_t flat = kTRUE;  ///< Amplitude within 2 sigma of zero: no decay.
    Double_t p0 = 0.0;    ///< Level at large dt, ADC.
    Double_t p1 = 0.0;    ///< Tail amplitude at dt = 0, ADC.
    Double_t p1err = 0.0; ///< Its error, ADC.
    Double_t p2 = 0.0;    ///< Decay time, us.
    Double_t p2err = 0.0; ///< Its error, us.
    Double_t p3 = 0.0;    ///< Linear lift (the baseline term), ADC/us.
    Double_t chi2 = 0.0;  ///< Chi2 of the fit over the window.
    Int_t ndf = 0;        ///< Its degrees of freedom.
  };
  TauFit tau;                ///< Summary of the dt profile, all hits.
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
  /// The group fits, one per Group: each the fit over the sum of its
  /// channels, the fallback for a channel without a kernel of its own.
  Kernel kernel[kNGroups];
  /// The applied kernels, one per (board, channel) index as in the group map
  /// (empty before Measure()): the channel's own fit, or a flagged copy of
  /// its group's.
  std::vector<Kernel> kernel_ch;
  std::vector<Int_t> group_of;  ///< The group map Measure() was given.
  std::vector<TString> name_ch; ///< Channel-map name per index.
  Bool_t
      enabled[kNGroups]; ///< Whether the group was configured for correction.
  Long64_t n_hits = 0;   ///< Hits examined.
  Long64_t n_seeds = 0;  ///< Hits usable as a previous pulse.
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
  /// Per channel index, deviation against time to previous pulse before and
  /// after the channel's own kernel; drawn per group on one canvas.
  std::vector<TProfile *> dtprev_before_ch;
  std::vector<TProfile *> dtprev_after_ch;
  /// @}

  Result();  ///< @brief Construct with zeroed counters and null diagnostics.
  ~Result(); ///< @brief Frees any diagnostic histograms still held.
  /// @brief Delete every diagnostic histogram and null the pointers.
  void FreeDiagnostics();
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
 * @return `kFALSE` when no channel could be fitted.
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
 * @param[in,out] res  Kernels to apply, each channel's own from
 *                     `Result::kernel_ch`; its clamp and shift counters are
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
