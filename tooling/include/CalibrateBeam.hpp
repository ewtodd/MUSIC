#ifndef CALIBRATE_BEAM_HPP
#define CALIBRATE_BEAM_HPP

#include "BeamFit2D.hpp"
#include "Constants.hpp"
#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Normalization.hpp"
#include "Paths.hpp"
#include "PlottingUtils.hpp"
#include <Rtypes.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TKey.h>
#include <TLegend.h>
#include <TLine.h>
#include <TList.h>
#include <TMath.h>
#include <TROOT.h>
#include <TSpectrum.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <queue>
#include <set>
#include <thread>
#include <vector>

/**
 * @brief One readout channel's beam-peak calibration.
 *
 * Produced per subfile by CalibrateBeam, and consumed by EnergyView to turn raw
 * ADC into calibrated units.
 */
struct ChannelCal {
  TString name; ///< Channel name from the active channel map.
  /// Readout end: `'L'` or `'R'` for the split strips 1-16, `'S'` for the
  /// unsplit end strips 0 and 17, `'C'` for the cathode, `'G'` for the grid.
  Char_t side;
  Int_t strip; ///< Anode strip index this channel reads; `-1` for the cathode.
  /// Modal beam peak, in ADC: the modal bin of the channel's beam-peak
  /// histogram (100 bins over the 5th-95th percentile of the samples), or the
  /// median for the cathode, which has no clean peak. No fit and no width.
  /// `0` when uncalibrated (too few samples).
  ///
  /// For the long end of a segmented strip this is its anchor; for the short
  /// end it only seeds the ridge fit: see #gain.
  Double_t fit_adc = 0.0;
  Long64_t n_samples = 0; ///< Events the anchor was taken from.
  /// Left/right gain-match override. When >= 0 this is used instead of
  /// `1 / fit_adc`; when negative the reciprocal fit is used.
  ///
  /// Set by the gain-match step. The long end is anchored on its own modal
  /// beam peak, which is the full deposit for the majority of beam events
  /// that leave the short end silent. A robust line through the beam-gated
  /// (short, long) pairs, the charge-sharing ridge, has slope `-k` and
  /// intercept `a`; the short end's anchor is `C_long / k` and its offset
  /// (#offset_adc) `(a - C_long) / k`, the sideways shift of the ridge
  /// relative to the long-only peak. With gain x (ADC - offset) on a short
  /// end that fired, events with and without a short hit sum to the same
  /// total. Where the ridge is not measurable the short side falls back to
  /// the parity's median `C_short / C_long` and median offset, or the global
  /// median ratio and no offset.
  /// This puts left and right on the same charge scale, so reaction events —
  /// which share charge between ends differently than beam events do — no
  /// longer sawtooth between even and odd strips.
  Double_t gain = -1.0;
  /// ADC subtracted from this channel before its gain, only when the channel
  /// fired (a reading of 0 stays 0). Short ends of strips 1-16 only, from the
  /// ridge; 0 elsewhere. On 37Cl the R-chain short ends carry 100-230 ADC of
  /// it, the L-chain ones a few tens.
  Double_t offset_adc = 0.0;
  /// Short channels only: `C_short / C_long = -1 / slope` from this subfile's
  /// own ridge fit, or `0` when the ridge was not measurable and a fallback
  /// was used.
  ///
  /// This is a preamplifier gain ratio and therefore fixed per channel. That
  /// is what lets AggregateRidgeRatiosForRun() take the run-level median over
  /// the subfiles that measured it (three or more) and rewrite each subfile's
  /// stored short-side gain from that median and its own long anchor; it
  /// edits the calibration tree, not this record.
  Double_t ridge_ratio = 0.0;
  /// Short channels only: #offset_adc in units of `C_long`, from this
  /// subfile's own ridge fit, aggregated over a run like #ridge_ratio.
  Double_t ridge_offset = 0.0;
};

/**
 * @brief Output of the strip centroid alignment, the post-gain calibration
 * step.
 *
 * Decodes events with the per-channel gains already applied and locates two
 * peaks in each strip's smoothed total (0.005 a.u. bins, parabolic sub-bin
 * refinement, no fit): the beam and the two-particle pile-up, which is
 * exactly twice the beam. The strip is then put on the scale `total' =
 * factor * total + offset` that has the beam at 1.0 and the pile-up at 2.0.
 * The second point is what separates the gain from a baseline offset, which
 * a beam-only anchor cannot do. A strip whose pile-up peak is not found keeps
 * a factor-only alignment, `1 / beam`.
 *
 * @note Runs **after** the anchor reduction and after the initial calibration
 *       tree is written — it needs that tree on disk so EnergyView can decode.
 *       The factors are then stored back into the calibration tree and applied
 *       by EnergyView to every end of a strip after the per-channel gain.
 */
struct StripAlignmentResult {
  Bool_t ok = kFALSE; ///< Whether the step produced usable results.
  /// Per-strip multiplicative alignment factors, applied to every end.
  /// Identity (1.0) when no peak was found, so an unaligned dataset decodes
  /// unchanged rather than being scaled by zero.
  Double_t factors[18] = {};
  /// Per-strip additive offsets in calibrated units, applied after the
  /// factor to the strip total; EnergyView carries each on the strip's long
  /// end. Zero when the pile-up peak was not found.
  Double_t offsets[18] = {};
  /// Measured beam peak per strip before alignment, in calibrated units.
  Double_t centroids[18] = {};
  /// Measured two-particle pile-up peak per strip before alignment, in
  /// calibrated units; 0 where none was found.
  Double_t pileups[18] = {};
};

/**
 * @brief Per-channel beam-peak calibration, gain matching and strip alignment.
 *
 * Runs after the pipeline over its events files. For each subfile it locates
 * every channel's modal beam peak, matches the left and right gains of each
 * segmented strip onto a common charge scale by anchoring both ends on the
 * charge-sharing ridge, then aligns every strip's total so that the beam
 * peak sits at 1.0 and the two-particle pile-up peak at 2.0.
 *
 * The ridge ratios are per-channel constants rather than per-subfile
 * measurements, and are aggregated to the run level once every subfile has
 * been processed.
 *
 * @todo Reintroduce the per-run energy-resolution TOML for Remix-MUSIC-Sim.
 *       It was removed on 2026-09-16 together with the per-channel peak
 *       width it was built from, when the beam anchor became the modal bin
 *       with no fit. It took each channel's relative resolution in percent
 *       FWHM from the calibration tree, medianed it over a run's subfiles,
 *       and wrote `sim_control/Calibration_Run<N>_eres.toml` as a
 *       `[detector.eres]` table keyed Cathode, S0, S17, L1..L16, R1..R16.
 *       Bringing it back needs a width per channel again, measured from the
 *       beam-peak histogram rather than a fit, and a run-level aggregation
 *       after the ridge ratios.
 */
class CalibrateBeam {
public:
  /**
   * @brief One ChannelCal per readout channel in the active channel map.
   * @return The templates, with names, sides and strips filled and the fitted
   *         quantities left at their defaults.
   */
  static std::vector<ChannelCal> BuildChannels();

  /**
   * @brief Calibrate one subfile end to end.
   * @param spec            Subfile to process.
   * @param chans_template  Channel templates from BuildChannels(); copied, not
   *                        modified.
   */
  static void
  CalibrateBeamOneSubfile(const FileSpec &spec,
                          const std::vector<ChannelCal> &chans_template);

  /**
   * @brief Replace each subfile's ridge ratio and short-end offset with the
   * run-level medians.
   *
   * The ratio is a preamplifier gain ratio and the offset a channel property,
   * so both are fixed per channel. Taking the median across the subfiles that
   * measured the ridge makes the values robust for subfiles where it was not
   * measurable.
   *
   * @param run   Run number.
   * @param specs That run's subfiles.
   */
  static void AggregateRidgeRatiosForRun(Int_t run,
                                         const std::vector<FileSpec> &specs);

  /**
   * @brief Calibrate every configured subfile, then aggregate per run.
   * @param file_label Restrict to one subfile by its FileSet::FileLabel().
   *                   Empty, the default, processes all of them.
   */
  static void Run(const TString &file_label = "");
};

#endif
