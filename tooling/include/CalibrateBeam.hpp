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
#include <toml++/toml.hpp>
#include <vector>

/**
 * @brief One readout channel's beam-peak calibration.
 *
 * Produced per subfile by CalibrateBeam, and consumed by EnergyView to turn raw
 * ADC into calibrated units.
 */
struct ChannelCal {
  TString name;                 ///< Channel name from the active channel map.
  Char_t side;                  ///< Readout end: `'L'` or `'R'`.
  Int_t strip;                  ///< Anode strip index this channel reads.
  Double_t fit_adc = 0.0;       ///< Fitted beam-peak centroid, in ADC.
  Double_t fit_sigma_adc = 0.0; ///< Fitted beam-peak width, in ADC.
  Long64_t n_samples = 0.0;     ///< Events entering the fit.
  /// Left/right gain-match override. When >= 0 this is used instead of
  /// `1 / fit_adc`; when negative the reciprocal fit is used.
  ///
  /// Set by the gain-match step: the long side is anchored on its own beam
  /// peak, the short side on `C_long / |slope|` of the charge-sharing ridge.
  /// This puts left and right on the same charge scale, so reaction events —
  /// which share charge between ends differently than beam events do — no
  /// longer sawtooth between even and odd strips.
  Double_t gain = -1.0;
  /// Short channels only: `C_short / C_long` from this subfile's own ridge fit,
  /// or `0` when the ridge was not measurable and a fallback was used.
  ///
  /// This is a preamplifier gain ratio and therefore fixed per channel, which
  /// is what lets AggregateRidgeRatiosForRun() replace each subfile's value
  /// with the run-level median over the subfiles that did measure it.
  Double_t ridge_ratio = 0.0;
};

/**
 * @brief Output of the two post-gain calibration steps.
 *
 * **Beam energy window.** Gaussian-fits strip 0 (or strip 17) and reports
 * `mu ± 3*sigma` in calibrated units, so downstream beam selection has one
 * data-driven window rather than a hand-set threshold.
 *
 * **Strip centroid alignment.** Decodes events with the per-channel gains
 * already applied, finds each strip's beam-peak centroid from the eSum
 * histogram, fits a robust degree-3 polynomial trend through the centroids with
 * iterative outlier rejection, and derives a multiplicative per-strip
 * `factor = reference / centroid` that pulls every strip onto that smooth
 * trend.
 *
 * @note Both steps run **after** the anchor reduction and after the initial
 *       calibration tree is written — the alignment step needs that tree on
 *       disk so EnergyView can decode. The factors are then stored back into
 *       the calibration tree and applied by EnergyView after the per-channel
 *       gain.
 */
struct StripAlignmentResult {
  Bool_t ok = kFALSE; ///< Whether the steps produced usable results.
  Double_t beam_e_min =
      0.0; ///< Lower edge of the beam window, calibrated units.
  Double_t beam_e_max = 0.0; ///< Upper edge of the beam window.
  /// Per-strip multiplicative alignment factors. Identity (1.0) when no
  /// centroid measurements were available, so an unaligned dataset decodes
  /// unchanged rather than being scaled by zero.
  Double_t factors[18] = {};
  /// Measured beam-peak centroid per strip, in calibrated units.
  Double_t centroids[18] = {};
};

/**
 * @brief Per-channel beam-peak calibration, gain matching and strip alignment.
 *
 * Runs after the pipeline over its events files. For each subfile it fits every
 * channel's beam peak, matches the left and right gains onto a common charge
 * scale, then derives the per-strip alignment that pulls all strips onto one
 * smooth trend.
 *
 * Some results are per-channel constants rather than per-subfile measurements
 * — the ridge ratios and the energy-resolution TOML — and those are aggregated
 * to the run level once every subfile has been processed.
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
   * @brief Save the per-channel calibration overlay figure.
   * @param chans        Calibrated channels.
   * @param samples      Sampled spectra, one vector per channel.
   * @param plot_subdir  Plot subdirectory under the plots base.
   * @param file_label   Subfile label, used in the plot name.
   */
  static void
  SaveCalibSampleOverlay(const std::vector<ChannelCal> &chans,
                         const std::vector<std::vector<Float_t>> &samples,
                         const TString &plot_subdir, const TString &file_label);

  /**
   * @brief Replace each subfile's ridge ratio with the run-level median.
   *
   * The ratio is a preamplifier gain ratio and so is fixed per channel. Taking
   * the median across the subfiles that measured it makes the value robust for
   * subfiles where the ridge was not measurable.
   *
   * @param run   Run number.
   * @param specs That run's subfiles.
   */
  static void AggregateRidgeRatiosForRun(Int_t run,
                                         const std::vector<FileSpec> &specs);
  /**
   * @brief Aggregate the run's energy-resolution measurements into its TOML.
   * @param run   Run number.
   * @param specs That run's subfiles.
   */
  static void AggregateEresTomlForRun(Int_t run,
                                      const std::vector<FileSpec> &specs);

  /**
   * @brief Calibrate every configured subfile, then aggregate per run.
   * @param file_label Restrict to one subfile by its FileSet::FileLabel().
   *                   Empty, the default, processes all of them.
   */
  static void Run(const TString &file_label = "");
};

#endif
