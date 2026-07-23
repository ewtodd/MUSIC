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

struct ChannelCal {
  TString name;
  Char_t side;
  Int_t strip;
  Double_t fit_adc = 0.0;
  Double_t fit_sigma_adc = 0.0;
  Long64_t n_samples = 0.0;
  // L/R gain-match override. When >= 0, Gain() returns this instead of
  // 1/fit_adc. Set by ComputeLRGainMatch (the check_LR notebook's two-pass
  // recipe): long side anchored on its beam peak, short side anchored on its
  // "shoulder" (its response in events where the long side reads low, i.e.
  // the charge went mostly to the short end), then a per-strip eSum median
  // alignment applied to the short side only. Puts L and R in the same
  // charge scale so reaction events (different L/R sharing than beam) don't
  // sawtooth between even/odd strips. When < 0, falls back to 1/fit_adc.
  Double_t gain = -1.0;
  // Effective gain for decoding the strip total from the LONG side alone when
  // the short side did not fire (its trigger threshold sits near the
  // short-side beam deposit, so it fires probabilistically — adding zero in
  // those events makes the total bimodal and the traces sawtooth). Set only
  // on the LONG channel of each strip:
  //   impute_gain = gain_long + gain_short * (short_anchor / long_anchor)
  // so a beam event decodes to the same total whether or not the short end
  // fired, and reaction events get a proportional first-order estimate.
  // <= 0 disables imputation for the strip.
  Double_t impute_gain = 0.0;
};

// Result of the two post-gain calibration steps:
//   1. DeriveBeamEnergyWindow — Gaussian-fit Strip0 (or Strip17) and report
//      mu ± 3*sigma in a.u. so downstream beam-selection code has a single,
//      data-driven beam-energy window.
//   2. FindStripCentroidAlignment — decode events with the per-channel gains
//      already applied, find each strip's beam-peak AND pileup-peak centroids,
//      fit a smooth polynomial reference trend (for outlier rejection only),
//      and derive a per-strip LINEAR correction (slope + intercept) that maps
//      the beam peak to 1.0 a.u. and the pileup peak to 2.0 a.u. This
//      flattens the Bragg curve to 1.0 and corrects per-strip nonlinearity
//      that a single-point (gain-only) calibration can't.
//
// Both steps run AFTER ReduceToAnchors and AFTER the initial
// WriteCalibrationToEvents (the alignment step needs the calibration tree on
// disk so EnergyView can decode). The slope/intercept are stored in the
// calibration tree and applied by EnergyView after the per-channel gain.
struct StripAlignmentResult {
  Bool_t ok = kFALSE;
  Double_t beam_e_min = 0.0;
  Double_t beam_e_max = 0.0;
  // Identity by default (slope=1, intercept=0) so that an uncomputed result
  // is a no-op when written to the calibration tree and applied by EnergyView.
  // If these were 0, EnergyView would zero every strip total.
  Double_t slopes[18] = {};     // set to 1.0 in FindStripCentroidAlignment
  Double_t intercepts[18] = {}; // 0.0 = no offset
  Double_t centroids[18] = {};  // measured beam-peak centroid per strip (a.u.)
  Double_t pileup_centroids[18] =
      {}; // measured pileup-peak centroid per strip (a.u.)
};

class CalibrateBeam {
public:
  static std::vector<ChannelCal> BuildChannels();

  static void
  CalibrateBeamOneSubfile(const FileSpec &spec,
                          const std::vector<ChannelCal> &chans_template);

  static void
  SaveCalibSampleOverlay(const std::vector<ChannelCal> &chans,
                         const std::vector<std::vector<Float_t>> &samples,
                         const TString &plot_subdir, const TString &file_label);

  static void AggregateEresTomlForRun(Int_t run,
                                      const std::vector<FileSpec> &specs);

  static void Run(const TString &file_label = "");
};

#endif
