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
  // recipe): long side anchored on its beam peak, short side anchored on
  // the correlation-ridge peak, then a per-strip eSum median alignment
  // applied to the short side only. Puts L and R in the same charge scale
  // so reaction events (different L/R sharing than beam) don't sawtooth
  // between even/odd strips. When < 0, falls back to 1/fit_adc.
  Double_t gain = -1.0;
};

// Result of the two post-gain calibration steps:
//   1. DeriveBeamEnergyWindow — Gaussian-fit Strip0 (or Strip17) and report
//      mu ± 3*sigma in a.u. so downstream beam-selection code has a single,
//      data-driven beam-energy window.
//   2. FindStripCentroidAlignment — decode events with the per-channel gains
//      already applied, find each strip's beam-peak centroid (B) and pileup-
//      peak centroid (P) from the eSum histogram (the pileup population is
//      isolated by gating on the double-beam total energy), and derive a
//      two-point LINEAR normalization per strip: total_corrected =
//      slope * total + intercept, with the line passing through (B, 1.0) and
//      (P, 2.0). The two-point line corrects the ADC sublinearity — the
//      pileup reads ~1.88x instead of 2x — that a single multiplicative
//      factor cannot (a factor maps beam onto 1.0 and leaves the pileup at
//      1.88). This matches the notebook's approach (37Cl_an.ipynb cell 5)
//      extended with the second point.
//
// Both steps run AFTER ReduceToAnchors and AFTER the initial
// WriteCalibrationToEvents (the alignment step needs the calibration tree on
// disk so EnergyView can decode). The slopes/intercepts are stored in the
// calibration tree and applied by EnergyView after the per-channel gain.
struct StripAlignmentResult {
  Bool_t ok = kFALSE;
  Double_t beam_e_min = 0.0;
  Double_t beam_e_max = 0.0;
  // Default slope = 1.0, intercept = 0.0 (identity — no correction). Set by
  // FindStripCentroidAlignment when centroid measurements are available.
  Double_t slopes[18] = {};
  Double_t intercepts[18] = {};
  Double_t centroids[18] = {};        // measured beam-peak centroid per strip
  Double_t pileup_centroids[18] = {}; // measured pileup-peak centroid per strip
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
