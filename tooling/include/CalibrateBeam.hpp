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
//      already applied, find each strip's beam-peak centroid from the eSum
//      histogram, fit a robust degree-3 polynomial reference trend through
//      the centroids with iterative outlier rejection, and derive a
//      multiplicative per-strip factor = reference / centroid that pulls
//      every strip onto the smooth trend. This matches the notebook's
//      approach (37Cl_an.ipynb cell 5).
//
// Both steps run AFTER ReduceToAnchors and AFTER the initial
// WriteCalibrationToEvents (the alignment step needs the calibration tree on
// disk so EnergyView can decode). The factors are stored in the calibration
// tree and applied by EnergyView after the per-channel gain.
struct StripAlignmentResult {
  Bool_t ok = kFALSE;
  Double_t beam_e_min = 0.0;
  Double_t beam_e_max = 0.0;
  // Default factor = 1.0 (identity — no scaling). Set by
  // FindStripCentroidAlignment when centroid measurements are available.
  Double_t factors[18] = {};
  Double_t centroids[18] = {}; // measured beam-peak centroid per strip (a.u.)
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
