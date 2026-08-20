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
  // L/R gain-match override: when >= 0, Gain() uses this instead of 1/fit_adc.
  // Set by ComputeLRGainMatch (the check_LR two-pass recipe), putting L and R
  // in the same charge scale.
  Double_t gain = -1.0;
};

// Post-gain results: the beam energy window (mu ± 3*sigma of the Gaussian fit
// to Strip0/17, a.u.) and the per-strip two-point line (B,1.0)-(P,2.0) that
// fixes the ADC sublinearity no single factor can; runs after ReduceToAnchors +
// the initial calibration write, and is consumed by EnergyView.
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
