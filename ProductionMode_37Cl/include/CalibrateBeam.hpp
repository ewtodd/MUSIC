#ifndef CALIBRATE_BEAM_HPP
#define CALIBRATE_BEAM_HPP

#include "BeamFit2D.hpp"
#include "Constants.hpp"
#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Normalization.hpp"
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
#include <TList.h>
#include <TMath.h>
#include <TROOT.h>
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

const Int_t kNPeaks = 2;

// Per-channel calibration: sim anchor (MeV, from beam sim file) plus the
// per-subfile ADC anchors and resulting pol2 fit.
struct ChannelCal {
  TString name;
  Char_t side;
  Int_t strip;
  Double_t sim_mu_mev = 0.0;
  Double_t sim_sigma_mev = 0.0;
  Double_t anchor_adc[kNPeaks] = {0, 0};
  Double_t anchor_sigma_adc[kNPeaks] = {0, 0};
  Long64_t n_samples[kNPeaks] = {0, 0};
  Double_t lin_a = 0, lin_b = 0;
  Double_t lin_chi2_ndf = -1;
  Bool_t lin_ok = kFALSE;
};

class CalibrateBeam {
public:
  static TString DefaultSimBeamPath(const TString &project_root);

  static std::vector<ChannelCal> LoadSimChans(const TString &sim_path);

  static void CalibrateBeamOneSubfile(const FileSpec &spec,
                                      const std::vector<ChannelCal> &sim_chans);

  // Overlay (one color per channel, log-y, MeV) of ONLY the events used for
  // calibration: the k=1 + k=2 anchor samples. Each beam-dE channel shows its
  // beam and single-pileup peak. Same axes as the dynamic_range_check plot but
  // restricted to calibration events.
  static void SaveCalibSampleOverlay(
      const std::vector<ChannelCal> &chans,
      const std::vector<std::vector<std::vector<Float_t>>> &samples,
      const TString &plot_subdir, const TString &file_label);

  static void AggregateEresTomlForRun(Int_t run,
                                      const std::vector<FileSpec> &specs);

  static void Run(const TString &file_label = "");
};

#endif
