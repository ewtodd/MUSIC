#ifndef CALIBRATE_BEAM_INTERNAL_HPP
#define CALIBRATE_BEAM_INTERNAL_HPP

#include "CalibrateBeam.hpp"

// Beam-gate ellipse half-widths (in sigma) and the minimum sample count a
// channel needs before its beam-peak fit is trusted.
const Double_t kEllipseNSigmaX = 3;
const Double_t kEllipseNSigmaY = 3;
const Long64_t kMinSamples = 200;

inline Char_t LongSide(Int_t strip) { return (strip % 2 == 0) ? 'R' : 'L'; }

inline Bool_t IsBeamdEChannel(const ChannelCal &c) {
  if (c.side == 'S')
    return kTRUE;
  if (c.side != 'L' && c.side != 'R')
    return kFALSE;
  if (c.strip < 1 || c.strip > 16)
    return kFALSE;
  return c.side == LongSide(c.strip);
}

inline Bool_t IsCalibrated(const ChannelCal &c) { return c.fit_adc > 0; }

inline Double_t Gain(const ChannelCal &c) {
  return c.gain >= 0.0 ? c.gain : 1.0 / c.fit_adc;
}

inline Double_t ResolutionFWHMPercent(const ChannelCal &c) {
  if (c.fit_adc <= 0)
    return 0.0;
  const Double_t kFwhmPerSigma = 2.0 * TMath::Sqrt(2.0 * TMath::Log(2.0));
  return 100.0 * kFwhmPerSigma * c.fit_sigma_adc / c.fit_adc;
}

inline Double_t ApplyCal(const ChannelCal &c, Double_t adc) {
  return Gain(c) * adc;
}

// Beam-peak histogram + fit helpers (CalibrateBeamFits.cpp).
Double_t Median(const std::vector<Float_t> &v);
Double_t InterquartileRange(const std::vector<Float_t> &v);
void RobustPeakSeed(const std::vector<Float_t> &v, Double_t &mode,
                    Double_t &sigma);
TH1F *MakeBeamPeakHist(const TString &name, const TString &title,
                       const std::vector<Float_t> &v, Double_t mode,
                       Double_t sigma);
Bool_t FitBeamPeakGaussian(const std::vector<Float_t> &v, const TString &fname,
                           Double_t &peak_adc, Double_t &sigma_adc,
                           TF1 *&fit_out);

// Beam gate fits (CalibrateBeamGates.cpp).
BeamFit2D FindBeamGateStp2VsStp1(const FileSpec &spec, const TString &run_label,
                                 const TString &plot_subdir,
                                 Bool_t save_plot = kTRUE);
BeamFit2D FindBeamGateStp0VsGrid(const FileSpec &spec, const TString &run_label,
                                 const TString &plot_subdir,
                                 const BeamFit2D &beam,
                                 Bool_t save_plot = kTRUE);

// L/R gain matching (CalibrateBeamLRGain.cpp).
void CollectAnchorSamplesOneSubfile(const FileSpec &spec,
                                    const std::vector<ChannelCal> &chans,
                                    const BeamFit2D &beam,
                                    const BeamFit2D &beam0vGrid,
                                    std::vector<std::vector<Float_t>> &samples);
void ReduceToAnchors(std::vector<ChannelCal> &chans,
                     std::vector<std::vector<Float_t>> &samples,
                     std::vector<TF1 *> &fits_out, const TString &run_label);

// Calibrated-spectrum overlays (CalibrateBeamOverlay.cpp).
void SaveDynamicRangeOverlay(const FileSpec &spec,
                             const std::vector<ChannelCal> &chans,
                             const TString &plot_subdir,
                             const TString &file_label);

#endif
