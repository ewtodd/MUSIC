#ifndef BEAM_FIT_2D_HPP
#define BEAM_FIT_2D_HPP

#include <Rtypes.h>
#include <TChain.h>
#include <TEllipse.h>
#include <TF2.h>
#include <TFitResult.h>
#include <TH2F.h>
#include <TMath.h>
#include <TString.h>
#include <algorithm>
#include <cmath>
#include <iostream>

struct BeamFit2D {
  Double_t amp = 0;
  Double_t mu_x = 0, mu_y = 0;
  Double_t sigma_x = 0, sigma_y = 0;
  Double_t rho = 0;
  Bool_t ok = kFALSE;
};

struct Moments2D {
  Double_t mu_x = 0, mu_y = 0;
  Double_t sigma_x = 0, sigma_y = 0;
  Double_t rho = 0;
  Double_t weight = 0;
};

class BeamFitUtils {
public:
  static void DiagonalizeCov(const BeamFit2D &b, Double_t &major,
                             Double_t &minor, Double_t &theta_deg);

  static Bool_t InEllipse(const BeamFit2D &b, Double_t x, Double_t y,
                          Double_t n_sigma);

  // Axis-aligned ellipse gate with independent horizontal/vertical
  // half-widths (in sigma_x / sigma_y units); ignores rho.
  static Bool_t InEllipseXY(const BeamFit2D &b, Double_t x, Double_t y,
                            Double_t nx, Double_t ny);

  // Weighted mean+covariance over bins above `thresh` in the bin box.
  static Moments2D ComputeMoments(TH2F *h, Int_t lo_bx, Int_t hi_bx,
                                  Int_t lo_by, Int_t hi_by, Double_t thresh,
                                  Double_t bw_x, Double_t bw_y);

  // bigaus fit inside an explicit window with explicit seed. Falls back to the
  // seed on a failed fit so callers can still gate.
  static BeamFit2D FitBigausInWindow(TH2F *h, Double_t x_lo, Double_t x_hi,
                                     Double_t y_lo, Double_t y_hi,
                                     const Moments2D &seed, Double_t seed_amp,
                                     const TString &fname, const TString &tag);

  // Shared Strip0-vs-Strip1-long-side beam gate (delta-e-scatter, gate-cache):
  // GateNSigma() is the ellipse half-width in sigma; BuildGateHist fills the 2D
  // MeV histogram; FitBigPeak moments-seeds a bigaus on its largest peak.
  static Double_t GateNSigma();
  static TH2F *BuildGateHist(TChain *chain, const TString &name);
  static BeamFit2D FitBigPeak(TH2F *h, const TString &tag);

  // Draw the diagonalized n_sigma gate ellipse on the current pad (no-op if the
  // fit failed). Caller has already drawn the gate histogram.
  static void DrawGateEllipse(const BeamFit2D &gate, Double_t n_sigma);
};

#endif
