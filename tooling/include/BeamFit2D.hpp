#ifndef BEAM_FIT_2D_HPP
#define BEAM_FIT_2D_HPP

#include <Rtypes.h>
#include <TH2F.h>
#include <TString.h>
#include <algorithm>
#include <cmath>

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
  static Bool_t InEllipseXY(const BeamFit2D &b, Double_t x, Double_t y,
                            Double_t nx, Double_t ny);

  static Moments2D ComputeMoments(TH2F *h, Int_t lo_bx, Int_t hi_bx,
                                  Int_t lo_by, Int_t hi_by, Double_t thresh,
                                  Double_t bw_x, Double_t bw_y);
};

#endif
