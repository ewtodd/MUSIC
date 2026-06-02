#ifndef CALIBRATION_HPP
#define CALIBRATION_HPP

#include <Rtypes.h>
#include <TMath.h>

//  Mulgin, S I, Okolovich, V N, and Zhdanov, S V. Two-parametric method for
//  silicon detector calibration in heavy ion and fission fragment spectrometry.
//  Netherlands: N. p., 1997. Web. doi:10.1016/S0168-9002(96)01211-9.
//
//   P = B * E_det(E) + C
//   E_det(E) = E - 0.55*E/(1 + 0.37568*E)
//
// Dataset-independent physics; the per-dataset gain-group boundary lives in the
// SiCalib config (SiCalib::GAIN_GROUP_BOUNDARY_RUN), not here.
class Calibration {
public:
  static Double_t MulginEffectiveDetected(Double_t E) {
    return E - 0.55 * E / (1.0 + 0.37568 * E);
  }

  static Double_t MulginPredictChannel(Double_t E, Double_t B, Double_t C) {
    return B * MulginEffectiveDetected(E) + C;
  }

  // Numerical inverse: ADC channel -> true energy (MeV).
  static Double_t MulginInverseEnergy(Double_t adc, Double_t B, Double_t C) {
    Double_t E_det_target = (adc - C) / B;
    Double_t lo = 0.0;
    Double_t hi = 500.0;
    for (Int_t k = 0; k < 200; k++) {
      Double_t mid = 0.5 * (lo + hi);
      Double_t f = MulginEffectiveDetected(mid) - E_det_target;
      if (f > 0.0) {
        hi = mid;
      } else {
        lo = mid;
      }
    }
    return 0.5 * (lo + hi);
  }

  // Propagate adc -> E uncertainty via dE/d(adc) = 1 / (B * dE_det/dE).
  static Double_t MulginInverseEnergyError(Double_t adc, Double_t B, Double_t C,
                                           Double_t adc_err) {
    Double_t E = MulginInverseEnergy(adc, B, C);
    Double_t denom = 1.0 + 0.37568 * E;
    Double_t dE_det_dE = 1.0 - 0.55 / (denom * denom);
    Double_t dadc_dE = B * dE_det_dE;
    if (dadc_dE == 0.0)
      return 0.0;
    return TMath::Abs(adc_err / dadc_dE);
  }
};

#endif
