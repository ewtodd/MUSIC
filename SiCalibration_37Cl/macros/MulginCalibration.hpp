#ifndef MULGIN_CALIBRATION_HPP
#define MULGIN_CALIBRATION_HPP

#include <TMath.h>
#include <TROOT.h>
#include <TString.h>

//  Mulgin, S I, Okolovich, V N, and Zhdanov, S V. Two-parametric method for
//  silicon detector calibration in heavy ion and fission fragment spectrometry.
//  Netherlands: N. p., 1997. Web. doi:10.1016/S0168-9002(96)01211-9.

//   P = B * E_det(E, gamma) + C
//   E_det(E, gamma) = E - 0.55*E/(1 + 0.37568*E) - gamma*E
// gamma is the effective linear PHD term (= 37*alpha + beta).

inline Double_t MulginEffectiveDetected(Double_t E, Double_t gamma) {
  return E - 0.55 * E / (1.0 + 0.37568 * E) - gamma * E;
}

inline Double_t MulginPredictChannel(Double_t E, Double_t B, Double_t C,
                                     Double_t gamma) {
  return B * MulginEffectiveDetected(E, gamma) + C;
}

// Numerical inverse: ADC channel → true energy (MeV).
inline Double_t MulginInverseEnergy(Double_t adc, Double_t B, Double_t C,
                                    Double_t gamma) {
  Double_t E_det_target = (adc - C) / B;
  Double_t lo = 0.0;
  Double_t hi = 500.0;
  for (Int_t k = 0; k < 200; k++) {
    Double_t mid = 0.5 * (lo + hi);
    Double_t f = MulginEffectiveDetected(mid, gamma) - E_det_target;
    if (f > 0.0) {
      hi = mid;
    } else {
      lo = mid;
    }
  }
  return 0.5 * (lo + hi);
}

// Propagate adc → E uncertainty via dE/d(adc) = 1 / (B * dE_det/dE).
inline Double_t MulginInverseEnergyError(Double_t adc, Double_t B, Double_t C,
                                         Double_t gamma, Double_t adc_err) {
  Double_t E = MulginInverseEnergy(adc, B, C, gamma);
  Double_t denom = 1.0 + 0.37568 * E;
  Double_t dE_det_dE = 1.0 - 0.55 / (denom * denom) - gamma;
  Double_t dadc_dE = B * dE_det_dE;
  if (dadc_dE == 0.0)
    return 0.0;
  return TMath::Abs(adc_err / dadc_dE);
}

// Gain group assignment: gain shift occurred between runs 22 and 23.
inline Bool_t MulginIsGroupA(Int_t run) { return run <= 22; }

#endif
