#include "Normalization.hpp"

Bool_t EnergyView::Attach(TTree *t) {
  is_mev = (t->FindBranch("TotaldEMeV") != nullptr);
  if (is_mev) {
    t->SetBranchAddress("LeftdEMeV", leftdE_mev);
    t->SetBranchAddress("RightdEMeV", rightdE_mev);
    t->SetBranchAddress("TotaldEMeV", totaldE_mev);
    t->SetBranchAddress("CathodeMeV", &cathode_mev);
  } else {
    t->SetBranchAddress("LeftdE", leftdE_adc);
    t->SetBranchAddress("RightdE", rightdE_adc);
    t->SetBranchAddress("TotaldE", totaldE_adc);
    t->SetBranchAddress("Cathode", &cathode_adc);
  }
  return is_mev;
}

void EnergyView::Decode() {
  if (is_mev) {
    for (Int_t s = 0; s < 18; s++) {
      left[s] = Double_t(leftdE_mev[s]);
      right[s] = Double_t(rightdE_mev[s]);
      total[s] = Double_t(totaldE_mev[s]);
    }
    cathode = Double_t(cathode_mev);
  } else {
    for (Int_t s = 0; s < 18; s++) {
      left[s] = Double_t(leftdE_adc[s]);
      right[s] = Double_t(rightdE_adc[s]);
      total[s] = Double_t(totaldE_adc[s]);
    }
    cathode = Double_t(cathode_adc);
  }
  if (Constants::IGNORE_SHORT_STRIPS && !is_mev) {
    for (Int_t s = 1; s <= 16; s++) {
      if ((s % 2) != 0) {
        total[s] = left[s];
        right[s] = 0.0;
      } else {
        total[s] = right[s];
        left[s] = 0.0;
      }
    }
  }
  if (Constants::IGNORE_STRIP_0) {
    left[0] = 0.0;
    right[0] = 0.0;
    total[0] = 0.0;
  }
  if (Constants::IGNORE_STRIP_17) {
    left[17] = 0.0;
    right[17] = 0.0;
    total[17] = 0.0;
  }
}

const char *EnergyView::Unit() const { return is_mev ? "MeV" : "ADC"; }
