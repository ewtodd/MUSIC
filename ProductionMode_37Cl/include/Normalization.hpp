#ifndef NORMALIZATION_HPP
#define NORMALIZATION_HPP

#include "Constants.hpp"
#include <Rtypes.h>
#include <TTree.h>

struct EnergyView {
  Bool_t is_mev;

  Int_t leftdE_adc[18], rightdE_adc[18], totaldE_adc[18];
  Int_t cathode_adc;

  Float_t leftdE_mev[18], rightdE_mev[18], totaldE_mev[18];
  Float_t cathode_mev;

  Double_t left[18], right[18], total[18];
  Double_t cathode;

  Bool_t Attach(TTree *t);
  void Decode();
  const char *Unit() const;
};

#endif
