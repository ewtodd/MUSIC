#ifndef NORMALIZATION_HPP
#define NORMALIZATION_HPP

#include "Constants.hpp"
#include <Rtypes.h>
#include <TTree.h>

struct EnergyView {
  UShort_t left_0_17_adc[18], rightdE_adc[18];
  UShort_t hits_adc[36];
  Short_t cathode_adc;
  Short_t grid_adc;

  Float_t gain_left[18], gain_right[18], gain_cathode;
  Bool_t is_normed; // a calibration tree was found -> Decode() yields a.u.

  // decoded per-event values (a.u. if is_normed, else raw ADC)
  Double_t left[18], right[18], total[18];
  Double_t cathode;
  Double_t grid;

  TTree *tree_;       // events tree/chain we are bound to
  Int_t loaded_tree_; // tree number whose gains are currently loaded (-1 none)

  EnergyView()
      : cathode_adc(0), grid_adc(0), gain_cathode(0.0f), is_normed(kFALSE),
        cathode(0.0), grid(0.0), tree_(nullptr), loaded_tree_(-1) {}

  Bool_t Attach(TTree *t);
  void Decode();
  void LoadGains();
  const char *Unit() const;
};

#endif
