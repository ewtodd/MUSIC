#include "Normalization.hpp"
#include <TFile.h>

Bool_t EnergyView::Attach(TTree *t) {
  tree_ = t;
  t->SetBranchAddress("Left_0_17_dE", left_0_17_adc);
  t->SetBranchAddress("RightdE", rightdE_adc);
  t->SetBranchAddress("Cathode", &cathode_adc);
  // Materialize the first tree so GetCurrentFile() resolves for a TChain (it is
  // null until a tree is loaded); this makes is_mev correct right after Attach.
  t->LoadTree(0);
  LoadGains();
  loaded_tree_ = t->GetTreeNumber();
  return is_mev;
}

void EnergyView::LoadGains() {
  for (Int_t s = 0; s < 18; s++) {
    gain_left[s] = 0.0f;
    gain_right[s] = 0.0f;
  }
  gain_cathode = 0.0f;
  is_mev = kFALSE;
  if (!tree_)
    return;
  TFile *f = tree_->GetCurrentFile();
  if (!f)
    return;
  TTree *cal = static_cast<TTree *>(f->Get("calibration"));
  if (!cal || cal->GetEntries() < 1)
    return;
  Float_t gl[18] = {0}, gr[18] = {0}, gc = 0.0f;
  cal->SetBranchAddress("GainLeft", gl);
  cal->SetBranchAddress("GainRight", gr);
  cal->SetBranchAddress("GainCathode", &gc);
  cal->GetEntry(0);
  for (Int_t s = 0; s < 18; s++) {
    gain_left[s] = gl[s];
    gain_right[s] = gr[s];
  }
  gain_cathode = gc;
  is_mev = kTRUE;
}

void EnergyView::Decode() {
  // A TChain advances fTreeNumber when GetEntry crosses into a new subfile;
  // reload that file's gains when it does. Plain TTrees report -1 forever, so
  // this never re-fires for them (gains already loaded in Attach).
  if (tree_) {
    Int_t tn = tree_->GetTreeNumber();
    if (tn != loaded_tree_) {
      LoadGains();
      loaded_tree_ = tn;
    }
  }
  if (is_mev) {
    for (Int_t s = 0; s < 18; s++) {
      left[s] = Double_t(gain_left[s]) * Double_t(left_0_17_adc[s]);
      right[s] = Double_t(gain_right[s]) * Double_t(rightdE_adc[s]);
      total[s] = left[s] + right[s];
    }
    // Guard the -1 "no cathode" sentinel: uncalibrated/absent -> 0 MeV (matches
    // the old per-event cal, which only applied the gain when cathode_adc > 0).
    cathode = (cathode_adc > 0) ? Double_t(gain_cathode) * Double_t(cathode_adc)
                                : 0.0;
  } else {
    for (Int_t s = 0; s < 18; s++) {
      left[s] = Double_t(left_0_17_adc[s]);
      right[s] = Double_t(rightdE_adc[s]);
      total[s] = left[s] + right[s];
    }
    cathode = Double_t(cathode_adc);
  }
  if (Constants::IGNORE_SHORT_STRIPS) {
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
