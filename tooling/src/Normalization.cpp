#include "Normalization.hpp"
#include <TFile.h>

Bool_t EnergyView::Attach(TTree *t) {
  tree_ = t;
  t->SetBranchAddress("Left_0_17_dE", left_0_17_adc);
  t->SetBranchAddress("RightdE", rightdE_adc);
  t->SetBranchAddress("Hits", hits_adc);
  t->SetBranchAddress("Cathode", &cathode_adc);
  t->SetBranchAddress("Grid", &grid_adc);
  // Materialize the first tree so GetCurrentFile() resolves for a TChain (it is
  // null until a tree is loaded); this makes is_normed correct right after
  // Attach.
  t->LoadTree(0);
  LoadGains();
  loaded_tree_ = t->GetTreeNumber();
  return is_normed;
}

void EnergyView::LoadGains() {
  for (Int_t s = 0; s < 18; s++) {
    gain_left[s] = 0.0f;
    gain_right[s] = 0.0f;
    slope[s] = 1.0f;
    intercept[s] = 0.0f;
  }
  gain_cathode = 0.0f;
  is_normed = kFALSE;
  if (!tree_)
    return;
  TFile *f = tree_->GetCurrentFile();
  if (!f)
    return;
  TTree *cal = static_cast<TTree *>(f->Get("calibration"));
  if (!cal || cal->GetEntries() < 1)
    return;
  Float_t gl[18] = {0}, gr[18] = {0}, gc = 0.0f;
  Float_t sl[18] = {0}, ic[18] = {0};
  cal->SetBranchAddress("GainLeft", gl);
  cal->SetBranchAddress("GainRight", gr);
  cal->SetBranchAddress("GainCathode", &gc);
  // StripSlope/StripIntercept are optional (written by the second
  // calibration write; older files carry StripFactor instead). When absent
  // the default slope = 1 / intercept = 0 (identity) is used.
  Bool_t has_linear = cal->GetBranch("StripSlope") != nullptr;
  if (has_linear) {
    cal->SetBranchAddress("StripSlope", sl);
    cal->SetBranchAddress("StripIntercept", ic);
  } else if (cal->GetBranch("StripFactor")) {
    // Legacy multiplicative factors -> linear form (factor * total + 0).
    cal->SetBranchAddress("StripFactor", sl);
  }
  cal->GetEntry(0);
  for (Int_t s = 0; s < 18; s++) {
    gain_left[s] = gl[s];
    gain_right[s] = gr[s];
    if (has_linear) {
      slope[s] = sl[s];
      intercept[s] = ic[s];
    } else if (sl[s] > 0.0f) {
      slope[s] = sl[s]; // legacy factor
      intercept[s] = 0.0f;
    }
  }
  gain_cathode = gc;
  is_normed = kTRUE;
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
  if (is_normed) {
    for (Int_t s = 0; s < 18; s++) {
      left[s] = Double_t(gain_left[s]) * Double_t(left_0_17_adc[s]);
      right[s] = Double_t(gain_right[s]) * Double_t(rightdE_adc[s]);
      total[s] = left[s] + right[s];
    }
    // Guard the -1 "no cathode" sentinel: uncalibrated/absent -> 0 a.u.
    cathode = (cathode_adc > 0) ? Double_t(gain_cathode) * Double_t(cathode_adc)
                                : 0.0;
    // Grid has no calibration gain; normalize to [0, 1] by dividing by max ADC.
    grid = Double_t(grid_adc) / 16384.0;
  } else {
    for (Int_t s = 0; s < 18; s++) {
      left[s] = Double_t(left_0_17_adc[s]);
      right[s] = Double_t(rightdE_adc[s]);
      total[s] = left[s] + right[s];
    }
    cathode = Double_t(cathode_adc);
    grid = Double_t(grid_adc);
  }
  if (Constants::cfg.IGNORE_SHORT_STRIPS) {
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
  // Apply the per-strip two-point linear correction (beam centroid -> 1.0,
  // pileup centroid -> 2.0). Applied after per-channel gain and
  // IGNORE_SHORT_STRIPS. This corrects the ADC sublinearity that a single
  // multiplicative factor cannot (with beam at 1.0 the pileup would read
  // ~1.88 instead of 2.0).
  if (is_normed) {
    for (Int_t s = 0; s <= 17; s++) {
      total[s] = Double_t(slope[s]) * total[s] + Double_t(intercept[s]);
    }
  }
}

const char *EnergyView::Unit() const { return is_normed ? "a.u." : "ADC"; }
