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
    // Identity transformation by default — no alignment applied.
    strip_slope[s] = 1.0f;
    strip_intercept[s] = 0.0f;
    gain_long_impute[s] = 0.0f;
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
  Float_t ss[18] = {0}, si[18] = {0}, gi[18] = {0};
  cal->SetBranchAddress("GainLeft", gl);
  cal->SetBranchAddress("GainRight", gr);
  cal->SetBranchAddress("GainCathode", &gc);
  // StripSlope/StripIntercept/GainLongImpute are optional (added by the
  // calibration passes). If absent, the defaults remain: identity alignment
  // and no short-side imputation.
  Bool_t has_slope = cal->GetBranch("StripSlope") != nullptr;
  Bool_t has_intercept = cal->GetBranch("StripIntercept") != nullptr;
  Bool_t has_impute = cal->GetBranch("GainLongImpute") != nullptr;
  if (has_slope)
    cal->SetBranchAddress("StripSlope", ss);
  if (has_intercept)
    cal->SetBranchAddress("StripIntercept", si);
  if (has_impute)
    cal->SetBranchAddress("GainLongImpute", gi);
  cal->GetEntry(0);
  for (Int_t s = 0; s < 18; s++) {
    gain_left[s] = gl[s];
    gain_right[s] = gr[s];
    if (has_slope)
      strip_slope[s] = ss[s];
    if (has_intercept)
      strip_intercept[s] = si[s];
    if (has_impute)
      gain_long_impute[s] = gi[s];
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
    // Short-side imputation: the short end's trigger threshold sits near its
    // beam deposit, so it fires probabilistically (strongly parity-dependent).
    // When it reads zero, decode the total from the long side alone with the
    // effective gain that restores the expected short contribution
    // proportionally — otherwise the total is bimodal per strip and traces
    // sawtooth between even/odd strips. Long side: L on odd strips, R on
    // even (same parity rule as IGNORE_SHORT_STRIPS below).
    for (Int_t s = 1; s <= 16; s++) {
      if (gain_long_impute[s] <= 0.0f)
        continue;
      Bool_t l_is_long = ((s % 2) != 0);
      UShort_t short_adc = l_is_long ? rightdE_adc[s] : left_0_17_adc[s];
      if (short_adc == 0) {
        UShort_t long_adc = l_is_long ? left_0_17_adc[s] : rightdE_adc[s];
        total[s] = Double_t(gain_long_impute[s]) * Double_t(long_adc);
      }
    }
    // Guard the -1 "no cathode" sentinel: uncalibrated/absent -> 0 a.u.
    // (matches the old per-event cal, which only applied the gain when
    // cathode_adc > 0).
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
  // Apply per-strip linear alignment (slope + intercept) to the strip total.
  // This flattens the Bragg curve to 1.0 a.u. using the beam + pileup peaks
  // (computed in FindStripCentroidAlignment). Applied after per-channel gain
  // and IGNORE_SHORT_STRIPS so it operates on the final total. Strips 0 and 17
  // are single-ended anode strips (not guards) and are aligned identically.
  if (is_normed) {
    for (Int_t s = 0; s <= 17; s++) {
      total[s] = strip_slope[s] * total[s] + strip_intercept[s];
    }
  }
}

const char *EnergyView::Unit() const { return is_normed ? "a.u." : "ADC"; }
