#include "Normalization.hpp"
#include <TFile.h>
#include <iostream>

Bool_t EnergyView::Attach(TTree *t) {
  tree_ = t;
  if (!t->GetBranch("LeftdE") || !t->GetBranch("Strip0dE"))
    std::cerr << "EnergyView: events tree has no LeftdE/Strip0dE branches; it "
                 "predates the per-end layout and must be rebuilt through the "
                 "pipeline"
              << std::endl;
  t->SetBranchAddress("LeftdE", leftdE_adc);
  t->SetBranchAddress("RightdE", rightdE_adc);
  t->SetBranchAddress("Strip0dE", &strip0_adc);
  t->SetBranchAddress("Strip17dE", &strip17_adc);
  t->SetBranchAddress("Hits", hits_adc);
  t->SetBranchAddress("Cathode", &cathode_adc);
  t->SetBranchAddress("Grid", &grid_adc);
  // Materialize the first tree so GetCurrentFile() resolves for a TChain (null
  // until a tree is loaded); this makes is_normed correct after Attach.
  t->LoadTree(0);
  LoadGains();
  loaded_tree_ = t->GetTreeNumber();
  return is_normed;
}

void EnergyView::LoadGains() {
  for (Int_t k = 0; k < 16; k++) {
    gain_left[k] = 0.0f;
    gain_right[k] = 0.0f;
  }
  for (Int_t s = 0; s < 18; s++)
    strip_factor[s] = 1.0f;
  gain_strip0 = 0.0f;
  gain_strip17 = 0.0f;
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
  Float_t gl[16] = {0}, gr[16] = {0}, g0 = 0.0f, g17 = 0.0f, gc = 0.0f;
  Float_t sf[18] = {0};
  cal->SetBranchAddress("GainLeft", gl);
  cal->SetBranchAddress("GainRight", gr);
  cal->SetBranchAddress("GainStrip0", &g0);
  cal->SetBranchAddress("GainStrip17", &g17);
  cal->SetBranchAddress("GainCathode", &gc);
  // StripFactor is optional (added after the initial gain write, filled in by
  // FindStripCentroidAlignment on the second); absent -> default 1.0 (identity)
  Bool_t has_factor = cal->GetBranch("StripFactor") != nullptr;
  if (has_factor)
    cal->SetBranchAddress("StripFactor", sf);
  cal->GetEntry(0);
  for (Int_t k = 0; k < 16; k++) {
    gain_left[k] = gl[k];
    gain_right[k] = gr[k];
  }
  if (has_factor)
    for (Int_t s = 0; s < 18; s++)
      strip_factor[s] = sf[s];
  gain_strip0 = g0;
  gain_strip17 = g17;
  gain_cathode = gc;
  is_normed = kTRUE;
}

void EnergyView::Decode() {
  // A TChain advances the tree number when GetEntry crosses into a new file,
  // so reload its gains; plain TTrees report -1 forever and never re-fire.
  if (tree_) {
    Int_t tn = tree_->GetTreeNumber();
    if (tn != loaded_tree_) {
      LoadGains();
      loaded_tree_ = tn;
    }
  }
  if (is_normed) {
    for (Int_t k = 0; k < 16; k++) {
      left[k] = Double_t(gain_left[k]) * Double_t(leftdE_adc[k]);
      right[k] = Double_t(gain_right[k]) * Double_t(rightdE_adc[k]);
    }
    strip0 = Double_t(gain_strip0) * Double_t(strip0_adc);
    strip17 = Double_t(gain_strip17) * Double_t(strip17_adc);
    // The -1 "no cathode" sentinel: uncalibrated/absent -> 0 a.u.
    cathode = (cathode_adc > 0) ? Double_t(gain_cathode) * Double_t(cathode_adc)
                                : 0.0;
    // Grid has no calibration gain; normalize to [0, 1] by dividing by max ADC.
    grid = Double_t(grid_adc) / 16384.0;
  } else {
    for (Int_t k = 0; k < 16; k++) {
      left[k] = Double_t(leftdE_adc[k]);
      right[k] = Double_t(rightdE_adc[k]);
    }
    strip0 = Double_t(strip0_adc);
    strip17 = Double_t(strip17_adc);
    cathode = Double_t(cathode_adc);
    grid = Double_t(grid_adc);
  }
  // Long end only: L on odd strips, R on even. Zeroing the short end keeps
  // Total() the plain sum of the ends.
  if (Constants::cfg.IGNORE_SHORT_STRIPS) {
    for (Int_t s = 1; s <= 16; s++) {
      if ((s % 2) != 0)
        right[s - 1] = 0.0;
      else
        left[s - 1] = 0.0;
    }
  }
  // Per-strip multiplicative alignment: pulls each strip's peak centroid onto
  // the pol3 trend. On every end of the strip, so the sum carries it.
  if (is_normed) {
    for (Int_t s = 1; s <= 16; s++) {
      left[s - 1] *= Double_t(strip_factor[s]);
      right[s - 1] *= Double_t(strip_factor[s]);
    }
    strip0 *= Double_t(strip_factor[0]);
    strip17 *= Double_t(strip_factor[17]);
  }
}

const char *EnergyView::Unit() const { return is_normed ? "a.u." : "ADC"; }
