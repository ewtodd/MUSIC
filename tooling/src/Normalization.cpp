#include "Normalization.hpp"
#include <TFile.h>
#include <TMath.h>
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
    offset_left[k] = 0.0f;
    offset_right[k] = 0.0f;
  }
  for (Int_t s = 0; s < 18; s++) {
    strip_factor[s] = 1.0f;
    strip_offset[s] = 0.0f;
  }
  gain_strip0 = 0.0f;
  gain_strip17 = 0.0f;
  gain_cathode = 0.0f;
  gain_grid = 0.0f;
  is_normed = kFALSE;
  if (!tree_)
    return;
  TFile *f = tree_->GetCurrentFile();
  if (!f)
    return;
  TTree *cal = static_cast<TTree *>(f->Get("calibration"));
  if (!cal || cal->GetEntries() < 1)
    return;
  Float_t gl[16] = {0}, gr[16] = {0}, g0 = 0.0f, g17 = 0.0f, gc = 0.0f,
          gg = 0.0f;
  Float_t ol[16] = {0}, orr[16] = {0};
  Float_t sf[18] = {0}, so[18] = {0};
  cal->SetBranchAddress("GainLeft", gl);
  cal->SetBranchAddress("GainRight", gr);
  cal->SetBranchAddress("GainStrip0", &g0);
  cal->SetBranchAddress("GainStrip17", &g17);
  cal->SetBranchAddress("GainCathode", &gc);
  // Absent on files calibrated before the grid had an anchor.
  if (cal->GetBranch("GainGrid"))
    cal->SetBranchAddress("GainGrid", &gg);
  // Per-end offsets are absent on files calibrated before they were
  // measured; absent -> 0 (no offset).
  Bool_t has_off = cal->GetBranch("OffsetLeft") != nullptr;
  if (has_off) {
    cal->SetBranchAddress("OffsetLeft", ol);
    cal->SetBranchAddress("OffsetRight", orr);
  }
  // StripFactor is optional (added after the initial gain write, filled in by
  // FindStripCentroidAlignment on the second); absent -> default 1.0 (identity)
  Bool_t has_factor = cal->GetBranch("StripFactor") != nullptr;
  if (has_factor)
    cal->SetBranchAddress("StripFactor", sf);
  Bool_t has_offset = cal->GetBranch("StripOffset") != nullptr;
  if (has_offset)
    cal->SetBranchAddress("StripOffset", so);
  cal->GetEntry(0);
  for (Int_t k = 0; k < 16; k++) {
    gain_left[k] = gl[k];
    gain_right[k] = gr[k];
    offset_left[k] = has_off ? ol[k] : 0.0f;
    offset_right[k] = has_off ? orr[k] : 0.0f;
  }
  if (has_factor)
    for (Int_t s = 0; s < 18; s++)
      strip_factor[s] = sf[s];
  if (has_offset)
    for (Int_t s = 0; s < 18; s++)
      strip_offset[s] = so[s];
  gain_strip0 = g0;
  gain_strip17 = g17;
  gain_cathode = gc;
  gain_grid = gg;
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
    // An end that fired carries its offset; one that did not reads 0 and
    // stays 0. Clamped at 0 so a reading below the offset is not negative.
    for (Int_t k = 0; k < 16; k++) {
      left[k] = leftdE_adc[k] > 0
                    ? Double_t(gain_left[k]) *
                          TMath::Max(0.0, Double_t(leftdE_adc[k]) -
                                              Double_t(offset_left[k]))
                    : 0.0;
      right[k] = rightdE_adc[k] > 0
                     ? Double_t(gain_right[k]) *
                           TMath::Max(0.0, Double_t(rightdE_adc[k]) -
                                               Double_t(offset_right[k]))
                     : 0.0;
    }
    strip0 = Double_t(gain_strip0) * Double_t(strip0_adc);
    strip17 = Double_t(gain_strip17) * Double_t(strip17_adc);
    // The -1 "no cathode" sentinel: uncalibrated/absent -> 0 a.u.
    cathode = (cathode_adc > 0) ? Double_t(gain_cathode) * Double_t(cathode_adc)
                                : 0.0;
    // The grid over its modal beam peak (1.0 for a beam event); files
    // calibrated before it had an anchor fall back to a full-scale fraction.
    grid = gain_grid > 0.0f ? Double_t(gain_grid) * Double_t(grid_adc)
                            : Double_t(grid_adc) / 16384.0;
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
  if (Constants::ActiveIgnoreShortStrips()) {
    for (Int_t s = 1; s <= 16; s++) {
      if ((s % 2) != 0)
        right[s - 1] = 0.0;
      else
        left[s - 1] = 0.0;
    }
  }
  // The factor on every end, the offset on the fired long end only:
  // unfired ends stay 0; on CoMPASS offset on nothing outnumbered the beam.
  if (is_normed) {
    for (Int_t s = 1; s <= 16; s++) {
      left[s - 1] *= Double_t(strip_factor[s]);
      right[s - 1] *= Double_t(strip_factor[s]);
      Double_t &long_end = ((s % 2) != 0) ? left[s - 1] : right[s - 1];
      if (long_end > 0.0)
        long_end += Double_t(strip_offset[s]);
    }
    if (strip0 > 0.0)
      strip0 = strip0 * Double_t(strip_factor[0]) + Double_t(strip_offset[0]);
    if (strip17 > 0.0)
      strip17 =
          strip17 * Double_t(strip_factor[17]) + Double_t(strip_offset[17]);
  }
}
