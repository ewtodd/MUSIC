#ifndef NORMALIZATION_HPP
#define NORMALIZATION_HPP

#include "Constants.hpp"
#include <Rtypes.h>
#include <TTree.h>

// Reads the raw-ADC events tree and calibrates on the fly: MeV = gain x ADC.
// Per-strip gains come from the one-row "calibration" tree that lives in the
// same file as the events tree (written by CalibrateBeam). No per-event
// calibrated values are stored on disk -- the raw tree plus the gain row fully
// determine them.
//
// GainLeft[s] multiplies Left_0_17_dE[s] and GainRight[s] multiplies RightdE[s]
// for all 18 indices, so total[s] = left[s] + right[s] holds uniformly: at the
// guard indices 0/17 the energy sits in the left slot and RightdE is 0.
//
// For a TChain spanning subfiles -- each with its own gains -- Decode() notices
// when the chain crosses into a new file (the tree number changes) and reloads
// that file's gains, so per-subfile calibration is exact (no per-run gain
// approximation). This is a non-owning view: it mutates nothing on the tree and
// holds no resource, so it is safe to copy and to destroy in any order relative
// to the tree/file it reads.
struct EnergyView {
  // raw ADC, bound to the events tree/chain. cathode_adc keeps the -1 "no
  // cathode hit" sentinel, so it is Short_t (matching the Cathode branch).
  UShort_t left_0_17_adc[18], rightdE_adc[18];
  Short_t cathode_adc;

  // per-strip gains (MeV/ADC) from the active file's calibration tree; 0 where
  // uncalibrated.
  Float_t gain_left[18], gain_right[18], gain_cathode;
  Bool_t is_mev; // a calibration tree was found -> Decode() yields MeV

  // decoded per-event values (MeV if is_mev, else raw ADC promoted to Double_t)
  Double_t left[18], right[18], total[18];
  Double_t cathode;

  TTree *tree_;       // events tree/chain we are bound to
  Int_t loaded_tree_; // tree number whose gains are currently loaded (-1 none)

  EnergyView()
      : cathode_adc(0), gain_cathode(0.0f), is_mev(kFALSE), cathode(0.0),
        tree_(nullptr), loaded_tree_(-1) {}

  // Bind raw branches and load the gains for the current file. Returns is_mev.
  Bool_t Attach(TTree *t);
  // Fill left/right/total/cathode for the current entry, reloading gains first
  // if the chain has advanced into a new subfile.
  void Decode();
  // (Re)load per-strip gains from the active file's calibration tree.
  void LoadGains();
  const char *Unit() const;
};

#endif
