#ifndef NORMALIZATION_HPP
#define NORMALIZATION_HPP

#include "Constants.hpp"
#include <TROOT.h>
#include <TString.h>
#include <TTree.h>
#include <cstdio>

struct EnergyView {
  Bool_t is_mev;

  // ADC source (from the event builder)
  Int_t leftdE_adc[18], rightdE_adc[18], totaldE_adc[18];
  Int_t cathode_adc;

  // MeV source (from BeamCalibration)
  Float_t leftdE_mev[18], rightdE_mev[18], totaldE_mev[18];
  Float_t cathode_mev;

  // Per-event Double_t view, populated by Decode().
  Double_t left[18], right[18], total[18];
  Double_t cathode;

  Bool_t Attach(TTree *t) {
    is_mev = (t->GetBranch("TotaldEMeV") != nullptr);
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

  void Decode() {
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
        if (Constants::LongAnodeSide(s) == 'L') {
          total[s] = left[s];
          right[s] = 0.0;
        } else {
          total[s] = right[s];
          left[s] = 0.0;
        }
      }
    }
  }

  const char *Unit() const { return is_mev ? "MeV" : "ADC"; }
};

struct CutXY {
  Double_t x;
  Double_t y;
};

inline CutXY ComputeCutXY(const TString &hist_name,
                          const Double_t normalized[18], Double_t event_total) {
  CutXY out = {0, 0};
  if (hist_name.BeginsWith("h2_totalE_vs_stripE_s")) {
    TString num = hist_name(TString("h2_totalE_vs_stripE_s").Length(),
                            hist_name.Length());
    out.x = normalized[num.Atoi()];
    out.y = event_total;
  } else if (hist_name.BeginsWith("h2_sum_s")) {
    Int_t sX, sY, sI, sJ;
    if (std::sscanf(hist_name.Data(), "h2_sum_s%d_s%d_vs_s%d_s%d", &sX, &sY,
                    &sI, &sJ) == 4) {
      out.x = normalized[sI] + normalized[sJ];
      out.y = normalized[sX] + normalized[sY];
    }
  }
  return out;
}

#endif
