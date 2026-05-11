#ifndef NORMALIZATION_HPP
#define NORMALIZATION_HPP

#include <TFile.h>
#include <TROOT.h>
#include <TString.h>
#include <TTree.h>
#include <cstdio>

struct Baseline {
  Float_t scale[18];
};

inline Bool_t LoadBaseline(TFile *file, Baseline &b) {
  for (Int_t s = 0; s < 18; s++)
    b.scale[s] = 1.0f;
  TTree *cal = static_cast<TTree *>(file->Get("baseline"));
  if (!cal || cal->GetEntries() == 0)
    return kFALSE;
  cal->SetBranchAddress("Scale", b.scale);
  cal->GetEntry(0);
  return kTRUE;
}

inline void ComputeNormalized(const Baseline &b, const Int_t leftdE[18],
                              const Int_t rightdE[18],
                              const Int_t totaldE[18],
                              Double_t normalized[18]) {
  for (Int_t s = 0; s < 18; s++) {
    if (s == 0 || s == 17)
      normalized[s] = b.scale[s] * Double_t(totaldE[s]);
    else
      normalized[s] = b.scale[s] * Double_t(leftdE[s] + rightdE[s]);
  }
}

struct CutXY {
  Double_t x;
  Double_t y;
};

inline CutXY ComputeCutXY(const TString &hist_name,
                          const Double_t normalized[18],
                          Double_t event_total) {
  CutXY out = {0, 0};
  if (hist_name.BeginsWith("h2_StripE_vs_TotalE_s")) {
    TString num = hist_name(TString("h2_StripE_vs_TotalE_s").Length(),
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
