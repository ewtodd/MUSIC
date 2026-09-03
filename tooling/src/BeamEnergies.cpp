#include "BeamEnergies.hpp"
#include "Constants.hpp"
#include "Paths.hpp"
#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>

namespace BeamEnergies {

TString SimPath() {
  return Paths::DatasetDir() + "/sim_root_files/" +
         Constants::cfg.CROSS_SECTION_CONFIG.BEAM_SIM_FILE;
}

Bool_t Profile(const TString &path, Double_t *dE, Double_t &e_strip0) {
  TFile f(path, "READ");
  if (f.IsZombie())
    return kFALSE;
  TTree *t = static_cast<TTree *>(f.Get("events_MeV"));
  TTree *mc = static_cast<TTree *>(f.Get("MC"));
  if (!t || !mc)
    return kFALSE;
  Float_t left[18], right[18];
  t->SetBranchAddress("Left_0_17_dE", left);
  t->SetBranchAddress("RightdE", right);
  for (Int_t s = 0; s < 18; s++)
    dE[s] = 0.0;
  const Long64_t n = t->GetEntries();
  if (n == 0)
    return kFALSE;
  for (Long64_t i = 0; i < n; i++) {
    t->GetEntry(i);
    for (Int_t s = 0; s < 18; s++)
      dE[s] += left[s] + right[s];
  }
  for (Int_t s = 0; s < 18; s++)
    dE[s] /= Double_t(n);

  TLeaf *gas = mc->GetLeaf("beam_energy_gas");
  TLeaf *dead = mc->GetLeaf("DeadUS_dE");
  if (!gas || !dead)
    return kFALSE;
  const Long64_t nm = mc->GetEntries();
  if (nm == 0)
    return kFALSE;
  Double_t sum = 0.0;
  for (Long64_t i = 0; i < nm; i++) {
    mc->GetEntry(i);
    sum += gas->GetValue() - dead->GetValue();
  }
  e_strip0 = sum / Double_t(nm);
  return kTRUE;
}

Double_t LabAtStrip(const Double_t *dE, Double_t e_strip0, Int_t reac) {
  Double_t e = e_strip0;
  for (Int_t s = 0; s < reac; s++)
    e -= dE[s];
  return e;
}

Double_t CmFraction() {
  const CrossSectionConfig &X = Constants::cfg.CROSS_SECTION_CONFIG;
  const Int_t target_a = TargetGasA(X.TARGET_GAS);
  return Double_t(target_a) / Double_t(X.BEAM_A + target_a);
}

} // namespace BeamEnergies
