#ifndef GATE_CACHE_HPP
#define GATE_CACHE_HPP

#include "BeamFit2D.hpp"
#include "Constants.hpp"
#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <Rtypes.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TEllipse.h>
#include <TFile.h>
#include <TH2F.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVectorD.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <mutex>
#include <vector>

// Per-run cache of gate-passing events. The bigaus ellipse gate (Strip1-L vs
// Strip2-R) is fit once per run; passing events are copied into
// root_files/Gated_Run{N}.root (tree "events_gated") alongside the gate
// histogram and fit params, so downstream tools (DeltaEScatter, StripScatter)
// reuse the selection instead of re-fitting and re-scanning every time.
class GateCache {
public:
  static Double_t NSigma();
  static TString FileSubpath(Int_t run);
  static Bool_t Exists(Int_t run);

  // Build the cache from the events_cal source chain if it is missing (or if
  // SKIP_EXISTING is off). No-op when an up-to-date file is already present.
  static void EnsureForRun(Int_t run, TChain *src_chain);

  // Read the stored bigaus gate params from an open cache file.
  static Bool_t LoadGate(TFile *f, BeamFit2D &gate);

  // Standalone entry: build caches for every configured run.
  static void Run();

private:
  static TString FileName(Int_t run);
  static void SaveGatePlot(TH2F *gate_hist, const BeamFit2D &gate,
                           const TString &subdir);
  static void BuildForRun(Int_t run, TChain *src_chain);
};

#endif
