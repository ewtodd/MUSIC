#ifndef STRIP_SUM_SCATTER_HPP
#define STRIP_SUM_SCATTER_HPP

#include "BeamFit2D.hpp"
#include "Constants.hpp"
#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Normalization.hpp"
#include "PlottingUtils.hpp"
#include "TraceCreator.hpp"
#include <Rtypes.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TCutG.h>
#include <TEllipse.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

// Aggregate 2D scatter (all runs combined) of (sum of TotaldE over the full
// split-anode range) vs (sum over a downstream reaction window), after a
// per-run beam gate + (a,n) selection. Strip windows are named constants in the
// .cpp. Reads the calibrated events_cal sidecars (MeV).
class StripSumScatter {
public:
  static void Run();

private:
  // Fit this run's beam gate and add its (a,n)-selected events into the shared
  // aggregate histogram h; gate_subdir receives the per-run gate plot. Returns
  // the fitted gate (ok=false on failure) so the interactive trace-sampling
  // pass can reuse it.
  static BeamFit2D FillRun(Int_t run, TChain *chain, const TString &gate_subdir,
                           TH2F *h);
};

#endif
