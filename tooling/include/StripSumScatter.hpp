#ifndef STRIP_SUM_SCATTER_HPP
#define STRIP_SUM_SCATTER_HPP

#include "BeamFit2D.hpp"
#include "Constants.hpp"
#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Normalization.hpp"
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
#include <algorithm>
#include <iostream>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

// Per-run 2D scatter of (sum of TotaldE over a downstream strip window) vs
// (sum of TotaldE over the full split-anode range). Strip windows are named
// constants in the .cpp. Reads the calibrated events_cal sidecars (MeV).
class StripSumScatter {
public:
  static void Run();

private:
  static void ProcessRun(Int_t run, TChain *chain, const TString &subdir);
};

#endif
