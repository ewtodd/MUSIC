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

// Per-reaction-strip 2D scatters (all runs combined) of (sum of dE over the
// full split-anode range) vs (sum over the downstream reaction window), after a
// per-run beam gate + (a,n) selection. One scatter is produced for every
// candidate reaction strip [REACTION_STRIP_MIN, REACTION_STRIP_MAX] in a single
// pass and cached (fingerprinted) under root_files, so re-runs and reaction-
// strip changes are near-instant. Energies are calibrated MeV via EnergyView.
class StripSumScatter {
public:
  static void Run();
};

#endif
