#ifndef SIM_TRACE_CREATOR_HPP
#define SIM_TRACE_CREATOR_HPP

#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Paths.hpp"
#include "PlottingUtils.hpp"
#include "RemixConstants.hpp"
#include "RemixSim.hpp"
#include <Rtypes.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <mutex>
#include <vector>

// Sim-side per-event trace builder. Distinct from the experimental tooling
// TraceCreator (which uses EnergyView +
// cal sidecars + FileSpec runs); this reads the calibrated sim "events_MeV"
// tree (Float_t TotaldE/LeftdE/RightdE/Cathode branches), writes per-strip H2
// summaries back into the sim file (UPDATE), and saves per-event trace plots.
class SimTraceCreator {
public:
  static void Run();

private:
  static void BuildTraces(std::vector<TString> input_output_filepaths,
                          std::vector<TString> file_labels,
                          Bool_t save_plots = kTRUE, Bool_t reprocess = kFALSE);
};

#endif
