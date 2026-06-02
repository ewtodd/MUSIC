#ifndef SIM_DELTA_E_SCATTER_HPP
#define SIM_DELTA_E_SCATTER_HPP

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
#include <TLegend.h>
#include <TString.h>
#include <TTree.h>
#include <iostream>
#include <mutex>
#include <vector>

// Sim-side overlay scatter plots. Distinct from the experimental tooling
// DeltaEScatter (which does per-run
// bigaus-gated cal-sidecar scatters); this overlays multiple sim datasets
// (truth vs eres, classified by tag) with no gating.
class SimDeltaEScatter {
public:
  static void Run();

private:
  typedef TGraph *(*ScatterBuilderFn)(TTree *, Int_t, Long64_t);

  static TGraph *BuildScatter(TTree *tree, Int_t color, Long64_t max_points);
  static TGraph *BuildTotalSumScatter(TTree *tree, Int_t color,
                                      Long64_t max_points);
  static void MakeOverlayScatter(
      const std::vector<TString> &filepaths, const std::vector<TString> &labels,
      const TString &output_name, const TString &output_subdir,
      ScatterBuilderFn builder, const TString &frame_title, Double_t xmin,
      Double_t xmax, Double_t ymin, Double_t ymax, Long64_t max_points = 25000);
  static std::vector<TGraph *> BuildTracesForDataset(TTree *tree, Int_t color,
                                                     Int_t n_traces);
  static void MakeTraceOverlay(const std::vector<TString> &filepaths,
                               const std::vector<TString> &labels,
                               const TString &output_name,
                               const TString &output_subdir,
                               Int_t n_per_dataset = 1000);
  static TString PrettyLabel(const TString &tag);
};

#endif
