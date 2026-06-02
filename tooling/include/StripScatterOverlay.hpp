#ifndef STRIP_SCATTER_OVERLAY_HPP
#define STRIP_SCATTER_OVERLAY_HPP

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

// Sim-side per-reaction-strip overlay scatter plots. Groups eres
// sim datasets by reaction strip (decoded from the "_s<N>" tag token), overlays
// each group's strip/window scatters, and folds in the unreacted beam
// reference.
class StripScatterOverlay {
public:
  static void Run();

private:
  typedef TGraph *(*ScatterBuilderFn)(TTree *, Int_t, Long64_t, Int_t);

  static void PrepareTree(TTree *tree, Float_t *totaldE, Long64_t max_points,
                          Long64_t &n_total, Long64_t &n_use, Long64_t &stride);
  static void StyleScatter(TGraph *g, Int_t color);
  static TGraph *BuildStripScatter(TTree *tree, Int_t color,
                                   Long64_t max_points, Int_t strip);
  static TGraph *BuildReactionToEndScatter(TTree *tree, Int_t color,
                                           Long64_t max_points, Int_t center);
  static void MakeOverlay(
      const std::vector<TString> &filepaths, const std::vector<TString> &labels,
      ScatterBuilderFn builder, Int_t builder_arg, const TString &output_name,
      const TString &output_subdir, const TString &frame_title, Double_t xmin,
      Double_t xmax, Double_t ymin, Double_t ymax, Long64_t max_points = 25000);
  static TString PrettyLabel(const TString &tag);
};

#endif
