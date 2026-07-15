#ifndef TRACE_CREATOR_HPP
#define TRACE_CREATOR_HPP

#include "Constants.hpp"
#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "Normalization.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TString.h>
#include <TTree.h>
#include <mutex>
#include <vector>

class TraceCreator {
public:
  // Same trace from a raw per-strip total[18] array (e.g. a cached event).
  // Caller owns the returned graph.
  static TGraph *BuildTraceFromTotals(const Double_t *total);

  // Draw and save a set of sample traces as an overlay plot. Traces are drawn
  // in a single color on a TH2F frame; the file is written and the canvas is
  // saved to plots if SAVE_PLOTS is set. The traces vector is NOT deleted.
  static void SaveSampleTraces(const std::vector<TGraph *> &traces,
                               const TString &save_name, const TString &subdir,
                               Double_t y_min, Double_t y_max,
                               const char *y_title);

  static void BuildNormedSummaryHistograms(const TString &input_filename,
                                           const TString &file_label);
};

#endif
