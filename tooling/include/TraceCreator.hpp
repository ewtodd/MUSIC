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
  static void BuildNormedSummaryHistograms(const TString &input_filename,
                                           const TString &file_label);
};

#endif
