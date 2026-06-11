#ifndef TRACE_CREATOR_HPP
#define TRACE_CREATOR_HPP

#include "Constants.hpp"
#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "Normalization.hpp"
#include "PlottingUtils.hpp"
#include <Rtypes.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2.h>
#include <TH2F.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TTree.h>
#include <iostream>
#include <mutex>
#include <vector>

class TraceCreator {
public:
  static void WriteH2CanvasToFile(TFile *file, TH2F *h,
                                  const TString &save_name,
                                  const TString &subdir);
  static TGraph *BuildEventTrace(const EnergyView &ev);
  // Same trace from a raw per-strip total[18] array (e.g. a cached event).
  // Caller owns the returned graph.
  static TGraph *BuildTraceFromTotals(const Double_t *total);
  static void BuildNormedSummaryHistograms(const TString &input_filename,
                                           const TString &file_label,
                                           const FileSpec &spec);
  static void BuildTraces(std::vector<TString> input_filenames,
                          std::vector<TString> file_labels,
                          Bool_t save_plots = kTRUE, Bool_t reprocess = kFALSE);
};

#endif
