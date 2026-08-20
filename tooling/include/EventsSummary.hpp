#ifndef EVENTS_SUMMARY_HPP
#define EVENTS_SUMMARY_HPP

#include "Constants.hpp"
#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "Normalization.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TTree.h>
#include <mutex>
#include <vector>

// Shared summary histogram set used by both EventBuilder (raw ADC) and
// EventsSummary (calibrated a.u.); the helpers eliminate duplicated build/plot
// logic.
struct SummaryHistograms {
  TH2F *h_music;
  TH1F *h_mult;
  TH2F *h2_R_vs_L[18];
  TH1F *h1_cathode;
  TH1F *h1_strip17;
  TH2F *h2_strip0_vs_grid;
  TH1F *h1_strip0;
  TH1F *h1_grid;

  SummaryHistograms()
      : h_music(nullptr), h_mult(nullptr), h1_cathode(nullptr),
        h1_strip17(nullptr), h2_strip0_vs_grid(nullptr), h1_strip0(nullptr),
        h1_grid(nullptr) {
    for (Int_t s = 0; s < 18; s++)
      h2_R_vs_L[s] = nullptr;
  }
};

// Bin range and label configuration for creating summary histograms.
struct SummaryHistConfig {
  TString unit_label; // e.g. "ADC" or "a.u."
  Double_t strip_e_min;
  Double_t strip_e_max;
  Bool_t odd_even_split; // if true, per-strip max from config
  Double_t left_odd_max;
  Double_t left_even_max;
  Double_t right_odd_max;
  Double_t right_even_max;
  Double_t cathode_max;
  Double_t strip17_max;
  Double_t grid_max;
  Double_t strip0_max;
  Int_t music_energy_bins; // bins for Y axis of h_music
};

void CreateSummaryHistograms(SummaryHistograms &h,
                             const SummaryHistConfig &cfg);

void SaveAndDeleteSummaryHistograms(SummaryHistograms &h, TFile *out_file,
                                    const TString &subdir,
                                    const TString &plot_suffix);

class EventsSummary {
public:
  // Same trace from a raw per-strip total[18] array (e.g. a cached event).
  // Caller owns the returned graph.
  static TGraph *BuildTraceFromTotals(const Double_t *total);

  // Draw and save a set of sample traces as an overlay plot (single color on a
  // TH2F frame); written to plots if SAVE_PLOTS. The traces vector is NOT
  // deleted.
  static void SaveSampleTraces(const std::vector<TGraph *> &traces,
                               const TString &save_name, const TString &subdir,
                               Double_t y_min, Double_t y_max,
                               const char *y_title);

  static void BuildNormedSummaryHistograms(const TString &input_filename,
                                           const TString &file_label);
};

#endif
