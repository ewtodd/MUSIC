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

/**
 * @brief The standard diagnostic histogram set for a subfile.
 *
 * Filled at two points with the same shapes but different units: by
 * EventBuilder in raw ADC, and by EventsSummary in calibrated arbitrary units.
 * Sharing the struct and its create/save helpers is what keeps the two views
 * comparable rather than drifting apart.
 *
 * @note Members are raw pointers, owned by whoever called
 *       CreateSummaryHistograms(). SaveAndDeleteSummaryHistograms() both writes
 *       and frees them.
 */
struct SummaryHistograms {
  TH2F *h_music;              ///< Energy against strip index: the MUSIC plot.
  TH1F *h_mult;               ///< Hit multiplicity per event.
  TH2F *h2_long_vs_short[18]; ///< Long against short gate, per strip.
  TH1F *h1_cathode;           ///< Cathode energy.
  TH1F *h1_strip17;           ///< Strip 17 energy.
  TH2F *h2_strip0_vs_grid;    ///< Strip 0 against grid energy.
  TH1F *h1_strip0;            ///< Strip 0 energy.
  TH1F *h1_grid;              ///< Grid energy.

  /// @brief Construct with every pointer null; call CreateSummaryHistograms().
  SummaryHistograms()
      : h_music(nullptr), h_mult(nullptr), h1_cathode(nullptr),
        h1_strip17(nullptr), h2_strip0_vs_grid(nullptr), h1_strip0(nullptr),
        h1_grid(nullptr) {
    for (Int_t s = 0; s < 18; s++)
      h2_long_vs_short[s] = nullptr;
  }
};

/**
 * @brief Axis ranges and labels for building a SummaryHistograms set.
 *
 * Supplied by the caller so the same shapes can be built in ADC or in
 * calibrated units without the histogram code knowing which.
 */
struct SummaryHistConfig {
  TString unit_label;   ///< Axis unit, e.g. `"ADC"` or `"a.u."`.
  Double_t strip_e_min; ///< Lower bound of the per-strip energy axis.
  Double_t strip_e_max; ///< Upper bound of the per-strip energy axis.
  /// When true, odd and even strips take separate ceilings from the four
  /// `*_odd_max` / `*_even_max` fields rather than sharing one.
  Bool_t odd_even_split;
  Double_t left_odd_max;   ///< Left-side ceiling, odd strips.
  Double_t left_even_max;  ///< Left-side ceiling, even strips.
  Double_t right_odd_max;  ///< Right-side ceiling, odd strips.
  Double_t right_even_max; ///< Right-side ceiling, even strips.
  Double_t cathode_max;    ///< Cathode axis maximum.
  Double_t strip17_max;    ///< Strip 17 axis maximum.
  Double_t grid_max;       ///< Grid axis maximum.
  Double_t strip0_max;     ///< Strip 0 axis maximum.
  Int_t music_energy_bins; ///< Bins on the energy axis of
                           ///< SummaryHistograms::h_music.
};

/**
 * @brief Allocate every histogram in the set.
 * @param[out] h   Set to populate; existing pointers are overwritten, not
 * freed.
 * @param cfg      Axis ranges and labels.
 */
void CreateSummaryHistograms(SummaryHistograms &h,
                             const SummaryHistConfig &cfg);

/**
 * @brief Write the set to file, save the plots, then free every histogram.
 * @param[in,out] h   Set to write; every pointer is deleted and left dangling.
 * @param out_file    Destination file. Must be open for writing.
 * @param subdir      Plot subdirectory under the plots base.
 * @param plot_suffix Suffix appended to each plot filename.
 * @warning Deletes the histograms. Do not touch @p h afterwards without
 *          calling CreateSummaryHistograms() again.
 */
void SaveAndDeleteSummaryHistograms(SummaryHistograms &h, TFile *out_file,
                                    const TString &subdir,
                                    const TString &plot_suffix);

/**
 * @brief Calibrated summary plots and sample traces, built from events files.
 *
 * Runs after the pipeline, over its output, using EnergyView to decode raw ADC
 * into calibrated units.
 */
class EventsSummary {
public:
  /**
   * @brief Build a trace graph from one event's per-strip totals.
   * @param total 18 per-strip energies, as from EnergyView::total or a cached
   *              event.
   * @return A newly allocated graph of energy against strip index. **The caller
   *         owns it.**
   */
  static TGraph *BuildTraceFromTotals(const Double_t *total);

  /**
   * @brief Draw a set of traces overlaid on one frame and save it.
   *
   * All traces are drawn in a single colour on a `TH2F` frame. The canvas is
   * saved to the plots directory when `SAVE_PLOTS` is set.
   *
   * @param traces    Traces to overlay. **Not deleted** — the caller keeps
   *                  ownership.
   * @param save_name Output basename, without extension.
   * @param subdir    Plot subdirectory under the plots base.
   * @param y_min     Lower bound of the energy axis.
   * @param y_max     Upper bound of the energy axis.
   * @param y_title   Energy axis title, carrying the unit.
   */
  static void SaveSampleTraces(const std::vector<TGraph *> &traces,
                               const TString &save_name, const TString &subdir,
                               Double_t y_min, Double_t y_max,
                               const char *y_title);

  /**
   * @brief Build the calibrated summary histograms for one events file.
   * @param input_filename Events ROOT file to read.
   * @param file_label     Label from FileSet::FileLabel(), used in plot names.
   */
  static void BuildNormedSummaryHistograms(const TString &input_filename,
                                           const TString &file_label);
};

#endif
