#ifndef FILE_SET_HPP
#define FILE_SET_HPP

#include "Constants.hpp"
#include "IOUtils.hpp"
#include <Rtypes.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <mutex>
#include <vector>

/**
 * @brief Serialises all plotting and canvas work.
 *
 * ROOT's global graphics state is not thread-safe and the pipeline runs file
 * workers concurrently, so every canvas operation must hold this.
 */
extern std::mutex g_plot_mutex;

/**
 * @brief Serialises multi-line progress logging from worker threads.
 *
 * A streaming insertion chain is not atomic. Without this, per-run lines from
 * concurrent workers interleave mid-line and the log stops being readable.
 */
extern std::mutex g_log_mutex;

/**
 * @brief One input file: a run number and the subfile suffix within it.
 *
 * CoMPASS splits a long run across numbered subfiles, so a run alone does not
 * identify a file.
 */
struct FileSpec {
  Int_t run;      ///< Run number.
  TString suffix; ///< Subfile suffix, empty for the first subfile.
};

/**
 * @brief Discovering input files and naming the artefacts derived from them.
 *
 * All-static. Every name here is derived rather than stored, so a stage can
 * reconstruct the path to any other stage's output from a FileSpec alone.
 */
class FileSet {
public:
  /// @brief Path to a CoMPASS binary subfile.
  /// @param s Run and subfile.
  static TString CompassBinPath(const FileSpec &s);

  /// @brief Subfile suffixes present on disk for a CoMPASS run.
  /// @param run Run number.
  /// @return The suffixes found, empty if the run has no files.
  static std::vector<TString> DiscoverRunSuffixes(Int_t run);

  /// @brief Path to a SOLARIS `.sol` subfile.
  /// @param s Run and subfile.
  static TString SolBinPath(const FileSpec &s);

  /// @brief Subfile suffixes present on disk for a SOLARIS run.
  /// @param run Run number.
  static std::vector<TString> DiscoverSolRunSuffixes(Int_t run);

  /// @brief Subfile suffixes for which processed output already exists.
  /// @param run Run number.
  static std::vector<TString> DiscoverProcessedRunSuffixes(Int_t run);

  /**
   * @brief Every raw input subfile for the configured runs.
   * @return One FileSpec per subfile, over the runs selected by `RUNS` and
   *         `N_CHUNKS` in the dataset configuration.
   */
  static std::vector<FileSpec> BuildFileSpecs();

  /// @brief Every subfile that already has processed output.
  static std::vector<FileSpec> BuildProcessedFileSpecs();

  /**
   * @brief The union of raw and processed subfiles, without duplicates.
   *
   * Lets a stage run over everything available whether or not a given subfile
   * has been preprocessed yet.
   *
   * @return Raw specs first, then any processed spec not already among them.
   */
  static std::vector<FileSpec> BuildRawOrProcessedFileSpecs();

  /// @brief Filename of the raw converted ROOT file for a subfile.
  static TString RawRootName(const FileSpec &s);

  /// @brief Filename of the timing-shift friend tree for a subfile.
  static TString ShiftFriendName(const FileSpec &s);

  /// @brief Filename of the built-events ROOT file for a subfile.
  static TString EventsName(const FileSpec &s);

  /**
   * @brief Human-readable label identifying a subfile.
   *
   * `run<N><suffix>`, prefixed with the active file tag when the dataset sets
   * one. Used in plot names, log lines, and as the key ResolveFileSpec()
   * reverses.
   */
  static TString FileLabel(const FileSpec &s);

  /**
   * @brief Recover the FileSpec behind a label from FileLabel().
   *
   * @param file_label Label to resolve.
   * @return The matching spec, or one with `run == -1` and an empty suffix if
   *         no configured subfile produces that label. Check `run` before use.
   */
  static FileSpec ResolveFileSpec(const TString &file_label);

  /**
   * @brief Chain every run's events files, grouped by run.
   *
   * @param[out] run_order Run numbers in the order they should be processed.
   *
   * @return One `TChain` per run, keyed by run number. **The caller owns every
   *         chain** and must delete them.
   */
  static std::map<Int_t, TChain *>
  GroupEventsByRun(std::vector<Int_t> &run_order);

  /**
   * @brief Stride that visits at most @p max_points of @p n_total entries.
   *
   * For scatter plots, where drawing every entry costs time and renders as a
   * solid block anyway.
   *
   * @param n_total    Entries available.
   * @param max_points Cap on entries to visit; `0` or less means no cap.
   *
   * @return The stride, never below 1.
   */
  static Long64_t SampleStride(Long64_t n_total, Long64_t max_points);
};

#endif
