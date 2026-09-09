#ifndef PATHS_HPP
#define PATHS_HPP

#include "Constants.hpp"
#include <Rtypes.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>

/**
 * @brief Where the active dataset lives, and where its output goes.
 *
 * The two are deliberately different kinds of value. The dataset directory is
 * baked in at build time, because a binary belongs to exactly one dataset. The
 * results directory is read at run time, because where processed output lands
 * is a per-machine deployment choice that should be redirectable without a
 * rebuild.
 */
class Paths {
public:
  /**
   * @brief Absolute path to the active dataset directory, `analysis/<iso>`.
   *
   * Read from `MUSIC_DATASET_DIR`, baked in at build time by the Makefile.
   *
   * @return The dataset directory, without a trailing slash.
   *
   * @warning Fatally exits the process if the variable is unset. A binary with
   *          no dataset has nothing to operate on, and continuing would write
   *          output to an arbitrary location.
   *
   * @note Prints the tooling banner on first call, once per process.
   */
  static TString DatasetDir();

  /// @brief The dataset's isotope name, e.g. `"37Cl"`.
  /// @return The name from the build-time `MUSIC_DATASET_NAME`.
  static TString DatasetName();

  /**
   * @brief Absolute path to the directory receiving generated output.
   *
   * Covers `root_files` and `plots`. Read at run time from
   * `MUSIC_RESULTS_DIR`, falling back to DatasetDir() when unset — so the
   * default writes output in-repo.
   *
   * @return The results directory, without a trailing slash.
   *
   * @note Set `MUSIC_RESULTS_DIR` to send output to a scratch drive without
   *       rebuilding. Running from the dev shell points it back at the in-repo
   *       dataset directory, which is writable.
   */
  static TString ResultsDir();

  /**
   * @brief Print the project wordmark.
   *
   * Public so it can run ahead of the lazy banner that DatasetDir() prints, so
   * that a log opens with the project name rather than with a GPU probe.
   */
  static void PrintLogo();

private:
  static void PrintBanner(const TString &dataset_dir);
};

#endif
