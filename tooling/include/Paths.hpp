#ifndef PATHS_HPP
#define PATHS_HPP

#include "Constants.hpp"
#include <Rtypes.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>

class Paths {
public:
  // Absolute path to the active dataset directory (analysis/<iso>), read from
  // the MUSIC_DATASET_DIR baked in at build by the Makefile. Fatal-exits if
  // unset. Prints the tooling banner once, on first call.
  static TString DatasetDir();

  // Dataset isotope name (e.g. "37Cl"), from the build-time MUSIC_DATASET_NAME.
  static TString DatasetName();

  // Absolute path to the directory that holds GENERATED outputs (root_files,
  // plots) for the active dataset. Read at runtime from the MUSIC_RESULTS_DIR
  // env var; falls back to DatasetDir() when unset, so default behaviour writes
  // outputs in-repo exactly as before. Unlike DatasetDir() this is a runtime
  // (not build-time) value on purpose: where processed output lands is a
  // per-machine deployment choice, redirectable without a rebuild.
  static TString ResultsDir();

private:
  static void PrintBanner(const TString &dataset_dir);
};

#endif
