#ifndef PATHS_HPP
#define PATHS_HPP

#include "Constants.hpp"
#include <Rtypes.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>

class Paths {
public:
  // Absolute path to the active dataset directory (analysis/<iso>), from the
  // build-time MUSIC_DATASET_DIR (fatal-exits if unset); banner prints once.
  static TString DatasetDir();

  // Dataset isotope name (e.g. "37Cl"), from the build-time MUSIC_DATASET_NAME.
  static TString DatasetName();

  // Directory that holds GENERATED outputs (root_files, plots). Read at RUNTIME
  // from MUSIC_RESULTS_DIR, falling back to DatasetDir() (in-repo, as before);
  // unlike DatasetDir() it is runtime because where output lands is a
  // per-machine deployment choice.
  static TString ResultsDir();

private:
  static void PrintBanner(const TString &dataset_dir);
};

#endif
