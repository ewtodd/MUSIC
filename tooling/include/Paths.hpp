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

  // Legacy: derive a repo-relative root from a source __FILE__. Retained for
  // any non-dataset use; dataset code should use DatasetDir().
  static TString ProjectRootOf(const char *file);

private:
  static void PrintBanner(const TString &dataset_dir);
};

#endif
