#include "Paths.hpp"

#ifndef MUSIC_DATASET_NAME
#define MUSIC_DATASET_NAME "unknown"
#endif
#ifndef MUSIC_GIT_HASH
#define MUSIC_GIT_HASH "unknown"
#endif
// Absolute dataset dir, baked in at build by the Makefile (-DMUSIC_DATASET_DIR)
// from the DATASET the binary was built for. The binary is self-locating; there
// is no runtime env-var dependency.
#ifndef MUSIC_DATASET_DIR
#define MUSIC_DATASET_DIR ""
#endif

TString Paths::ProjectRootOf(const char *file) {
  TString path = file;
  path.ReplaceAll("/./", "/");
  TString src_dir = gSystem->DirName(path);
  return gSystem->DirName(src_dir);
}

TString Paths::DatasetName() { return TString(MUSIC_DATASET_NAME); }

void Paths::PrintBanner(const TString &dataset_dir) {
  std::cout << "============================================================"
            << std::endl;
  std::cout << " MUSIC tooling | dataset=" << MUSIC_DATASET_NAME
            << " | git=" << MUSIC_GIT_HASH << std::endl;
  std::cout << " dataset dir : " << dataset_dir << std::endl;
  std::cout << " CoMPASS base: " << Constants::COMPASS_BASE_DIR << std::endl;
  std::cout << " runs        :";
  for (Int_t i = 0; i < Int_t(Constants::RUN_NUMBERS.size()); i++)
    std::cout << " " << Constants::RUN_NUMBERS[i];
  std::cout << std::endl;
  std::cout << " event mode  : "
            << (Constants::USE_TIME_WINDOW_EVENTS
                    ? Form("time-window (%.1f us)",
                           Constants::EVENT_TIME_WINDOW_US)
                    : "grid-triggered")
            << " | calibration "
            << (Constants::SKIP_CALIBRATION ? "SKIPPED" : "on") << " | plots "
            << (Constants::SAVE_PLOTS ? "on" : "SKIPPED") << std::endl;
  std::cout << "============================================================"
            << std::endl;
}

TString Paths::DatasetDir() {
  TString d = MUSIC_DATASET_DIR;
  if (d.Length() == 0) {
    std::cerr << "FATAL: this binary was built without a dataset dir baked in "
                 "(MUSIC_DATASET_DIR). Rebuild via the Makefile."
              << std::endl;
    gSystem->Exit(1);
  }
  static Bool_t printed = kFALSE;
  if (!printed) {
    PrintBanner(d);
    printed = kTRUE;
  }
  return d;
}
