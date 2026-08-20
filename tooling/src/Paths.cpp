#include "Paths.hpp"

#ifndef MUSIC_DATASET_NAME
#define MUSIC_DATASET_NAME "unknown"
#endif
#ifndef MUSIC_GIT_HASH
#define MUSIC_GIT_HASH "unknown"
#endif
// Absolute dataset dir, baked in at build time (-DMUSIC_DATASET_DIR); the
// binary is self-locating with no runtime env-var dependency.
#ifndef MUSIC_DATASET_DIR
#define MUSIC_DATASET_DIR ""
#endif

TString Paths::DatasetName() { return TString(MUSIC_DATASET_NAME); }

TString Paths::ResultsDir() {
  const Char_t *env = gSystem->Getenv("MUSIC_RESULTS_DIR");
  if (env && env[0] != '\0')
    return TString(env);
  return DatasetDir();
}

void Paths::PrintBanner(const TString &dataset_dir) {
  std::cout << "============================================================"
            << std::endl;
  std::cout << " MUSIC tooling | dataset=" << MUSIC_DATASET_NAME
            << " | git=" << MUSIC_GIT_HASH << std::endl;
  std::cout << " dataset dir : " << dataset_dir << std::endl;
  std::cout << " CoMPASS base: " << Constants::cfg.COMPASS_BASE_DIR
            << std::endl;
  std::cout << " runs        :";
  for (Int_t i = 0; i < Int_t(Constants::cfg.RUN_NUMBERS.size()); i++)
    std::cout << " " << Constants::cfg.RUN_NUMBERS[i];
  std::cout << std::endl;
  std::cout << " event mode  : " << Constants::cfg.REFERENCE_CHANNEL << " ("
            << Constants::cfg.EVENT_TIME_WINDOW_US << " us window)"
            << " | calibration "
            << (Constants::cfg.SKIP_CALIBRATION ? "SKIPPED" : "on")
            << " | plots " << (Constants::cfg.SAVE_PLOTS ? "on" : "SKIPPED")
            << std::endl;
  std::cout << "============================================================"
            << std::endl;
}

TString Paths::DatasetDir() {
  const Char_t *env = gSystem->Getenv("MUSIC_DATASET_DIR");
  TString d;
  if (env && env[0] != '\0') {
    d = TString(env);
    // Reject a MUSIC_DATASET_DIR from a different dataset: it would silently
    // mix one dataset's config/data into another's output tree.
    TString name = TString(MUSIC_DATASET_NAME);
    if (name.Length() > 0 && name != "unknown" && !d.Contains("/" + name)) {
      std::cerr << "FATAL: this binary was built for dataset '" << name
                << "' but MUSIC_DATASET_DIR points to '" << d << "'."
                << std::endl;
      std::cerr << "You are probably running the wrong ./result symlink or "
                   "are in the wrong dataset dev shell."
                << std::endl;
      std::cerr << "Rebuild with `nix build .#" << name
                << "` or enter the matching shell (`nix develop .#" << name
                << "`)." << std::endl;
      gSystem->Exit(1);
    }
  } else {
    d = MUSIC_DATASET_DIR;
    if (d.Length() == 0) {
      std::cerr
          << "FATAL: this binary was built without a dataset dir baked in "
             "(MUSIC_DATASET_DIR). Rebuild via the Makefile."
          << std::endl;
      gSystem->Exit(1);
    }
  }
  static Bool_t printed = kFALSE;
  if (!printed) {
    PrintBanner(d);
    printed = kTRUE;
  }
  return d;
}
