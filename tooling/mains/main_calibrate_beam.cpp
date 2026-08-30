#include "CalibrateBeam.hpp"
#include "Constants.hpp"
#include <TString.h>
#include <cstdlib>
#include <iostream>

// Resolve the epoch from the run number in a file label ("run100_8_c000") and
// make it active for the run. Without this the tool would calibrate on the flat
// DatasetConfig while the pipeline uses the epoch's hardware settings, and the
// two would silently disagree.
static Bool_t SetEpochFromLabel(const TString &file_label) {
  if (Constants::cfg.EPOCHS.empty())
    return kTRUE;
  // A tagged label ("compass_run37_1") names its epoch outright. Run numbers
  // repeat across acquisition systems, so the tag is the only thing that
  // disambiguates them; fall back to the run number only for untagged eras.
  for (Int_t e = 0; e < Int_t(Constants::cfg.EPOCHS.size()); e++) {
    const RunEpoch &ep = Constants::cfg.EPOCHS[e];
    if (ep.file_tag.Length() > 0 && file_label.BeginsWith(ep.file_tag + "_")) {
      Constants::SetActiveEpoch(&ep);
      std::cout << "label '" << file_label << "' -> epoch " << ep.name
                << std::endl;
      return kTRUE;
    }
  }
  Int_t start = file_label.Index("run");
  if (start < 0)
    return kFALSE;
  Int_t i = start + 3, run = 0, ndig = 0;
  while (i < file_label.Length() && isdigit(file_label[i])) {
    run = run * 10 + (file_label[i] - '0');
    i++;
    ndig++;
  }
  if (ndig == 0)
    return kFALSE;
  const RunEpoch *ep = Constants::EpochForRun(run);
  if (!ep) {
    std::cerr << "run " << run << " belongs to no declared epoch" << std::endl;
    return kFALSE;
  }
  Constants::SetActiveEpoch(ep);
  std::cout << "run " << run << " -> epoch " << ep->name << std::endl;
  return kTRUE;
}

int main(int argc, char **argv) {
  TString file_label = (argc > 1) ? TString(argv[1]) : TString("");
  if (!Constants::cfg.EPOCHS.empty() && file_label.Length() == 0) {
    std::cerr << "this dataset declares epochs, so a file label is required "
                 "(e.g. run100_8_c000): the epoch cannot be resolved otherwise"
              << std::endl;
    return 1;
  }
  if (!SetEpochFromLabel(file_label)) {
    std::cerr << "cannot resolve an epoch for '" << file_label << "'"
              << std::endl;
    return 1;
  }
  CalibrateBeam::Run(file_label);
  Constants::SetActiveEpoch(nullptr);
  return 0;
}
