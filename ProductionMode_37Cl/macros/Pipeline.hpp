#ifndef PIPELINE_HPP
#define PIPELINE_HPP

#include "Constants.hpp"
#include <TString.h>
#include <mutex>
#include <vector>

inline std::mutex g_plot_mutex;

// Identifies one input file by run number and suffix.
// suffix is "_1", "_2", ... for indexed files, or "" for the unindexed
// file that COMPASS also writes (DataR_run_<run>.BIN).
struct FileSpec {
  Int_t run;
  TString suffix;
};

// Returns one FileSpec per (run in RUN_NUMBERS) x (suffix in
// {_0, _1, ..., _<N_FILES-1>, ""}).
inline std::vector<FileSpec> BuildFileSpecs() {
  std::vector<FileSpec> specs;
  for (Int_t r = 0; r < Int_t(Constants::RUN_NUMBERS.size()); r++) {
    Int_t run = Constants::RUN_NUMBERS[r];
    for (Int_t i = 1; i < Constants::N_FILES; i++) {
      FileSpec s;
      s.run = run;
      s.suffix = Form("_%d", i);
      specs.push_back(s);
    }
    FileSpec s;
    s.run = run;
    s.suffix = "";
    specs.push_back(s);
  }
  return specs;
}

// File-name helpers. The lowercase "run" / capitalized "Run" split
// matches the existing on-disk convention.
inline TString RawRootName(const FileSpec &s) {
  return Form("DataR_run_%d%s.root", s.run, s.suffix.Data());
}
inline TString SortedName(const FileSpec &s) {
  return Form("Sorted_Run%d%s", s.run, s.suffix.Data());
}
inline TString EventsClassicName(const FileSpec &s) {
  return Form("EventsClassic_Run%d%s", s.run, s.suffix.Data());
}
inline TString EventsNearestName(const FileSpec &s) {
  return Form("EventsNearest_Run%d%s", s.run, s.suffix.Data());
}
inline TString EventsName(const FileSpec &s) {
  if (Constants::EVENT_BUILDER_MODE == Constants::EVENT_BUILDER_NEAREST_GRID)
    return EventsNearestName(s);
  return EventsClassicName(s);
}
inline TString FileLabel(const FileSpec &s) {
  return Form("run%d%s", s.run, s.suffix.Data());
}

// Resolves a user-supplied label ("16", "37", "run16", "Run_16",
// "EventsNearest_Run16.root", ...) to a FileSpec. Returns spec with
// run = -1 if no integer can be parsed out.
inline FileSpec ResolveFileSpec(const TString &file_label) {
  for (Int_t k = 0; k < Int_t(BuildFileSpecs().size()); k++) {
    FileSpec s = BuildFileSpecs()[k];
    if (FileLabel(s) == file_label)
      return s;
  }
  FileSpec s;
  s.run = -1;
  s.suffix = "";
  Int_t i = 0;
  while (i < file_label.Length() &&
         !(file_label[i] >= '0' && file_label[i] <= '9'))
    i++;
  Int_t j = i;
  while (j < file_label.Length() &&
         (file_label[j] >= '0' && file_label[j] <= '9'))
    j++;
  if (j > i)
    s.run = TString(file_label(i, j - i)).Atoi();
  return s;
}

#endif
