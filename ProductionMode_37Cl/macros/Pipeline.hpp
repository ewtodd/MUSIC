#ifndef PIPELINE_HPP
#define PIPELINE_HPP

#include "Constants.hpp"
#include <TString.h>
#include <TSystem.h>
#include <algorithm>
#include <iostream>
#include <mutex>
#include <vector>

inline std::mutex g_plot_mutex;

struct FileSpec {
  Int_t run;
  TString suffix;
};

inline TString CompassBinPath(const FileSpec &s) {
  return Constants::COMPASS_BASE_DIR +
         Form("run_%d/RAW/DataR_run_%d%s.BIN", s.run, s.run, s.suffix.Data());
}

inline std::vector<TString> DiscoverRunSuffixes(Int_t run) {
  std::vector<TString> suffixes;
  TString dir = Constants::COMPASS_BASE_DIR + Form("run_%d/RAW/", run);
  void *dirp = gSystem->OpenDirectory(dir);
  if (!dirp) {
    std::cerr << "DiscoverRunSuffixes: cannot open " << dir << std::endl;
    return suffixes;
  }
  TString prefix = Form("DataR_run_%d", run);
  TString ext = ".BIN";
  const Char_t *name;
  while ((name = gSystem->GetDirEntry(dirp))) {
    TString fname(name);
    if (!fname.BeginsWith(prefix))
      continue;
    if (!fname.EndsWith(ext))
      continue;
    TString rest =
        fname(prefix.Length(), fname.Length() - prefix.Length() - ext.Length());
    if (rest == "") {
      suffixes.push_back("");
      continue;
    }
    if (rest.Length() < 2 || rest[0] != '_')
      continue;
    TString num = rest(1, rest.Length() - 1);
    if (!num.IsDigit())
      continue;
    suffixes.push_back(rest);
  }
  gSystem->FreeDirectory(dirp);
  std::sort(suffixes.begin(), suffixes.end(),
            [](const TString &a, const TString &b) {
              if (a == "")
                return true;
              if (b == "")
                return false;
              return TString(a(1, a.Length() - 1)).Atoi() <
                     TString(b(1, b.Length() - 1)).Atoi();
            });
  return suffixes;
}

inline std::vector<FileSpec> BuildFileSpecs() {
  std::vector<FileSpec> specs;
  for (Int_t r = 0; r < Int_t(Constants::RUN_NUMBERS.size()); r++) {
    Int_t run = Constants::RUN_NUMBERS[r];
    if (Constants::N_FILES < 0) {
      std::vector<TString> suffixes = DiscoverRunSuffixes(run);
      for (Int_t k = 0; k < Int_t(suffixes.size()); k++) {
        FileSpec s;
        s.run = run;
        s.suffix = suffixes[k];
        specs.push_back(s);
      }
    } else {
      FileSpec s0;
      s0.run = run;
      s0.suffix = "";
      specs.push_back(s0);
      for (Int_t i = 1; i < Constants::N_FILES; i++) {
        FileSpec s;
        s.run = run;
        s.suffix = Form("_%d", i);
        specs.push_back(s);
      }
    }
  }
  return specs;
}

inline TString RawRootName(const FileSpec &s) {
  return Form("DataR_run_%d%s.root", s.run, s.suffix.Data());
}
inline TString ShiftFriendName(const FileSpec &s) {
  return Form("DataR_run_%d%s.shift.root", s.run, s.suffix.Data());
}
inline TString EventsName(const FileSpec &s) {
  return Form("Events_Run%d%s", s.run, s.suffix.Data());
}
inline TString TracesSidecarName(const FileSpec &s) {
  return Form("Events_Run%d%s.traces.root", s.run, s.suffix.Data());
}
inline TString FileLabel(const FileSpec &s) {
  return Form("run%d%s", s.run, s.suffix.Data());
}

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
