#include "FileSet.hpp"

std::mutex g_plot_mutex;

TString FileSet::CompassBinPath(const FileSpec &s) {
  return Constants::COMPASS_BASE_DIR +
         Form("run_%d/RAW/DataR_run_%d%s.BIN", s.run, s.run, s.suffix.Data());
}

std::vector<TString> FileSet::DiscoverRunSuffixes(Int_t run) {
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

std::vector<FileSpec> FileSet::BuildFileSpecs() {
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

TString FileSet::RawRootName(const FileSpec &s) {
  return Form("DataR_run_%d%s.root", s.run, s.suffix.Data());
}

TString FileSet::ShiftFriendName(const FileSpec &s) {
  return Form("DataR_run_%d%s.shift.root", s.run, s.suffix.Data());
}

TString FileSet::EventsName(const FileSpec &s) {
  return Form("Events_Run%d%s", s.run, s.suffix.Data());
}

TString FileSet::CalSidecarName(const FileSpec &s) {
  return Form("Events_Run%d%s.cal.root", s.run, s.suffix.Data());
}

TString FileSet::FileLabel(const FileSpec &s) {
  return Form("run%d%s", s.run, s.suffix.Data());
}

TFile *FileSet::AttachCalSidecar(TTree *events, const FileSpec &spec) {
  TString cal_subpath = CalSidecarName(spec);
  TString full = IO::GetRootFilesBaseDir() + "/" + cal_subpath;
  if (gSystem->AccessPathName(full))
    return nullptr;
  TFile *cal = IO::OpenForReading(cal_subpath);
  if (!cal || cal->IsZombie()) {
    if (cal)
      delete cal;
    return nullptr;
  }
  if (!cal->Get("events_cal")) {
    cal->Close();
    delete cal;
    return nullptr;
  }
  events->AddFriend("events_cal", full);
  return cal;
}

FileSpec FileSet::ResolveFileSpec(const TString &file_label) {
  std::vector<FileSpec> specs = BuildFileSpecs();
  for (Int_t k = 0; k < Int_t(specs.size()); k++) {
    if (FileLabel(specs[k]) == file_label)
      return specs[k];
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
