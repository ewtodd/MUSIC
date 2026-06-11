#include "FileSet.hpp"

std::mutex g_plot_mutex;

TString FileSet::CompassBinPath(const FileSpec &s) {
  return Constants::COMPASS_BASE_DIR +
         Form("run_%d/RAW/DataR_run_%d%s.BIN", s.run, s.run, s.suffix.Data());
}

namespace {
std::vector<TString> DiscoverSuffixesIn(const TString &dir,
                                        const TString &prefix,
                                        const TString &ext) {
  std::vector<TString> suffixes;
  void *dirp = gSystem->OpenDirectory(dir);
  if (!dirp) {
    std::cerr << "DiscoverSuffixesIn: cannot open " << dir << std::endl;
    return suffixes;
  }
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

std::vector<FileSpec> BuildSpecsImpl(Bool_t processed) {
  std::vector<FileSpec> specs;
  for (Int_t r = 0; r < Int_t(Constants::RUN_NUMBERS.size()); r++) {
    Int_t run = Constants::RUN_NUMBERS[r];
    if (Constants::N_FILES < 0) {
      std::vector<TString> suffixes =
          processed ? FileSet::DiscoverProcessedRunSuffixes(run)
                    : FileSet::DiscoverRunSuffixes(run);
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
} // namespace

std::vector<TString> FileSet::DiscoverRunSuffixes(Int_t run) {
  return DiscoverSuffixesIn(Constants::COMPASS_BASE_DIR +
                                Form("run_%d/RAW/", run),
                            Form("DataR_run_%d", run), ".BIN");
}

std::vector<TString> FileSet::DiscoverProcessedRunSuffixes(Int_t run) {
  return DiscoverSuffixesIn(IO::GetRootFilesBaseDir(),
                            Form("Events_Run%d", run), ".root");
}

std::vector<FileSpec> FileSet::BuildFileSpecs() {
  return BuildSpecsImpl(kFALSE);
}

std::vector<FileSpec> FileSet::BuildProcessedFileSpecs() {
  return BuildSpecsImpl(kTRUE);
}

std::vector<FileSpec> FileSet::BuildRawOrProcessedFileSpecs() {
  std::vector<FileSpec> specs = BuildFileSpecs();
  std::vector<FileSpec> processed = BuildProcessedFileSpecs();
  for (Int_t k = 0; k < Int_t(processed.size()); k++) {
    Bool_t already = kFALSE;
    for (Int_t j = 0; j < Int_t(specs.size()); j++) {
      if (specs[j].run == processed[k].run &&
          specs[j].suffix == processed[k].suffix) {
        already = kTRUE;
        break;
      }
    }
    if (!already)
      specs.push_back(processed[k]);
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

TString FileSet::FileLabel(const FileSpec &s) {
  return Form("run%d%s", s.run, s.suffix.Data());
}

std::map<Int_t, TChain *>
FileSet::GroupEventsByRun(std::vector<Int_t> &run_order) {
  std::map<Int_t, TChain *> chain_by_run;
  std::vector<FileSpec> all_specs = BuildProcessedFileSpecs();
  for (Int_t i = 0; i < Int_t(all_specs.size()); i++) {
    const FileSpec &s = all_specs[i];
    TString full = IO::GetRootFilesBaseDir() + "/" + EventsName(s) + ".root";
    if (gSystem->AccessPathName(full)) {
      std::cerr << "Missing events file: " << full << std::endl;
      continue;
    }
    if (chain_by_run.find(s.run) == chain_by_run.end()) {
      chain_by_run[s.run] = new TChain("events");
      run_order.push_back(s.run);
    }
    chain_by_run[s.run]->Add(full);
  }
  return chain_by_run;
}

Long64_t FileSet::SampleStride(Long64_t n_total, Long64_t max_points) {
  Long64_t n_visit =
      (max_points > 0 && n_total > max_points) ? max_points : n_total;
  Long64_t stride = (n_visit > 0) ? (n_total / n_visit) : 1;
  if (stride < 1)
    stride = 1;
  return stride;
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
