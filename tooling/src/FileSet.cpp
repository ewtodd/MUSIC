#include "FileSet.hpp"

std::mutex g_plot_mutex;

TString FileSet::CompassBinPath(const FileSpec &s) {
  return Constants::cfg.COMPASS_BASE_DIR +
         Form("run_%d/RAW/DataR_run_%d%s.BIN", s.run, s.run, s.suffix.Data());
}

TString FileSet::SolBinPath(const FileSpec &s) {
  // SOLARIS naming: music_exp1915_<RUN3DIGITS>_00_66222_<SEQ3DIGITS>.sol
  // suffix is empty for _000, or "_1" -> _001, etc.
  // Chunk suffix "_cNNN" resolves to _chunkNNN.sol in split dir.
  if (s.suffix.BeginsWith("_c")) {
    TString chunkIdx = s.suffix(2, s.suffix.Length() - 2);
    Int_t seq = 0;
    if (s.suffix.Length() > 4) {
      TString seqStr = s.suffix(4, s.suffix.Length() - 4);
      seq = seqStr.Atoi();
    }
    return Constants::cfg.SOL_SPLIT_DIR +
           Form("music_exp1915_%03d_00_66222_%03d_chunk%s.sol", s.run, seq,
                chunkIdx.Data());
  }

  Int_t seq = 0;
  if (s.suffix != "" && s.suffix[0] == '_') {
    TString seq_str = s.suffix(1, s.suffix.Length() - 1);
    seq = seq_str.Atoi();
  }
  return Constants::cfg.SOL_BASE_DIR +
         Form("music_exp1915_%03d_00_66222_%03d.sol", s.run, seq);
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
  for (Int_t r = 0; r < Int_t(Constants::cfg.RUN_NUMBERS.size()); r++) {
    Int_t run = Constants::cfg.RUN_NUMBERS[r];
    std::vector<TString> suffixes;
    if (processed) {
      suffixes = FileSet::DiscoverProcessedRunSuffixes(run);
    } else if (Constants::cfg.USE_SOLARIS_DATA) {
      suffixes = FileSet::DiscoverSolRunSuffixes(run);
    } else {
      suffixes = FileSet::DiscoverRunSuffixes(run);
    }
    Int_t limit = suffixes.size();
    if (Constants::cfg.N_CHUNKS > 0 && Constants::cfg.N_CHUNKS < limit)
      limit = Constants::cfg.N_CHUNKS;
    for (Int_t k = 0; k < limit; k++) {
      FileSpec s;
      s.run = run;
      s.suffix = suffixes[k];
      specs.push_back(s);
    }
  }
  return specs;
}
} // namespace

std::vector<TString> FileSet::DiscoverRunSuffixes(Int_t run) {
  return DiscoverSuffixesIn(Constants::cfg.COMPASS_BASE_DIR +
                                Form("run_%d/RAW/", run),
                            Form("DataR_run_%d", run), ".BIN");
}

std::vector<TString> FileSet::DiscoverSolRunSuffixes(Int_t run) {
  std::vector<TString> suffixes;
  TString prefix = Form("music_exp1915_%03d_00_66222_", run);

  // Check for split chunks first
  TString split_dir = Constants::cfg.SOL_SPLIT_DIR;
  void *dirp = gSystem->OpenDirectory(split_dir);
  if (dirp) {
    const Char_t *name;
    while ((name = gSystem->GetDirEntry(dirp))) {
      TString fname(name);
      if (!fname.BeginsWith(prefix))
        continue;
      if (!fname.EndsWith(".sol"))
        continue;
      if (!fname.Contains("_chunk"))
        continue;

      Int_t chunkPos = fname.Index("_chunk");
      if (chunkPos < 0)
        continue;

      TString seqStr = fname(prefix.Length(), chunkPos - prefix.Length());
      if (!seqStr.IsDigit())
        continue;
      Int_t seq = seqStr.Atoi();

      TString chunkStr = fname(chunkPos + 6, fname.Length() - chunkPos - 6 - 4);
      if (!chunkStr.IsDigit())
        continue;

      if (seq == 0) {
        suffixes.push_back(Form("_c%s", chunkStr.Data()));
      } else {
        suffixes.push_back(Form("_%d_c%s", seq, chunkStr.Data()));
      }
    }
    gSystem->FreeDirectory(dirp);

    if (!suffixes.empty()) {
      std::sort(suffixes.begin(), suffixes.end(),
                [](const TString &a, const TString &b) { return a < b; });
      return suffixes;
    }
  }

  // Fall back to original files
  TString sol_dir = Constants::cfg.SOL_BASE_DIR;
  dirp = gSystem->OpenDirectory(sol_dir);
  if (!dirp) {
    std::cerr << "DiscoverSolRunSuffixes: cannot open " << sol_dir << std::endl;
    return suffixes;
  }

  const Char_t *name;
  while ((name = gSystem->GetDirEntry(dirp))) {
    TString fname(name);
    if (!fname.BeginsWith(prefix))
      continue;
    if (!fname.EndsWith(".sol"))
      continue;

    TString rest = fname(prefix.Length(), fname.Length() - prefix.Length() - 4);
    if (!rest.IsDigit())
      continue;

    Int_t seq = rest.Atoi();
    if (seq == 0) {
      suffixes.push_back("");
    } else {
      suffixes.push_back(Form("_%d", seq));
    }
  }
  gSystem->FreeDirectory(dirp);

  // Sort by sequence number
  std::sort(suffixes.begin(), suffixes.end(),
            [](const TString &a, const TString &b) {
              if (a == "")
                return true;
              if (b == "")
                return false;
              return a.Atoi() < b.Atoi();
            });

  return suffixes;
}

std::vector<TString> FileSet::DiscoverProcessedRunSuffixes(Int_t run) {
  std::vector<TString> suffixes;
  TString dir = IO::GetRootFilesBaseDir();
  TString prefix = Form("Events_Run%d", run);
  void *dirp = gSystem->OpenDirectory(dir);
  if (!dirp)
    return suffixes;
  const Char_t *name;
  while ((name = gSystem->GetDirEntry(dirp))) {
    TString fname(name);
    if (!fname.BeginsWith(prefix))
      continue;
    if (!fname.EndsWith(".root"))
      continue;
    TString rest = fname(prefix.Length(), fname.Length() - prefix.Length() - 5);
    suffixes.push_back(rest);
  }
  gSystem->FreeDirectory(dirp);
  std::sort(suffixes.begin(), suffixes.end(),
            [](const TString &a, const TString &b) {
              if (a == "")
                return true;
              if (b == "")
                return false;
              // Extract leading numeric part for sorting: _1_c000 -> 1, _c000
              // -> 0
              Int_t na = 0, nb = 0;
              if (a.Length() > 1 && a[0] == '_') {
                TString numA = a(1, a.Length() - 1);
                Int_t dash = numA.Index('_');
                if (dash > 0)
                  numA = numA(0, dash);
                if (numA.IsDigit())
                  na = numA.Atoi();
              }
              if (b.Length() > 1 && b[0] == '_') {
                TString numB = b(1, b.Length() - 1);
                Int_t dash = numB.Index('_');
                if (dash > 0)
                  numB = numB(0, dash);
                if (numB.IsDigit())
                  nb = numB.Atoi();
              }
              if (na != nb)
                return na < nb;
              return a < b;
            });
  return suffixes;
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
