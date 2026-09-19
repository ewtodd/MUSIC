#include "BinaryUtils.hpp"
#include "Constants.hpp"
#include "FileSet.hpp"
#include <condition_variable>
#include <iostream>
#include <mutex>
#include <queue>
#include <set>
#include <sstream>
#include <thread>
#include <vector>

// Check if a SOL file uses Minimum format by reading the first block header.
Bool_t IsMinimumFormat(const char *filePath) {
  std::ifstream file(filePath, std::ios::binary);
  if (!file.is_open()) {
    return kFALSE;
  }

  UShort_t blockHeader;
  file.read(reinterpret_cast<char *>(&blockHeader), sizeof(UShort_t));
  if (file.fail()) {
    file.close();
    return kFALSE;
  }
  file.close();

  if ((blockHeader & 0xAA00) != 0xAA00) {
    return kFALSE;
  }

  UChar_t dataType = blockHeader & 0xF;
  return dataType == SOLData::Minimum;
}

// Discover original SOL files for a run from the base dir (not split dir).
std::vector<TString> DiscoverSolRunSuffixesFromBase(Int_t run) {
  std::vector<TString> suffixes;
  TString sol_dir = Constants::cfg.SOL_BASE_DIR;
  void *dirp = gSystem->OpenDirectory(sol_dir);
  if (!dirp) {
    std::cerr << "DiscoverSolRunSuffixesFromBase: cannot open " << sol_dir
              << std::endl;
    return suffixes;
  }

  TString prefix = Form("music_exp1915_%03d_00_66222_", run);
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

struct WorkItem {
  Int_t run;
  TString suffix;
  Double_t chunk_seconds; // the run's epoch's chunk length
};

struct SplitResult {
  Int_t nSplit;
  Int_t nSkipped;
  Int_t nMissing;
  Int_t nAlreadySplit;
};

SplitResult SplitWorker(std::queue<WorkItem> &work, std::mutex &work_mutex,
                        const char *outputDir) {
  SplitResult result;
  result.nSplit = 0;
  result.nSkipped = 0;
  result.nMissing = 0;
  result.nAlreadySplit = 0;

  std::mutex log_mutex;

  while (true) {
    WorkItem item;
    {
      std::lock_guard<std::mutex> lk(work_mutex);
      if (work.empty())
        break;
      item = work.front();
      work.pop();
    }

    FileSpec spec;
    spec.run = item.run;
    spec.suffix = item.suffix;
    TString solPath = FileSet::SolBinPath(spec);

    if (gSystem->AccessPathName(solPath)) {
      {
        std::lock_guard<std::mutex> lk(log_mutex);
        std::cerr << "  [missing] " << solPath.Data() << std::endl;
      }
      result.nMissing++;
      continue;
    }

    TString baseName = solPath;
    Int_t lastSlash = baseName.Last('/');
    if (lastSlash >= 0) {
      baseName = baseName(lastSlash + 1, baseName.Length() - lastSlash - 1);
    }

    // Check if already split
    TString chunk0Path = TString(outputDir) + "/" + baseName + "_chunk000.sol";
    if (!gSystem->AccessPathName(chunk0Path)) {
      // Counted in the closing summary; one line per file is noise on a
      // re-run over hundreds of subfiles.
      result.nAlreadySplit++;
      continue;
    }

    // Check format
    if (!IsMinimumFormat(solPath.Data())) {
      {
        std::lock_guard<std::mutex> lk(log_mutex);
        std::cerr << "  [skip] Not Minimum format: " << baseName << std::endl;
      }
      result.nSkipped++;
      continue;
    }

    Int_t totalBlocks = 0;
    Int_t totalChunks = 0;
    std::vector<TString> outputFiles = SOLReader::SplitSolFileByTime(
        solPath.Data(), outputDir, item.chunk_seconds, totalBlocks,
        totalChunks);

    if (outputFiles.empty()) {
      {
        std::lock_guard<std::mutex> lk(log_mutex);
        std::cerr << "  [error] Failed to split: " << solPath.Data()
                  << std::endl;
      }
      result.nSkipped++;
      continue;
    }

    {
      std::lock_guard<std::mutex> lk(log_mutex);
      std::cout << "  [split] " << baseName << ": " << totalBlocks
                << " blocks -> " << totalChunks << " chunks" << std::endl;
    }
    result.nSplit++;
  }

  return result;
}

int main(int argc, char *argv[]) {
  // A chunk length on the command line overrides every epoch's own.
  Double_t forcedSeconds = 0.0;
  Bool_t forced = kFALSE;
  Int_t nWorkers = Constants::cfg.SOL_N_SPLIT_WORKERS;

  if (argc >= 2) {
    forcedSeconds = std::stod(argv[1]);
    forced = kTRUE;
  }
  if (argc >= 3) {
    nWorkers = std::stoi(argv[2]);
  }

  // The runs to split, each with its chunk length: every enabled SOLARIS
  // epoch's runs at the epoch's `split_chunk_seconds`, or the flat list at
  // SOL_SPLIT_CHUNK_SECONDS. CoMPASS has no .sol and is skipped. A
  // non-positive length means no splitting for those runs: the pipeline then
  // reads their whole files from SOL_BASE_DIR (FileSet::DiscoverSolRunSuffixes
  // ignores the split directory for that epoch).
  std::vector<std::pair<Int_t, Double_t>> runs; // run, chunk seconds
  if (Constants::cfg.EPOCHS.empty()) {
    const Double_t sec =
        forced ? forcedSeconds : Constants::cfg.SOL_SPLIT_CHUNK_SECONDS;
    for (Int_t r = 0; r < Int_t(Constants::cfg.RUN_NUMBERS.size()); r++)
      runs.push_back(std::make_pair(Constants::cfg.RUN_NUMBERS[r], sec));
  } else {
    for (Int_t e = 0; e < Int_t(Constants::cfg.EPOCHS.size()); e++) {
      const RunEpoch &ep = Constants::cfg.EPOCHS[e];
      if (!ep.enabled || ep.source != kSolaris)
        continue;
      const Double_t sec = forced ? forcedSeconds : ep.split_chunk_seconds;
      std::cout << "epoch " << ep.name << ": " << ep.runs.size() << " run(s), ";
      if (sec > 0)
        std::cout << sec << "s chunks" << std::endl;
      else
        std::cout << "no splitting (whole files)" << std::endl;
      for (Int_t r = 0; r < Int_t(ep.runs.size()); r++)
        runs.push_back(std::make_pair(ep.runs[r], sec));
    }
  }
  std::queue<WorkItem> work;
  for (Int_t r = 0; r < Int_t(runs.size()); r++) {
    if (!(runs[r].second > 0))
      continue;
    // A run declared twice takes its first epoch's length.
    Bool_t seen = kFALSE;
    for (Int_t q = 0; q < r && !seen; q++)
      seen = runs[q].first == runs[r].first;
    if (seen)
      continue;
    std::vector<TString> suffixes =
        DiscoverSolRunSuffixesFromBase(runs[r].first);
    for (Int_t k = 0; k < Int_t(suffixes.size()); k++) {
      WorkItem item;
      item.run = runs[r].first;
      item.suffix = suffixes[k];
      item.chunk_seconds = runs[r].second;
      work.push(item);
    }
  }
  if (work.empty()) {
    std::cout << "SOLARIS preprocessing: nothing to split; the pipeline reads "
                 "the whole files from "
              << Constants::cfg.SOL_BASE_DIR.Data() << std::endl;
    return 0;
  }
  std::cout << "SOLARIS preprocessing: splitting Minimum files (" << nWorkers
            << " workers)" << std::endl;
  std::cout << "Input dir:  " << Constants::cfg.SOL_BASE_DIR.Data()
            << std::endl;
  std::cout << "Output dir: " << Constants::cfg.SOL_SPLIT_DIR.Data()
            << std::endl;
  std::cout << std::endl;

  gSystem->mkdir(Constants::cfg.SOL_SPLIT_DIR, kTRUE);

  std::cout << "Total files to process: " << work.size() << std::endl;

  // Launch workers
  std::mutex work_mutex;
  std::vector<std::thread> workers;
  std::vector<SplitResult> results(nWorkers);

  for (Int_t w = 0; w < nWorkers; w++) {
    workers.emplace_back([&work, &work_mutex, &results, w]() {
      results[w] =
          SplitWorker(work, work_mutex, Constants::cfg.SOL_SPLIT_DIR.Data());
    });
  }

  for (Int_t w = 0; w < nWorkers; w++) {
    workers[w].join();
  }

  Int_t totalSplit = 0;
  Int_t totalSkipped = 0;
  Int_t totalMissing = 0;
  Int_t totalAlreadySplit = 0;
  for (Int_t w = 0; w < nWorkers; w++) {
    totalSplit += results[w].nSplit;
    totalSkipped += results[w].nSkipped;
    totalMissing += results[w].nMissing;
    totalAlreadySplit += results[w].nAlreadySplit;
  }

  std::cout << std::endl;
  std::cout << "Preprocessing complete:" << std::endl;
  std::cout << "  Split:        " << totalSplit << std::endl;
  std::cout << "  Already split:" << totalAlreadySplit << std::endl;
  std::cout << "  Skipped:      " << totalSkipped << std::endl;
  std::cout << "  Missing:      " << totalMissing << std::endl;

  return (totalMissing > 0) ? 1 : 0;
}
