#ifndef FILE_SET_HPP
#define FILE_SET_HPP

#include "Constants.hpp"
#include "IOUtils.hpp"
#include <Rtypes.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <mutex>
#include <vector>

// Serializes all plotting/canvas work; ROOT's global graphics state is not
// thread-safe and the pipeline runs file workers concurrently.
extern std::mutex g_plot_mutex;

struct FileSpec {
  Int_t run;
  TString suffix;
};

class FileSet {
public:
  static TString CompassBinPath(const FileSpec &s);
  static std::vector<TString> DiscoverRunSuffixes(Int_t run);
  static std::vector<FileSpec> BuildFileSpecs();
  static TString RawRootName(const FileSpec &s);
  static TString ShiftFriendName(const FileSpec &s);
  static TString EventsName(const FileSpec &s);
  static TString FileLabel(const FileSpec &s);

  static FileSpec ResolveFileSpec(const TString &file_label);

  // Group every configured subfile's events tree by run number into one TChain
  // per run. Each chain carries the per-subfile "calibration" gain table in its
  // files, so an EnergyView attached to the chain calibrates on the fly and
  // reloads gains at each file boundary. Subfiles missing on disk are skipped
  // with a warning. run_order receives the run numbers in first-seen order; the
  // caller owns and deletes the returned chains.
  static std::map<Int_t, TChain *>
  GroupEventsByRun(std::vector<Int_t> &run_order);

  // Stride for visiting at most max_points entries spread across n_total (so a
  // scan covers the whole chain rather than just its head). <=0 max_points
  // means "visit all". Always >= 1.
  static Long64_t SampleStride(Long64_t n_total, Long64_t max_points);
};

#endif
