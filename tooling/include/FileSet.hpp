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
  static std::vector<TString> DiscoverProcessedRunSuffixes(Int_t run);
  static std::vector<FileSpec> BuildFileSpecs();
  static std::vector<FileSpec> BuildProcessedFileSpecs();
  static std::vector<FileSpec> BuildRawOrProcessedFileSpecs();
  static TString RawRootName(const FileSpec &s);
  static TString ShiftFriendName(const FileSpec &s);
  static TString EventsName(const FileSpec &s);
  static TString FileLabel(const FileSpec &s);

  static FileSpec ResolveFileSpec(const TString &file_label);

  static std::map<Int_t, TChain *>
  GroupEventsByRun(std::vector<Int_t> &run_order);

  static Long64_t SampleStride(Long64_t n_total, Long64_t max_points);
};

#endif
