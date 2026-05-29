#ifndef FILE_SET_HPP
#define FILE_SET_HPP

#include "Constants.hpp"
#include "IOUtils.hpp"
#include <Rtypes.h>
#include <TFile.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <iostream>
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
  static TString CalSidecarName(const FileSpec &s);
  static TString FileLabel(const FileSpec &s);

  // Attach the per-subfile calibration sidecar as a TTree friend so EnergyView
  // (via FindBranch) sees TotaldEMeV/etc. Returns the opened sidecar TFile
  // (caller closes after the events file is done) or nullptr if no sidecar.
  static TFile *AttachCalSidecar(TTree *events, const FileSpec &spec);

  static FileSpec ResolveFileSpec(const TString &file_label);
};

#endif
