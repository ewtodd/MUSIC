#ifndef REMIX_SIM_HPP
#define REMIX_SIM_HPP

#include "Paths.hpp"
#include <Rtypes.h>
#include <TString.h>
#include <TSystem.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Sim-side ("Remix") pipeline support. Enumerates one sim FileSpec per control
// file under <dataset>/control that writes a traces_<iso>_ output, and resolves
// the on-disk sim ROOT paths under <dataset>/sim_root_files. The shared
// g_plot_mutex lives in FileSet.hpp; this header does NOT redeclare it.
//
// Sim-side data-blob constants (strip ranges etc.) live in namespace Remix in
// RemixConstants.hpp; this class holds only the sim-pipeline behaviour.
class RemixSim {
public:
  // One sim file. `tag` becomes both the on-disk basename
  // (sim_root_files/traces_<iso>_<tag>.root) and the plot subdirectory leaf.
  struct SimFileSpec {
    TString tag;
  };

  // Absolute path to the control directory: <dataset>/control.
  static TString ControlDir();

  // Pull the traces tag out of a control file's
  // `output = ".../traces_<iso>_<tag>.root"` line. Returns "" for control files
  // that write something other than a traces_<iso>_ file.
  static TString TagFromControlFile(const TString &filepath);

  // Enumerate one SimFileSpec per control file in <dataset>/control that writes
  // a traces_<iso>_ output, sorted by tag for stable colour/legend ordering.
  static std::vector<SimFileSpec> BuildFileSpecs();

  // On-disk basename (no directory, no extension): traces_<iso>_<tag>.
  static TString TracesName(const SimFileSpec &s);

  // Absolute path to the sim ROOT file for this spec:
  // <dataset>/sim_root_files/traces_<iso>_<tag>.root.
  static TString SimRootPath(const SimFileSpec &s);

  // Plot/label tag: <iso>_<tag>.
  static TString FileLabel(const SimFileSpec &s);

  // Reaction strip encoded as a trailing "_s<N>" token on the tag. Returns -1
  // for unreacted-beam tags (no suffix).
  static Int_t ReactionStripOf(const TString &tag);

  // Tag with the trailing "_s<N>" reaction-strip token removed. Tags without a
  // suffix pass through unchanged.
  static TString TagWithoutStrip(const TString &tag);

  static Bool_t IsEresTag(const TString &tag);

private:
  static Bool_t SimFileSpecTagLess(const SimFileSpec &a, const SimFileSpec &b);
};

#endif
