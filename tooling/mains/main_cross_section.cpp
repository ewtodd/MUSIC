// cross-section: absolute (a,xn) cross section per reaction strip. The work
// is CrossSection (tooling/src/CrossSection.cpp); this only sets the output
// locations and runs it.
#include "CrossSection.hpp"
#include "InitUtils.hpp"
#include "Paths.hpp"
#include <TROOT.h>

int main() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                Paths::ResultsDir() + "/plots",
                                Paths::ResultsDir() + "/root_files");
  gROOT->SetBatch(kTRUE);
  CrossSection xs;
  return xs.Run() ? 0 : 1;
}
