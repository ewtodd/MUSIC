// tag-efficiency: how often a real (a,n) event at each reaction strip is left
// out of that strip's count by the tag and the region cut, and how often it
// is counted one strip late instead. The method is TagEfficiency
// (tooling/src/TagEfficiency.cpp); the record it writes is what cross-section
// unfolds by. This only sets the output locations and runs it.
#include "InitUtils.hpp"
#include "Paths.hpp"
#include "TagEfficiency.hpp"
#include <TROOT.h>

int main() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                Paths::ResultsDir() + "/plots",
                                Paths::ResultsDir() + "/root_files");
  gROOT->SetBatch(kTRUE);
  TagEfficiency eff;
  return eff.Run() ? 0 : 1;
}
