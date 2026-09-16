// template-match: an (a,n) count per strip by chi-square template matching
// against the mean tagged trace, from a second pass over the data, set beside
// the tagged count and the bootstrap efficiencies of both. The method is
// TemplateMatch (tooling/src/TemplateMatch.cpp); this only sets the output
// locations and runs it. Needs the scatter cache from strip-sum-scatter and,
// for the tag side of the comparison, tag-efficiency's store.
#include "InitUtils.hpp"
#include "Paths.hpp"
#include "TemplateMatch.hpp"
#include <TROOT.h>

int main() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                Paths::ResultsDir() + "/plots",
                                Paths::ResultsDir() + "/root_files");
  gROOT->SetBatch(kTRUE);
  TemplateMatch tm;
  return tm.Run() ? 0 : 1;
}
