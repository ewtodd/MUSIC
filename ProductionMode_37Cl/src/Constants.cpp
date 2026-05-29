#include "Constants.hpp"

namespace Paths {

TString ProjectRootOf(const char *file) {
  TString path = file;
  path.ReplaceAll("/./", "/");
  TString src_dir = gSystem->DirName(path);
  return gSystem->DirName(src_dir);
}

} // namespace Paths

namespace Constants {

Long64_t LookupTTFOffsetPs(Int_t board, Int_t channel) {
  std::map<std::pair<Int_t, Int_t>, Long64_t>::const_iterator it =
      ttfOffsetPs.find(std::pair<Int_t, Int_t>(board, channel));
  return (it == ttfOffsetPs.end()) ? 0 : it->second;
}

} // namespace Constants
