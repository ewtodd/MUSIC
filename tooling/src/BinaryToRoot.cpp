#include "BinaryToRoot.hpp"

TString BinaryToRoot::HeaderSidecarName(Int_t run) {
  return Form("DataR_run_%d.header", run);
}

void BinaryToRoot::WriteHeaderSidecar(Int_t run, UShort_t header) {
  TString path = IO::GetRootFilesBaseDir() + "/" + HeaderSidecarName(run);
  gSystem->mkdir(gSystem->DirName(path), kTRUE);
  std::ofstream f(path.Data());
  f << "0x" << std::hex << header << std::endl;
}

Bool_t BinaryToRoot::ReadHeaderSidecar(Int_t run, UShort_t &header) {
  TString path = IO::GetRootFilesBaseDir() + "/" + HeaderSidecarName(run);
  if (gSystem->AccessPathName(path))
    return kFALSE;
  std::ifstream f(path.Data());
  std::string s;
  f >> s;
  if (s.empty())
    return kFALSE;
  header = static_cast<UShort_t>(std::stoul(s, nullptr, 0));
  return kTRUE;
}
