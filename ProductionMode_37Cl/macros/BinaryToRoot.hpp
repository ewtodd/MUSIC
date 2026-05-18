#ifndef BINARY_TO_ROOT_HPP
#define BINARY_TO_ROOT_HPP

#include "IOUtils.hpp"
#include "Pipeline.hpp"
#include <TString.h>
#include <TSystem.h>
#include <fstream>
#include <ios>
#include <string>

inline TString DataRBaseName(const FileSpec &s) {
  return Form("DataR_run_%d%s", s.run, s.suffix.Data());
}

inline TString HeaderSidecarName(Int_t run) {
  return Form("DataR_run_%d.header", run);
}

inline void WriteHeaderSidecar(Int_t run, UShort_t header) {
  TString path = IO::GetRootFilesBaseDir() + "/" + HeaderSidecarName(run);
  gSystem->mkdir(gSystem->DirName(path), kTRUE);
  std::ofstream f(path.Data());
  f << "0x" << std::hex << header << std::endl;
}

inline Bool_t ReadHeaderSidecar(Int_t run, UShort_t &header) {
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

#endif
