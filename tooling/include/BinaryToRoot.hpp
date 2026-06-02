#ifndef BINARY_TO_ROOT_HPP
#define BINARY_TO_ROOT_HPP

#include "FileSet.hpp"
#include "IOUtils.hpp"
#include <Rtypes.h>
#include <TString.h>
#include <TSystem.h>
#include <fstream>
#include <ios>
#include <string>

class BinaryToRoot {
public:
  static TString DataRBaseName(const FileSpec &s);
  static TString HeaderSidecarName(Int_t run);
  static void WriteHeaderSidecar(Int_t run, UShort_t header);
  static Bool_t ReadHeaderSidecar(Int_t run, UShort_t &header);
};

#endif
