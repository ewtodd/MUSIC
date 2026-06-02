#include "CalibrateBeam.hpp"
#include <TString.h>

int main(int argc, char **argv) {
  TString file_label = (argc > 1) ? TString(argv[1]) : TString("");
  CalibrateBeam::Run(file_label);
  return 0;
}
