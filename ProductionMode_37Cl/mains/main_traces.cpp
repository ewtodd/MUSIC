#include "Traces.hpp"
#include <TString.h>

int main(int argc, char **argv) {
  TString file_label = (argc > 1) ? TString(argv[1]) : TString("");
  Traces::Run(file_label);
  return 0;
}
