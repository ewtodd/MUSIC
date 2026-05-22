#include <TString.h>

void CalibrateBeam(TString file_label = "");

int main(int argc, char **argv) {
  TString file_label = (argc > 1) ? TString(argv[1]) : TString("");
  CalibrateBeam(file_label);
  return 0;
}
