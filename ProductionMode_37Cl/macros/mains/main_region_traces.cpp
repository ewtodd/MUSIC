#include <TString.h>

void RegionTraces(TString cut_name, TString file_label);

int main(int argc, char **argv) {
  TString cut_name = (argc > 1) ? TString(argv[1]) : TString("");
  TString file_label = (argc > 2) ? TString(argv[2]) : TString("");
  RegionTraces(cut_name, file_label);
  return 0;
}
