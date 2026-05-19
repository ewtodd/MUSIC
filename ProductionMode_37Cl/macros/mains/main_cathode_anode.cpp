#include <TString.h>

void CathodeAnodeRatio(TString file_label);

int main(int argc, char **argv) {
  TString file_label = (argc > 1) ? TString(argv[1]) : TString("");
  CathodeAnodeRatio(file_label);
  return 0;
}
