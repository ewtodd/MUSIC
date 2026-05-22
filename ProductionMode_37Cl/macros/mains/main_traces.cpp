#include <TString.h>

void Traces(TString file_label = "");

int main(int argc, char **argv) {
  TString file_label = (argc > 1) ? TString(argv[1]) : TString("");
  Traces(file_label);
  return 0;
}
