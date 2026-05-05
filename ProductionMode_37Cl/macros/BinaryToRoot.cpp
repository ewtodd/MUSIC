#include "Constants.hpp"
#include "InitUtils.hpp"
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <ostream>
#include <vector>

void BinaryToRoot() {
  Bool_t reprocess_initial = kTRUE;

  std::vector<TString> filepaths, output_names;

  TString path_prefix = "/home/e-work/LabData/MUSIC/37Cl/ProductionMode/RAW/";
  for (Int_t i = 0; i < Constants::N_FILES; i++) {
    filepaths.push_back(Form("%sDataR_run_37_%d.BIN", path_prefix.Data(), i));
    std::cout << "Processing file: " << std::endl;
    std::cout << filepaths[i] << std::endl;
    output_names.push_back(Form("DataR_run_37_%d", i));
  }

  InitUtils::SetROOTPreferences();
  UShort_t saved_global_header = 0;

  for (Int_t i = 0; i < Constants::N_FILES; i++) {
    UShort_t header_to_use = (i == 0) ? 0 : saved_global_header;
    UShort_t returned_header = InitUtils::ConvertCoMPASSBinToROOT(
        filepaths[i], output_names[i], header_to_use);

    if (i == 0) {
      saved_global_header = returned_header;
      std::cout << "Saved global header: 0x" << std::hex << saved_global_header
                << std::dec << std::endl;
    }
  }
}
