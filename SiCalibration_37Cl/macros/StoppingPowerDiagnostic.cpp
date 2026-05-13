#include "Constants.hpp"
#include "StoppingPowerLISE.hpp"
#include <iostream>
#include <TROOT.h>
#include <TString.h>

void StoppingPowerDiagnostic() {

  const Int_t A_Cl37 = 37;
  const Double_t titanium_front_mg_cm2 = 0.9;
  const Double_t titanium_back_mg_cm2  = 1.3;
  const TString ti_lise_filename = "lise/37Cl_in_Ti_NOT_MICRON.lise";

  const Double_t ppac_layer_micron = 1.5; // single mylar foil
  const Int_t    ppac_n_layers     = 4;   // four foils per PPAC
  const TString  mylar_lise_filename = "lise/37Cl_in_Mylar.lise";

  const Double_t threshold_MeV = 71.35;

  const char *model_names[7] = {"(unused)", "Ziegler", "ATIMA12LS",
                                 "ATIMA12NoLS", "ATIMA14Weick",
                                 "Electrical", "(unused)"};

  struct Config {
    Int_t         run;  // 0 = no single run maps to this config
    const char   *label;
    Double_t      beam_energy_MeV;
    Bool_t        ppac;
    Bool_t        front_window;
    Bool_t        both_windows;
  };

  Config configs[] = {
    {21, "91.87 MeV |    PPAC | no windows",      91.87, kTRUE,  kFALSE, kFALSE},
    {22, "71.35 MeV |    PPAC | no windows",      71.35, kTRUE,  kFALSE, kFALSE},
    {23, "71.35 MeV | no PPAC | front Ti only",   71.35, kFALSE, kTRUE,  kFALSE},
    {23, "71.35 MeV |    PPAC | front Ti only",   71.35, kTRUE,  kTRUE,  kFALSE},
    { 0, "92.00 MeV |    PPAC | no windows",      92.00, kTRUE,  kFALSE, kFALSE},
    {24, "92.00 MeV | no PPAC | front Ti only",   92.00, kFALSE, kTRUE,  kFALSE},
    {24, "92.00 MeV |    PPAC | front Ti only",   92.00, kTRUE,  kTRUE,  kFALSE},
    {25, "92.00 MeV | no PPAC | both Ti windows", 92.00, kFALSE, kTRUE,  kTRUE },
    {25, "92.00 MeV |    PPAC | both Ti windows", 92.00, kTRUE,  kTRUE,  kTRUE },
  };

  for (auto &cfg : configs) {
    TString header = (cfg.run > 0) ? Form("Run %2d  -  %s", cfg.run, cfg.label)
                                   : Form("(no run)  %s", cfg.label);
    std::cout << "\n=== " << header << " ===" << std::endl;
    std::cout << Form("  %-16s  %12s  %12s  %12s",
                      "Model", "dE_PPAC(MeV)", "dE_Ti(MeV)", "E_final(MeV)")
              << std::endl;

    Double_t sum_dE_ppac = 0, sum_dE_ti = 0, sum_E_final = 0;

    for (Int_t model = 1; model < 5; ++model) {
      Double_t energy = cfg.beam_energy_MeV;
      Double_t dE_ppac = 0;

      if (cfg.ppac) {
        for (Int_t layer = 0; layer < ppac_n_layers; ++layer) {
          Double_t dedx = GetStoppingPowerFromLISE(mylar_lise_filename, model,
                                                   energy / A_Cl37);
          Double_t dE = dedx * ppac_layer_micron;
          dE_ppac += dE;
          energy  -= dE;
        }
      }

      Double_t dE_ti = 0;

      if (cfg.front_window) {
        Double_t dedx_front = GetStoppingPowerFromLISE(ti_lise_filename, model,
                                                       energy / A_Cl37);
        Double_t dE_front = dedx_front * titanium_front_mg_cm2;
        energy -= dE_front;
        dE_ti  += dE_front;
      }

      if (cfg.both_windows) {
        Double_t dedx_back = GetStoppingPowerFromLISE(ti_lise_filename, model,
                                                      energy / A_Cl37);
        Double_t dE_back = dedx_back * titanium_back_mg_cm2;
        energy -= dE_back;
        dE_ti  += dE_back;
      }

      std::cout << Form("  %-16s  %12.4f  %12.4f  %12.4f",
                        model_names[model], dE_ppac, dE_ti, energy)
                << std::endl;

      sum_dE_ppac  += dE_ppac;
      sum_dE_ti    += dE_ti;
      sum_E_final  += energy;
    }

    std::cout << Form("  %-16s  %12.4f  %12.4f  %12.4f",
                      "--- average ---",
                      sum_dE_ppac / 4, sum_dE_ti / 4, sum_E_final / 4)
              << std::endl;
  }

  // Crucial question: can 92 MeV through both Ti windows survive above 71.35 MeV?
  std::cout << "\n=== CRUCIAL QUESTION (Run 25 geometry): 92 MeV through both"
               " Ti windows — can E_final > " << threshold_MeV << " MeV? ==="
            << std::endl;

  for (Bool_t ppac : {kFALSE, kTRUE}) {
    std::cout << (ppac ? "  With PPAC:" : "  No PPAC:  ") << std::endl;
    Bool_t any_yes = kFALSE;
    Double_t sum_E_final = 0;
    for (Int_t model = 1; model < 5; ++model) {
      Double_t energy = 92.00;

      if (ppac) {
        for (Int_t layer = 0; layer < ppac_n_layers; ++layer) {
          Double_t dedx = GetStoppingPowerFromLISE(mylar_lise_filename, model,
                                                   energy / A_Cl37);
          energy -= dedx * ppac_layer_micron;
        }
      }

      Double_t dedx_front = GetStoppingPowerFromLISE(ti_lise_filename, model,
                                                     energy / A_Cl37);
      energy -= dedx_front * titanium_front_mg_cm2;

      Double_t dedx_back = GetStoppingPowerFromLISE(ti_lise_filename, model,
                                                    energy / A_Cl37);
      energy -= dedx_back * titanium_back_mg_cm2;

      Bool_t above = energy > threshold_MeV;
      if (above) any_yes = kTRUE;
      sum_E_final += energy;
      std::cout << Form("    %-16s  E_final = %8.4f MeV  ->  %s",
                        model_names[model], energy,
                        above ? "YES (> 71.35)" : "NO  (<= 71.35)")
                << std::endl;
    }
    Double_t avg = sum_E_final / 4;
    std::cout << Form("    %-16s  E_final = %8.4f MeV  ->  %s",
                      "--- average ---", avg,
                      avg > threshold_MeV ? "YES (> 71.35)" : "NO  (<= 71.35)")
              << std::endl;
    std::cout << "    VERDICT: " << (any_yes ? "YES for at least one model"
                                             : "NO for all models")
              << std::endl;
  }
}
