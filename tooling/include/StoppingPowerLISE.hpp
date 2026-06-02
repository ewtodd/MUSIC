#ifndef STOPPING_POWER_LISE_HPP
#define STOPPING_POWER_LISE_HPP

#include <Rtypes.h>
#include <TString.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Reader for LISE++ stopping-power tables (tab-delimited, 7 models per row).
// Dataset-independent: the caller picks the isotope-specific .lise file.
class Lise {
public:
  // Linearly interpolated dE/dx for `modelIndex` at `energy_MeV_per_u`.
  static Double_t StoppingPower(const TString &filename, Int_t modelIndex,
                                Double_t energy_MeV_per_u) {

    std::ifstream file(filename.Data());
    if (!file.is_open()) {
      std::cerr << "ERROR: Could not open file " << filename << std::endl;
      return -1.0;
    }

    std::string line;
    for (Int_t i = 0; i < 3; i++) {
      std::getline(file, line);
    }

    std::vector<Double_t> energies;
    std::vector<std::vector<Double_t>> dedx_values;

    while (std::getline(file, line)) {
      if (line.empty() || line[0] == '!')
        continue;

      std::istringstream iss(line);
      std::vector<std::string> tokens;
      std::string token;

      while (std::getline(iss, token, '\t')) {
        tokens.push_back(token);
      }

      if (tokens.size() < 2)
        continue;

      try {
        Double_t E = std::stod(tokens[0]);

        if (energies.empty() || energies.back() != E) {
          energies.push_back(E);
          std::vector<Double_t> dedx_for_models(7);

          for (Int_t m = 0; m < 7; m++) {
            Int_t col_index = 1 + m * 2;
            if (col_index < (Int_t)tokens.size()) {
              dedx_for_models[m] = std::stod(tokens[col_index]);
            }
          }

          dedx_values.push_back(dedx_for_models);
        }
      } catch (const std::exception &e) {
        continue;
      }
    }

    file.close();

    if (energies.size() < 2) {
      std::cerr << "ERROR: Not enough data points in " << filename << std::endl;
      return -1.0;
    }

    Int_t idx = 0;
    if (energy_MeV_per_u <= energies.front()) {
      idx = 0;
    } else if (energy_MeV_per_u >= energies.back()) {
      idx = (Int_t)energies.size() - 2;
    } else {
      for (idx = 0; idx < (Int_t)energies.size() - 1; idx++) {
        if (energies[idx] <= energy_MeV_per_u &&
            energy_MeV_per_u <= energies[idx + 1]) {
          break;
        }
      }
    }

    Double_t E1 = energies[idx];
    Double_t E2 = energies[idx + 1];
    Double_t dedx1 = dedx_values[idx][modelIndex];
    Double_t dedx2 = dedx_values[idx + 1][modelIndex];

    return dedx1 + (energy_MeV_per_u - E1) / (E2 - E1) * (dedx2 - dedx1);
  }
};

#endif
