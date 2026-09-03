#ifndef BEAM_ENERGIES_HPP
#define BEAM_ENERGIES_HPP

#include <TString.h>

// The beam's energy at each strip, from the simulated unreacted beam.
//
// The simulation is calibrated against the measured per-strip energy loss, so
// the centre-of-mass energy of each strip is the simulation's rather than a
// nominal dE/dx table's. The energy entering strip 0 is what came through the
// entrance window minus what the upstream dead layer of gas took -- MUSIC has
// 3.6 cm of gas before the first anode strip, worth ~25 MeV for 87Rb, and
// starting the bookkeeping at the window instead puts every strip a full
// strip and a half too high. Both numbers come from the simulation's own
// truth record, so every energy is on one calibrated footing. Shared by the
// cross section, which reports at these energies, and talys-xs, which
// computes over them.
namespace BeamEnergies {
// The dataset's beam sim file (CROSS_SECTION_CONFIG.BEAM_SIM_FILE, in
// sim_root_files).
TString SimPath();
// Mean energy deposited in each of the 18 strips, and the energy entering
// strip 0. False when the file or its trees are missing.
Bool_t Profile(const TString &path, Double_t *dE, Double_t &e_strip0);
// Lab energy entering strip `reac`: what entered strip 0, less each earlier
// strip's mean loss.
Double_t LabAtStrip(const Double_t *dE, Double_t e_strip0, Int_t reac);
// Centre-of-mass fraction A_t / (A_b + A_t) from the dataset's config.
Double_t CmFraction();
} // namespace BeamEnergies

#endif
