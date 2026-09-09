#ifndef BEAM_ENERGIES_HPP
#define BEAM_ENERGIES_HPP

#include <TString.h>

/**
 * @brief The beam's energy at each anode strip, from the simulated unreacted
 * beam.
 *
 * The simulation is calibrated against the measured per-strip energy loss, so
 * each strip's centre-of-mass energy is the simulation's own rather than a
 * nominal dE/dx table's.
 *
 * The energy entering strip 0 is what came through the entrance window **minus
 * what the upstream dead layer of gas took**. MUSIC has 3.6 cm of gas ahead of
 * the first anode strip — worth roughly 25 MeV for 87Rb — so starting the
 * bookkeeping at the window instead puts every strip a full strip and a half
 * too high. Both numbers come from the simulation's own truth record, which is
 * what puts every energy on one calibrated footing.
 *
 * Shared by the cross-section calculation, which reports at these energies, and
 * by talys-xs, which computes over them.
 */
namespace BeamEnergies {

/// @brief Path to the dataset's beam simulation file.
/// @return `CROSS_SECTION_CONFIG.BEAM_SIM_FILE`, resolved under
///         `sim_root_files`.
TString SimPath();

/**
 * @brief Read the per-strip energy loss profile from a beam simulation.
 *
 * @param path Beam simulation file, normally from SimPath().
 * @param[out] dE Mean energy deposited in each of the 18 strips, in MeV. Must
 *                point at 18 elements.
 * @param[out] e_strip0 Lab energy entering strip 0, in MeV — already net of the
 *                      upstream gas dead layer.
 *
 * @return `kFALSE` if the file or the trees inside it are missing, in which
 *         case the outputs are not meaningful.
 */
Bool_t Profile(const TString &path, Double_t *dE, Double_t &e_strip0);

/**
 * @brief Lab energy entering a given strip.
 *
 * What entered strip 0, less the mean loss of every strip before @p reac.
 *
 * @param dE       Per-strip losses from Profile().
 * @param e_strip0 Energy entering strip 0, from Profile().
 * @param reac     Strip index, 0 to 17.
 *
 * @return The lab-frame energy entering that strip, in MeV.
 */
Double_t LabAtStrip(const Double_t *dE, Double_t e_strip0, Int_t reac);

/// @brief Centre-of-mass energy fraction for this dataset's reaction.
/// @return `A_t / (A_b + A_t)`, from the dataset's configuration. Multiply a
///         lab energy by this to get the centre-of-mass energy.
Double_t CmFraction();

} // namespace BeamEnergies

#endif
