#ifndef REMIX_SIM_HPP
#define REMIX_SIM_HPP

#include "Paths.hpp"
#include <Rtypes.h>
#include <TString.h>
#include <TSystem.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

/**
 * @brief Locating and naming the simulated data for a dataset's reactions.
 *
 * Simulations come from
 * [Remix-MUSIC-Sim](https://github.com/ewtodd/Remix-MUSIC-Sim), driven by the
 * TOML control files in `<dataset>/control`. Each control file that writes a
 * `traces_<iso>_` output contributes one simulation to compare against.
 */
class RemixSim {
public:
  /**
   * @brief One simulated dataset, identified by its tag.
   *
   * The tag becomes both the on-disk basename
   * (`sim_root_files/traces_<iso>_<tag>.root`) and the plot subdirectory leaf,
   * so one string names the simulation everywhere it appears.
   */
  struct SimFileSpec {
    TString tag; ///< Simulation tag, from the control file's output line.
  };

  /// @brief Absolute path to the control directory, `<dataset>/control`.
  static TString ControlDir();

  /**
   * @brief Extract the simulation tag from a control file.
   *
   * Reads the `output = ".../traces_<iso>_<tag>.root"` line.
   *
   * @param filepath Control file to read.
   * @return The tag, or an empty string for a control file that writes
   *         something other than a `traces_<iso>_` output.
   */
  static TString TagFromControlFile(const TString &filepath);

  /**
   * @brief Every simulation this dataset defines.
   * @return One spec per control file that writes a `traces_<iso>_` output,
   *         sorted by tag so colour and legend ordering stay stable between
   *         runs.
   */
  static std::vector<SimFileSpec> BuildFileSpecs();

  /// @brief On-disk basename, `traces_<iso>_<tag>`, with no directory or
  /// extension.
  /// @param s Simulation to name.
  static TString TracesName(const SimFileSpec &s);

  /// @brief Absolute path to a simulation's ROOT file.
  /// @param s Simulation to locate.
  /// @return `<dataset>/sim_root_files/traces_<iso>_<tag>.root`.
  static TString SimRootPath(const SimFileSpec &s);

  /**
   * @brief Which strip a simulated reaction occurs on.
   * @param tag Simulation tag.
   * @return The strip index from a trailing `_s<N>` token, or `-1` for an
   *         unreacted-beam tag, which carries no such suffix.
   */
  static Int_t ReactionStripOf(const TString &tag);

  /// @brief The tag with any trailing `_s<N>` reaction-strip token removed.
  /// @param tag Simulation tag; one without a suffix passes through unchanged.
  static TString TagWithoutStrip(const TString &tag);

  /// @brief Whether a tag denotes an energy-resolution simulation.
  /// @param tag Simulation tag.
  static Bool_t IsEresTag(const TString &tag);

private:
  static Bool_t SimFileSpecTagLess(const SimFileSpec &a, const SimFileSpec &b);
};

#endif
