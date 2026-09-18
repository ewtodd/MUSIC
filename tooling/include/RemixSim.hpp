#ifndef REMIX_SIM_HPP
#define REMIX_SIM_HPP

#include "Paths.hpp"
#include <Rtypes.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
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

  /**
   * @brief One entry of a simulation's `events_MeV` tree.
   *
   * The layout is the experimental events tree's: strips 1-16 are read at a
   * left and a right end, held in arrays of 16 indexed by `strip - 1`; strips
   * 0 and 17 are unsegmented and each a single value. Energies are MeV.
   * Left() and Right() give the 18-strip view the analysis works in, with the
   * unsegmented strips counted wholly on the left and zero on the right.
   */
  struct Event {
    Float_t left[16];  ///< Left end of strips 1-16, index `strip - 1`.
    Float_t right[16]; ///< Right end of strips 1-16, index `strip - 1`.
    Float_t strip0;    ///< Strip 0, unsegmented.
    Float_t strip17;   ///< Strip 17, unsegmented.

    Event();

    /**
     * @brief Bind to an `events_MeV` tree and set up the branch addresses.
     * @param t Tree to read. Borrowed; must outlive this event.
     * @return `kTRUE` if the expected branches were found.
     */
    Bool_t Attach(TTree *t);

    /// @brief Left end of strip `s` (1-16), or the whole of strip 0 or 17.
    Double_t Left(Int_t s) const;
    /// @brief Right end of strip `s` (1-16); zero for strips 0 and 17.
    Double_t Right(Int_t s) const;
    /// @brief The strip's deposit, `Left(s) + Right(s)`.
    Double_t Total(Int_t s) const { return Left(s) + Right(s); }
  };

private:
  static Bool_t SimFileSpecTagLess(const SimFileSpec &a, const SimFileSpec &b);
};

#endif
