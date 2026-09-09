#ifndef BINARY_TO_ROOT_HPP
#define BINARY_TO_ROOT_HPP

#include "FileSet.hpp"
#include "IOUtils.hpp"
#include <Rtypes.h>
#include <TString.h>
#include <TSystem.h>
#include <fstream>
#include <ios>
#include <string>

/**
 * @brief Persistence for a run's CoMPASS global header.
 *
 * The global header says which fields each record in a CoMPASS binary carries,
 * and it is read once when the binary is first converted. Later stages work
 * from the ROOT output and no longer have the binary to hand, so the header is
 * written beside the ROOT files as a small sidecar and read back when needed.
 */
class BinaryToRoot {
public:
  /// @brief Sidecar filename for a run.
  /// @param run Run number.
  /// @return The bare filename, `DataR_run_<run>.header`, with no directory.
  static TString HeaderSidecarName(Int_t run);

  /**
   * @brief Write a run's global header to its sidecar.
   *
   * The file lands beside the run's ROOT output, under the ROOT files base
   * directory; parent directories are created as needed. The value is written
   * as `0x`-prefixed hexadecimal.
   *
   * @param run    Run number.
   * @param header CoMPASS global header word.
   *
   * @note Failures to open the file are silent — the sidecar is a cache, and a
   *       missing one is recoverable by re-reading the binary.
   */
  static void WriteHeaderSidecar(Int_t run, UShort_t header);

  /**
   * @brief Read a run's global header back from its sidecar.
   *
   * @param run    Run number.
   * @param[out] header Set to the stored header on success; untouched
   *                    otherwise. Parsed with a base of `0`, so the `0x` prefix
   *                    is honoured.
   *
   * @return `kTRUE` if the sidecar exists and held a value; `kFALSE` if it is
   *         missing or empty.
   */
  static Bool_t ReadHeaderSidecar(Int_t run, UShort_t &header);
};

#endif
