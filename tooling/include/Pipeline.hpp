#ifndef PIPELINE_HPP
#define PIPELINE_HPP

#include "BinaryToRoot.hpp"
#include "BinaryUtils.hpp"
#include "CalibrateBeam.hpp"
#include "Constants.hpp"
#include "EventBuilder.hpp"
#include "EventsSummary.hpp"
#include "FileSet.hpp"
#include "GpuAccel.hpp"
#include "HitAdapter.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Timing.hpp"
#include <TError.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <queue>
#include <set>
#include <streambuf>
#include <thread>
#include <utility>
#include <vector>

/**
 * @brief The main entry point: raw binaries through to built events.
 *
 * For each CoMPASS binary subfile, Run() performs: binary to raw hits,
 * multi-board timing alignment, GPU timestamp sort, event building,
 * beam-energy calibration, and trace creation.
 *
 * Subfiles are processed in parallel, over the runs selected by `RUNS` and
 * `N_CHUNKS` in the dataset configuration.
 */
class Pipeline {
public:
  /**
   * @brief Run the whole pipeline over the configured subfiles.
   *
   * Takes no arguments — every knob lives in `Constants.hpp` /
   * `Constants.cpp` and the `control/` TOMLs. Output goes to
   * Paths::ResultsDir(), and progress is logged to
   * `<dataset dir>/pipeline_fused.log`.
   *
   * @note Concurrency is bounded by the dataset configuration, and GPU sorts
   *       additionally by GpuAccel::TryAcquireSortSlot(); a worker that cannot
   *       get a GPU slot takes the CPU path rather than waiting.
   */
  static void Run();
};

#endif
