#ifndef RUNEPOCH_HPP
#define RUNEPOCH_HPP

#include "DedupStrategy.hpp"
#include <Rtypes.h>
#include <RtypesCore.h>
#include <TString.h>
#include <map>
#include <utility>
#include <vector>

/// @brief Which acquisition system produced an epoch's data.
enum RunSource {
  kCoMPASS, ///< CAEN CoMPASS.
  kSolaris  ///< SOLARIS DAQ.
};

/**
 * @brief One acquisition period of a dataset.
 *
 * The 37Cl data spans a CoMPASS era and two SOLARIS eras, which differ in
 * digitiser, channel map, board count, trigger window and ADC scale. Those
 * settings therefore cannot live as one flat block per dataset.
 *
 * @warning **Every epoch must declare every field.** `Constants::Active*()`
 *          reads the active epoch when one is set and the flat `DatasetConfig`
 *          otherwise, with no "declared / inherit" heuristic in between. A
 *          half-populated epoch is a configuration error, not a silent fallback
 *          to another era's value — and it will not announce itself.
 *
 * @note Deliberately absent: any selector for the left/right gain-matching
 *       method. There is one method — the ridge/shoulder anchor plus the eSum
 *       pass — it is the one that reproduces the reference traces, and making
 *       it configurable is how the `kPerEndBeamPeak` regression reached the
 *       late runs.
 */
struct RunEpoch {
  TString name;     ///< Epoch name, used in logs and plot paths.
  RunSource source; ///< Acquisition system for this era.
  Bool_t enabled;   ///< Whether this epoch participates in the analysis.

  /// Prefix for this epoch's output files and plot directories.
  ///
  /// Run numbers are only unique **within** an acquisition system: CoMPASS run
  /// 37 and SOLARIS run 37 are unrelated experiments years apart that happen to
  /// share a number, and without a tag both would write `Events_Run37_*.root`
  /// over each other. Leave empty for the era whose filenames are already on
  /// disk (the SOLARIS runs); set it on any era that would otherwise collide.
  TString file_tag;

  std::vector<Int_t> runs; ///< Run numbers belonging to this epoch.
  /// Cap on subfiles processed per run; `-1` for all of them.
  Int_t max_files;

  /// @name Digitiser layout
  /// @{
  Int_t n_boards;   ///< Boards in this era's setup.
  Int_t n_channels; ///< Channels per board.
  /// Board and channel to detector-element name, e.g. `"Grid"`, `"de_l_03"`.
  std::map<std::pair<Int_t, Int_t>, TString> channel_map;
  /// @}

  /// @name Timing and event building
  /// @{
  UShort_t timing_ref_board;
  /// Channels on the reference board used to align the others against it.
  std::vector<UShort_t> timing_ref_board_channels;
  Bool_t do_board_sync; ///< Run the multi-board timing alignment.
  Bool_t do_sort;       ///< Time-sort hits before event building.
  /// Coincidence window, in microseconds, for grouping hits into one event.
  Double_t event_time_window_us;
  TString reference_channel; ///< Channel whose hits seed events, e.g. `"Grid"`.
  Double_t
      reference_channel_min_adc; ///< Lower energy gate on the seed channel.
  Double_t
      reference_channel_max_adc; ///< Upper energy gate on the seed channel.
  DedupStrategy dedup_strategy;  ///< How repeated hits on a slot resolve.
  /// @}

  Bool_t has_cathode; ///< Whether this era instrumented the cathode.

  /// @name ADC scale
  /// The CoMPASS era digitised to a different full scale than SOLARIS, so the
  /// per-end ceilings and the strip range differ between epochs.
  /// @{
  Double_t strip_e_min_adc;    ///< Lower bound of the per-strip energy range.
  Double_t strip_e_max_adc;    ///< Upper bound of the per-strip energy range.
  Double_t cathode_max_adc;    ///< Cathode full scale.
  Double_t grid_max_adc;       ///< Grid full scale.
  Double_t strip0_max_adc;     ///< Strip 0 full scale.
  Double_t strip17_max_adc;    ///< Strip 17 full scale.
  Double_t left_even_max_adc;  ///< Left-side ceiling, even strips.
  Double_t left_odd_max_adc;   ///< Left-side ceiling, odd strips.
  Double_t right_even_max_adc; ///< Right-side ceiling, even strips.
  Double_t right_odd_max_adc;
  /// @}

  /// @brief Construct with SOLARIS-era defaults; override per epoch.
  RunEpoch()
      : name(""), source(kSolaris), enabled(kTRUE), file_tag(""), max_files(-1),
        n_boards(1), n_channels(64), timing_ref_board(0), do_board_sync(kFALSE),
        do_sort(kFALSE), event_time_window_us(8.0), reference_channel("Grid"),
        reference_channel_min_adc(0.0), reference_channel_max_adc(16384.0),
        dedup_strategy(kLARGEST_ENERGY), has_cathode(kFALSE),
        strip_e_min_adc(0.0), strip_e_max_adc(4096.0), cathode_max_adc(16384.0),
        grid_max_adc(16384.0), strip0_max_adc(16384.0),
        strip17_max_adc(16384.0), left_even_max_adc(16384.0),
        left_odd_max_adc(16384.0), right_even_max_adc(16384.0),
        right_odd_max_adc(16384.0) {}
};

#endif
