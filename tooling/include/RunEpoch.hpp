#ifndef RUNEPOCH_HPP
#define RUNEPOCH_HPP

#include "DedupStrategy.hpp"
#include <Rtypes.h>
#include <RtypesCore.h>
#include <TString.h>
#include <map>
#include <utility>
#include <vector>

// One acquisition period of a dataset. The 37Cl data spans a CoMPASS era and
// two SOLARIS eras, which differ in digitiser, channel map, board count,
// trigger window and ADC scale, so those settings cannot live as one flat block
// per dataset.
//
// Every epoch declares every field: Constants::Active*() reads the active
// epoch when one is set and the flat DatasetConfig otherwise, with no
// "declared / inherit" heuristic in between. A half-populated epoch is
// therefore a config error rather than a silent fallback to another era's
// value.
//
// Deliberately absent: any selector for the L/R gain-matching method. There is
// one method (the ridge/shoulder anchor plus the eSum pass), it is the one that
// reproduces the reference traces, and making it configurable is how the
// kPerEndBeamPeak regression reached the late runs.
enum RunSource { kCoMPASS, kSolaris };

struct RunEpoch {
  TString name;
  RunSource source;
  Bool_t enabled;

  // Prefix for this epoch's output files and plot directories. Run numbers are
  // only unique within an acquisition system: CoMPASS run 37 and SOLARIS run 37
  // are unrelated experiments years apart that happen to share a number, and
  // without a tag both would write Events_Run37_*.root over each other. Leave
  // empty for the era whose filenames are already on disk (the SOLARIS runs);
  // set it on any era that would otherwise collide.
  TString file_tag;

  std::vector<Int_t> runs;
  // Cap on subfiles processed per run; -1 for all.
  Int_t max_files;

  // Digitiser layout.
  Int_t n_boards;
  Int_t n_channels;
  std::map<std::pair<Int_t, Int_t>, TString> channel_map;

  // Timing / event building.
  UShort_t timing_ref_board;
  std::vector<UShort_t> timing_ref_board_channels;
  Bool_t do_board_sync;
  Bool_t do_sort;
  Double_t event_time_window_us;
  TString reference_channel;
  Double_t reference_channel_min_adc;
  Double_t reference_channel_max_adc;
  DedupStrategy dedup_strategy;

  Bool_t has_cathode;

  // ADC scale. The CoMPASS era digitised to a different full scale than
  // SOLARIS, so the per-end ceilings and the strip range differ.
  Double_t strip_e_min_adc;
  Double_t strip_e_max_adc;
  Double_t cathode_max_adc;
  Double_t grid_max_adc;
  Double_t strip0_max_adc;
  Double_t strip17_max_adc;
  Double_t left_even_max_adc;
  Double_t left_odd_max_adc;
  Double_t right_even_max_adc;
  Double_t right_odd_max_adc;

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
