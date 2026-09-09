#ifndef EVENT_BUILDER_HPP
#define EVENT_BUILDER_HPP

#include "BinaryUtils.hpp"
#include "DedupStrategy.hpp"
#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "PlottingUtils.hpp"
#include "SlotLayout.hpp"
#include <Rtypes.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TObject.h>
#include <TParameter.h>
#include <TString.h>
#include <TTree.h>
#include <array>
#include <iostream>
#include <map>
#include <mutex>
#include <utility>
#include <vector>

/**
 * @brief One built event: the per-slot energies and the flags seen while
 * building it.
 *
 * Array indices follow the fixed layout in SlotLayout.hpp, so slot
 * `Constants::ARR_SLOT_CATHODE` is the cathode in every array here.
 */
struct EventState {
  Int_t leftdE[18];  ///< Left-side energy per anode strip, in ADC.
  Int_t rightdE[18]; ///< Right-side energy per anode strip, in ADC.
  Int_t totaldE[18]; ///< Summed energy per anode strip, in ADC.
  Int_t
      hits[Constants::N_ARR_SLOTS]; ///< Hit count per slot; see SlotLayout.hpp.
  Int_t cathode;                    ///< Cathode energy, in ADC.
  Int_t grid;                       ///< Frisch grid energy, in ADC.
  /// Bitwise OR of the CoMPASS status flags of every hit in the event, so a
  /// single test covers pileup or saturation anywhere in it.
  UInt_t flags_or;
  Bool_t had_cathode; ///< Whether a cathode hit was seen at all.
  /// Timestamp of the reference grid hit that seeded this event. Unique per
  /// event and stable across re-chunking, so it joins a cached event back to
  /// its source without depending on file or entry numbering.
  ULong64_t ref_ts;
};

/**
 * @brief Raw per-slot detail for one event, kept alongside EventState.
 *
 * Where EventState holds the reduced quantities, this keeps what each slot
 * actually contributed, for deduplication decisions and diagnostics.
 */
struct PerChannelData {
  ULong64_t timestamps[Constants::N_ARR_SLOTS]; ///< Hit timestamp per slot.
  UShort_t energies[Constants::N_ARR_SLOTS];    ///< Hit energy per slot, ADC.
  UInt_t flags[Constants::N_ARR_SLOTS];         ///< Status flags per slot.
};

/**
 * @brief Groups time-sorted hits into events.
 *
 * A grid hit seeds an event; hits within the coincidence window are assigned to
 * its slots. A second hit on an already-filled slot is resolved by the
 * configured DedupStrategy.
 *
 * All-static: the builder holds no state of its own, so the same functions
 * serve the parallel file workers without synchronisation.
 */
class EventBuilder {
public:
  /// @brief Board/channel to array slot lookup, indexed as BuildSlotMap()
  /// defines.
  typedef std::vector<Int_t> SlotMap;

  /// @brief Clear an event back to its empty state, ready to be seeded.
  /// @param[out] e Event to reset.
  static void ResetEventState(EventState &e);
  /// @brief Clear the per-slot detail alongside a reset event.
  /// @param[out] p Data to reset.
  static void ResetPerChannelData(PerChannelData &p);
  /**
   * @brief Decide which of two hits on the same slot survives.
   *
   * @param cand_ts     Candidate hit's timestamp.
   * @param prev_ts     Timestamp of the hit already assigned.
   * @param cand_energy Candidate hit's energy, in ADC.
   * @param prev_energy Energy of the hit already assigned.
   * @param ref_ts      Timestamp of the grid hit that seeded the event, which
   *                    DedupStrategy::kCLOSEST_TO_FIRST_GRID measures against.
   * @param strategy    Resolution rule.
   *
   * @return `kTRUE` to replace the existing hit with the candidate.
   *
   * @note DedupStrategy::kDISCARD is not decided here — it rejects the whole
   *       event, which is handled by the caller.
   */
  static Bool_t ShouldKeepHit(ULong64_t cand_ts, ULong64_t prev_ts,
                              UShort_t cand_energy, UShort_t prev_energy,
                              ULong64_t ref_ts, DedupStrategy strategy);
  /**
   * @brief Build the board/channel to slot lookup for this dataset.
   * @return The map, sized and ordered so a hit's board and channel index into
   *         it directly.
   */
  static SlotMap BuildSlotMap();
  /**
   * @brief Place one hit into an event, resolving any collision.
   *
   * @param[in,out] e  Event to fill.
   * @param[in,out] pc Per-slot detail, or null to skip recording it.
   * @param ref_ts     Seeding grid hit's timestamp.
   * @param slot       Target slot; see SlotLayout.hpp.
   * @param energy     Hit energy, in ADC.
   * @param timestamp  Hit timestamp.
   * @param flags      Hit status flags, OR-ed into EventState::flags_or.
   * @param strategy   Rule for resolving a slot that is already filled.
   */
  static void AssignHit(EventState &e, PerChannelData *pc, ULong64_t ref_ts,
                        Int_t slot, UShort_t energy, ULong64_t timestamp,
                        UInt_t flags, DedupStrategy strategy);
  /**
   * @brief Whether an event has everything required to be written out.
   * @param e Event to test.
   * @return `kTRUE` if it is complete by the dataset's criteria.
   */
  static Bool_t CheckEventComplete(const EventState &e);
  /**
   * @brief Unpack the OR-ed status flags into the conditions that matter.
   * @param e Event to inspect.
   * @param[out] has_fake       A synthetic event was inserted by the
   * acquisition.
   * @param[out] has_saturation Some channel saturated within its gate.
   * @param[out] has_pileup     Pileup was flagged on some channel.
   */
  static void GetFlagSummary(const EventState &e, Bool_t &has_fake,
                             Bool_t &has_saturation, Bool_t &has_pileup);
  /**
   * @brief Build every event in a subfile and write them to a ROOT tree.
   *
   * @param hits        Hits for one subfile, **already sorted by timestamp**.
   * @param slot_map    Lookup from BuildSlotMap().
   * @param output_name Output basename for the events file.
   * @param file_label  Label from FileSet::FileLabel(), used in logs and plots.
   *
   * @return `kTRUE` on success.
   *
   * @warning The hits must be time-sorted. The builder walks them in order and
   *          closes an event once the coincidence window passes, so unsorted
   *          input silently produces wrong events rather than an error.
   */
  static Bool_t BuildEventsFromSortedHits(const std::vector<RawHit> &hits,
                                          const SlotMap &slot_map,
                                          const TString &output_name,
                                          const TString &file_label);
};

#endif
