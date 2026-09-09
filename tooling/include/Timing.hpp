#ifndef TIMING_HPP
#define TIMING_HPP

#include "BinaryUtils.hpp"
#include "ChannelTiming.hpp"
#include "Constants.hpp"
#include "FileSet.hpp"
#include "GpuAccel.hpp"
#include "PlottingUtils.hpp"
#include <Rtypes.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TLine.h>
#include <TMath.h>
#include <TString.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <mutex>
#include <utility>
#include <vector>

/**
 * @file Timing.hpp
 * @brief Aligning the clocks of independent digitiser boards.
 *
 * Boards free-run on their own clocks, so hits from one board carry a constant
 * offset against another's. Events cannot be built across boards until that
 * offset is measured and removed.
 *
 * The beam method measures it from the data itself: the beam's own time
 * structure appears on every board, so the shift that best superimposes one
 * board's structure on the reference board's is the offset between them.
 */

/// @brief Per-board timing offsets, indexed by board number.
struct TimeShiftResult {
  /// Offset to add to each board's timestamps, in picoseconds. The reference
  /// board's entry is zero by construction.
  std::vector<Long64_t> board_shifts;
};

/// @brief One long-end channel, used as a timing reference.
struct LongChan {
  UShort_t board;   ///< Board id.
  UShort_t channel; ///< Channel on that board.
  TString name;     ///< Channel name from the active map.
};

/**
 * @brief Multi-board timing alignment and hit sorting.
 *
 * All-static; the parallel file workers share these without synchronisation.
 */
class Timing {
public:
  /**
   * @brief Start of the busiest window of a given width.
   *
   * Beam is not delivered uniformly across a subfile, and the alignment wants
   * the stretch with the most hits rather than an arbitrary one.
   *
   * @param times           Hit times in seconds. Taken by value and sorted.
   * @param window_width_s  Window width, in seconds.
   * @param fallback_start_s Returned when @p times is empty or the width is
   *                         not positive.
   * @return The window start, in seconds.
   */
  static Double_t FindDensestTimeWindowStartS(std::vector<Double_t> times,
                                              Double_t window_width_s,
                                              Double_t fallback_start_s);
  /**
   * @brief Normalised squared deviation between two timing structures at a
   * shift.
   *
   * The cost function the shift search minimises: how badly a candidate shift
   * superimposes one board's beam structure on the reference's.
   *
   * @param ref_x       Reference times.
   * @param ref_y       Reference values at those times.
   * @param gr_x        Candidate board's times.
   * @param gr_y        Candidate board's values.
   * @param shift       Trial shift applied to the candidate, in seconds.
   * @param thresh_dt_us Pairing tolerance, in microseconds.
   * @param[out] npts   Points that paired up under @p thresh_dt_us.
   * @param[out] nsd2   The normalised squared deviation; smaller is better.
   *
   * @note A shift is only trusted when @p npts exceeds
   *       `Constants::cfg.TIMING_SHIFT_MIN_NPTS`; too few pairs make a small
   *       deviation meaningless.
   */
  static void ComputeNSD2(const std::vector<Double_t> &ref_x,
                          const std::vector<Double_t> &ref_y,
                          const std::vector<Double_t> &gr_x,
                          const std::vector<Double_t> &gr_y, Double_t shift,
                          Double_t thresh_dt_us, Int_t &npts, Double_t &nsd2);

  /// @brief Whether a channel name denotes a long end.
  /// @param name Channel name from the active map.
  static Bool_t IsLongChannel(const TString &name);
  /// @brief Ordering predicate giving a stable channel sequence.
  /// @param a First channel.
  /// @param b Second channel.
  /// @return `kTRUE` when @p a sorts before @p b.
  static Bool_t LongChanOrder(const LongChan &a, const LongChan &b);
  /// @brief Every long-end channel in the active map, in stable order.
  static std::vector<LongChan> BuildLongChannelList();

  /**
   * @brief Save the before-and-after alignment diagnostic figures.
   *
   * Full-range and zoomed views of the same data, so both the gross alignment
   * and the residual structure are visible.
   *
   * @param h_before          Timing structure before alignment.
   * @param h_after           And after.
   * @param h_before_zoom     Zoomed view before alignment.
   * @param h_after_zoom      And after.
   * @param file_label        Subfile label, used in the plot names.
   * @param before_zoom_t0_s  Zoom window start before alignment, in seconds.
   * @param after_zoom_t0_s   And after.
   */
  static void PlotExtremeEvents2D(TH2F *h_before, TH2F *h_after,
                                  TH2F *h_before_zoom, TH2F *h_after_zoom,
                                  const TString &file_label,
                                  Double_t before_zoom_t0_s,
                                  Double_t after_zoom_t0_s);
  /**
   * @brief Save the shift-search cost landscape for one board pair.
   *
   * Shows whether the chosen shift sits at a clear maximum or on a plateau,
   * which is the difference between a measurement and a coincidence.
   *
   * @param shifts          Trial shifts searched, in seconds.
   * @param inv_nsd2_values Inverse deviation at each, so higher is better.
   * @param best_shift      The shift selected.
   * @param ref_board       Reference board id.
   * @param board           Board being aligned.
   * @param file_label      Subfile label.
   * @param tag             Extra tag distinguishing passes.
   */
  static void PlotCostLandscape(const std::vector<Double_t> &shifts,
                                const std::vector<Double_t> &inv_nsd2_values,
                                Double_t best_shift, UShort_t ref_board,
                                UShort_t board, const TString &file_label,
                                const TString &tag);

  /**
   * @brief Best shift aligning one board's beam structure to the reference's.
   *
   * Scans candidate shifts and keeps the one maximising the inverse normalised
   * squared deviation, over the overlap window where both boards have data.
   *
   * @param ref            Reference board's timing structure.
   * @param gr             Board being aligned.
   * @param overlap_tmin_s Overlap window start, in seconds.
   * @param overlap_tmax_s Overlap window end.
   * @param thresh_dt_us   Pairing tolerance, in microseconds.
   * @param ref_board      Reference board id, for the diagnostic plot.
   * @param board          Board being aligned, likewise.
   * @param file_label     Subfile label.
   *
   * @return The best shift, in seconds.
   */
  static Double_t FindShiftBeam(TGraph *ref, TGraph *gr,
                                Double_t overlap_tmin_s,
                                Double_t overlap_tmax_s, Double_t thresh_dt_us,
                                UShort_t ref_board, UShort_t board,
                                const TString &file_label);

  /**
   * @brief Build the beam timing structure for each channel.
   *
   * @param hits          Raw hits for the subfile.
   * @param channels      Channels to extract, from BuildLongChannelList().
   * @param min_energy    Lower energy gate, selecting beam-like hits.
   * @param max_energy    Upper energy gate.
   * @param tmin_s        Window start, in seconds.
   * @param tmax_s        Window end.
   * @param thresh_dt_us  Pairing tolerance, in microseconds.
   *
   * @return One graph per channel, in @p channels order. **The caller owns
   *         every graph.**
   */
  static std::vector<TGraph *> ExtractAllChannelsTimingStructureFromHits(
      const std::vector<RawHit> &hits, const std::vector<LongChan> &channels,
      Double_t min_energy, Double_t max_energy, Double_t tmin_s,
      Double_t tmax_s, Double_t thresh_dt_us);

  /**
   * @brief Measure every board's offset against a reference board.
   *
   * @param hits             Raw hits for the subfile.
   * @param file_label       Subfile label, for diagnostics.
   * @param ref_board        Board every other is aligned to.
   * @param board_channels   Reference channel per board.
   * @param min_energy       Lower energy gate for beam-like hits.
   * @param max_energy       Upper energy gate.
   * @param overlap_margin_s Margin trimmed from the overlap window, in seconds.
   * @param thresh_dt_us     Pairing tolerance, in microseconds.
   *
   * @return The per-board shifts, with the reference board at zero.
   */
  static TimeShiftResult CalcTimeShiftsBeamMethodFromHits(
      const std::vector<RawHit> &hits, const TString &file_label,
      UShort_t ref_board, const std::vector<UShort_t> &board_channels,
      Double_t min_energy, Double_t max_energy, Double_t overlap_margin_s,
      Double_t thresh_dt_us);

  /**
   * @brief Add the measured offsets to the hits, in place.
   * @param[in,out] hits  Hits to shift.
   * @param board_shifts  Offsets from CalcTimeShiftsBeamMethodFromHits().
   * @note Leaves the hits out of time order; sort afterwards with
   *       SortHitsByTimestamp() before building events.
   */
  static void ApplyShiftsInPlace(std::vector<RawHit> &hits,
                                 const std::vector<Long64_t> &board_shifts);

  /**
   * @brief Sort hits into ascending timestamp order, in place.
   * @param[in,out] hits Hits to sort.
   * @note Dispatches to the CUDA kernel when GpuAccel reports it available and
   *       a sort slot is free; otherwise sorts on the host. The result is
   *       identical either way.
   */
  static void SortHitsByTimestamp(std::vector<RawHit> &hits);
};

#endif
