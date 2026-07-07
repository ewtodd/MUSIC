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

struct TimeShiftResult {
  std::vector<Long64_t> board_shifts;
};

struct LongChan {
  UShort_t board;
  UShort_t channel;
  TString name;
};

class Timing {
public:
  static Double_t FindDensestTimeWindowStartS(std::vector<Double_t> times,
                                              Double_t window_width_s,
                                              Double_t fallback_start_s);
  static void ComputeNSD2(const std::vector<Double_t> &ref_x,
                          const std::vector<Double_t> &ref_y,
                          const std::vector<Double_t> &gr_x,
                          const std::vector<Double_t> &gr_y, Double_t shift,
                          Double_t thresh_dt_us, Int_t &npts, Double_t &nsd2);

  static Bool_t IsLongChannel(const TString &name);
  static Bool_t LongChanOrder(const LongChan &a, const LongChan &b);
  static std::vector<LongChan> BuildLongChannelList();

  static void PlotExtremeEvents2D(TH2F *h_before, TH2F *h_after,
                                  TH2F *h_before_zoom, TH2F *h_after_zoom,
                                  const TString &file_label,
                                  Double_t before_zoom_t0_s,
                                  Double_t after_zoom_t0_s);
  static void PlotCostLandscape(const std::vector<Double_t> &shifts,
                                const std::vector<Double_t> &inv_nsd2_values,
                                Double_t best_shift, UShort_t ref_board,
                                UShort_t board, const TString &file_label,
                                const TString &tag);

  static Double_t FindShiftBeam(TGraph *ref, TGraph *gr,
                                Double_t overlap_tmin_s,
                                Double_t overlap_tmax_s, Double_t thresh_dt_us,
                                UShort_t ref_board, UShort_t board,
                                const TString &file_label);

  static std::vector<TGraph *> ExtractAllChannelsTimingStructureFromHits(
      const std::vector<RawHit> &hits, const std::vector<LongChan> &channels,
      Double_t min_energy, Double_t max_energy, Double_t tmin_s,
      Double_t tmax_s, Double_t thresh_dt_us);

  static TimeShiftResult CalcTimeShiftsBeamMethodFromHits(
      const std::vector<RawHit> &hits, const TString &file_label,
      UShort_t ref_board, const std::vector<UShort_t> &board_channels,
      Double_t min_energy, Double_t max_energy, Double_t overlap_margin_s,
      Double_t thresh_dt_us);

  static void ApplyShiftsInPlace(std::vector<RawHit> &hits,
                                 const std::vector<Long64_t> &board_shifts);

  static void SortHitsByTimestamp(std::vector<RawHit> &hits);
};

#endif
