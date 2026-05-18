#ifndef TIMING_HPP
#define TIMING_HPP

#include "BinaryUtils.hpp"
#include "Constants.hpp"
#include "Pipeline.hpp"
#include "PlottingUtils.hpp"
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

inline void ComputeNSD2(const std::vector<Double_t> &ref_x,
                        const std::vector<Double_t> &ref_y,
                        const std::vector<Double_t> &gr_x,
                        const std::vector<Double_t> &gr_y, Double_t shift,
                        Double_t thresh_dt_us, Int_t &npts, Double_t &nsd2) {
  nsd2 = 0;
  npts = 0;

  Double_t tmin_gr = gr_x.empty() ? 0 : gr_x.front() - shift;
  Double_t tmax_gr = gr_x.empty() ? 0 : gr_x.back() - shift;
  Double_t tmin_ref = ref_x.empty() ? 0 : ref_x.front();
  Double_t tmax_ref = ref_x.empty() ? 0 : ref_x.back();

  Double_t tmin = TMath::Max(tmin_ref, tmin_gr);
  Double_t tmax = TMath::Min(tmax_ref, tmax_gr);

  Int_t gr_idx = 0;

  for (Int_t p = 0; p < Int_t(ref_x.size()); p++) {
    Double_t tref = ref_x[p];
    Double_t dtref = ref_y[p];

    if (tref < tmin || tref > tmax || dtref <= thresh_dt_us)
      continue;

    Double_t tref_shifted = tref + shift;

    while (gr_idx < Int_t(gr_x.size()) - 1 && gr_x[gr_idx + 1] < tref_shifted)
      gr_idx++;

    if (gr_idx >= Int_t(gr_x.size()) - 1)
      continue;

    Double_t x0 = gr_x[gr_idx];
    Double_t x1 = gr_x[gr_idx + 1];
    Double_t y0 = gr_y[gr_idx];
    Double_t y1 = gr_y[gr_idx + 1];

    Double_t dt = y0 + (y1 - y0) * (tref_shifted - x0) / (x1 - x0);
    nsd2 += pow(dt - dtref, 2);
    npts++;
  }

  if (nsd2 > 0 && npts > 0) {
    nsd2 = sqrt(nsd2) / npts;
  } else {
    nsd2 = 1e12;
  }
}

struct LongChan {
  UShort_t board;
  UShort_t channel;
  TString name;
};

inline Bool_t IsLongChannel(const TString &name) {
  if (name == "Strip0" || name == "Strip17")
    return kTRUE;
  if (name.Length() < 2)
    return kFALSE;
  Char_t side = name[0];
  if (side != 'L' && side != 'R')
    return kFALSE;
  TString num = name;
  num.Remove(0, 1);
  if (!num.IsDigit())
    return kFALSE;
  Int_t n = num.Atoi();
  if (side == 'L' && (n % 2 == 1))
    return kTRUE;
  if (side == 'R' && (n % 2 == 0))
    return kTRUE;
  return kFALSE;
}

inline Bool_t LongChanOrder(const LongChan &a, const LongChan &b) {
  if (a.board != b.board)
    return a.board < b.board;
  return a.channel < b.channel;
}

inline std::vector<LongChan> BuildLongChannelList() {
  std::vector<LongChan> list;
  std::map<std::pair<Int_t, Int_t>, TString>::const_iterator it;
  for (it = Constants::channelMap.begin(); it != Constants::channelMap.end();
       ++it) {
    if (!IsLongChannel(it->second))
      continue;
    LongChan lc;
    lc.board = static_cast<UShort_t>(it->first.first);
    lc.channel = static_cast<UShort_t>(it->first.second);
    lc.name = it->second;
    list.push_back(lc);
  }
  std::sort(list.begin(), list.end(), LongChanOrder);
  return list;
}

inline void PlotExtremeEvents2D(TH2F *h_before, TH2F *h_after,
                                const TString &file_label) {
  std::lock_guard<std::mutex> lock(g_plot_mutex);
  TString subdir = "timing/" + file_label;

  TCanvas *c_before = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureAndDraw2DHistogram(h_before, c_before);
  c_before->SetLogz(kFALSE);
  PlottingUtils::SaveFigure(c_before, "extreme_events_before", subdir,
                            PlotSaveOptions::kLINEAR);
  delete c_before;

  TCanvas *c_after = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureAndDraw2DHistogram(h_after, c_after);
  c_after->SetLogz(kFALSE);
  PlottingUtils::SaveFigure(c_after, "extreme_events_after", subdir,
                            PlotSaveOptions::kLINEAR);
  delete c_after;
}

inline void PlotCostLandscape(const std::vector<Double_t> &shifts,
                              const std::vector<Double_t> &inv_nsd2_values,
                              Double_t best_shift, UShort_t ref_board,
                              UShort_t board, const TString &file_label) {
  std::lock_guard<std::mutex> lock(g_plot_mutex);

  TGraph *g = new TGraph();
  Double_t y_max = 0;
  for (Int_t i = 0; i < Int_t(shifts.size()); i++) {
    g->SetPoint(g->GetN(), shifts[i], inv_nsd2_values[i]);
    if (inv_nsd2_values[i] > y_max)
      y_max = inv_nsd2_values[i];
  }

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureGraph(
      g, kBlue + 1,
      Form("Cost landscape Board %d-%d;Candidate shift [s];1/NSD^{2}",
           ref_board, board));
  g->SetLineWidth(PlottingUtils::GetLineWidth());
  g->Draw("AL");

  TLine *best = new TLine(best_shift, 0, best_shift, 1.05 * y_max);
  best->SetLineColor(kRed + 1);
  best->SetLineStyle(2);
  best->SetLineWidth(PlottingUtils::GetLineWidth());
  best->Draw();

  PlottingUtils::AddText(Form("best shift = %.6f s", best_shift), 0.85, 0.85);

  PlottingUtils::SaveFigure(canvas, Form("cost_landscape_board_%d", board),
                            "timing/" + file_label, PlotSaveOptions::kLINEAR);
  delete canvas;
}

inline Double_t FindShiftBeam(TGraph *ref, TGraph *gr, Double_t overlap_tmin_s,
                              Double_t overlap_tmax_s, Double_t thresh_dt_us,
                              UShort_t ref_board, UShort_t board,
                              const TString &file_label) {

  std::vector<Double_t> ref_x(ref->GetN()), ref_y(ref->GetN());
  std::vector<Double_t> gr_x(gr->GetN()), gr_y(gr->GetN());

  std::cout << "Caching graph data..." << std::endl;
  for (Int_t i = 0; i < ref->GetN(); i++) {
    ref->GetPoint(i, ref_x[i], ref_y[i]);
  }
  for (Int_t i = 0; i < gr->GetN(); i++) {
    gr->GetPoint(i, gr_x[i], gr_y[i]);
  }

  Double_t tini = 0;
  for (Int_t p = 0; p < Int_t(ref_x.size()); p++) {
    if (ref_x[p] >= overlap_tmin_s) {
      tini = ref_x[p];
      break;
    }
  }

  std::cout << "Looking for timeshift between " << overlap_tmin_s << " - "
            << overlap_tmax_s << " sec" << std::endl;
  std::cout << "tini = " << tini << " s" << std::endl;
  std::cout << "Scanning " << gr->GetN() << " candidate shifts..." << std::endl;

  Double_t best_shift = 0;
  Double_t max_inv_nsd2 = 0;

  std::vector<Double_t> scan_shifts;
  std::vector<Double_t> scan_inv_nsd2;
  scan_shifts.reserve(gr->GetN());
  scan_inv_nsd2.reserve(gr->GetN());

  for (Int_t p = 0; p < gr->GetN(); p++) {
    Int_t npts = 0;
    Double_t nsd2 = 0;
    Double_t shift = gr_x[p] - tini;

    if (TMath::Abs(shift) > Constants::TIMING_MAX_ABS_SHIFT_S)
      continue;

    if (p % 1000 == 0)
      std::cout << "  Testing shift: " << shift << " s (" << p << "/"
                << gr->GetN() << ")" << std::endl;

    ComputeNSD2(ref_x, ref_y, gr_x, gr_y, shift, thresh_dt_us, npts, nsd2);

    if (npts > 10) {
      Double_t inv_nsd2 = 1.0 / nsd2;
      scan_shifts.push_back(shift);
      scan_inv_nsd2.push_back(inv_nsd2);

      if (inv_nsd2 > max_inv_nsd2) {
        best_shift = shift;
        max_inv_nsd2 = inv_nsd2;
      }
    }
  }

  Int_t npts0 = 0;
  Double_t nsd2_0 = 0;
  ComputeNSD2(ref_x, ref_y, gr_x, gr_y, 0.0, thresh_dt_us, npts0, nsd2_0);
  Double_t inv_nsd2_zero = (npts0 > 10) ? 1.0 / nsd2_0 : 0.0;

  std::cout << "Best shift: " << best_shift << " s (1/NSD2 = " << max_inv_nsd2
            << ")" << std::endl;
  std::cout << "No-shift   : 0 s         (1/NSD2 = " << inv_nsd2_zero
            << ", npts=" << npts0 << ")" << std::endl;
  if (inv_nsd2_zero > 0)
    std::cout << "  Improvement vs no-shift: " << max_inv_nsd2 / inv_nsd2_zero
              << "x" << std::endl;

  PlotCostLandscape(scan_shifts, scan_inv_nsd2, best_shift, ref_board, board,
                    file_label);

  return best_shift;
}

std::vector<TGraph *> ExtractAllChannelsTimingStructureFromHits(
    const std::vector<RawHit> &hits, const std::vector<LongChan> &channels,
    Double_t min_energy, Double_t max_energy, Double_t tmin_s, Double_t tmax_s,
    Double_t thresh_dt_us);

TimeShiftResult CalcTimeShiftsBeamMethodFromHits(
    const std::vector<RawHit> &hits, const TString &file_label,
    UShort_t ref_board, const std::vector<UShort_t> &board_channels,
    Double_t min_energy, Double_t max_energy, Double_t overlap_margin_s,
    Double_t thresh_dt_us);

void ApplyShiftsInPlace(std::vector<RawHit> &hits,
                        const std::vector<Long64_t> &board_shifts);

void SortHitsByTimestamp(std::vector<RawHit> &hits);

#endif
