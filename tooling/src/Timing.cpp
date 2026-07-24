#include "Timing.hpp"
#include <algorithm>
#include <cmath>

struct ShiftScanResult {
  Double_t best_shift = 0.0;
  Double_t best_inv_nsd2 = 0.0;
  Int_t best_npts = 0;
  Bool_t found = kFALSE;
  std::vector<Double_t> shifts;
  std::vector<Double_t> inv_nsd2_values;
};

ShiftScanResult ScanShiftRange(const std::vector<Double_t> &ref_x,
                               const std::vector<Double_t> &ref_y,
                               const std::vector<Double_t> &gr_x,
                               const std::vector<Double_t> &gr_y,
                               Double_t shift_min_s, Double_t shift_max_s,
                               Double_t shift_step_s, Double_t thresh_dt_us,
                               const TString &label) {

  ShiftScanResult result;

  if (shift_step_s <= 0) {
    std::cerr << label << ": shift_step_s <= 0 (" << shift_step_s
              << "), skipping scan." << std::endl;
    return result;
  }

  if (shift_min_s > shift_max_s) {
    std::cerr << label << ": shift_min_s > shift_max_s (" << shift_min_s
              << " > " << shift_max_s << "), skipping scan." << std::endl;
    return result;
  }

  Int_t candidate_count = static_cast<Int_t>(std::floor(
                              (shift_max_s - shift_min_s) / shift_step_s)) +
                          1;

  if (candidate_count > Constants::cfg.TIMING_SHIFT_MAX_SCAN_CANDIDATES) {
    std::cerr << label << ": candidate count " << candidate_count
              << " exceeds limit "
              << Constants::cfg.TIMING_SHIFT_MAX_SCAN_CANDIDATES
              << ", skipping scan." << std::endl;
    return result;
  }

  std::cout << "  " << label << " range: [" << shift_min_s << ", "
            << shift_max_s << "] s ([" << shift_min_s * 1e6 << ", "
            << shift_max_s * 1e6 << "] us)" << std::endl;
  std::cout << "  " << label << " step: " << shift_step_s << " s ("
            << shift_step_s * 1e6 << " us)" << std::endl;
  std::cout << "  " << label << " candidates: " << candidate_count << std::endl;

  result.shifts.reserve(candidate_count);
  result.inv_nsd2_values.reserve(candidate_count);

  for (Int_t i = 0; i < candidate_count; i++) {
    Double_t shift = shift_min_s + i * shift_step_s;
    Int_t npts = 0;
    Double_t nsd2 = 0;

    Timing::ComputeNSD2(ref_x, ref_y, gr_x, gr_y, shift, thresh_dt_us, npts,
                        nsd2);

    if (npts > Constants::cfg.TIMING_SHIFT_MIN_NPTS) {
      Double_t inv_nsd2 = 1.0 / nsd2;
      result.shifts.push_back(shift);
      result.inv_nsd2_values.push_back(inv_nsd2);

      if (!result.found || inv_nsd2 > result.best_inv_nsd2) {
        result.best_shift = shift;
        result.best_inv_nsd2 = inv_nsd2;
        result.best_npts = npts;
        result.found = kTRUE;
      }
    }

    if (i % 50000 == 0 && i > 0) {
      std::cout << "    Progress: " << i << "/" << candidate_count << std::endl;
    }
  }

  return result;
}

void Timing::ComputeNSD2(const std::vector<Double_t> &ref_x,
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

Bool_t Timing::IsLongChannel(const TString &name) {
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

Bool_t Timing::LongChanOrder(const LongChan &a, const LongChan &b) {
  if (a.board != b.board)
    return a.board < b.board;
  return a.channel < b.channel;
}

std::vector<LongChan> Timing::BuildLongChannelList() {
  std::vector<LongChan> list;
  std::map<std::pair<Int_t, Int_t>, TString>::const_iterator it;
  for (it = Constants::ActiveChannelMap().begin();
       it != Constants::ActiveChannelMap().end(); ++it) {
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

void Timing::PlotExtremeEvents2D(TH2F *h_before, TH2F *h_after,
                                 TH2F *h_before_zoom, TH2F *h_after_zoom,
                                 const TString &file_label,
                                 Double_t before_zoom_t0_s,
                                 Double_t after_zoom_t0_s) {
  std::lock_guard<std::mutex> lock(g_plot_mutex);
  TString subdir = "timing/" + file_label;

  TCanvas *c_before = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureAndDraw2DHistogram(h_before, c_before);
  c_before->SetLogz(kFALSE);
  if (Constants::cfg.SAVE_PLOTS)
    PlottingUtils::SaveFigure(c_before, "extreme_events_before", subdir,
                              PlotSaveOptions::kLINEAR);
  delete c_before;

  TCanvas *c_after = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureAndDraw2DHistogram(h_after, c_after);
  c_after->SetLogz(kFALSE);
  if (Constants::cfg.SAVE_PLOTS)
    PlottingUtils::SaveFigure(c_after, "extreme_events_after", subdir,
                              PlotSaveOptions::kLINEAR);
  delete c_after;

  Int_t n = gStyle->GetNumberOfColors();
  std::vector<Int_t> old_palette;
  old_palette.reserve(n);

  for (Int_t i = 0; i < n; i++) {
    old_palette.push_back(gStyle->GetColorPalette(i));
  }

  Int_t palette[2];
  palette[0] = TColor::GetColor("#FFFFFF"); // background / zero
  palette[1] = TColor::GetColor("#D62728"); // occupied bin / event

  if (h_before_zoom) {
    TCanvas *c_before_zoom = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::Configure2DHistogram(h_before_zoom, c_before_zoom);

    gStyle->SetPalette(2, palette);
    h_before_zoom->SetMinimum(0);
    h_before_zoom->SetMaximum(1);
    c_before_zoom->SetLogz(kFALSE);
    h_before_zoom->Draw();
    c_before_zoom->SetRightMargin(0.07);
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(c_before_zoom, "extreme_events_before_zoom_us",
                                subdir, PlotSaveOptions::kLINEAR);
    delete c_before_zoom;
  }

  if (h_after_zoom) {
    TCanvas *c_after_zoom = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::Configure2DHistogram(h_after_zoom, c_after_zoom);

    gStyle->SetPalette(2, palette);
    h_after_zoom->SetMinimum(0);
    h_after_zoom->SetMaximum(1);
    c_after_zoom->SetLogz(kFALSE);
    h_after_zoom->Draw();
    c_after_zoom->SetRightMargin(0.07);
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(c_after_zoom, "extreme_events_after_zoom_us",
                                subdir, PlotSaveOptions::kLINEAR);
    delete c_after_zoom;
  }
  if (!old_palette.empty()) {
    gStyle->SetPalette(static_cast<Int_t>(old_palette.size()),
                       old_palette.data());
  }
}

void Timing::PlotCostLandscape(const std::vector<Double_t> &shifts,
                               const std::vector<Double_t> &inv_nsd2_values,
                               Double_t best_shift, UShort_t ref_board,
                               UShort_t board, const TString &file_label,
                               const TString &tag) {
  std::lock_guard<std::mutex> lock(g_plot_mutex);

  std::vector<std::pair<Double_t, Double_t>> points;
  points.reserve(shifts.size());
  for (Int_t i = 0; i < Int_t(shifts.size()); i++) {
    points.push_back(std::make_pair(shifts[i], inv_nsd2_values[i]));
  }
  std::sort(points.begin(), points.end());

  TGraph *g = new TGraph();
  Double_t y_max = 0;
  for (Int_t i = 0; i < Int_t(points.size()); i++) {
    g->SetPoint(g->GetN(), points[i].first, points[i].second);
    if (points[i].second > y_max)
      y_max = points[i].second;
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

  if (Constants::cfg.SAVE_PLOTS)
    PlottingUtils::SaveFigure(
        canvas, Form("cost_landscape_board_%d%s", board, tag.Data()),
        "timing/" + file_label, PlotSaveOptions::kLINEAR);
  delete canvas;
}

Double_t Timing::FindShiftBeam(TGraph *ref, TGraph *gr, Double_t overlap_tmin_s,
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

  const Double_t coarse_step_s =
      Constants::cfg.TIMING_SHIFT_COARSE_STEP_US * 1e-6;
  const Double_t fine_step_s = Constants::cfg.TIMING_SHIFT_FINE_STEP_US * 1e-6;
  const Double_t fine_half_width_s =
      Constants::cfg.TIMING_SHIFT_FINE_HALF_WIDTH_US * 1e-6;

  const Double_t scan_min_s = -Constants::cfg.TIMING_MAX_ABS_SHIFT_S;
  const Double_t scan_max_s = Constants::cfg.TIMING_MAX_ABS_SHIFT_S;

  std::cout << "Looking for timeshift between " << overlap_tmin_s << " - "
            << overlap_tmax_s << " sec" << std::endl;
  std::cout << "Max absolute shift range: [" << scan_min_s << ", " << scan_max_s
            << "] s ([" << scan_min_s * 1e6 << ", " << scan_max_s * 1e6
            << "] us)" << std::endl;
  std::cout << "Coarse step: " << coarse_step_s << " s ("
            << Constants::cfg.TIMING_SHIFT_COARSE_STEP_US << " us)"
            << std::endl;
  std::cout << "Fine step: " << fine_step_s << " s ("
            << Constants::cfg.TIMING_SHIFT_FINE_STEP_US << " us)" << std::endl;
  std::cout << "Fine half-width: " << fine_half_width_s << " s ("
            << Constants::cfg.TIMING_SHIFT_FINE_HALF_WIDTH_US << " us)"
            << std::endl;

  ShiftScanResult coarse =
      ScanShiftRange(ref_x, ref_y, gr_x, gr_y, scan_min_s, scan_max_s,
                     coarse_step_s, thresh_dt_us, "Coarse");

  Double_t final_shift = 0.0;
  Double_t final_inv_nsd2 = 0.0;

  ShiftScanResult fine;

  if (coarse.found) {
    std::cout << "Coarse best: shift = " << coarse.best_shift << " s ("
              << coarse.best_shift * 1e6
              << " us), 1/NSD2 = " << coarse.best_inv_nsd2
              << ", npts = " << coarse.best_npts << std::endl;

    Double_t fine_min_s = coarse.best_shift - fine_half_width_s;
    Double_t fine_max_s = coarse.best_shift + fine_half_width_s;

    if (fine_min_s < scan_min_s)
      fine_min_s = scan_min_s;
    if (fine_max_s > scan_max_s)
      fine_max_s = scan_max_s;

    fine = ScanShiftRange(ref_x, ref_y, gr_x, gr_y, fine_min_s, fine_max_s,
                          fine_step_s, thresh_dt_us, "Fine");

    if (fine.found) {
      final_shift = fine.best_shift;
      final_inv_nsd2 = fine.best_inv_nsd2;
      std::cout << "Fine best: shift = " << fine.best_shift << " s ("
                << fine.best_shift * 1e6
                << " us), 1/NSD2 = " << fine.best_inv_nsd2
                << ", npts = " << fine.best_npts << std::endl;
    } else {
      final_shift = coarse.best_shift;
      final_inv_nsd2 = coarse.best_inv_nsd2;
      std::cout << "Fine scan found no valid candidate; falling back to coarse "
                   "best."
                << std::endl;
    }
  } else {
    std::cerr << "WARNING: Coarse scan found no valid candidate. Returning 0.0 "
                 "shift."
              << std::endl;
  }

  Int_t npts0 = 0;
  Double_t nsd2_0 = 0;
  ComputeNSD2(ref_x, ref_y, gr_x, gr_y, 0.0, thresh_dt_us, npts0, nsd2_0);
  Double_t inv_nsd2_zero =
      (npts0 > Constants::cfg.TIMING_SHIFT_MIN_NPTS) ? 1.0 / nsd2_0 : 0.0;

  std::cout << "No-shift   : 0 s         (1/NSD2 = " << inv_nsd2_zero
            << ", npts=" << npts0 << ")" << std::endl;
  if (inv_nsd2_zero > 0 && final_inv_nsd2 > 0)
    std::cout << "  Improvement vs no-shift: " << final_inv_nsd2 / inv_nsd2_zero
              << "x" << std::endl;

  std::cout << "Final selected shift: " << final_shift << " s ("
            << final_shift * 1e6 << " us)" << std::endl;

  if (fine.found) {
    PlotCostLandscape(fine.shifts, fine.inv_nsd2_values, final_shift, ref_board,
                      board, file_label, "_fine");
    if (!coarse.shifts.empty())
      PlotCostLandscape(coarse.shifts, coarse.inv_nsd2_values,
                        coarse.best_shift, ref_board, board, file_label,
                        "_coarse");
  } else if (!coarse.shifts.empty()) {
    PlotCostLandscape(coarse.shifts, coarse.inv_nsd2_values, final_shift,
                      ref_board, board, file_label, "_coarse");
  }

  return final_shift;
}

std::vector<TGraph *> Timing::ExtractAllChannelsTimingStructureFromHits(
    const std::vector<RawHit> &hits, const std::vector<LongChan> &channels,
    Double_t min_energy, Double_t max_energy, Double_t tmin_s, Double_t tmax_s,
    Double_t thresh_dt_us) {

  std::map<std::pair<Int_t, Int_t>, Int_t> chan_to_idx;
  for (Int_t i = 0; i < Int_t(channels.size()); i++) {
    chan_to_idx[std::pair<Int_t, Int_t>(channels[i].board,
                                        channels[i].channel)] = i;
  }

  std::vector<std::vector<ULong64_t>> per_chan_ts(channels.size());
  for (Int_t i = 0; i < Int_t(channels.size()); i++) {
    per_chan_ts[i].reserve(10000);
  }

  Long64_t n_entries = Long64_t(hits.size());
  std::cout << "Single-pass extraction across " << channels.size()
            << " long channels (in-memory, " << n_entries << " hits)..."
            << std::endl;

  for (Long64_t i = 0; i < n_entries; i++) {
    const RawHit &h = hits[i];
    if (h.energy < min_energy || h.energy > max_energy)
      continue;
    Double_t time_s = h.timestamp / 1e12;
    if (time_s < tmin_s || time_s > tmax_s)
      continue;

    std::map<std::pair<Int_t, Int_t>, Int_t>::const_iterator it =
        chan_to_idx.find(std::pair<Int_t, Int_t>(h.board, h.channel));
    if (it == chan_to_idx.end())
      continue;

    per_chan_ts[it->second].push_back(h.timestamp);

    if (i % 10000000 == 0)
      std::cout << "  Progress: " << i << "/" << n_entries << std::endl;
  }

  std::vector<TGraph *> graphs;
  graphs.reserve(channels.size());
  for (Int_t c = 0; c < Int_t(channels.size()); c++) {
    TGraph *g = new TGraph();
    const std::vector<ULong64_t> &ts = per_chan_ts[c];
    for (Int_t i = 0; i + 1 < Int_t(ts.size()); i++) {
      Double_t dt_us = (ts[i + 1] - ts[i]) / 1e6;
      if (dt_us > thresh_dt_us) {
        Double_t time_s = ts[i] / 1e12;
        g->SetPoint(g->GetN(), time_s, dt_us);
      }
    }
    std::cout << "  " << channels[c].name << " (B" << channels[c].board << "C"
              << channels[c].channel << "): " << ts.size() << " hits, "
              << g->GetN() << " extreme events" << std::endl;
    graphs.push_back(g);
  }

  return graphs;
}

TimeShiftResult Timing::CalcTimeShiftsBeamMethodFromHits(
    const std::vector<RawHit> &hits, const TString &file_label,
    UShort_t ref_board, const std::vector<UShort_t> &board_channels,
    Double_t min_energy, Double_t max_energy, Double_t overlap_margin_s,
    Double_t thresh_dt_us) {

  TimeShiftResult result;
  result.board_shifts.assign(Constants::cfg.N_BOARDS, 0);

  // Board sync disabled for this dataset (e.g. 87Rb): nothing to compute, so
  // skip the whole extract/scan/extreme-events pipeline and leave every board
  // at zero shift. Per-channel TTF correction still happens in ApplyShifts.
  if (!Constants::cfg.TIMING_DO_BOARD_SYNC) {
    std::cout << "Board sync disabled for this dataset; skipping timeshift "
                 "calculation (all board shifts = 0)."
              << std::endl;
    return result;
  }

  if (hits.empty()) {
    std::cerr << "CalcTimeShiftsBeamMethodFromHits: empty hits vector"
              << std::endl;
    return result;
  }

  ULong64_t ts_min = hits[0].timestamp;
  ULong64_t ts_max = hits[0].timestamp;
  for (Int_t i = 1; i < Int_t(hits.size()); i++) {
    if (hits[i].timestamp < ts_min)
      ts_min = hits[i].timestamp;
    if (hits[i].timestamp > ts_max)
      ts_max = hits[i].timestamp;
  }
  Double_t file_tmin_s = ts_min / 1e12;
  Double_t file_tmax_s = ts_max / 1e12;
  Double_t overlap_tmin_s = file_tmin_s + overlap_margin_s;
  Double_t overlap_tmax_s = file_tmax_s - overlap_margin_s;
  std::cout << "File span [" << file_tmin_s << ", " << file_tmax_s
            << "] s, using overlap [" << overlap_tmin_s << ", "
            << overlap_tmax_s << "] s" << std::endl;

  std::vector<LongChan> long_channels = BuildLongChannelList();
  std::vector<TGraph *> long_graphs = ExtractAllChannelsTimingStructureFromHits(
      hits, long_channels, min_energy, max_energy, overlap_tmin_s,
      overlap_tmax_s, thresh_dt_us);

  std::map<UShort_t, Int_t> board_to_ref_idx;
  for (Int_t i = 0; i < Int_t(long_channels.size()); i++) {
    if (long_channels[i].channel == board_channels[long_channels[i].board])
      board_to_ref_idx[long_channels[i].board] = i;
  }

  std::map<UShort_t, Int_t>::const_iterator ref_it =
      board_to_ref_idx.find(ref_board);
  if (ref_it == board_to_ref_idx.end()) {
    std::cerr << "Reference board " << ref_board
              << " ref channel is not in long-channel list" << std::endl;
    for (Int_t k = 0; k < Int_t(long_graphs.size()); k++)
      delete long_graphs[k];
    return result;
  }
  TGraph *ref_graph = long_graphs[ref_it->second];

  std::vector<Double_t> board_shifts_s(Constants::cfg.N_BOARDS, 0.0);

  for (UShort_t board = 0; board < Constants::cfg.N_BOARDS; board++) {
    if (board == ref_board)
      continue;

    std::map<UShort_t, Int_t>::const_iterator it = board_to_ref_idx.find(board);
    if (it == board_to_ref_idx.end()) {
      std::cout << "WARNING: Board " << board
                << " ref channel not in long-channel list" << std::endl;
      continue;
    }

    TGraph *board_graph = long_graphs[it->second];
    if (board_graph->GetN() < 10) {
      std::cout << "WARNING: Not enough data points for Board " << board
                << std::endl;
      continue;
    }

    std::cout << "Processing Board " << board << " Channel "
              << board_channels[board] << std::endl;

    Double_t shift_s =
        FindShiftBeam(ref_graph, board_graph,
                      overlap_tmin_s + 0.1 * (overlap_tmax_s - overlap_tmin_s),
                      overlap_tmax_s - 0.1 * (overlap_tmax_s - overlap_tmin_s),
                      thresh_dt_us, ref_board, board, file_label);

    board_shifts_s[board] = shift_s;
    Long64_t shift_ps = static_cast<Long64_t>(shift_s * 1e12);
    result.board_shifts[board] = -shift_ps;

    std::cout << "Board " << ref_board << "-" << board << " shift: " << shift_s
              << " s (" << shift_ps << " ps)" << std::endl;
  }

  Int_t time_bins = TMath::Min(
      500,
      TMath::Max(100, static_cast<Int_t>((file_tmax_s - file_tmin_s) * 50.0)));
  Int_t n_y = static_cast<Int_t>(long_channels.size());

  TH2F *h_extreme_before =
      new TH2F("hExtremeBefore", ";Time [s];", time_bins, file_tmin_s,
               file_tmax_s, n_y, -0.5, n_y - 0.5);

  TH2F *h_extreme_after =
      new TH2F("hExtremeAfter", ";Time [s];", time_bins, file_tmin_s,
               file_tmax_s, n_y, -0.5, n_y - 0.5);

  std::vector<Double_t> extreme_times_before;
  std::vector<Double_t> extreme_times_after;

  extreme_times_before.reserve(10000);
  extreme_times_after.reserve(10000);

  for (Int_t c = 0; c < Int_t(long_channels.size()); c++) {
    TGraph *g = long_graphs[c];
    Double_t shift_s = board_shifts_s[long_channels[c].board];

    for (Int_t k = 0; k < g->GetN(); k++) {
      Double_t x, y;
      g->GetPoint(k, x, y);

      Double_t x_before_s = x;
      Double_t x_after_s = x - shift_s;

      h_extreme_before->Fill(x_before_s, Double_t(c));
      h_extreme_after->Fill(x_after_s, Double_t(c));

      extreme_times_before.push_back(x_before_s);
      extreme_times_after.push_back(x_after_s);
    }
  }

  const Double_t zoom_width_us = 1000.0;
  const Double_t zoom_width_s = zoom_width_us * 1e-6;

  // About 1 us/bin. Increase this if you want sub-us binning.
  Int_t zoom_bins =
      TMath::Min(2000, TMath::Max(100, static_cast<Int_t>(zoom_width_us)));

  // Choose a dense after-correction window, then use the same absolute window
  // for before and after. This makes before/after directly comparable.
  Double_t after_zoom_t0_s = FindDensestTimeWindowStartS(
      extreme_times_after, zoom_width_s, file_tmin_s);

  Double_t before_zoom_t0_s = after_zoom_t0_s;

  std::cout << "Zoom window before: [" << before_zoom_t0_s << ", "
            << before_zoom_t0_s + zoom_width_s << "] s" << std::endl;

  std::cout << "Zoom window after : [" << after_zoom_t0_s << ", "
            << after_zoom_t0_s + zoom_width_s << "] s" << std::endl;

  TH2F *h_extreme_before_zoom =
      new TH2F("hExtremeBeforeZoom", ";Time [#mus];", zoom_bins, 0.0,
               zoom_width_us, n_y, -0.5, n_y - 0.5);

  TH2F *h_extreme_after_zoom =
      new TH2F("hExtremeAfterZoom", ";Time [#mus];", zoom_bins, 0.0,
               zoom_width_us, n_y, -0.5, n_y - 0.5);

  // Restore Y-axis labels for all four plots.
  for (Int_t b = 0; b < n_y; b++) {
    TString label =
        Form("B%d %s", long_channels[b].board, long_channels[b].name.Data());

    h_extreme_before->GetYaxis()->SetBinLabel(b + 1, label);
    h_extreme_after->GetYaxis()->SetBinLabel(b + 1, label);

    h_extreme_before_zoom->GetYaxis()->SetBinLabel(b + 1, label);
    h_extreme_after_zoom->GetYaxis()->SetBinLabel(b + 1, label);
  }

  // Second pass:
  // Fill only the zoom-window histograms, in microseconds relative to t0.
  Int_t n_zoom_before = 0;
  Int_t n_zoom_after = 0;

  for (Int_t c = 0; c < Int_t(long_channels.size()); c++) {
    TGraph *g = long_graphs[c];
    Double_t shift_s = board_shifts_s[long_channels[c].board];

    for (Int_t k = 0; k < g->GetN(); k++) {
      Double_t x, y;
      g->GetPoint(k, x, y);

      Double_t x_before_s = x;
      Double_t x_after_s = x - shift_s;

      Double_t x_before_zoom_us = (x_before_s - before_zoom_t0_s) * 1e6;
      Double_t x_after_zoom_us = (x_after_s - after_zoom_t0_s) * 1e6;

      if (x_before_zoom_us >= 0.0 && x_before_zoom_us <= zoom_width_us) {
        h_extreme_before_zoom->Fill(x_before_zoom_us, Double_t(c));
        n_zoom_before++;
      }

      if (x_after_zoom_us >= 0.0 && x_after_zoom_us <= zoom_width_us) {
        h_extreme_after_zoom->Fill(x_after_zoom_us, Double_t(c));
        n_zoom_after++;
      }
    }
  }

  std::cout << "Zoom plot entries before correction: " << n_zoom_before
            << std::endl;
  std::cout << "Zoom plot entries after correction : " << n_zoom_after
            << std::endl;

  PlotExtremeEvents2D(h_extreme_before, h_extreme_after, h_extreme_before_zoom,
                      h_extreme_after_zoom, file_label, before_zoom_t0_s,
                      after_zoom_t0_s);
  std::cout << "RESULTS (ref board = " << ref_board << ")" << std::endl;
  for (UShort_t board = 0; board < Constants::cfg.N_BOARDS; board++) {
    if (board == ref_board)
      continue;
    std::cout << "  Board " << ref_board << "-" << board << ": "
              << result.board_shifts[board] * 1e-12 << " s" << std::endl;
  }

  delete h_extreme_before;
  delete h_extreme_after;
  delete h_extreme_before_zoom;
  delete h_extreme_after_zoom;

  for (Int_t k = 0; k < Int_t(long_graphs.size()); k++)
    delete long_graphs[k];

  return result;
}

void Timing::ApplyShiftsInPlace(std::vector<RawHit> &hits,
                                const std::vector<Long64_t> &board_shifts) {
  for (Int_t i = 0; i < Int_t(hits.size()); i++) {
    Long64_t board_shift = (hits[i].board < UShort_t(board_shifts.size()))
                               ? board_shifts[hits[i].board]
                               : 0;
    // Per-channel TTF-delay correction (per-dataset; 0 when the map is empty,
    // e.g. 87Rb). Independent of the second-scale board-pattern shift.
    Long64_t ttf_offset =
        Constants::LookupTTFOffsetPs(hits[i].board, hits[i].channel);
    hits[i].timestamp = hits[i].timestamp + board_shift - ttf_offset;
  }
}

void Timing::SortHitsByTimestamp(std::vector<RawHit> &hits) {
  if (GpuAccel::Available() && GpuAccel::TryAcquireSortSlot()) {
    Int_t rc = GpuAccel::GetSort()(hits.data(), Long64_t(hits.size()));
    GpuAccel::ReleaseSortSlot();
    if (rc == 0)
      return;
    std::cerr << "[GPU] Sort failed (rc=" << rc << "), falling back to CPU."
              << std::endl;
  }
  std::sort(hits.begin(), hits.end(), [](const RawHit &a, const RawHit &b) {
    return a.timestamp < b.timestamp;
  });
}

Double_t Timing::FindDensestTimeWindowStartS(std::vector<Double_t> times,
                                             Double_t window_width_s,
                                             Double_t fallback_start_s) {
  if (times.empty() || window_width_s <= 0)
    return fallback_start_s;

  std::sort(times.begin(), times.end());

  Long64_t best_i = 0;
  Long64_t best_count = 0;
  Long64_t j = 0;

  for (Long64_t i = 0; i < Long64_t(times.size()); i++) {
    if (j < i)
      j = i;

    while (j < Long64_t(times.size()) && times[j] <= times[i] + window_width_s)
      j++;

    Long64_t count = j - i;
    if (count > best_count) {
      best_count = count;
      best_i = i;
    }
  }

  // Add a small left padding so the first event is not exactly on the axis
  // edge.
  return times[best_i] - 0.05 * window_width_s;
}
