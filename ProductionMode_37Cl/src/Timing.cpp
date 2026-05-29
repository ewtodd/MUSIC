#include "Timing.hpp"

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

void Timing::PlotExtremeEvents2D(TH2F *h_before, TH2F *h_after,
                                 const TString &file_label) {
  std::lock_guard<std::mutex> lock(g_plot_mutex);
  TString subdir = "timing/" + file_label;

  TCanvas *c_before = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureAndDraw2DHistogram(h_before, c_before);
  c_before->SetLogz(kFALSE);
  if (Constants::SAVE_PLOTS)
    PlottingUtils::SaveFigure(c_before, "extreme_events_before", subdir,
                              PlotSaveOptions::kLINEAR);
  delete c_before;

  TCanvas *c_after = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureAndDraw2DHistogram(h_after, c_after);
  c_after->SetLogz(kFALSE);
  if (Constants::SAVE_PLOTS)
    PlottingUtils::SaveFigure(c_after, "extreme_events_after", subdir,
                              PlotSaveOptions::kLINEAR);
  delete c_after;
}

void Timing::PlotCostLandscape(const std::vector<Double_t> &shifts,
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

  if (Constants::SAVE_PLOTS)
    PlottingUtils::SaveFigure(canvas, Form("cost_landscape_board_%d", board),
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
  result.board_shifts.assign(Constants::N_BOARDS, 0);

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

  std::vector<Double_t> board_shifts_s(Constants::N_BOARDS, 0.0);

  for (UShort_t board = 0; board < Constants::N_BOARDS; board++) {
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

  for (Int_t b = 0; b < n_y; b++) {
    TString label =
        Form("B%d %s", long_channels[b].board, long_channels[b].name.Data());
    h_extreme_before->GetYaxis()->SetBinLabel(b + 1, label);
    h_extreme_after->GetYaxis()->SetBinLabel(b + 1, label);
  }

  for (Int_t c = 0; c < Int_t(long_channels.size()); c++) {
    TGraph *g = long_graphs[c];
    Double_t shift_s = board_shifts_s[long_channels[c].board];
    for (Int_t k = 0; k < g->GetN(); k++) {
      Double_t x, y;
      g->GetPoint(k, x, y);
      h_extreme_before->Fill(x, Double_t(c));
      h_extreme_after->Fill(x - shift_s, Double_t(c));
    }
  }

  PlotExtremeEvents2D(h_extreme_before, h_extreme_after, file_label);

  std::cout << "RESULTS (ref board = " << ref_board << ")" << std::endl;
  for (UShort_t board = 0; board < Constants::N_BOARDS; board++) {
    if (board == ref_board)
      continue;
    std::cout << "  Board " << ref_board << "-" << board << ": "
              << result.board_shifts[board] * 1e-12 << " s" << std::endl;
  }

  delete h_extreme_before;
  delete h_extreme_after;
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
