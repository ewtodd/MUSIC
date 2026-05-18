#include "Timing.hpp"
#include "BinaryUtils.hpp"
#include "Constants.hpp"
#include <TGraph.h>
#include <TH2F.h>
#include <TMath.h>
#include <TString.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

std::vector<TGraph *> ExtractAllChannelsTimingStructureFromHits(
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

TimeShiftResult CalcTimeShiftsBeamMethodFromHits(
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

void ApplyShiftsInPlace(std::vector<RawHit> &hits,
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

void SortHitsByTimestamp(std::vector<RawHit> &hits) {
  std::sort(hits.begin(), hits.end(), [](const RawHit &a, const RawHit &b) {
    return a.timestamp < b.timestamp;
  });
}
