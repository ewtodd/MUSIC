#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PipelineMutex.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <cstddef>
#include <mutex>
#include <utility>
#include <vector>

struct TimeShiftResult {
  std::vector<Long64_t> board_shifts;
};

void ComputeNSD2(const std::vector<Double_t> &ref_x,
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

  std::size_t gr_idx = 0;

  for (std::size_t p = 0; p < ref_x.size(); p++) {
    Double_t tref = ref_x[p];
    Double_t dtref = ref_y[p];

    if (tref < tmin || tref > tmax || dtref <= thresh_dt_us)
      continue;

    Double_t tref_shifted = tref + shift;

    while (gr_idx < gr_x.size() - 1 && gr_x[gr_idx + 1] < tref_shifted)
      gr_idx++;

    if (gr_idx >= gr_x.size() - 1)
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

TGraph *ExtractBeamTimingStructure(TTree *tree, UShort_t board,
                                   UShort_t channel, Double_t min_energy,
                                   Double_t max_energy, Double_t tmin_s,
                                   Double_t tmax_s, Double_t thresh_dt_us) {
  UShort_t tree_board, tree_channel, tree_energy;
  ULong64_t tree_timestamp;

  tree->SetBranchAddress("Board", &tree_board);
  tree->SetBranchAddress("Channel", &tree_channel);
  tree->SetBranchAddress("Energy", &tree_energy);
  tree->SetBranchAddress("Timestamp", &tree_timestamp);

  std::vector<ULong64_t> timestamps;
  timestamps.reserve(10000);

  Long64_t n_entries = tree->GetEntries();
  std::cout << "Extracting timing structure for Board " << board << " Channel "
            << channel << std::endl;

  for (Long64_t i = 0; i < n_entries; i++) {
    tree->GetEntry(i);

    Double_t time_s = tree_timestamp / 1e12;

    if (tree_board == board && tree_channel == channel &&
        tree_energy >= min_energy && tree_energy <= max_energy &&
        time_s >= tmin_s && time_s <= tmax_s) {
      timestamps.push_back(tree_timestamp);
    }

    if (i % 10000000 == 0)
      std::cout << "  Progress: " << i << "/" << n_entries << std::endl;
  }

  std::cout << "  Found " << timestamps.size() << " events in window"
            << std::endl;

  TGraph *delta_graph = new TGraph();

  for (std::size_t i = 0; i + 1 < timestamps.size(); i++) {
    Double_t dt_us = (timestamps[i + 1] - timestamps[i]) / 1e6;

    if (dt_us > thresh_dt_us) {
      Double_t time_s = timestamps[i] / 1e12;
      delta_graph->SetPoint(delta_graph->GetN(), time_s, dt_us);
    }
  }

  std::cout << "  Built graph with " << delta_graph->GetN()
            << " large delta-T points" << std::endl;

  return delta_graph;
}

struct BinnedShift {
  Double_t bin_center_s;
  Double_t shift_s;
  Int_t n_ref;
  Int_t n_gr;
  Bool_t valid;
};

std::vector<BinnedShift> FindShiftBeamBinned(TGraph *ref, TGraph *gr,
                                             Double_t overlap_tmin_s,
                                             Double_t overlap_tmax_s,
                                             Double_t thresh_dt_us,
                                             Int_t n_bins,
                                             Double_t global_shift_s,
                                             Double_t max_drift_s) {
  std::vector<BinnedShift> results;
  Double_t bin_width_s = (overlap_tmax_s - overlap_tmin_s) / n_bins;
  Double_t step_s = 1e-6;
  Int_t n_steps = 2 * static_cast<Int_t>(max_drift_s / step_s) + 1;

  for (Int_t b = 0; b < n_bins; b++) {
    Double_t bin_start = overlap_tmin_s + b * bin_width_s;
    Double_t bin_end = bin_start + bin_width_s;
    Double_t bin_center = 0.5 * (bin_start + bin_end);

    std::vector<Double_t> ref_x, ref_y, gr_x, gr_y;
    for (Int_t i = 0; i < ref->GetN(); i++) {
      Double_t x, y;
      ref->GetPoint(i, x, y);
      if (x >= bin_start && x < bin_end) {
        ref_x.push_back(x);
        ref_y.push_back(y);
      }
    }
    for (Int_t i = 0; i < gr->GetN(); i++) {
      Double_t x, y;
      gr->GetPoint(i, x, y);
      if (x >= bin_start && x < bin_end) {
        gr_x.push_back(x);
        gr_y.push_back(y);
      }
    }

    BinnedShift r;
    r.bin_center_s = bin_center;
    r.shift_s = global_shift_s;
    r.n_ref = static_cast<Int_t>(ref_x.size());
    r.n_gr = static_cast<Int_t>(gr_x.size());
    r.valid = kFALSE;

    if (r.n_ref < 20 || r.n_gr < 20) {
      results.push_back(r);
      continue;
    }

    Double_t best_shift = global_shift_s;
    Double_t max_inv_nsd2 = 0;

    for (Int_t k = 0; k < n_steps; k++) {
      Double_t shift = global_shift_s - max_drift_s + k * step_s;
      Int_t npts = 0;
      Double_t nsd2 = 0;
      ComputeNSD2(ref_x, ref_y, gr_x, gr_y, shift, thresh_dt_us, npts, nsd2);
      if (npts > 5) {
        Double_t inv_nsd2 = 1.0 / nsd2;
        if (inv_nsd2 > max_inv_nsd2) {
          best_shift = shift;
          max_inv_nsd2 = inv_nsd2;
        }
      }
    }

    r.shift_s = best_shift;
    r.valid = (max_inv_nsd2 > 0);
    results.push_back(r);
  }

  return results;
}

void PlotBinnedShifts(const std::vector<BinnedShift> &binned, UShort_t board,
                      const TString &file_label) {
  std::lock_guard<std::mutex> lock(g_plot_mutex);

  TGraph *gr_valid = new TGraph();
  TGraph *gr_invalid = new TGraph();
  for (std::size_t i = 0; i < binned.size(); i++) {
    if (binned[i].valid)
      gr_valid->SetPoint(gr_valid->GetN(), binned[i].bin_center_s,
                         binned[i].shift_s * 1e3);
    else
      gr_invalid->SetPoint(gr_invalid->GetN(), binned[i].bin_center_s, 0.0);
  }

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureGraph(
      gr_valid, kBlue + 1,
      Form("Board 0-%d binned shift;Time bin center [s];Shift [ms]", board));
  gr_valid->SetMarkerStyle(20);
  gr_valid->SetMarkerSize(0.9);
  gr_valid->SetMarkerColor(kBlue + 1);
  gr_valid->Draw("AP");

  if (gr_invalid->GetN() > 0) {
    gr_invalid->SetMarkerStyle(24);
    gr_invalid->SetMarkerSize(0.9);
    gr_invalid->SetMarkerColor(kRed + 1);
    gr_invalid->Draw("P SAME");
  }

  PlottingUtils::SaveFigure(canvas, Form("binned_shift_board_%d", board),
                            "timing/" + file_label, PlotSaveOptions::kLINEAR);
  delete canvas;
}

void PlotDeltaTTimeDistribution(TGraph *graph, UShort_t board,
                                Double_t overlap_tmin_s,
                                Double_t overlap_tmax_s,
                                const TString &file_label) {
  std::lock_guard<std::mutex> lock(g_plot_mutex);

  TH1F *h = new TH1F(
      PlottingUtils::GetRandomName(),
      Form("Board %d large-#Delta t event times;Time [s];Counts", board), 200,
      overlap_tmin_s, overlap_tmax_s);

  for (Int_t i = 0; i < graph->GetN(); i++) {
    Double_t x, y;
    graph->GetPoint(i, x, y);
    h->Fill(x);
  }

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureAndDrawHistogram(h, kBlue + 1);
  PlottingUtils::SaveFigure(canvas, Form("dt_time_dist_board_%d", board),
                            "timing/" + file_label, PlotSaveOptions::kLINEAR);
  delete canvas;
}

void PlotResiduals(TGraph *ref_graph, TGraph *board_graph, Double_t shift_s,
                   UShort_t board, const TString &file_label) {
  std::lock_guard<std::mutex> lock(g_plot_mutex);
  TGraph *board_shifted = static_cast<TGraph *>(board_graph->Clone());
  for (Int_t i = 0; i < board_shifted->GetN(); i++) {
    Double_t x, y;
    board_shifted->GetPoint(i, x, y);
    board_shifted->SetPoint(i, x - shift_s, y);
  }

  TH1D *h_residuals = new TH1D(
      PlottingUtils::GetRandomName(),
      Form("Residuals Board 0-%d;#Delta t_{ref} - #Delta t_{board} [#mus];"
           "Counts",
           board),
      200, -5, 5);

  for (Int_t i = 0; i < ref_graph->GetN(); i++) {
    Double_t t_ref, dt_ref;
    ref_graph->GetPoint(i, t_ref, dt_ref);

    Double_t dt_board = board_shifted->Eval(t_ref, 0, "S");
    if (dt_board > 0)
      h_residuals->Fill(dt_ref - dt_board);
  }

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureAndDrawHistogram(h_residuals, kBlue + 1);

  h_residuals->Fit("gaus", "Q");
  TF1 *fit = h_residuals->GetFunction("gaus");
  if (fit) {
    fit->SetLineColor(kRed + 1);
    fit->SetLineWidth(PlottingUtils::GetLineWidth());
    PlottingUtils::AddText(Form("#mu = %.2f #mus", fit->GetParameter(1)), 0.85,
                           0.85);
    PlottingUtils::AddText(Form("#sigma = %.2f #mus", fit->GetParameter(2)),
                           0.85, 0.78);
  }

  PlottingUtils::SaveFigure(canvas, Form("residuals_board_%d", board),
                            "timing/" + file_label, PlotSaveOptions::kLINEAR);
  delete canvas;
  delete board_shifted;
}

void PlotCorrelation(TGraph *ref_graph, TGraph *board_graph, Double_t shift_s,
                     UShort_t board, const TString &file_label) {
  std::lock_guard<std::mutex> lock(g_plot_mutex);
  TGraph *board_shifted = static_cast<TGraph *>(board_graph->Clone());
  for (Int_t i = 0; i < board_shifted->GetN(); i++) {
    Double_t x, y;
    board_shifted->GetPoint(i, x, y);
    board_shifted->SetPoint(i, x - shift_s, y);
  }

  Double_t corr_max = 2200;
  TH2D *h_corr = new TH2D(PlottingUtils::GetRandomName(),
                          Form("Correlation Board 0-%d;#Delta t_{ref} "
                               "[#mus];#Delta t_{board %d} [#mus]",
                               board, board),
                          100, 0, corr_max, 100, 0, corr_max);

  for (Int_t i = 0; i < ref_graph->GetN(); i++) {
    Double_t t_ref, dt_ref;
    ref_graph->GetPoint(i, t_ref, dt_ref);

    Double_t dt_board = board_shifted->Eval(t_ref, 0, "S");
    if (dt_board > 0 && dt_ref > 0)
      h_corr->Fill(dt_ref, dt_board);
  }

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  canvas->SetLeftMargin(0.18);
  PlottingUtils::ConfigureAndDraw2DHistogram(h_corr, canvas);
  h_corr->GetYaxis()->SetTitleOffset(1.4);

  TLine *diag = new TLine(0, 0, corr_max, corr_max);
  diag->SetLineColor(kRed + 1);
  diag->SetLineWidth(PlottingUtils::GetLineWidth());
  diag->SetLineStyle(2);
  diag->Draw();

  PlottingUtils::SaveFigure(canvas, Form("correlation_board_%d", board),
                            "timing/" + file_label, PlotSaveOptions::kLINEAR);
  delete canvas;
  delete board_shifted;
}

Double_t FindShiftBeam(TGraph *ref, TGraph *gr, Double_t overlap_tmin_s,
                       Double_t overlap_tmax_s, Double_t thresh_dt_us,
                       UShort_t board, const TString &file_label) {

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
  for (std::size_t p = 0; p < ref_x.size(); p++) {
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

  for (Int_t p = 0; p < gr->GetN(); p++) {
    Int_t npts = 0;
    Double_t nsd2 = 0;
    Double_t shift = gr_x[p] - tini;

    if (p % 1000 == 0)
      std::cout << "  Testing shift: " << shift << " s (" << p << "/"
                << gr->GetN() << ")" << std::endl;

    ComputeNSD2(ref_x, ref_y, gr_x, gr_y, shift, thresh_dt_us, npts, nsd2);

    if (npts > 10) {
      Double_t inv_nsd2 = 1.0 / nsd2;

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
    std::cout << "  Improvement vs no-shift: "
              << max_inv_nsd2 / inv_nsd2_zero << "x" << std::endl;

  return best_shift;
}

std::vector<TimeShiftResult> CalcTimeShiftsBeamMethod(
    std::vector<TString> input_filepaths, std::vector<TString> file_labels,
    UShort_t ref_board, UShort_t ref_channel,
    std::vector<UShort_t> board_channels, Double_t min_energy,
    Double_t max_energy, Double_t overlap_margin_s, Double_t thresh_dt_us,
    Bool_t reprocess) {

  std::vector<TimeShiftResult> results;

  if (!reprocess) {
    std::cout << "Skipping beam-based time shift calculation" << std::endl;
    return results;
  }

  Int_t n_files = input_filepaths.size();

  for (Int_t i = 0; i < n_files; i++) {
    TString input_filepath = input_filepaths[i];
    TString file_label = file_labels[i];
    std::cout << "Processing: " << input_filepath << std::endl;

    TFile *input_file = IO::OpenForReading(input_filepath);
    if (!input_file || input_file->IsZombie()) {
      std::cerr << "Error opening file" << std::endl;
      continue;
    }

    TTree *tree = static_cast<TTree *>(input_file->Get("Data_R"));
    if (!tree) {
      std::cerr << "Error getting tree" << std::endl;
      input_file->Close();
      continue;
    }
    tree->LoadBaskets();

    Double_t file_tmin_s = tree->GetMinimum("Timestamp") / 1e12;
    Double_t file_tmax_s = tree->GetMaximum("Timestamp") / 1e12;
    Double_t overlap_tmin_s = file_tmin_s + overlap_margin_s;
    Double_t overlap_tmax_s = file_tmax_s - overlap_margin_s;
    std::cout << "File span [" << file_tmin_s << ", " << file_tmax_s
              << "] s, using overlap [" << overlap_tmin_s << ", "
              << overlap_tmax_s << "] s" << std::endl;

    std::cout << "Extracting reference board" << std::endl;
    TGraph *ref_graph = ExtractBeamTimingStructure(
        tree, ref_board, ref_channel, min_energy, max_energy, overlap_tmin_s,
        overlap_tmax_s, thresh_dt_us);
    PlotDeltaTTimeDistribution(ref_graph, ref_board, overlap_tmin_s,
                               overlap_tmax_s, file_label);

    TimeShiftResult result;
    result.board_shifts.assign(Constants::N_BOARDS, 0);

    for (UShort_t board = 0; board < Constants::N_BOARDS; board++) {
      if (board == ref_board)
        continue;

      std::cout << "Processing Board " << board << " Channel "
                << board_channels[board] << std::endl;

      TGraph *board_graph = ExtractBeamTimingStructure(
          tree, board, board_channels[board], min_energy, max_energy,
          overlap_tmin_s, overlap_tmax_s, thresh_dt_us);
      PlotDeltaTTimeDistribution(board_graph, board, overlap_tmin_s,
                                 overlap_tmax_s, file_label);

      if (board_graph->GetN() < 10) {
        std::cout << "WARNING: Not enough data points for Board " << board
                  << std::endl;
        delete board_graph;
        continue;
      }

      Double_t shift_s = FindShiftBeam(
          ref_graph, board_graph,
          overlap_tmin_s + 0.1 * (overlap_tmax_s - overlap_tmin_s),
          overlap_tmax_s - 0.1 * (overlap_tmax_s - overlap_tmin_s),
          thresh_dt_us, board, file_label);

      PlotResiduals(ref_graph, board_graph, shift_s, board, file_label);
      PlotCorrelation(ref_graph, board_graph, shift_s, board, file_label);

      std::vector<BinnedShift> binned = FindShiftBeamBinned(
          ref_graph, board_graph, overlap_tmin_s, overlap_tmax_s, thresh_dt_us,
          Constants::N_TIME_BINS, shift_s, Constants::MAX_BIN_DRIFT_S);
      PlotBinnedShifts(binned, board, file_label);

      Long64_t shift_ps = static_cast<Long64_t>(shift_s * 1e12);
      result.board_shifts[board] = -shift_ps;

      std::cout << "Board " << ref_board << "-" << board
                << " shift: " << shift_s << " s (" << shift_ps << " ps)"
                << std::endl;

      delete board_graph;
    }

    results.push_back(result);

    std::cout << "RESULTS (ref board = " << ref_board << ")" << std::endl;
    for (UShort_t board = 0; board < Constants::N_BOARDS; board++) {
      if (board == ref_board)
        continue;
      std::cout << "  Board " << ref_board << "-" << board << ": "
                << result.board_shifts[board] * 1e-12 << " s" << std::endl;
    }

    delete ref_graph;
    input_file->Close();
    delete input_file;
  }

  return results;
}

void ApplyTimeShift(std::vector<TString> input_filepaths,
                    std::vector<TimeShiftResult> results, Bool_t reprocess) {
  if (!reprocess)
    return;

  Int_t n_files = input_filepaths.size();
  for (Int_t i = 0; i < n_files; i++) {
    TString input_filepath = input_filepaths[i];
    TimeShiftResult res = results[i];

    TFile *input_output_file = IO::OpenForWriting(input_filepath, "UPDATE");

    TTree *input_tree = static_cast<TTree *>(input_output_file->Get("Data_R"));
    if (!input_tree) {
      std::cerr << "Error getting tree for filepath " << input_filepath
                << std::endl;
      input_output_file->Close();
      continue;
    }

    UShort_t board;
    ULong64_t timestamp, shifted_timestamp;
    input_tree->SetBranchAddress("Board", &board);
    input_tree->SetBranchAddress("Timestamp", &timestamp);

    TBranch *existing_branch = input_tree->GetBranch("ShiftedTimestamp");
    if (existing_branch) {
      std::cout << "ShiftedTimestamp branch already exists, will overwrite"
                << std::endl;
      input_tree->SetBranchAddress("ShiftedTimestamp", &shifted_timestamp);
    } else {
      input_tree->Branch("ShiftedTimestamp", &shifted_timestamp,
                         "ShiftedTimestamp/l");
    }

    Int_t n_entries = input_tree->GetEntries();
    input_tree->LoadBaskets();
    for (Int_t j = 0; j < n_entries; j++) {
      input_tree->GetEntry(j);
      Long64_t shifted_timestamp_calc = timestamp;
      if (board < res.board_shifts.size())
        shifted_timestamp_calc = timestamp + res.board_shifts[board];
      shifted_timestamp = shifted_timestamp_calc;

      if (existing_branch) {
        input_tree->GetBranch("ShiftedTimestamp")->Fill();
      }

      if (j % 1000000 == 0) {
        std::cout << "Applying shifts: " << j << "/" << n_entries << std::endl;
      }
    }

    input_output_file->cd();
    input_tree->Write("", TObject::kOverwrite);
    input_output_file->Close();

    std::cout << "Applied time shifts to: " << input_filepath << std::endl;
  }
}

void TimesortData(std::vector<TString> input_filepaths,
                  std::vector<TString> output_names, Bool_t reprocess = kTRUE) {
  if (!reprocess)
    return;

  Int_t n_files = input_filepaths.size();
  for (Int_t i = 0; i < n_files; i++) {

    TString input_filepath = input_filepaths[i];

    TFile *input_file = IO::OpenForReading(input_filepath);
    TTree *input_tree = static_cast<TTree *>(input_file->Get("Data_R"));

    Long64_t n_entries = input_tree->GetEntries();
    std::cout << "Total entries: " << n_entries << std::endl;

    ULong64_t shifted_timestamp;
    input_tree->SetBranchAddress("ShiftedTimestamp", &shifted_timestamp);

    std::vector<std::pair<ULong64_t, Long64_t>> timestamp_index;
    timestamp_index.reserve(n_entries);

    std::cout << "Reading timestamps" << std::endl;
    for (Long64_t j = 0; j < n_entries; j++) {
      input_tree->GetEntry(j);
      timestamp_index.push_back(std::make_pair(shifted_timestamp, j));

      if (j % 10000000 == 0) {
        std::cout << "Progress: " << j << "/" << n_entries << std::endl;
      }
    }

    std::cout << "Sorting" << std::endl;
    std::sort(timestamp_index.begin(), timestamp_index.end());

    std::cout << "Writing sorted tree" << std::endl;

    TString output_name = output_names[i];
    TString output_filepath =
        IO::GetRootFilesBaseDir() + "/" + output_name + ".root";
    gSystem->mkdir(gSystem->DirName(output_filepath), kTRUE);
    TFile *output_file = new TFile(output_filepath, "RECREATE", "", 0);

    TTree *output_tree = input_tree->CloneTree(0);

    std::cout << "Loading baskets into memory" << std::endl;
    input_tree->LoadBaskets();

    std::cout << "Copying entries in sorted order" << std::endl;
    for (Long64_t j = 0; j < n_entries; j++) {
      input_tree->GetEntry(timestamp_index[j].second);
      output_tree->Fill();

      if (j % 10000000 == 0) {
        std::cout << "Progress: " << j << "/" << n_entries << std::endl;
      }
    }

    output_file->cd();
    output_tree->Write();
    output_file->Close();
    input_file->Close();
    std::cout << "Sorted file saved to: " << output_filepath << std::endl;
  }
}

void Timing() {
  Bool_t reprocess_calculation = kTRUE;
  Bool_t reprocess_sorting = kTRUE;

  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                TString(gSystem->pwd()) + "/plots",
                                TString(gSystem->pwd()) + "/root_files");

  std::vector<TString> filepaths, sorted_output_names, file_labels;

  std::vector<FileSpec> specs = BuildFileSpecs();
  for (std::size_t k = 0; k < specs.size(); k++) {
    filepaths.push_back(RawRootName(specs[k]));
    std::cout << "Will process: " << filepaths.back() << std::endl;
    sorted_output_names.push_back(SortedName(specs[k]));
    file_labels.push_back(FileLabel(specs[k]));
  }

  std::vector<TimeShiftResult> results = CalcTimeShiftsBeamMethod(
      filepaths, file_labels, Constants::REF_BOARD,
      Constants::BOARD_CHANNELS[Constants::REF_BOARD], Constants::BOARD_CHANNELS,
      Constants::MIN_ENERGY, Constants::MAX_ENERGY, Constants::OVERLAP_MARGIN_S,
      Constants::THRESH_DT_US, reprocess_calculation);
  ApplyTimeShift(filepaths, results, reprocess_calculation);

  TimesortData(filepaths, sorted_output_names, reprocess_sorting);

  std::cout << "TIMING ANALYSIS COMPLETE" << std::endl;
}
