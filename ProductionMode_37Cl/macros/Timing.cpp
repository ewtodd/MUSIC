#include "Constants.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TGraph.h>
#include <TMarker.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <vector>

struct TimeShiftResult {
  Long64_t board_0_1_average;
  Long64_t board_0_1_stddev;
  Long64_t board_0_2_average;
  Long64_t board_0_2_stddev;
  Long64_t board_0_3_average;
  Long64_t board_0_3_stddev;
};

void ComputeNSD2(const std::vector<Double_t> &ref_x,
                 const std::vector<Double_t> &ref_y,
                 const std::vector<Double_t> &gr_x,
                 const std::vector<Double_t> &gr_y, Double_t shift,
                 Double_t thresh_dt_us, Double_t overlap_tmin_s,
                 Double_t overlap_tmax_s, Int_t &npts, Double_t &nsd2) {
  nsd2 = 0;
  npts = 0;

  Double_t tmin_gr = gr_x.empty() ? 0 : gr_x.front() - shift;
  Double_t tmax_gr = gr_x.empty() ? 0 : gr_x.back() - shift;
  Double_t tmin_ref = ref_x.empty() ? 0 : ref_x.front();
  Double_t tmax_ref = ref_x.empty() ? 0 : ref_x.back();

  Double_t tmin = TMath::Max(tmin_ref, tmin_gr);
  Double_t tmax = TMath::Min(tmax_ref, tmax_gr);

  size_t gr_idx = 0;

  for (size_t p = 0; p < ref_x.size(); p++) {
    Double_t tref = ref_x[p];
    Double_t dtref = ref_y[p];

    if (tref >= tmin && tref <= tmax && dtref > thresh_dt_us) {
      Double_t tref_shifted = tref + shift;

      while (gr_idx < gr_x.size() - 1 && gr_x[gr_idx + 1] < tref_shifted) {
        gr_idx++;
      }

      if (gr_idx < gr_x.size() - 1) {
        Double_t x0 = gr_x[gr_idx];
        Double_t x1 = gr_x[gr_idx + 1];
        Double_t y0 = gr_y[gr_idx];
        Double_t y1 = gr_y[gr_idx + 1];

        Double_t dt = y0 + (y1 - y0) * (tref_shifted - x0) / (x1 - x0);
        nsd2 += pow(dt - dtref, 2);
        npts++;
      }
    }
  }

  if (nsd2 > 0 && npts > 0) {
    nsd2 = sqrt(nsd2) / npts;
  } else {
    nsd2 = 1e12;
  }
}

Double_t FindShiftBeam(TGraph *ref, TGraph *gr, Double_t overlap_tmin_s,
                       Double_t overlap_tmax_s, Double_t thresh_dt_us,
                       Int_t max_iterations = 10000) {

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
  for (size_t p = 0; p < ref_x.size(); p++) {
    if (ref_x[p] >= overlap_tmin_s) {
      tini = ref_x[p];
      break;
    }
  }

  std::cout << "Looking for timeshift between " << overlap_tmin_s << " - "
            << overlap_tmax_s << " sec" << std::endl;
  std::cout << "tini = " << tini << " s" << std::endl;

  Int_t step = TMath::Max(1, gr->GetN() / max_iterations);

  std::vector<Double_t> shifts;
  std::vector<Double_t> inv_nsd2_vals;
  shifts.reserve(max_iterations);
  inv_nsd2_vals.reserve(max_iterations);

  Double_t best_shift = 0;
  Double_t max_inv_nsd2 = 0;

  std::cout << "Testing " << gr->GetN() / step << " shift values..."
            << std::endl;

  for (Int_t p = 0; p < gr->GetN(); p += step) {
    Int_t npts = 0;
    Double_t nsd2 = 0;
    Double_t shift = gr_x[p] - tini;

    if (p % 1000 == 0)
      std::cout << "Testing shift: " << shift << " s (" << p / step << "/"
                << gr->GetN() / step << ")" << std::endl;

    ComputeNSD2(ref_x, ref_y, gr_x, gr_y, shift, thresh_dt_us, overlap_tmin_s,
                overlap_tmax_s, npts, nsd2);

    if (npts > 10) {
      Double_t inv_nsd2 = 1.0 / nsd2;
      shifts.push_back(shift);
      inv_nsd2_vals.push_back(inv_nsd2);

      if (inv_nsd2 > max_inv_nsd2) {
        best_shift = shift;
        max_inv_nsd2 = inv_nsd2;
      }
    }
  }

  std::cout << "Best shift: " << best_shift << " s (1/NSD2 = " << max_inv_nsd2
            << ")" << std::endl;

  return best_shift;
}

TGraph *ExtractBeamTimingStructure(TTree *tree, UShort_t board,
                                   UShort_t channel, Double_t min_energy,
                                   Double_t max_energy, Double_t tmin_s,
                                   Double_t tmax_s, Double_t thresh_dt_us,
                                   Double_t timeshift_s = 0.0) {
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

    if (i % 1000000 == 0) {
      std::cout << "Progress: " << i << "/" << n_entries << std::endl;
    }
  }

  std::cout << "Found " << timestamps.size() << " events" << std::endl;

  TGraph *delta_graph = new TGraph();

  for (size_t i = 0; i < timestamps.size() - 1; i++) {
    Double_t dt_us = (timestamps[i + 1] - timestamps[i]) / 1e6;

    if (dt_us > thresh_dt_us) {
      Double_t time_s = (timestamps[i] / 1e12) - timeshift_s;
      delta_graph->SetPoint(delta_graph->GetN(), time_s, dt_us);
    }
  }

  std::cout << "Created graph with " << delta_graph->GetN()
            << " large delta-T points" << std::endl;

  return delta_graph;
}

void VisualizeTimeShiftQuality(TGraph *ref_graph, TGraph *shifted_graph,
                               Double_t shift_s, UShort_t board,
                               const TString &output_name) {
  TCanvas *c =
      new TCanvas(Form("c_board_%d", board),
                  Form("Time Shift Quality - Board %d", board), 1200, 800);
  c->Divide(2, 2);

  c->cd(1);
  gPad->SetLogy();
  ref_graph->SetMarkerStyle(20);
  ref_graph->SetMarkerSize(0.5);
  ref_graph->SetMarkerColor(kBlue);
  ref_graph->SetTitle(Form(
      "Board %d Timing Structure Overlay;Time [s];#Delta t [#mus]", board));
  ref_graph->Draw("AP");

  TGraph *shifted_clone = (TGraph *)shifted_graph->Clone();
  for (Int_t i = 0; i < shifted_clone->GetN(); i++) {
    Double_t x, y;
    shifted_clone->GetPoint(i, x, y);
    shifted_clone->SetPoint(i, x - shift_s, y);
  }
  shifted_clone->SetMarkerStyle(24);
  shifted_clone->SetMarkerSize(0.5);
  shifted_clone->SetMarkerColor(kRed);
  shifted_clone->Draw("P SAME");

  TLegend *leg1 = new TLegend(0.65, 0.75, 0.89, 0.89);
  leg1->AddEntry(ref_graph, "Reference (Board 0)", "p");
  leg1->AddEntry(shifted_clone, Form("Board %d (shifted)", board), "p");
  leg1->Draw();

  c->cd(2);
  TH1D *h_residuals = new TH1D(
      Form("h_res_board_%d", board),
      Form("Residuals Board %d;#Delta t_{ref} - #Delta t_{board} [#mus];Counts",
           board),
      200, -50, 50);

  for (Int_t i = 0; i < ref_graph->GetN(); i++) {
    Double_t t_ref, dt_ref;
    ref_graph->GetPoint(i, t_ref, dt_ref);

    Double_t dt_board = shifted_clone->Eval(t_ref, 0, "S");
    if (dt_board > 0) {
      h_residuals->Fill(dt_ref - dt_board);
    }
  }

  h_residuals->SetLineColor(kBlue);
  h_residuals->SetFillColor(kBlue);
  h_residuals->SetFillStyle(3004);
  h_residuals->Draw();

  h_residuals->Fit("gaus", "Q");
  TF1 *fit = h_residuals->GetFunction("gaus");
  if (fit) {
    fit->SetLineColor(kRed);
    fit->SetLineWidth(2);

    TPaveText *pt = new TPaveText(0.55, 0.65, 0.89, 0.89, "NDC");
    pt->AddText(Form("Mean: %.2f #mus", fit->GetParameter(1)));
    pt->AddText(Form("Sigma: %.2f #mus", fit->GetParameter(2)));
    pt->SetFillColor(kWhite);
    pt->Draw();
  }

  c->cd(3);
  gPad->SetLogy();
  TGraph *gr_zoom = (TGraph *)ref_graph->Clone();
  TGraph *shifted_zoom = (TGraph *)shifted_clone->Clone();

  gr_zoom->SetTitle("Zoomed Overlap Region;Time [s];#Delta t [#mus]");
  gr_zoom->GetXaxis()->SetRangeUser(1.0, 5.0);
  gr_zoom->Draw("AP");
  shifted_zoom->Draw("P SAME");

  c->cd(4);
  TH2D *h_corr = new TH2D(Form("h_corr_board_%d", board),
                          Form("Correlation Board %d;#Delta t_{ref} "
                               "[#mus];#Delta t_{board %d} [#mus]",
                               board, board),
                          100, 0, 2000, 100, 0, 2000);

  for (Int_t i = 0; i < ref_graph->GetN(); i++) {
    Double_t t_ref, dt_ref;
    ref_graph->GetPoint(i, t_ref, dt_ref);

    Double_t dt_board = shifted_clone->Eval(t_ref, 0, "S");
    if (dt_board > 0 && dt_ref > 0) {
      h_corr->Fill(dt_ref, dt_board);
    }
  }

  h_corr->Draw("COLZ");

  TLine *diag = new TLine(0, 0, 2000, 2000);
  diag->SetLineColor(kRed);
  diag->SetLineWidth(2);
  diag->SetLineStyle(2);
  diag->Draw();

  PlottingUtils::SaveFigure(
      c, Form("timeshift_quality_board_%d_%s", board, output_name.Data()),
      "", PlotSaveOptions::kLINEAR);

  delete shifted_clone;
  delete shifted_zoom;
  delete c;
}

void VisualizeShiftScan(const std::vector<Double_t> &shifts,
                        const std::vector<Double_t> &inv_nsd2_vals,
                        Double_t best_shift, UShort_t board,
                        const TString &output_name) {
  TCanvas *c = new TCanvas(Form("c_scan_board_%d", board),
                           Form("Shift Scan - Board %d", board), 800, 600);

  TGraph *gr_scan = new TGraph(shifts.size());
  for (size_t i = 0; i < shifts.size(); i++) {
    gr_scan->SetPoint(i, shifts[i], inv_nsd2_vals[i]);
  }

  gr_scan->SetTitle(
      Form("Time Shift Scan Board %d;Shift [s];1/NSD^{2}", board));
  gr_scan->SetMarkerStyle(20);
  gr_scan->SetMarkerSize(0.8);
  gr_scan->SetMarkerColor(kBlue);
  gr_scan->Draw("AP");

  TMarker *best_marker = new TMarker(best_shift, gr_scan->Eval(best_shift), 29);
  best_marker->SetMarkerColor(kRed);
  best_marker->SetMarkerSize(3);
  best_marker->Draw();

  TLegend *leg = new TLegend(0.15, 0.75, 0.45, 0.89);
  leg->AddEntry(gr_scan, "Scan points", "p");
  leg->AddEntry(best_marker, Form("Best: %.3f s", best_shift), "p");
  leg->Draw();

  PlottingUtils::SaveFigure(
      c, Form("shift_scan_board_%d_%s", board, output_name.Data()),
      "", PlotSaveOptions::kLINEAR);
  delete c;
}

Double_t FindShiftBeam(TGraph *ref, TGraph *gr, Double_t overlap_tmin_s,
                       Double_t overlap_tmax_s, Double_t thresh_dt_us,
                       UShort_t board, const TString &output_name,
                       Int_t max_iterations = 10000) {

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
  for (size_t p = 0; p < ref_x.size(); p++) {
    if (ref_x[p] >= overlap_tmin_s) {
      tini = ref_x[p];
      break;
    }
  }

  std::cout << "Looking for timeshift between " << overlap_tmin_s << " - "
            << overlap_tmax_s << " sec" << std::endl;
  std::cout << "tini = " << tini << " s" << std::endl;

  Int_t step = TMath::Max(1, gr->GetN() / max_iterations);

  std::vector<Double_t> shifts;
  std::vector<Double_t> inv_nsd2_vals;
  shifts.reserve(max_iterations);
  inv_nsd2_vals.reserve(max_iterations);

  Double_t best_shift = 0;
  Double_t max_inv_nsd2 = 0;

  std::cout << "Testing " << gr->GetN() / step << " shift values..."
            << std::endl;

  for (Int_t p = 0; p < gr->GetN(); p += step) {
    Int_t npts = 0;
    Double_t nsd2 = 0;
    Double_t shift = gr_x[p] - tini;

    if (p % 1000 == 0)
      std::cout << "Testing shift: " << shift << " s (" << p / step << "/"
                << gr->GetN() / step << ")" << std::endl;

    ComputeNSD2(ref_x, ref_y, gr_x, gr_y, shift, thresh_dt_us, overlap_tmin_s,
                overlap_tmax_s, npts, nsd2);

    if (npts > 10) {
      Double_t inv_nsd2 = 1.0 / nsd2;
      shifts.push_back(shift);
      inv_nsd2_vals.push_back(inv_nsd2);

      if (inv_nsd2 > max_inv_nsd2) {
        best_shift = shift;
        max_inv_nsd2 = inv_nsd2;
      }
    }
  }

  std::cout << "Best shift: " << best_shift << " s (1/NSD2 = " << max_inv_nsd2
            << ")" << std::endl;

  VisualizeShiftScan(shifts, inv_nsd2_vals, best_shift, board, output_name);

  return best_shift;
}

std::vector<TimeShiftResult> CalcTimeShiftsBeamMethod(
    std::vector<TString> input_filepaths, UShort_t ref_board = 0,
    UShort_t ref_channel = 0,
    std::vector<UShort_t> board_channels = {0, 0, 0, 0},
    Double_t min_energy = 2000, Double_t max_energy = 20000,
    Double_t overlap_tmin_s = 1.0, Double_t overlap_tmax_s = 10.0,
    Double_t thresh_dt_us = 200.0, Bool_t reprocess = kTRUE) {

  std::vector<TimeShiftResult> results;

  if (!reprocess) {
    std::cout << "Skipping beam-based time shift calculation" << std::endl;
    return results;
  }

  Int_t n_files = input_filepaths.size();

  for (Int_t i = 0; i < n_files; i++) {
    TString input_filepath = input_filepaths[i];
    TString file_label = Form("run37_%d", i);
    std::cout << "Processing: " << input_filepath << std::endl;

    TFile *input_file = new TFile(input_filepath, "READ");
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

    std::cout << "Extracting reference board" << std::endl;
    TGraph *ref_graph = ExtractBeamTimingStructure(
        tree, ref_board, ref_channel, min_energy, max_energy, overlap_tmin_s,
        overlap_tmax_s, thresh_dt_us, 0.0);

    TimeShiftResult result;
    result.board_0_1_average = 0;
    result.board_0_1_stddev = 0;
    result.board_0_2_average = 0;
    result.board_0_2_stddev = 0;
    result.board_0_3_average = 0;
    result.board_0_3_stddev = 0;

    for (UShort_t board = 1; board < 4; board++) {
      std::cout << "Processing Board " << board << " Channel "
                << board_channels[board] << std::endl;

      TGraph *board_graph = ExtractBeamTimingStructure(
          tree, board, board_channels[board], min_energy, max_energy,
          overlap_tmin_s, overlap_tmax_s, thresh_dt_us, 0.0);

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

      VisualizeTimeShiftQuality(ref_graph, board_graph, shift_s, board,
                                file_label);

      Long64_t shift_ps = static_cast<Long64_t>(shift_s * 1e12);

      switch (board) {
      case 1:
        result.board_0_1_average = -shift_ps;
        break;
      case 2:
        result.board_0_2_average = -shift_ps;
        break;
      case 3:
        result.board_0_3_average = -shift_ps;
        break;
      }

      std::cout << "Board 0-" << board << " shift: " << shift_s << " s ("
                << shift_ps << " ps)" << std::endl;

      delete board_graph;
    }

    results.push_back(result);

    std::cout << "RESULTS" << std::endl;
    std::cout << "Board 0-1: " << result.board_0_1_average * 1e-12 << " s"
              << std::endl;
    std::cout << "Board 0-2: " << result.board_0_2_average * 1e-12 << " s"
              << std::endl;
    std::cout << "Board 0-3: " << result.board_0_3_average * 1e-12 << " s"
              << std::endl;

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

    Long64_t mean_0_1, mean_0_2, mean_0_3;
    mean_0_1 = res.board_0_1_average;
    mean_0_2 = res.board_0_2_average;
    mean_0_3 = res.board_0_3_average;

    TFile *input_output_file = new TFile(input_filepath, "UPDATE");

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

    Long64_t shifted_timestamp_calc;

    Int_t n_entries = input_tree->GetEntries();
    input_tree->LoadBaskets();
    for (Int_t j = 0; j < n_entries; j++) {
      input_tree->GetEntry(j);
      switch (board) {
      case 0:
        shifted_timestamp_calc = timestamp;
        break;
      case 1:
        shifted_timestamp_calc = timestamp + mean_0_1;
        break;
      case 2:
        shifted_timestamp_calc = timestamp + mean_0_2;
        break;
      case 3:
        shifted_timestamp_calc = timestamp + mean_0_3;
        break;
      }
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

    TFile *input_file = new TFile(input_filepath, "READ");
    TTree *input_tree = (TTree *)input_file->Get("Data_R");

    Long64_t n_entries = input_tree->GetEntries();
    std::cout << "Total entries: " << n_entries << std::endl;

    ULong64_t shifted_timestamp;
    input_tree->SetBranchAddress("ShiftedTimestamp", &shifted_timestamp);

    std::vector<std::pair<ULong64_t, Long64_t>> timestamp_index;
    timestamp_index.reserve(n_entries);

    std::cout << "Reading timestamps" << std::endl;
    for (Long64_t j = 0; j < n_entries; j++) {
      input_tree->GetEntry(j);
      timestamp_index.push_back({shifted_timestamp, j});

      if (j % 10000000 == 0) {
        std::cout << "Progress: " << j << "/" << n_entries << std::endl;
      }
    }

    std::cout << "Sorting" << std::endl;
    std::sort(timestamp_index.begin(), timestamp_index.end());

    std::cout << "Writing sorted tree" << std::endl;

    TString output_name = output_names[i];
    TString output_filepath = "root_files/" + output_name + ".root";
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

  InitUtils::SetROOTPreferences();

  std::vector<TString> filepaths, sorted_output_names;

  TString path_prefix = "./root_files/";
  for (Int_t i = 0; i < Constants::N_FILES; i++) {
    filepaths.push_back(Form("%sDataR_run_37_%d.root", path_prefix.Data(), i));
    std::cout << "Will process: " << filepaths[i] << std::endl;
    sorted_output_names.push_back(Form("Sorted_Run37_%d", i));
  }

  std::vector<UShort_t> channels_to_use = {
      8, // Board 0: L1
      0, // Board 1: L3
      0, // Board 2: R4
      0  // Board 3: L5
  };
  UShort_t ref_board = 0;
  UShort_t ref_channel = channels_to_use[0];

  Double_t min_energy = 0;
  Double_t max_energy = 20000;
  Double_t overlap_tmin_s = 1;
  Double_t overlap_tmax_s = 75.0;
  Double_t thresh_dt_us = 150.0;

  std::vector<TimeShiftResult> results = CalcTimeShiftsBeamMethod(
      filepaths, ref_board, ref_channel, channels_to_use, min_energy,
      max_energy, overlap_tmin_s, overlap_tmax_s, thresh_dt_us,
      reprocess_calculation);
  ApplyTimeShift(filepaths, results, reprocess_calculation);

  TimesortData(filepaths, sorted_output_names, reprocess_sorting);

  std::cout << "TIMING ANALYSIS COMPLETE" << std::endl;
}
