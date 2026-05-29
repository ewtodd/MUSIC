#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Pipeline.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TTree.h>
#include <iostream>
#include <mutex>
#include <vector>

namespace {

const Int_t kStripA1 = 10;
const Int_t kStripA2 = 11;
const Int_t kStripB1 = 7;
const Int_t kStripB2 = 8;
const Int_t kStripB3 = 9;
const Double_t kAxisMin = 3.0;
const Double_t kAxisMax = 10.0;

const Double_t kTraceYMin = 1.0;
const Double_t kTraceYMax = 4.0;

const Double_t kTotalSumXMin = 3.0;
const Double_t kTotalSumXMax = 10.0;
const Double_t kTotalSumYMin = 0.0;
const Double_t kTotalSumYMax = 50.0;

typedef TGraph *(*ScatterBuilderFn)(TTree *, Int_t, Long64_t);

TGraph *BuildScatter(TTree *tree, Int_t color, Long64_t max_points) {
  // Accept either Remix sim ("TotaldE", MeV) or ProductionMode calibrated
  // ("TotaldEMeV") so this macro is drop-in usable in both projects.
  const char *total_name =
      tree->GetBranch("TotaldE") ? "TotaldE" : "TotaldEMeV";

  Float_t totaldE[18] = {0};
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus(total_name, 1);
  tree->SetBranchAddress(total_name, totaldE);

  Long64_t n_total = tree->GetEntries();
  Long64_t n_use =
      (max_points > 0 && n_total > max_points) ? max_points : n_total;
  // Stride-sample so we cover the full file rather than just its head, which
  // would bias when files were written in time order.
  Long64_t stride = (n_use > 0) ? (n_total / n_use) : 1;
  if (stride < 1)
    stride = 1;

  TGraph *g = new TGraph(n_use);
  Long64_t k = 0;
  for (Long64_t j = 0; j < n_total && k < n_use; j += stride) {
    tree->GetEntry(j);
    Double_t x = Double_t(totaldE[kStripA1]) + Double_t(totaldE[kStripA2]);
    Double_t y = Double_t(totaldE[kStripB1]) + Double_t(totaldE[kStripB2]) +
                 Double_t(totaldE[kStripB3]);
    g->SetPoint(k++, x, y);
  }
  g->Set(k);

  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.3);
  g->SetMarkerColorAlpha(color, 0.35);
  g->SetLineColor(color);
  return g;
}

TGraph *BuildTotalSumScatter(TTree *tree, Int_t color, Long64_t max_points) {
  const char *total_name =
      tree->GetBranch("TotaldE") ? "TotaldE" : "TotaldEMeV";

  Float_t totaldE[18] = {0};
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus(total_name, 1);
  tree->SetBranchAddress(total_name, totaldE);

  Long64_t n_total = tree->GetEntries();
  Long64_t n_use =
      (max_points > 0 && n_total > max_points) ? max_points : n_total;
  Long64_t stride = (n_use > 0) ? (n_total / n_use) : 1;
  if (stride < 1)
    stride = 1;

  const Int_t s_lo = Constants::FirstStrip();
  const Int_t s_hi = Constants::LastStrip();

  TGraph *g = new TGraph(n_use);
  Long64_t k = 0;
  for (Long64_t j = 0; j < n_total && k < n_use; j += stride) {
    tree->GetEntry(j);
    Double_t x = Double_t(totaldE[kStripB1]) + Double_t(totaldE[kStripB2]) +
                 Double_t(totaldE[kStripB3]);
    Double_t y = 0.0;
    for (Int_t s = s_lo; s <= s_hi; s++)
      y += Double_t(totaldE[s]);
    g->SetPoint(k++, x, y);
  }
  g->Set(k);

  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.3);
  g->SetMarkerColorAlpha(color, 0.35);
  g->SetLineColor(color);
  return g;
}

void MakeOverlayScatter(const std::vector<TString> &filepaths,
                        const std::vector<TString> &labels,
                        const TString &output_name,
                        const TString &output_subdir, ScatterBuilderFn builder,
                        const TString &frame_title, Double_t xmin,
                        Double_t xmax, Double_t ymin, Double_t ymax,
                        Long64_t max_points = 25000) {
  if (filepaths.size() != labels.size() || filepaths.empty()) {
    std::cerr << "MakeOverlayScatter: bad inputs" << std::endl;
    return;
  }

  std::vector<Int_t> palette = PlottingUtils::GetDefaultColors();
  std::vector<TFile *> files;
  std::vector<TGraph *> graphs;
  std::vector<TString> kept_labels;

  for (std::size_t i = 0; i < filepaths.size(); i++) {
    TFile *f = IO::OpenForReading(filepaths[i]);
    if (!f || f->IsZombie()) {
      std::cerr << "Error opening: " << filepaths[i] << std::endl;
      continue;
    }
    TTree *tree = static_cast<TTree *>(f->Get("events_MeV"));
    if (!tree)
      tree = static_cast<TTree *>(f->Get("events_cal"));
    if (!tree) {
      std::cerr << "No events tree in: " << filepaths[i] << std::endl;
      f->Close();
      continue;
    }

    Int_t color = palette[i % palette.size()];
    std::cout << "Building scatter for " << labels[i] << " ("
              << tree->GetEntries() << " entries)..." << std::endl;
    graphs.push_back(builder(tree, color, max_points));
    kept_labels.push_back(labels[i]);
    files.push_back(f);
  }
  if (graphs.empty())
    return;

  // Legend uses solid-color proxy markers (the real graphs are alpha-blended
  // so points overlap legibly, which makes the legend swatches too faint).
  std::vector<TGraph *> legend_proxies;
  for (std::size_t i = 0; i < graphs.size(); i++) {
    Int_t color = palette[i % palette.size()];
    TGraph *proxy = new TGraph(1);
    proxy->SetPoint(0, -1e9, -1e9);
    proxy->SetMarkerStyle(20);
    proxy->SetMarkerSize(1.4);
    proxy->SetMarkerColor(color);
    legend_proxies.push_back(proxy);
  }

  std::lock_guard<std::mutex> lock(g_plot_mutex);

  TString frame_name = Form("frame_%s", output_name.Data());
  TH2D *frame =
      new TH2D(frame_name, frame_title, 10, xmin, xmax, 10, ymin, ymax);
  frame->SetStats(0);
  frame->GetXaxis()->SetTitleSize(0.06);
  frame->GetYaxis()->SetTitleSize(0.06);
  frame->GetXaxis()->SetLabelSize(0.06);
  frame->GetYaxis()->SetLabelSize(0.06);
  frame->GetXaxis()->SetTitleOffset(1.2);
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->GetXaxis()->SetNdivisions(506);
  frame->GetYaxis()->SetNdivisions(506);

  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  frame->Draw();
  for (std::size_t i = 0; i < graphs.size(); i++)
    graphs[i]->Draw("P SAME");

  TLegend *legend = PlottingUtils::AddLegend(0.16, 0.40, 0.16, 0.32);
  for (std::size_t i = 0; i < graphs.size(); i++)
    legend->AddEntry(legend_proxies[i], kept_labels[i], "p");
  legend->Draw();

  PlottingUtils::SaveFigure(c, output_name, output_subdir,
                            PlotSaveOptions::kLINEAR);

  delete legend;
  delete c;
  delete frame;
  for (std::size_t i = 0; i < graphs.size(); i++)
    delete graphs[i];
  for (std::size_t i = 0; i < legend_proxies.size(); i++)
    delete legend_proxies[i];
  for (std::size_t i = 0; i < files.size(); i++)
    files[i]->Close();
}

std::vector<TGraph *> BuildTracesForDataset(TTree *tree, Int_t color,
                                            Int_t n_traces) {
  const char *total_name =
      tree->GetBranch("TotaldE") ? "TotaldE" : "TotaldEMeV";

  Float_t totaldE[18] = {0};
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus(total_name, 1);
  tree->SetBranchAddress(total_name, totaldE);

  Long64_t n_total = tree->GetEntries();
  Long64_t stride =
      (n_traces > 0 && n_total > Long64_t(n_traces)) ? (n_total / n_traces) : 1;
  if (stride < 1)
    stride = 1;

  const Int_t s_lo = Constants::FirstStrip();
  const Int_t n_active = Constants::NumActiveStrips();

  std::vector<TGraph *> traces;
  for (Int_t i = 0; i < n_traces; i++) {
    Long64_t j = Long64_t(i) * stride;
    if (j >= n_total)
      break;
    tree->GetEntry(j);
    TGraph *g = new TGraph(n_active);
    for (Int_t k = 0; k < n_active; k++) {
      Int_t s = s_lo + k;
      g->SetPoint(k, s, Double_t(totaldE[s]));
    }
    g->SetLineColorAlpha(color, 0.5);
    g->SetMarkerColorAlpha(color, 0.5);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.4);
    g->SetLineWidth(2);
    traces.push_back(g);
  }
  return traces;
}

void MakeTraceOverlay(const std::vector<TString> &filepaths,
                      const std::vector<TString> &labels,
                      const TString &output_name, const TString &output_subdir,
                      Int_t n_per_dataset = 1000) {
  if (filepaths.size() != labels.size() || filepaths.empty()) {
    std::cerr << "MakeTraceOverlay: bad inputs" << std::endl;
    return;
  }

  std::vector<Int_t> palette = PlottingUtils::GetDefaultColors();
  std::vector<TFile *> files;
  std::vector<std::vector<TGraph *>> traces_per_dataset;
  std::vector<TString> kept_labels;
  std::vector<Int_t> kept_colors;

  for (std::size_t i = 0; i < filepaths.size(); i++) {
    TFile *f = IO::OpenForReading(filepaths[i]);
    if (!f || f->IsZombie()) {
      std::cerr << "Error opening: " << filepaths[i] << std::endl;
      continue;
    }
    TTree *tree = static_cast<TTree *>(f->Get("events_MeV"));
    if (!tree)
      tree = static_cast<TTree *>(f->Get("events_cal"));
    if (!tree) {
      std::cerr << "No events tree in: " << filepaths[i] << std::endl;
      f->Close();
      continue;
    }

    Int_t color = palette[i % palette.size()];
    std::cout << "Building " << n_per_dataset << " traces for " << labels[i]
              << "..." << std::endl;
    traces_per_dataset.push_back(
        BuildTracesForDataset(tree, color, n_per_dataset));
    kept_labels.push_back(labels[i]);
    kept_colors.push_back(color);
    files.push_back(f);
  }
  if (traces_per_dataset.empty())
    return;

  std::vector<TGraph *> legend_proxies;
  for (std::size_t i = 0; i < kept_colors.size(); i++) {
    TGraph *proxy = new TGraph(1);
    proxy->SetPoint(0, -1e9, -1e9);
    proxy->SetMarkerStyle(20);
    proxy->SetMarkerSize(1.4);
    proxy->SetMarkerColor(kept_colors[i]);
    proxy->SetLineColor(kept_colors[i]);
    proxy->SetLineWidth(2);
    legend_proxies.push_back(proxy);
  }

  std::lock_guard<std::mutex> lock(g_plot_mutex);

  const Int_t s_lo = Constants::FirstStrip();
  const Int_t s_hi = Constants::LastStrip();

  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);

  Bool_t first = kTRUE;
  for (std::size_t i = 0; i < traces_per_dataset.size(); i++) {
    for (std::size_t j = 0; j < traces_per_dataset[i].size(); j++) {
      TGraph *g = traces_per_dataset[i][j];
      if (first) {
        g->SetTitle(";Strip Index;#DeltaE [MeV]");
        g->GetXaxis()->SetTitleSize(0.06);
        g->GetYaxis()->SetTitleSize(0.06);
        g->GetXaxis()->SetLabelSize(0.06);
        g->GetYaxis()->SetLabelSize(0.06);
        g->GetXaxis()->SetTitleOffset(1.2);
        g->GetYaxis()->SetTitleOffset(1.2);
        g->GetXaxis()->SetNdivisions(s_hi - s_lo + 1);
        g->GetYaxis()->SetNdivisions(506);
        g->GetXaxis()->SetLimits(s_lo - 0.5, s_hi + 0.5);
        g->GetYaxis()->SetRangeUser(kTraceYMin, kTraceYMax);
        g->Draw("ALP");
        first = kFALSE;
      } else {
        g->Draw("LP SAME");
      }
    }
  }

  TLegend *legend = PlottingUtils::AddLegend(0.70, 0.90, 0.74, 0.90);
  for (std::size_t i = 0; i < legend_proxies.size(); i++)
    legend->AddEntry(legend_proxies[i], kept_labels[i], "lp");
  legend->Draw();

  PlottingUtils::SaveFigure(c, output_name, output_subdir,
                            PlotSaveOptions::kLINEAR);

  delete legend;
  delete c;
  for (std::size_t i = 0; i < traces_per_dataset.size(); i++)
    for (std::size_t j = 0; j < traces_per_dataset[i].size(); j++)
      delete traces_per_dataset[i][j];
  for (std::size_t i = 0; i < legend_proxies.size(); i++)
    delete legend_proxies[i];
  for (std::size_t i = 0; i < files.size(); i++)
    files[i]->Close();
}

TString PrettyLabel(const TString &tag) {
  TString base = tag;
  base.ReplaceAll("_eres", "");
  if (base == "aa")
    return "#alpha + #alpha";
  if (base == "an")
    return "#alpha + n";
  if (base == "beam")
    return "Beam (^{37}Cl)";
  return base;
}

} // namespace

void DeltaEScatter() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");

  std::vector<FileSpec> all_specs = BuildFileSpecs();
  std::vector<TString> truth_files, truth_labels;
  std::vector<TString> eres_files, eres_labels;
  std::vector<TString> truth_trace_files, truth_trace_labels;
  std::vector<TString> eres_trace_files, eres_trace_labels;
  for (std::size_t i = 0; i < all_specs.size(); i++) {
    const FileSpec &s = all_specs[i];
    TString filepath = TracesName(s) + ".root";
    TString label = PrettyLabel(s.tag);
    TString base = s.tag;
    base.ReplaceAll("_eres", "");
    Bool_t is_beam = (base == "beam");
    if (s.tag.EndsWith("_eres")) {
      eres_files.push_back(filepath);
      eres_labels.push_back(label);
      if (!is_beam) {
        eres_trace_files.push_back(filepath);
        eres_trace_labels.push_back(label);
      }
    } else {
      truth_files.push_back(filepath);
      truth_labels.push_back(label);
      if (!is_beam) {
        truth_trace_files.push_back(filepath);
        truth_trace_labels.push_back(label);
      }
    }
  }

  TString strip_pair_title =
      Form(";#DeltaE strip %d+%d [MeV];#DeltaE strip %d+%d+%d [MeV]", kStripA1,
           kStripA2, kStripB1, kStripB2, kStripB3);
  MakeOverlayScatter(truth_files, truth_labels, "scatter_s89_vs_s1011_truth",
                     "scatter", BuildScatter, strip_pair_title, kAxisMin,
                     kAxisMax, kAxisMin, kAxisMax);
  MakeOverlayScatter(eres_files, eres_labels, "scatter_s89_vs_s1011_eres",
                     "scatter", BuildScatter, strip_pair_title, kAxisMin,
                     kAxisMax, kAxisMin, kAxisMax);

  TString total_sum_title =
      Form(";#DeltaE strip %d+%d+%d [MeV];Total #DeltaE [MeV]", kStripB1,
           kStripB2, kStripB3);
  MakeOverlayScatter(truth_files, truth_labels, "scatter_total_vs_s8910_truth",
                     "scatter", BuildTotalSumScatter, total_sum_title,
                     kTotalSumXMin, kTotalSumXMax, kTotalSumYMin,
                     kTotalSumYMax);
  MakeOverlayScatter(eres_files, eres_labels, "scatter_total_vs_s8910_eres",
                     "scatter", BuildTotalSumScatter, total_sum_title,
                     kTotalSumXMin, kTotalSumXMax, kTotalSumYMin,
                     kTotalSumYMax);

  MakeTraceOverlay(truth_trace_files, truth_trace_labels, "trace_overlay_truth",
                   "trace_overlay");
  MakeTraceOverlay(eres_trace_files, eres_trace_labels, "trace_overlay_eres",
                   "trace_overlay");
}
