#include "DeltaEScatter.hpp"

const Int_t kStripA1 = 10;
const Int_t kStripA2 = 11;
const Int_t kStripB1 = 7;
const Int_t kStripB2 = 8;
const Int_t kStripB3 = 9;
const Double_t kAxisMin = 3.0;
const Double_t kAxisMax = 10.0;

// Pass 2: build the scatter graph, gating each event on the bigaus ellipse.
// Stride-samples across the full chain so we cover the entire run rather than
// just its head; the stride is chosen so total visited entries land near
// max_points (gated count will be lower).
TGraph *BuildGatedScatter(TChain *chain, const BeamFit2D &gate, Int_t color,
                          Long64_t max_points) {
  Float_t leftdE[18] = {0};
  Float_t rightdE[18] = {0};
  Float_t totaldE[18] = {0};
  chain->SetBranchStatus("*", 0);
  chain->SetBranchStatus("LeftdEMeV", 1);
  chain->SetBranchStatus("RightdEMeV", 1);
  chain->SetBranchStatus("TotaldEMeV", 1);
  chain->SetBranchAddress("LeftdEMeV", leftdE);
  chain->SetBranchAddress("RightdEMeV", rightdE);
  chain->SetBranchAddress("TotaldEMeV", totaldE);

  Long64_t n_total = chain->GetEntries();
  Long64_t stride = FileSet::SampleStride(n_total, max_points);

  TGraph *g = new TGraph();
  Long64_t k = 0;
  for (Long64_t j = 0; j < n_total; j += stride) {
    chain->GetEntry(j);
    Double_t xg = Double_t(leftdE[1]);
    Double_t yg = Double_t(rightdE[2]);
    if (!(xg > 0 && yg > 0))
      continue;
    if (!BeamFitUtils::InEllipse(gate, xg, yg, BeamFitUtils::GateNSigma()))
      continue;
    Double_t x = Double_t(totaldE[kStripA1]) + Double_t(totaldE[kStripA2]);
    Double_t y = Double_t(totaldE[kStripB1]) + Double_t(totaldE[kStripB2]) +
                 Double_t(totaldE[kStripB3]);
    g->SetPoint(k++, x, y);
  }

  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.3);
  g->SetMarkerColorAlpha(color, 0.35);
  g->SetLineColor(color);
  return g;
}

// Same gating as BuildGatedScatter, but plots the strip B sum (7+8+9) against
// the calibrated cathode energy.
TGraph *BuildGatedCathodeScatter(TChain *chain, const BeamFit2D &gate,
                                 Int_t color, Long64_t max_points) {
  Float_t leftdE[18] = {0};
  Float_t rightdE[18] = {0};
  Float_t totaldE[18] = {0};
  Float_t cathode = 0;
  chain->SetBranchStatus("*", 0);
  chain->SetBranchStatus("LeftdEMeV", 1);
  chain->SetBranchStatus("RightdEMeV", 1);
  chain->SetBranchStatus("TotaldEMeV", 1);
  chain->SetBranchStatus("CathodeMeV", 1);
  chain->SetBranchAddress("LeftdEMeV", leftdE);
  chain->SetBranchAddress("RightdEMeV", rightdE);
  chain->SetBranchAddress("TotaldEMeV", totaldE);
  chain->SetBranchAddress("CathodeMeV", &cathode);

  Long64_t n_total = chain->GetEntries();
  Long64_t stride = FileSet::SampleStride(n_total, max_points);

  TGraph *g = new TGraph();
  Long64_t k = 0;
  for (Long64_t j = 0; j < n_total; j += stride) {
    chain->GetEntry(j);
    Double_t xg = Double_t(leftdE[1]);
    Double_t yg = Double_t(rightdE[2]);
    if (!(xg > 0 && yg > 0))
      continue;
    if (!BeamFitUtils::InEllipse(gate, xg, yg, BeamFitUtils::GateNSigma()))
      continue;
    Double_t x = Double_t(cathode);
    Double_t y = Double_t(totaldE[kStripB1]) + Double_t(totaldE[kStripB2]) +
                 Double_t(totaldE[kStripB3]);
    g->SetPoint(k++, x, y);
  }

  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.3);
  g->SetMarkerColorAlpha(color, 0.35);
  g->SetLineColor(color);
  return g;
}

void SaveScatterPlot(TGraph *graph, const TString &legend_text,
                     const TString &output_name, const TString &output_subdir) {
  std::lock_guard<std::mutex> lock(g_plot_mutex);

  TString frame_name = Form("frame_%s", output_name.Data());
  TString frame_title =
      Form(";#DeltaE strip %d+%d [MeV];#DeltaE strip %d+%d+%d [MeV]", kStripA1,
           kStripA2, kStripB1, kStripB2, kStripB3);
  TH2D *frame = new TH2D(frame_name, frame_title, 10, kAxisMin, kAxisMax, 10,
                         kAxisMin, kAxisMax);
  frame->SetStats(0);
  frame->GetXaxis()->SetTitleSize(0.06);
  frame->GetYaxis()->SetTitleSize(0.06);
  frame->GetXaxis()->SetLabelSize(0.06);
  frame->GetYaxis()->SetLabelSize(0.06);
  frame->GetXaxis()->SetTitleOffset(1.2);
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->GetXaxis()->SetNdivisions(506);
  frame->GetYaxis()->SetNdivisions(506);

  TGraph *proxy = new TGraph(1);
  proxy->SetPoint(0, -1e9, -1e9);
  proxy->SetMarkerStyle(20);
  proxy->SetMarkerSize(1.4);
  proxy->SetMarkerColor(graph->GetMarkerColor());

  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  frame->Draw();
  graph->Draw("P SAME");

  TLegend *legend = PlottingUtils::AddLegend(0.16, 0.40, 0.16, 0.24);
  legend->AddEntry(proxy, legend_text, "p");
  legend->Draw();

  PlottingUtils::SaveFigure(c, output_name, output_subdir,
                            PlotSaveOptions::kLINEAR);

  delete legend;
  delete c;
  delete frame;
  delete proxy;
}

// Cathode energy scale is run-dependent, so the frame range is taken from the
// graph extents (with a small margin) rather than the fixed kAxis* bounds.
void SaveCathodeScatterPlot(TGraph *graph, const TString &legend_text,
                            const TString &output_name,
                            const TString &output_subdir) {
  std::lock_guard<std::mutex> lock(g_plot_mutex);

  Double_t x_min = 0.0, x_max = 1.0, y_min = kAxisMin, y_max = kAxisMax;
  if (graph->GetN() > 0) {
    x_min = TMath::MinElement(graph->GetN(), graph->GetX());
    x_max = TMath::MaxElement(graph->GetN(), graph->GetX());
    y_min = TMath::MinElement(graph->GetN(), graph->GetY());
    y_max = TMath::MaxElement(graph->GetN(), graph->GetY());
    Double_t mx = 0.05 * (x_max - x_min);
    Double_t my = 0.05 * (y_max - y_min);
    x_min -= mx;
    x_max += mx;
    y_min -= my;
    y_max += my;
  }

  TString frame_name = Form("frame_%s", output_name.Data());
  TString frame_title =
      Form(";Cathode #DeltaE [MeV];#DeltaE strip %d+%d+%d [MeV]", kStripB1,
           kStripB2, kStripB3);
  TH2D *frame =
      new TH2D(frame_name, frame_title, 10, x_min, x_max, 10, y_min, y_max);
  frame->SetStats(0);
  frame->GetXaxis()->SetTitleSize(0.06);
  frame->GetYaxis()->SetTitleSize(0.06);
  frame->GetXaxis()->SetLabelSize(0.06);
  frame->GetYaxis()->SetLabelSize(0.06);
  frame->GetXaxis()->SetTitleOffset(1.2);
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->GetXaxis()->SetNdivisions(506);
  frame->GetYaxis()->SetNdivisions(506);

  TGraph *proxy = new TGraph(1);
  proxy->SetPoint(0, -1e9, -1e9);
  proxy->SetMarkerStyle(20);
  proxy->SetMarkerSize(1.4);
  proxy->SetMarkerColor(graph->GetMarkerColor());

  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  frame->Draw();
  graph->Draw("P SAME");

  TLegend *legend = PlottingUtils::AddLegend(0.16, 0.40, 0.16, 0.24);
  legend->AddEntry(proxy, legend_text, "p");
  legend->Draw();

  PlottingUtils::SaveFigure(c, output_name, output_subdir,
                            PlotSaveOptions::kLINEAR);

  delete legend;
  delete c;
  delete frame;
  delete proxy;
}

void SaveGatePlot(TH2F *gate_hist, const BeamFit2D &gate,
                  const TString &output_name, const TString &output_subdir) {
  std::lock_guard<std::mutex> lock(g_plot_mutex);
  TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureAndDraw2DHistogram(gate_hist, cv);
  BeamFitUtils::DrawGateEllipse(gate, BeamFitUtils::GateNSigma());
  PlottingUtils::SaveFigure(cv, output_name, output_subdir,
                            PlotSaveOptions::kLINEAR);
  delete cv;
}

void ProcessRun(Int_t run, TChain *chain, const TString &output_subdir) {
  if (!chain || chain->GetEntries() == 0) {
    std::cerr << "Run " << run << ": empty chain; skipping" << std::endl;
    return;
  }
  std::cout << "Run " << run << ": " << chain->GetEntries() << " events"
            << std::endl;

  TString tag = Form("run%d", run);
  TString gate_name = Form("h2_gate_S0_L1_%s", tag.Data());
  TH2F *gate_hist = BeamFitUtils::BuildGateHist(chain, gate_name);
  BeamFit2D gate = BeamFitUtils::FitBigPeak(gate_hist, tag);
  SaveGatePlot(gate_hist, gate, "gate_S0_vs_L1", output_subdir);

  if (!gate.ok) {
    std::cerr << "Run " << run << ": gate fit failed; skipping scatter"
              << std::endl;
    delete gate_hist;
    return;
  }

  TGraph *graph =
      BuildGatedScatter(chain, gate, kAzure + 2, /*max_points=*/25000);
  std::cout << "  gated scatter: " << graph->GetN() << " points" << std::endl;
  TString legend_text = Form("Run %d (S0 vs L1 gate)", run);
  SaveScatterPlot(graph, legend_text, "scatter_s1112_vs_s8910_cal",
                  output_subdir);

  TGraph *cathode_graph =
      BuildGatedCathodeScatter(chain, gate, kAzure + 2, /*max_points=*/25000);
  std::cout << "  gated cathode scatter: " << cathode_graph->GetN() << " points"
            << std::endl;
  SaveCathodeScatterPlot(cathode_graph, legend_text,
                         "scatter_s789_vs_cathode_cal", output_subdir);

  delete cathode_graph;
  delete graph;
  delete gate_hist;
}

void DeltaEScatter::Run() {
  const TString project_root = Paths::DatasetDir();
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");

  // Group cal sidecars by run; each run becomes one TChain that drives a
  // standalone gate fit and scatter plot.
  std::vector<Int_t> run_order;
  std::map<Int_t, TChain *> chain_by_run =
      FileSet::GroupCalSidecarsByRun(run_order);

  for (std::size_t i = 0; i < run_order.size(); i++) {
    Int_t run = run_order[i];
    TString subdir = Form("scatter/run%d", run);
    ProcessRun(run, chain_by_run[run], subdir);
  }

  for (std::size_t i = 0; i < run_order.size(); i++)
    delete chain_by_run[run_order[i]];
}
