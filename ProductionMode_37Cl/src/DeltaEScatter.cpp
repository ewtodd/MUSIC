#include "DeltaEScatter.hpp"

const Int_t kStripA1 = 10;
const Int_t kStripA2 = 11;
const Int_t kStripB1 = 7;
const Int_t kStripB2 = 8;
const Int_t kStripB3 = 9;
const Double_t kAxisMin = 3.0;
const Double_t kAxisMax = 10.0;

// Gate: 2D big-peak fit on Strip0 (S, TotaldEMeV[0]) vs Strip1 long-side
// (L, LeftdEMeV[1]). Mirrors the bigaus + 2-sigma ellipse pattern from
// CalibrateBeam's k=1 fit.
const Double_t kGateMin = 0.0;
const Double_t kGateMax = 12.0;
const Int_t kGateBins = 240;
const Double_t kGateNSigma = 2.0;
const Int_t kSeedHalfBins = 40;
const Double_t kSeedFrac = 0.30;

// Moments-seeded bigaus fit on the single largest peak. Same shape as the
// k=1 path in CalibrateBeam's FitBeamMultiplicities.
BeamFit2D FitBigPeak(TH2F *h, const TString &tag) {
  BeamFit2D out;
  if (!h || h->GetEntries() < 100)
    return out;

  Double_t bw_x = h->GetXaxis()->GetBinWidth(1);
  Double_t bw_y = h->GetYaxis()->GetBinWidth(1);

  Int_t bx = 0, by = 0, bz = 0;
  h->GetMaximumBin(bx, by, bz);
  Double_t peak_val = h->GetBinContent(bx, by);
  if (peak_val <= 0)
    return out;

  Int_t lo_bx = std::max(1, bx - kSeedHalfBins);
  Int_t hi_bx = std::min(h->GetNbinsX(), bx + kSeedHalfBins);
  Int_t lo_by = std::max(1, by - kSeedHalfBins);
  Int_t hi_by = std::min(h->GetNbinsY(), by + kSeedHalfBins);

  Moments2D m = BeamFitUtils::ComputeMoments(h, lo_bx, hi_bx, lo_by, hi_by,
                                             kSeedFrac * peak_val, bw_x, bw_y);
  if (m.weight <= 0) {
    std::cerr << "  FitBigPeak " << tag << ": no bins above seed threshold"
              << std::endl;
    return out;
  }

  Double_t x_lo = std::max(h->GetXaxis()->GetXmin(), m.mu_x - 3.0 * m.sigma_x);
  Double_t x_hi = std::min(h->GetXaxis()->GetXmax(), m.mu_x + 3.0 * m.sigma_x);
  Double_t y_lo = std::max(h->GetYaxis()->GetXmin(), m.mu_y - 3.0 * m.sigma_y);
  Double_t y_hi = std::min(h->GetYaxis()->GetXmax(), m.mu_y + 3.0 * m.sigma_y);

  out = BeamFitUtils::FitBigausInWindow(h, x_lo, x_hi, y_lo, y_hi, m, peak_val,
                                        Form("f_gate_%s", tag.Data()),
                                        TString("FitBigPeak ") + tag);

  std::cout << "  gate fit " << tag << ": mu=(" << out.mu_x << "," << out.mu_y
            << ") sigma=(" << out.sigma_x << "," << out.sigma_y
            << ") rho=" << out.rho << std::endl;
  return out;
}

// Pass 1: build the gate H2 (S0 vs L1) for the chain.
TH2F *BuildGateHist(TChain *chain, const TString &name) {
  Float_t leftdE[18] = {0};
  Float_t rightdE[18] = {0};
  chain->SetBranchStatus("*", 0);
  chain->SetBranchStatus("LeftdEMeV", 1);
  chain->SetBranchStatus("RightdEMeV", 1);
  chain->SetBranchAddress("LeftdEMeV", leftdE);
  chain->SetBranchAddress("RightdEMeV", rightdE);

  TH2F *h = new TH2F(name, ";#DeltaE S0 [MeV];#DeltaE L1 [MeV]", kGateBins,
                     kGateMin, kGateMax, kGateBins, kGateMin, kGateMax);
  h->SetDirectory(nullptr);

  Long64_t n = chain->GetEntries();
  for (Long64_t j = 0; j < n; j++) {
    chain->GetEntry(j);
    Double_t x = Double_t(leftdE[1]);
    Double_t y = Double_t(rightdE[2]);
    if (x > 0 && y > 0)
      h->Fill(x, y);
  }
  return h;
}

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
  Long64_t n_visit =
      (max_points > 0 && n_total > max_points) ? max_points : n_total;
  Long64_t stride = (n_visit > 0) ? (n_total / n_visit) : 1;
  if (stride < 1)
    stride = 1;

  TGraph *g = new TGraph();
  Long64_t k = 0;
  for (Long64_t j = 0; j < n_total; j += stride) {
    chain->GetEntry(j);
    Double_t xg = Double_t(leftdE[1]);
    Double_t yg = Double_t(rightdE[2]);
    if (!(xg > 0 && yg > 0))
      continue;
    if (!BeamFitUtils::InEllipse(gate, xg, yg, kGateNSigma))
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
// the raw cathode ADC. The cathode is kept in ADC (no calibration applied);
// it lives in the events tree, not the cal sidecar, so events_chain is the
// parallel events TChain aligned entry-for-entry with chain.
TGraph *BuildGatedCathodeScatter(TChain *chain, TChain *events_chain,
                                 const BeamFit2D &gate, Int_t color,
                                 Long64_t max_points) {
  Float_t leftdE[18] = {0};
  Float_t rightdE[18] = {0};
  Float_t totaldE[18] = {0};
  Int_t cathode_adc = 0;
  chain->SetBranchStatus("*", 0);
  chain->SetBranchStatus("LeftdEMeV", 1);
  chain->SetBranchStatus("RightdEMeV", 1);
  chain->SetBranchStatus("TotaldEMeV", 1);
  chain->SetBranchAddress("LeftdEMeV", leftdE);
  chain->SetBranchAddress("RightdEMeV", rightdE);
  chain->SetBranchAddress("TotaldEMeV", totaldE);
  events_chain->SetBranchStatus("*", 0);
  events_chain->SetBranchStatus("Cathode", 1);
  events_chain->SetBranchAddress("Cathode", &cathode_adc);

  Long64_t n_total = chain->GetEntries();
  Long64_t n_visit =
      (max_points > 0 && n_total > max_points) ? max_points : n_total;
  Long64_t stride = (n_visit > 0) ? (n_total / n_visit) : 1;
  if (stride < 1)
    stride = 1;

  TGraph *g = new TGraph();
  Long64_t k = 0;
  for (Long64_t j = 0; j < n_total; j += stride) {
    chain->GetEntry(j);
    Double_t xg = Double_t(leftdE[1]);
    Double_t yg = Double_t(rightdE[2]);
    if (!(xg > 0 && yg > 0))
      continue;
    if (!BeamFitUtils::InEllipse(gate, xg, yg, kGateNSigma))
      continue;
    events_chain->GetEntry(j);
    if (cathode_adc < 0)
      continue;
    Double_t x = Double_t(cathode_adc);
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

// Cathode ADC range is run-dependent, so the frame range is taken from the
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
      Form(";Cathode #DeltaE [ADC];#DeltaE strip %d+%d+%d [MeV]", kStripB1,
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
  if (gate.ok) {
    Double_t major, minor, theta_deg;
    BeamFitUtils::DiagonalizeCov(gate, major, minor, theta_deg);
    TEllipse *e = new TEllipse(gate.mu_x, gate.mu_y, kGateNSigma * major,
                               kGateNSigma * minor, 0.0, 360.0, theta_deg);
    e->SetFillStyle(0);
    e->SetLineColor(kRed + 1);
    e->SetLineWidth(2);
    e->Draw();
  }
  PlottingUtils::SaveFigure(cv, output_name, output_subdir,
                            PlotSaveOptions::kLINEAR);
  delete cv;
}

void ProcessRun(Int_t run, TChain *chain, TChain *events_chain,
                const TString &output_subdir) {
  if (!chain || chain->GetEntries() == 0) {
    std::cerr << "Run " << run << ": empty chain; skipping" << std::endl;
    return;
  }
  std::cout << "Run " << run << ": " << chain->GetEntries() << " events"
            << std::endl;

  TString tag = Form("run%d", run);
  TString gate_name = Form("h2_gate_S0_L1_%s", tag.Data());
  TH2F *gate_hist = BuildGateHist(chain, gate_name);
  BeamFit2D gate = FitBigPeak(gate_hist, tag);
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

  TGraph *cathode_graph = BuildGatedCathodeScatter(
      chain, events_chain, gate, kAzure + 2, /*max_points=*/25000);
  std::cout << "  gated cathode scatter: " << cathode_graph->GetN() << " points"
            << std::endl;
  SaveCathodeScatterPlot(cathode_graph, legend_text,
                         "scatter_s789_vs_cathode_adc", output_subdir);

  delete cathode_graph;
  delete graph;
  delete gate_hist;
}

void DeltaEScatter::Run() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");

  // Group cal sidecars by run; each run becomes one TChain that drives a
  // standalone gate fit and scatter plot.
  // The cal sidecar (events_cal) drives the gate and MeV strip sums; the raw
  // cathode ADC lives in the events tree, so build a parallel events TChain
  // per run with the same file add order. They align entry-for-entry because
  // each sidecar is built per-event from its events file.
  std::vector<FileSpec> all_specs = FileSet::BuildFileSpecs();
  std::map<Int_t, TChain *> chain_by_run;
  std::map<Int_t, TChain *> events_by_run;
  std::vector<Int_t> run_order;
  for (std::size_t i = 0; i < all_specs.size(); i++) {
    const FileSpec &s = all_specs[i];
    TString cal_full =
        IO::GetRootFilesBaseDir() + "/" + FileSet::CalSidecarName(s);
    TString events_full =
        IO::GetRootFilesBaseDir() + "/" + FileSet::EventsName(s) + ".root";
    if (gSystem->AccessPathName(cal_full)) {
      std::cerr << "Missing cal sidecar: " << cal_full << std::endl;
      continue;
    }
    if (gSystem->AccessPathName(events_full)) {
      std::cerr << "Missing events file: " << events_full << std::endl;
      continue;
    }
    if (chain_by_run.find(s.run) == chain_by_run.end()) {
      chain_by_run[s.run] = new TChain("events_cal");
      events_by_run[s.run] = new TChain("events");
      run_order.push_back(s.run);
    }
    chain_by_run[s.run]->Add(cal_full);
    events_by_run[s.run]->Add(events_full);
  }

  for (std::size_t i = 0; i < run_order.size(); i++) {
    Int_t run = run_order[i];
    TString subdir = Form("scatter/run%d", run);
    ProcessRun(run, chain_by_run[run], events_by_run[run], subdir);
  }

  for (std::size_t i = 0; i < run_order.size(); i++) {
    delete chain_by_run[run_order[i]];
    delete events_by_run[run_order[i]];
  }
}
