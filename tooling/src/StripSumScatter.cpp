#include "StripSumScatter.hpp"
#include "RemixSim.hpp"
#include <TKey.h>
#include <TNamed.h>
#include <limits>

void StripSumScatter::EnableEventBranches(TChain *chain) {
  chain->SetBranchStatus("*", 0);
  chain->SetBranchStatus("Left_0_17_dE", 1);
  chain->SetBranchStatus("RightdE", 1);
  chain->SetBranchStatus("Cathode", 1);
}

Bool_t StripSumScatter::AllStripsFired(const EnergyView &ev) {
  if (!(ev.total[0] > 0.0 && ev.total[17] > 0.0))
    return kFALSE;
  for (Int_t s = 1; s <= 16; s++)
    if (!(ev.total[s] > 0.0))
      return kFALSE;
  return kTRUE;
}

Bool_t StripSumScatter::PassesReaction(const EnergyView &ev, Int_t reac) {
  if (!AllStripsFired(ev))
    return kFALSE;
  for (Int_t s = 0; s < reac; s++)
    if (TMath::Abs(ev.total[s] - 1.0) > kBeamFlatTol)
      return kFALSE;
  Double_t reac_jump = ev.total[reac] - ev.total[reac - 1];
  if (!(reac_jump > kReacJumpMin && reac_jump < kReacJumpMax))
    return kFALSE;
  if (!(ev.total[reac] > 1.0 + kReacJumpMin &&
        ev.total[reac] < 1.0 + kReacJumpMax))
    return kFALSE;
  if (Constants::STRIP_SUM_APPLY_SMOOTHNESS)
    for (Int_t s = reac + 1; s <= kSmoothHiStrip; s++)
      if (TMath::Abs(ev.total[s] - ev.total[s - 1]) > kSmoothMaxStep)
        return kFALSE;
  return ev.total[17] < kStrip17Max;
}

// Clean un-reacted beam: all long anodes flat near the beam level (1 a.u.).
Bool_t StripSumScatter::IsBeamFlat(const EnergyView &ev) {
  if (!AllStripsFired(ev))
    return kFALSE;
  for (Int_t s = 1; s <= 16; s++)
    if (TMath::Abs(ev.total[s] - 1.0) > kBeamFlatTol)
      return kFALSE;
  return kTRUE;
}

Bool_t StripSumScatter::IsPileup(const EnergyView &ev) {
  const Double_t kThresh = 1.75;
  Int_t n = 0;
  for (Int_t s = 1; s <= 16; s++)
    if (ev.total[s] >= kThresh && ++n >= 2)
      return kTRUE;
  return kFALSE;
}

Bool_t StripSumScatter::IsNoise(const EnergyView &ev) {
  Int_t n = 0;
  for (Int_t s = 1; s <= 16; s++)
    if (ev.total[s] <= kLongStripMinValue && ++n >= 3)
      return kTRUE;
  return kFALSE;
}

Double_t StripSumScatter::SumRange(const Double_t *total, Int_t lo, Int_t hi) {
  Double_t sum = 0.0;
  for (Int_t s = lo; s <= hi; s++)
    sum += total[s];
  return sum;
}

// Mandatory strip1-vs-strip2 gate first, then any opt-in pairs (3/4, 5/6). An
// event must pass every active gate.
std::vector<StripSumScatter::GateSpec> StripSumScatter::ActiveGates() {
  std::vector<GateSpec> gates;
  GateSpec g;
  g.sx = kGateStripX;
  g.sy = kGateStripY;
  gates.push_back(g);
  if (Constants::STRIP_SUM_GATE_S3_S4) {
    g.sx = 3;
    g.sy = 4;
    gates.push_back(g);
  }
  if (Constants::STRIP_SUM_GATE_S5_S6) {
    g.sx = 5;
    g.sy = 6;
    gates.push_back(g);
  }
  return gates;
}

// Both extra gates off -> the original StripSumScatter_cache.root (existing
// caches stay valid); each enabled gate appends its own tag.
TString StripSumScatter::CacheName() {
  TString name = "StripSumScatter_cache";
  if (Constants::STRIP_SUM_GATE_S3_S4)
    name += "_g34";
  if (Constants::STRIP_SUM_GATE_S5_S6)
    name += "_g56";
  name += ".root";
  return name;
}

Bool_t StripSumScatter::PassesGate(const BeamFit2D &gate, const EnergyView &ev,
                                   Int_t sx, Int_t sy) {
  Double_t g0 = ev.total[sx];
  Double_t g1 = ev.total[sy];
  if (!(g0 > 0.0 && g1 > 0.0))
    return kFALSE;
  return BeamFitUtils::InEllipseXY(gate, g0, g1, kGateNSigmaX, kGateNSigmaY);
}

BeamFit2D
StripSumScatter::FindBeamGate(TChain *chain, Int_t sx, Int_t sy,
                              const std::vector<GateSpec> &prior_specs,
                              const std::vector<BeamFit2D> &prior_gates,
                              const TString &tag, const TString &subdir) {
  BeamFit2D out;
  EnergyView ev;
  ev.Attach(chain);
  EnableEventBranches(chain);
  TH2F *h =
      new TH2F(Form("h2_beamgate_s%d_s%d_%s", sx, sy, tag.Data()),
               Form(";#DeltaE strip %d [a.u.];#DeltaE strip %d [a.u.]", sx, sy),
               kGateBins, kGateMin, kGateMax, kGateBins, kGateMin, kGateMax);
  h->SetDirectory(nullptr);
  Long64_t n = chain->GetEntries();
  Long64_t stride = FileSet::SampleStride(n, kSampleMaxPoints);
  for (Long64_t j = 0; j < n; j += stride) {
    chain->GetEntry(j);
    ev.Decode();
    // Series gating: only events passing every prior gate feed this fit.
    Bool_t prior_ok = kTRUE;
    for (std::size_t gi = 0; gi < prior_specs.size(); gi++)
      if (!PassesGate(prior_gates[gi], ev, prior_specs[gi].sx,
                      prior_specs[gi].sy)) {
        prior_ok = kFALSE;
        break;
      }
    if (!prior_ok)
      continue;
    Double_t x = ev.total[sx];
    Double_t y = ev.total[sy];
    if (x > 0.0 && y > 0.0)
      h->Fill(x, y);
  }
  if (h->GetEntries() < 100) {
    delete h;
    return out;
  }
  Double_t bw_x = h->GetXaxis()->GetBinWidth(1);
  Double_t bw_y = h->GetYaxis()->GetBinWidth(1);
  Int_t bx = 0, by = 0, bz = 0;
  h->GetMaximumBin(bx, by, bz);
  Double_t peak_val = h->GetBinContent(bx, by);
  Int_t lo_bx = std::max(1, bx - kSeedHalfBins);
  Int_t hi_bx = std::min(h->GetNbinsX(), bx + kSeedHalfBins);
  Int_t lo_by = std::max(1, by - kSeedHalfBins);
  Int_t hi_by = std::min(h->GetNbinsY(), by + kSeedHalfBins);
  Moments2D m = BeamFitUtils::ComputeMoments(h, lo_bx, hi_bx, lo_by, hi_by,
                                             kSeedFrac * peak_val, bw_x, bw_y);
  if (m.weight <= 0) {
    delete h;
    return out;
  }
  out.amp = peak_val;
  out.mu_x = m.mu_x;
  out.mu_y = m.mu_y;
  out.sigma_x = m.sigma_x;
  out.sigma_y = m.sigma_y;
  out.rho = m.rho;
  out.ok = kTRUE;

  std::lock_guard<std::mutex> lock(g_plot_mutex);
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureAndDraw2DHistogram(h, c);
  TEllipse *e = new TEllipse(out.mu_x, out.mu_y, kGateNSigmaX * out.sigma_x,
                             kGateNSigmaY * out.sigma_y);
  e->SetFillStyle(0);
  e->SetLineColor(kRed + 1);
  e->SetLineWidth(2);
  e->Draw();
  PlottingUtils::SaveFigure(c, Form("beam_gate_s%d_s%d", sx, sy), subdir,
                            PlotSaveOptions::kLINEAR);
  delete c;

  delete h;
  return out;
}

void StripSumScatter::DrawTraceSet(const std::vector<TGraph *> &traces,
                                   Int_t color) {
  for (std::size_t i = 0; i < traces.size(); i++) {
    // Opaque lines, not a fractional alpha: the experimental overlay renders on
    // a live (non-batch) canvas whose backend ignores line alpha and draws it
    // opaque, while the sim overlay renders in batch mode where alpha IS
    // honored
    // -- so any alpha < 1 makes the two overlays look different. Full opacity
    // is the only setting that renders identically in both, keeping sim and exp
    // consistent.
    traces[i]->SetLineColor(color);
    traces[i]->SetLineWidth(1);
    traces[i]->Draw("L SAME");
  }
}

TGraph *StripSumScatter::TraceFromTotal(const Float_t *total) {
  Double_t td[18];
  for (Int_t s = 0; s < 18; s++)
    td[s] = Double_t(total[s]);
  return TraceCreator::BuildTraceFromTotals(td);
}

void StripSumScatter::DrawRegionTraces(const TString &save_name,
                                       const TString &subdir,
                                       const std::vector<TGraph *> &beam,
                                       const std::vector<TGraph *> &aa,
                                       const std::vector<TGraph *> &an,
                                       Double_t y_min, Double_t y_max,
                                       const char *y_title) {
  std::lock_guard<std::mutex> lock(g_plot_mutex);
  TH2F *frame = new TH2F("h_region_trace_frame", Form(";Strip;%s", y_title), 18,
                         -0.5, 17.5, 100, y_min, y_max);
  frame->SetStats(0);
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  frame->Draw();
  DrawTraceSet(beam, kGray + 2);
  DrawTraceSet(aa, kAzure + 2);
  DrawTraceSet(an, kRed + 1);

  TGraph *p_beam = new TGraph(1);
  TGraph *p_aa = new TGraph(1);
  TGraph *p_an = new TGraph(1);
  TGraph *proxies[3] = {p_beam, p_aa, p_an};
  Int_t pcol[3] = {kGray + 2, kAzure + 2, kRed + 1};
  for (Int_t i = 0; i < 3; i++) {
    proxies[i]->SetPoint(0, -1e9, -1e9);
    proxies[i]->SetLineColor(pcol[i]);
    proxies[i]->SetLineWidth(3);
  }
  TLegend *leg = PlottingUtils::AddLegend(0.725, 0.875, 0.70, 0.86);
  leg->AddEntry(p_beam, "Beam", "l");
  leg->AddEntry(p_aa, "(#alpha,#alpha')", "l");
  leg->AddEntry(p_an, "(#alpha,n)", "l");
  leg->Draw();

  PlottingUtils::SaveFigure(c, save_name, subdir, PlotSaveOptions::kLINEAR);
  delete leg;
  delete p_beam;
  delete p_aa;
  delete p_an;
  delete c;
  delete frame;
}

void StripSumScatter::TraceYRange(const std::vector<TGraph *> &beam,
                                  const std::vector<TGraph *> &aa,
                                  const std::vector<TGraph *> &an,
                                  Double_t &y_min, Double_t &y_max) {
  y_min = std::numeric_limits<Double_t>::max();
  y_max = -std::numeric_limits<Double_t>::max();
  const std::vector<TGraph *> *sets[3] = {&beam, &aa, &an};
  for (Int_t si = 0; si < 3; si++) {
    const std::vector<TGraph *> &v = *sets[si];
    for (std::size_t i = 0; i < v.size(); i++) {
      Double_t x = 0.0, y = 0.0;
      for (Int_t k = 0; k < v[i]->GetN(); k++) {
        v[i]->GetPoint(k, x, y);
        // Only the long anodes (1-16) frame the range; guard strips 0/17 are
        // single-ended and would stretch it.
        if (x < 0.5 || x > 16.5)
          continue;
        if (y < y_min)
          y_min = y;
        if (y > y_max)
          y_max = y;
      }
    }
  }
  if (y_min > y_max) { // no in-range points sampled
    y_min = 0.0;
    y_max = 1.0;
  }
  Double_t pad = 0.05 * (y_max - y_min);
  if (pad <= 0.0)
    pad = 1.0;
  y_min -= pad;
  y_max += pad;
}

TCutG *StripSumScatter::PromptCut(TCanvas *c, const char *name,
                                  const char *label) {
  std::cout << "  >>> draw the " << label
            << " region: left-click vertices, double-click to close"
            << std::endl;
  c->cd();
  TCutG *cut = dynamic_cast<TCutG *>(c->WaitPrimitive("CUTG", "CutG"));
  if (!cut) {
    std::cerr << "  no " << label << " cut drawn" << std::endl;
    return nullptr;
  }
  cut->SetName(name);
  cut->SetLineColor(kBlack);
  cut->SetLineWidth(2);
  return cut;
}

TString StripSumScatter::BuildFingerprint(const std::vector<Int_t> &run_order,
                                          std::map<Int_t, TChain *> &chains) {
  // v4: pileup rejection (>=3 long strips at >=2x beam); per-strip y-bounds
  // come from Constants::STRIP_SUM_Y_RANGE (tunable), so they're folded in and
  // the scatters no longer depend on the sim.
  // v5: reservoir gained a raw-ADC total_adc[18] branch; bump so v4 caches
  // (which lack it) rebuild instead of loading uninitialized raw traces.
  // v6: normalization hardcoded to 1 a.u. (NORM_MUSIC_MEV removed), so the
  // norm field is gone; bump so old caches built at other norms rebuild.
  TString s = Form(
      "v6 reac[%d,%d] flatTol=%.3f jump[%.3f,%.3f] smooth=%d,%d "
      "step=%.3f s17=%.3f gate[s%d,s%d,%.2f,%.2f,%d,%.3f,%.3f] x[%.3f,%.3f,%d] "
      "ybins=%d",
      kReacMin, kReacMax, kBeamFlatTol, kReacJumpMin, kReacJumpMax,
      Int_t(Constants::STRIP_SUM_APPLY_SMOOTHNESS), kSmoothHiStrip,
      kSmoothMaxStep, kStrip17Max, kGateStripX, kGateStripY, kGateNSigmaX,
      kGateNSigmaY, kGateBins, kGateMin, kGateMax, kXMin, kXMax, kXBins,
      kYBins);
  // Active beam gates (also keyed by cache filename, but folded in here too so
  // a mismatch never silently reuses a stale same-named cache).
  std::vector<GateSpec> gates = ActiveGates();
  for (std::size_t i = 0; i < gates.size(); i++)
    s += Form(" g[s%d,s%d]", gates[i].sx, gates[i].sy);
  Double_t y_lo[64], y_hi[64];
  YBounds(y_lo, y_hi);
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++)
    s += Form(" y%d[%.3f,%.3f]", reac, y_lo[ReacIndex(reac)],
              y_hi[ReacIndex(reac)]);
  for (std::size_t i = 0; i < run_order.size(); i++) {
    Int_t run = run_order[i];
    s += Form(" r%d:%lld", run, chains[run]->GetEntries());
  }
  return s;
}

// Per-reaction-strip y-axis bounds straight from Constants::STRIP_SUM_Y_RANGE
// (tune the map per dataset, per strip); strips absent from the map fall back
// to STRIP_SUM_YMIN/STRIP_SUM_YMAX. x stays fixed (strip-independent).
void StripSumScatter::YBounds(Double_t *y_lo, Double_t *y_hi) {
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    Int_t ri = ReacIndex(reac);
    std::map<Int_t, std::pair<Double_t, Double_t>>::const_iterator it =
        Constants::STRIP_SUM_Y_RANGE.find(reac);
    if (it != Constants::STRIP_SUM_Y_RANGE.end()) {
      y_lo[ri] = it->second.first;
      y_hi[ri] = it->second.second;
    } else {
      y_lo[ri] = Constants::STRIP_SUM_YMIN;
      y_hi[ri] = Constants::STRIP_SUM_YMAX;
    }
  }
}

TString StripSumScatter::PrettyLabel(const TString &tag) {
  TString base = RemixSim::TagWithoutStrip(tag);
  base.ReplaceAll("_eres", "");
  if (base == "aa")
    return "(#alpha,#alpha')";
  if (base == "an")
    return "(#alpha,n)";
  if (base == "beam")
    return "Beam";
  return base;
}

// Per-strip sim normalization gains: read the sim beam file and, for each
// strip, average the (unit-gain) per-strip beam total, then set gain[s] = 1 /
// mean[s] so every strip's sim beam lands on 1 a.u. This flattens the sim beam
// the SAME way the per-channel data normalization flattens the experimental
// beam
// -- a single global factor would not, since it preserves the sim's per-strip
// structure. Strips with no beam signal keep gain 0 (drop out like an
// uncalibrated channel).
Bool_t StripSumScatter::SimBeamGains(Double_t *gain) {
  for (Int_t s = 0; s < 18; s++)
    gain[s] = 0.0;
  // Reference the ERES beam file -- the same file type as the eres populations
  // plotted in SimOverlay/SimTraceOverlay -- so the normalized eres beam lands
  // exactly on 1. Falls back to the non-eres SIM_BEAM_FILE if no eres beam
  // control file is present.
  TString file;
  std::vector<RemixSim::SimFileSpec> specs = RemixSim::BuildFileSpecs();
  for (std::size_t i = 0; i < specs.size(); i++) {
    if (!RemixSim::IsEresTag(specs[i].tag))
      continue;
    TString base = RemixSim::TagWithoutStrip(specs[i].tag);
    base.ReplaceAll("_eres", "");
    if (base == "beam") {
      file = RemixSim::SimRootPath(specs[i]);
      break;
    }
  }
  if (file.Length() == 0)
    file = Paths::DatasetDir() + "/sim_root_files/" + Constants::SIM_BEAM_FILE;
  TFile *f = IO::OpenForReading(file);
  if (!f || f->IsZombie()) {
    std::cerr << "strip-sum-scatter: cannot open sim beam file " << file
              << "; sim overlay stays in raw sim units." << std::endl;
    if (f)
      delete f;
    return kFALSE;
  }
  TTree *t = static_cast<TTree *>(f->Get("events_MeV"));
  if (!t) {
    std::cerr << "strip-sum-scatter: no events_MeV tree in sim beam file "
              << file << "; sim overlay stays in raw sim units." << std::endl;
    f->Close();
    delete f;
    return kFALSE;
  }
  Float_t left[18] = {0}, right[18] = {0};
  t->SetBranchAddress("Left_0_17_dE", left);
  t->SetBranchAddress("RightdE", right);
  Long64_t n = t->GetEntries();
  Long64_t stride = FileSet::SampleStride(n, kSampleMaxPoints);
  Double_t sum[18] = {0};
  Long64_t cnt[18] = {0};
  // Unit gains so SimTotal yields the raw per-strip beam total (IGNORE_SHORT
  // aware), which is exactly the quantity these gains will later normalize.
  Double_t unit[18];
  for (Int_t s = 0; s < 18; s++)
    unit[s] = 1.0;
  for (Long64_t j = 0; j < n; j += stride) {
    t->GetEntry(j);
    Double_t total[18];
    SimTotal(left, right, unit, total);
    for (Int_t s = 0; s < 18; s++)
      if (total[s] > 0.0) {
        sum[s] += total[s];
        cnt[s]++;
      }
  }
  f->Close();
  delete f;
  Int_t n_set = 0;
  for (Int_t s = 0; s < 18; s++) {
    if (cnt[s] > 0 && sum[s] > 0.0) {
      gain[s] = 1.0 / (sum[s] / Double_t(cnt[s]));
      n_set++;
    }
  }
  if (n_set == 0)
    return kFALSE;
  std::cout << "strip-sum-scatter: sim per-strip beam normalization to 1 a.u. ("
            << n_set << " strips)." << std::endl;
  return kTRUE;
}

void StripSumScatter::SimTotal(const Float_t *left, const Float_t *right,
                               const Double_t *gain, Double_t *total) {
  for (Int_t s = 0; s < 18; s++)
    total[s] = gain[s] * (Double_t(left[s]) + Double_t(right[s]));
  if (Constants::IGNORE_SHORT_STRIPS)
    for (Int_t s = 1; s <= 16; s++)
      total[s] =
          gain[s] * ((s % 2) != 0 ? Double_t(left[s]) : Double_t(right[s]));
}

TGraph *StripSumScatter::SimPopScatter(const TString &file, Int_t reac,
                                       const Double_t *gain,
                                       Long64_t max_points) {
  TFile *f = IO::OpenForReading(file);
  if (!f || f->IsZombie()) {
    if (f)
      delete f;
    return nullptr;
  }
  TTree *t = static_cast<TTree *>(f->Get("events_MeV"));
  if (!t) {
    std::cerr << "  no events_MeV tree in " << file << std::endl;
    f->Close();
    delete f;
    return nullptr;
  }
  Float_t left[18] = {0}, right[18] = {0};
  t->SetBranchAddress("Left_0_17_dE", left);
  t->SetBranchAddress("RightdE", right);
  Long64_t n = t->GetEntries();
  Long64_t stride = FileSet::SampleStride(n, max_points);
  TGraph *g = new TGraph();
  Long64_t k = 0;
  for (Long64_t j = 0; j < n; j += stride) {
    t->GetEntry(j);
    Double_t total[18];
    SimTotal(left, right, gain, total);
    Double_t x = SumRange(total, kXLo, kXHi);
    Double_t y = SumRange(total, YLoOf(reac), YHiOf(reac));
    if (x > 0.0)
      g->SetPoint(k++, x, y);
  }
  g->Set(k);
  f->Close();
  delete f;
  return g;
}

// Up to max_traces per-strip dE-profile traces, stride-sampled across one sim
// file. Sim energies are arbitrary-unit floats; total[s] is the per-strip
// normalized gain[s]*(left[s] + right[s]) so the traces share the data's axis
// (and the sim beam is flat at NORM, like the data).
std::vector<TGraph *> StripSumScatter::SimPopTraces(const TString &file,
                                                    const Double_t *gain,
                                                    Long64_t max_traces) {
  std::vector<TGraph *> traces;
  TFile *f = IO::OpenForReading(file);
  if (!f || f->IsZombie()) {
    if (f)
      delete f;
    return traces;
  }
  TTree *t = static_cast<TTree *>(f->Get("events_MeV"));
  if (!t) {
    f->Close();
    delete f;
    return traces;
  }
  Float_t left[18] = {0}, right[18] = {0};
  t->SetBranchAddress("Left_0_17_dE", left);
  t->SetBranchAddress("RightdE", right);
  Long64_t n = t->GetEntries();
  Long64_t stride = FileSet::SampleStride(n, max_traces);
  for (Long64_t j = 0; j < n && Int_t(traces.size()) < max_traces;
       j += stride) {
    t->GetEntry(j);
    Double_t total[18];
    SimTotal(left, right, gain, total);
    traces.push_back(TraceCreator::BuildTraceFromTotals(total));
  }
  f->Close();
  delete f;
  return traces;
}

// Per reaction strip, overlay kTracesPerRegion sampled per-strip traces of each
// sim population in the experimental DrawRegionTraces style (beam grey, (a,a')
// azure, (a,n) red). The beam reference is the same for every strip. Sampled
// fresh each run (40 traces/file is trivial).
void StripSumScatter::SimTraceOverlay() {
  std::vector<RemixSim::SimFileSpec> specs = RemixSim::BuildFileSpecs();
  if (specs.empty())
    return;
  std::map<Int_t, TString> aa_file, an_file; // reaction strip -> sim file
  std::vector<TString> beam_files;
  for (std::size_t i = 0; i < specs.size(); i++) {
    if (!RemixSim::IsEresTag(specs[i].tag))
      continue;
    TString base = RemixSim::TagWithoutStrip(specs[i].tag);
    base.ReplaceAll("_eres", "");
    TString file = RemixSim::SimRootPath(specs[i]);
    Int_t strip = RemixSim::ReactionStripOf(specs[i].tag);
    if (base == "beam")
      beam_files.push_back(file);
    else if (strip >= kReacMin && strip <= kReacMax) {
      if (base == "aa")
        aa_file[strip] = file;
      else if (base == "an")
        an_file[strip] = file;
    }
  }

  Double_t gain[18];
  if (!SimBeamGains(gain))
    for (Int_t s = 0; s < 18; s++)
      gain[s] = 1.0;

  std::vector<TGraph *> beam_traces;
  for (std::size_t i = 0;
       i < beam_files.size() && Int_t(beam_traces.size()) < kTracesPerRegion;
       i++) {
    std::vector<TGraph *> t = SimPopTraces(
        beam_files[i], gain, kTracesPerRegion - Int_t(beam_traces.size()));
    for (std::size_t k = 0; k < t.size(); k++)
      beam_traces.push_back(t[k]);
  }

  for (Int_t r = kReacMin; r <= kReacMax; r++) {
    std::vector<TGraph *> aa_traces, an_traces;
    if (aa_file.find(r) != aa_file.end())
      aa_traces = SimPopTraces(aa_file[r], gain, kTracesPerRegion);
    if (an_file.find(r) != an_file.end())
      an_traces = SimPopTraces(an_file[r], gain, kTracesPerRegion);
    if (aa_traces.empty() && an_traces.empty())
      continue;
    DrawRegionTraces(Form("sim_region_traces_reac%d", r), "sim_scatter",
                     beam_traces, aa_traces, an_traces);
    for (std::size_t i = 0; i < aa_traces.size(); i++)
      delete aa_traces[i];
    for (std::size_t i = 0; i < an_traces.size(); i++)
      delete an_traces[i];
  }
  for (std::size_t i = 0; i < beam_traces.size(); i++)
    delete beam_traces[i];
}

// Fingerprint of the sim inputs + window geometry: each eres file's size+mtime
// (cheap, no open) plus the reaction-strip range and x window. Regenerating the
// sim (new mtimes) or changing the windows invalidates the cached overlay.
TString StripSumScatter::SimFingerprint(
    const std::vector<RemixSim::SimFileSpec> &specs) {
  // v2: sim is per-strip normalized via SimBeamGains(); the gain source (the
  // eres beam file) is stamped by the per-spec loop below (it covers every
  // eres file, including beam_eres).
  // v3: normalization hardcoded to 1 a.u. (NORM_MUSIC_MEV removed); bump so
  // overlays cached at other norms rebuild.
  TString s = Form("v3 reac[%d,%d] x[%d,%d]", kReacMin, kReacMax, kXLo, kXHi);
  for (std::size_t i = 0; i < specs.size(); i++) {
    if (!RemixSim::IsEresTag(specs[i].tag))
      continue;
    TString f = RemixSim::SimRootPath(specs[i]);
    Long_t id = 0, flags = 0, mtime = 0;
    Long64_t size = -1;
    if (gSystem->GetPathInfo(f, &id, &size, &flags, &mtime) != 0) {
      size = -1;
      mtime = 0;
    }
    s += Form(" %s:%lld:%ld", specs[i].tag.Data(), size, mtime);
  }
  return s;
}

// Reload cached sim scatter graphs (grouped by reaction strip; each graph's
// title holds its population label) if the fingerprint matches. Caller owns the
// returned graphs.
Bool_t StripSumScatter::LoadSimCache(
    const TString &fp, std::map<Int_t, std::vector<TGraph *>> &by_strip) {
  TString full = IO::GetRootFilesBaseDir() + TString("/") + kSimCacheName;
  if (gSystem->AccessPathName(full))
    return kFALSE;
  TFile *f = IO::OpenForReading(kSimCacheName);
  if (!f || f->IsZombie()) {
    if (f)
      delete f;
    return kFALSE;
  }
  TNamed *cfp = dynamic_cast<TNamed *>(f->Get("sim_fingerprint"));
  if (!cfp || fp != cfp->GetTitle()) {
    f->Close();
    delete f;
    return kFALSE;
  }
  TIter next(f->GetListOfKeys());
  TKey *key;
  while ((key = static_cast<TKey *>(next()))) {
    TString name = key->GetName();
    if (!name.BeginsWith("simg_r"))
      continue;
    TString rest = name(6, name.Length() - 6); // after "simg_r": <strip>_p<idx>
    Int_t us = rest.Index("_p");
    if (us < 0)
      continue;
    Int_t r = TString(rest(0, us)).Atoi();
    TGraph *g = dynamic_cast<TGraph *>(f->Get(name));
    if (!g)
      continue;
    by_strip[r].push_back(static_cast<TGraph *>(g->Clone()));
  }
  f->Close();
  delete f;
  return kTRUE;
}

void StripSumScatter::WriteSimCache(
    const TString &fp, const std::map<Int_t, std::vector<TGraph *>> &by_strip) {
  TFile *out = IO::OpenForWriting(kSimCacheName, "RECREATE");
  if (!out || out->IsZombie()) {
    if (out)
      delete out;
    return;
  }
  out->cd();
  TNamed cfp("sim_fingerprint", fp.Data());
  cfp.Write();
  std::map<Int_t, std::vector<TGraph *>>::const_iterator it;
  for (it = by_strip.begin(); it != by_strip.end(); ++it)
    for (std::size_t i = 0; i < it->second.size(); i++)
      it->second[i]->Write(Form("simg_r%d_p%d", it->first, Int_t(i)));
  out->Close();
  delete out;
}

// Sim-only comparison plots: one per reaction strip, each sim population a
// coloured+labelled point cloud on the same axes as that strip's data scatter,
// for side-by-side comparison with the data PID scatters. Beam (no reaction
// strip) overlays on every strip. The scatter graphs are fingerprint-cached
// (sim file sizes/mtimes + window geometry), so re-runs reload them instead of
// rescanning the sim files.
void StripSumScatter::SimOverlay(const std::map<Int_t, TH2F *> &scatter) {
  std::vector<RemixSim::SimFileSpec> specs = RemixSim::BuildFileSpecs();
  if (specs.empty()) {
    std::cerr
        << "strip-sum-scatter: no sim control files; skipping sim overlay."
        << std::endl;
    return;
  }
  TString fp = SimFingerprint(specs);

  std::map<Int_t, std::vector<TGraph *>>
      by_strip; // strip -> graphs (title=label)
  Bool_t loaded = LoadSimCache(fp, by_strip);

  if (!loaded) {
    std::map<Int_t, std::vector<SimPop>> reacted;
    std::vector<SimPop> refs;
    for (std::size_t i = 0; i < specs.size(); i++) {
      if (!RemixSim::IsEresTag(specs[i].tag))
        continue;
      SimPop p;
      p.file = RemixSim::SimRootPath(specs[i]);
      p.label = PrettyLabel(specs[i].tag);
      Int_t strip = RemixSim::ReactionStripOf(specs[i].tag);
      if (strip < 0)
        refs.push_back(p);
      else
        reacted[strip].push_back(p);
    }
    Double_t gain[18];
    if (!SimBeamGains(gain))
      for (Int_t s = 0; s < 18; s++)
        gain[s] = 1.0;
    const Long64_t kSimMaxPoints = 25000;
    for (Int_t r = kReacMin; r <= kReacMax; r++) {
      std::vector<SimPop> group = reacted[r];
      for (std::size_t i = 0; i < refs.size(); i++)
        group.push_back(refs[i]);
      for (std::size_t i = 0; i < group.size(); i++) {
        TGraph *g = SimPopScatter(group[i].file, r, gain, kSimMaxPoints);
        if (!g || g->GetN() == 0) {
          if (g)
            delete g;
          continue;
        }
        g->SetTitle(group[i].label);
        by_strip[r].push_back(g);
      }
    }
    Int_t n_graphs = 0;
    std::map<Int_t, std::vector<TGraph *>>::const_iterator cit;
    for (cit = by_strip.begin(); cit != by_strip.end(); ++cit)
      n_graphs += Int_t(cit->second.size());
    if (n_graphs == 0) {
      std::cerr << "strip-sum-scatter: no sim data found (regenerate "
                   "sim_root_files); skipping sim overlay."
                << std::endl;
      return;
    }
    WriteSimCache(fp, by_strip);
    std::cout << "strip-sum-scatter: built + cached sim overlay (" << n_graphs
              << " population graphs)." << std::endl;
  } else {
    std::cout
        << "strip-sum-scatter: loaded cached sim overlay (fingerprint match)."
        << std::endl;
  }

  std::map<Int_t, std::vector<TGraph *>>::iterator it;
  for (it = by_strip.begin(); it != by_strip.end(); ++it) {
    Int_t r = it->first;
    std::map<Int_t, TH2F *>::const_iterator sit = scatter.find(r);
    if (sit == scatter.end())
      continue;
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    TH2F *ref = sit->second;
    TH2F *frame =
        new TH2F(Form("sim_frame_r%d", r), "", 10, ref->GetXaxis()->GetXmin(),
                 ref->GetXaxis()->GetXmax(), 10, ref->GetYaxis()->GetXmin(),
                 ref->GetYaxis()->GetXmax());
    frame->SetStats(0);
    frame->GetXaxis()->SetTitle(ref->GetXaxis()->GetTitle());
    frame->GetYaxis()->SetTitle(ref->GetYaxis()->GetTitle());
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    c->SetLeftMargin(0.18);
    frame->Draw();
    // Match the experimental region-traces legend placement (top-right).
    TLegend *leg = PlottingUtils::AddLegend(0.725, 0.875, 0.70, 0.86);
    for (std::size_t i = 0; i < it->second.size(); i++) {
      TGraph *g = it->second[i];
      // Match the experimental region-trace colours (DrawRegionTraces): beam
      // grey, (a,a') azure, (a,n) red -- keyed off the population label.
      TString lab = g->GetTitle();
      Int_t color = kBlack;
      if (lab == "Beam")
        color = kGray + 2;
      else if (lab == "(#alpha,#alpha')")
        color = kAzure + 2;
      else if (lab == "(#alpha,n)")
        color = kRed + 1;
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.3);
      g->SetMarkerColorAlpha(color, 0.35);
      g->SetLineColor(color);
      g->Draw("P SAME");
      leg->AddEntry(g, g->GetTitle(), "p");
    }
    leg->Draw();
    PlottingUtils::SaveFigure(c,
                              Form("sim_normsumE_reac%d_s%d_%d_vs_s%d_%d", r,
                                   YLoOf(r), YHiOf(r), kXLo, kXHi),
                              "sim_scatter", PlotSaveOptions::kLINEAR);
    delete leg;
    delete c;
    delete frame;
  }

  std::map<Int_t, std::vector<TGraph *>>::iterator dit;
  for (dit = by_strip.begin(); dit != by_strip.end(); ++dit)
    for (std::size_t i = 0; i < dit->second.size(); i++)
      delete dit->second[i];
}

void StripSumScatter::Run() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                Paths::ResultsDir() + "/plots",
                                Paths::ResultsDir() + "/root_files");
  gROOT->SetBatch(kTRUE);

  std::vector<Int_t> run_order;
  std::map<Int_t, TChain *> chain_by_run = FileSet::GroupEventsByRun(run_order);
  if (run_order.empty()) {
    std::cerr << "strip-sum-scatter: no runs found" << std::endl;
    return;
  }

  TString fingerprint = BuildFingerprint(run_order, chain_by_run);
  // Per-gate-config cache file, so toggling a gate switches files instead of
  // clobbering another config's cached results.
  TString cache_name = CacheName();

  std::map<Int_t, TH2F *> scatter; // keyed by reaction strip
  std::vector<TraceEvt> reservoir;

  Bool_t loaded = kFALSE;
  TString cache_full = IO::GetRootFilesBaseDir() + TString("/") + cache_name;
  if (!gSystem->AccessPathName(cache_full)) {
    TFile *cf = IO::OpenForReading(cache_name);
    if (cf && !cf->IsZombie()) {
      TNamed *fp = dynamic_cast<TNamed *>(cf->Get("fingerprint"));
      if (fp && fingerprint == fp->GetTitle()) {
        Bool_t ok = kTRUE;
        for (Int_t reac = kReacMin; reac <= kReacMax && ok; reac++) {
          TH2F *h = dynamic_cast<TH2F *>(cf->Get(Form("scatter_r%d", reac)));
          if (!h) {
            ok = kFALSE;
            break;
          }
          TH2F *hc = static_cast<TH2F *>(h->Clone());
          hc->SetDirectory(nullptr);
          scatter[reac] = hc;
        }
        TTree *tt = dynamic_cast<TTree *>(cf->Get("traces"));
        if (ok && tt) {
          TraceEvt e;
          tt->SetBranchAddress("total", e.total);
          tt->SetBranchAddress("total_adc", e.total_adc);
          tt->SetBranchAddress("reac_mask", &e.reac_mask);
          tt->SetBranchAddress("beam_flat", &e.beam_flat);
          Long64_t nt = tt->GetEntries();
          reservoir.reserve(nt);
          for (Long64_t j = 0; j < nt; j++) {
            tt->GetEntry(j);
            reservoir.push_back(e);
          }
          loaded = ok;
        }
      }
      cf->Close();
      delete cf;
    }
    if (loaded)
      std::cout << "strip-sum-scatter: loaded cached scatters + "
                << reservoir.size() << " reservoir events (fingerprint match)."
                << std::endl;
    else
      std::cout << "strip-sum-scatter: cache present but stale; rebuilding."
                << std::endl;
  }

  if (!loaded) {
    std::vector<GateSpec> active_gates = ActiveGates();
    // One ellipse fit per active gate per run; a run is usable only if every
    // active gate fit succeeds.
    std::map<Int_t, std::vector<BeamFit2D>> gates;
    std::map<Int_t, Bool_t> gates_ok;
    for (std::size_t i = 0; i < run_order.size(); i++) {
      Int_t run = run_order[i];
      gates_ok[run] = kFALSE;
      if (!chain_by_run[run] || chain_by_run[run]->GetEntries() == 0)
        continue;
      std::vector<BeamFit2D> run_gates;
      std::vector<GateSpec> prior_specs;
      Bool_t all_ok = kTRUE;
      for (std::size_t gi = 0; gi < active_gates.size(); gi++) {
        // Fit each gate in series on data already passing the prior gates.
        BeamFit2D g = FindBeamGate(chain_by_run[run], active_gates[gi].sx,
                                   active_gates[gi].sy, prior_specs, run_gates,
                                   Form("run%d", run),
                                   Form("strip_sum_scatter/run%d", run));
        if (g.ok)
          std::cout << "  run " << run << " beam gate s" << active_gates[gi].sx
                    << "/s" << active_gates[gi].sy << ": mu=(" << g.mu_x << ","
                    << g.mu_y << ")" << std::endl;
        else {
          std::cerr << "  run " << run << " beam gate s" << active_gates[gi].sx
                    << "/s" << active_gates[gi].sy << " failed; skipping run"
                    << std::endl;
          all_ok = kFALSE;
        }
        run_gates.push_back(g);
        prior_specs.push_back(active_gates[gi]);
      }
      gates[run] = run_gates;
      gates_ok[run] = all_ok;
    }

    Double_t y_lo[64], y_hi[64];
    YBounds(y_lo, y_hi);

    for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
      Int_t ri = ReacIndex(reac);
      TH2F *h = new TH2F(
          Form("scatter_r%d", reac),
          Form(";norm. #DeltaE strips %d#rightarrow%d [a.u.];norm. #DeltaE "
               "strips %d#rightarrow%d [a.u.]",
               kXLo, kXHi, YLoOf(reac), YHiOf(reac)),
          kXBins, kXMin, kXMax, kYBins, y_lo[ri], y_hi[ri]);
      h->SetDirectory(nullptr);
      h->SetStats(0);
      scatter[reac] = h;
    }

    Long64_t total_gated = 0, total_seen = 0;
    Int_t n_beam_kept = 0;
    for (std::size_t i = 0; i < run_order.size(); i++) {
      Int_t run = run_order[i];
      TChain *chain = chain_by_run[run];
      if (!chain || !gates_ok[run])
        continue;
      const std::vector<BeamFit2D> &run_gates = gates[run];
      EnergyView ev;
      ev.Attach(chain);
      EnableEventBranches(chain);
      Long64_t n = chain->GetEntries();
      std::cout << "Run " << run << ": filling " << kNReac
                << " reaction-strip scatters over " << n << " events..."
                << std::endl;
      for (Long64_t j = 0; j < n; j++) {
        chain->GetEntry(j);
        ev.Decode();
        total_seen++;
        Bool_t passes_all = kTRUE;
        for (std::size_t gi = 0; gi < active_gates.size(); gi++)
          if (!PassesGate(run_gates[gi], ev, active_gates[gi].sx,
                          active_gates[gi].sy)) {
            passes_all = kFALSE;
            break;
          }
        if (!passes_all)
          continue;
        if (IsPileup(ev)) // reject overlapping-beam pileup
          continue;
        if (IsNoise(ev))
          continue;
        Double_t x = SumRange(ev.total, kXLo, kXHi);
        UInt_t mask = 0;
        for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
          if (!PassesReaction(ev, reac))
            continue;
          mask |= (1u << ReacIndex(reac));
          scatter[reac]->Fill(x, SumRange(ev.total, YLoOf(reac), YHiOf(reac)));
        }
        // Keep every reaction-passing event for traces; cap clean flat-beam
        // events (only ~kTracesPerRegion are ever drawn). The two are mutually
        // exclusive -- a flat-beam event has no reaction jump.
        Bool_t beam = (mask == 0) && IsBeamFlat(ev);
        if (mask == 0 && !(beam && n_beam_kept < kBeamReservoirCap))
          continue;
        if (beam)
          n_beam_kept++;
        TraceEvt e;
        for (Int_t s = 0; s < 18; s++) {
          e.total[s] = Float_t(ev.total[s]);
          e.total_adc[s] =
              Float_t(ev.left_0_17_adc[s]) + Float_t(ev.rightdE_adc[s]);
        }
        // Mirror IGNORE_SHORT_STRIPS: the normed total keeps only the long side
        // of a split strip, so the raw trace must drop the same side to stay
        // comparable.
        if (Constants::IGNORE_SHORT_STRIPS)
          for (Int_t s = 1; s <= 16; s++)
            e.total_adc[s] = ((s % 2) != 0) ? Float_t(ev.left_0_17_adc[s])
                                            : Float_t(ev.rightdE_adc[s]);
        e.reac_mask = mask;
        e.beam_flat = beam;
        reservoir.push_back(e);
        if (mask != 0)
          total_gated++;
      }
    }
    std::cout << "Built scatters: " << total_gated
              << " reaction events across strips " << kReacMin << "-"
              << kReacMax << " (" << total_seen << " seen), reservoir "
              << reservoir.size() << " events." << std::endl;

    TFile *out = IO::OpenForWriting(cache_name, "RECREATE");
    if (out && !out->IsZombie()) {
      out->cd();
      TNamed fp("fingerprint", fingerprint.Data());
      fp.Write();
      for (Int_t reac = kReacMin; reac <= kReacMax; reac++)
        scatter[reac]->Write(Form("scatter_r%d", reac));
      TTree *tt = new TTree("traces", "strip-sum trace reservoir");
      TraceEvt e;
      tt->Branch("total", e.total, "total[18]/F");
      tt->Branch("total_adc", e.total_adc, "total_adc[18]/F");
      tt->Branch("reac_mask", &e.reac_mask, "reac_mask/i");
      tt->Branch("beam_flat", &e.beam_flat, "beam_flat/O");
      for (std::size_t k = 0; k < reservoir.size(); k++) {
        e = reservoir[k];
        tt->Fill();
      }
      tt->Write();
      out->Close();
      delete out;
      std::cout << "strip-sum-scatter: wrote cache " << cache_name << std::endl;
    }
  }

  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
      TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDraw2DHistogram(scatter[reac], c);
      scatter[reac]->GetYaxis()->SetTitleOffset(1.3);
      c->SetLeftMargin(0.18);
      PlottingUtils::SaveFigure(c,
                                Form("normsumE_reac%d_s%d_%d_vs_s%d_%d", reac,
                                     YLoOf(reac), YHiOf(reac), kXLo, kXHi),
                                "strip_sum_scatter", PlotSaveOptions::kLINEAR);
      delete c;
    }
  }

  if (Constants::STRIP_SUM_RERUN_SIM) {
    SimOverlay(scatter);
    SimTraceOverlay();
  }

  Int_t reac = Constants::STRIP_SUM_CANDIDATE_REACTION_STRIP;
  if (reac < kReacMin || reac > kReacMax) {
    std::cerr << "strip-sum-scatter: candidate reaction strip " << reac
              << " outside [" << kReacMin << "," << kReacMax
              << "]; skipping interactive overlay." << std::endl;
  } else if (!gSystem->Getenv("DISPLAY")) {
    std::cerr << "strip-sum-scatter: no DISPLAY; skipping interactive "
                 "region-trace overlay (scatters already saved)."
              << std::endl;
  } else {
    Int_t app_argc = 1;
    char app_arg0[] = "strip-sum-scatter";
    char *app_argv[] = {app_arg0};
    TApplication app("strip-sum-scatter", &app_argc, app_argv);
    gROOT->SetBatch(kFALSE);

    TCanvas *cut_canvas = new TCanvas(
        "c_strip_sum_regions", "Draw (a,n) then (a,a') regions", 900, 700);
    cut_canvas->SetLogz(kTRUE); // match the saved scatter's z-scale
    scatter[reac]->Draw("COLZ");
    cut_canvas->Update();
    TCutG *cut_an = PromptCut(cut_canvas, "region_an", "(a,n)");
    TCutG *cut_aa = PromptCut(cut_canvas, "region_aa", "(a,a')");

    std::vector<TGraph *> tr_an, tr_aa, tr_beam;
    // Same selected events, raw (un-normalized) ADC -- one entry per normed
    // trace, kept in lock-step so the two overlays show the identical events.
    std::vector<TGraph *> tr_an_adc, tr_aa_adc, tr_beam_adc;
    UInt_t bit = (1u << ReacIndex(reac));
    for (std::size_t k = 0; k < reservoir.size(); k++) {
      if (tr_an.size() >= std::size_t(kTracesPerRegion) &&
          tr_aa.size() >= std::size_t(kTracesPerRegion) &&
          tr_beam.size() >= std::size_t(kTracesPerRegion))
        break;
      const TraceEvt &e = reservoir[k];
      if (e.beam_flat && tr_beam.size() < std::size_t(kTracesPerRegion)) {
        tr_beam.push_back(TraceFromTotal(e.total));
        tr_beam_adc.push_back(TraceFromTotal(e.total_adc));
        continue;
      }
      if (!(e.reac_mask & bit))
        continue;
      Double_t td[18];
      for (Int_t s = 0; s < 18; s++)
        td[s] = Double_t(e.total[s]);
      Double_t x = SumRange(td, kXLo, kXHi);
      Double_t y = SumRange(td, YLoOf(reac), YHiOf(reac));
      if (cut_an && tr_an.size() < std::size_t(kTracesPerRegion) &&
          cut_an->IsInside(x, y)) {
        tr_an.push_back(TraceFromTotal(e.total));
        tr_an_adc.push_back(TraceFromTotal(e.total_adc));
      } else if (cut_aa && tr_aa.size() < std::size_t(kTracesPerRegion) &&
                 cut_aa->IsInside(x, y)) {
        tr_aa.push_back(TraceFromTotal(e.total));
        tr_aa_adc.push_back(TraceFromTotal(e.total_adc));
      }
    }

    std::cout << "Sampled traces: beam=" << tr_beam.size()
              << " (a,a')=" << tr_aa.size() << " (a,n)=" << tr_an.size()
              << std::endl;
    DrawRegionTraces(Form("region_traces_reac%d", reac), "strip_sum_scatter",
                     tr_beam, tr_aa, tr_an);
    // Second overlay: identical events, raw ADC. y auto-ranges over the long
    // anodes (1-16) since the un-gained per-strip levels span the gain spread
    // the normalization removes -- a fixed normed band would hide it.
    Double_t adc_y_lo = 0.0, adc_y_hi = 0.0;
    TraceYRange(tr_beam_adc, tr_aa_adc, tr_an_adc, adc_y_lo, adc_y_hi);
    DrawRegionTraces(Form("region_traces_reac%d_adc", reac),
                     "strip_sum_scatter", tr_beam_adc, tr_aa_adc, tr_an_adc,
                     adc_y_lo, adc_y_hi, "#DeltaE [ADC]");

    for (std::size_t i = 0; i < tr_an.size(); i++)
      delete tr_an[i];
    for (std::size_t i = 0; i < tr_aa.size(); i++)
      delete tr_aa[i];
    for (std::size_t i = 0; i < tr_beam.size(); i++)
      delete tr_beam[i];
    for (std::size_t i = 0; i < tr_an_adc.size(); i++)
      delete tr_an_adc[i];
    for (std::size_t i = 0; i < tr_aa_adc.size(); i++)
      delete tr_aa_adc[i];
    for (std::size_t i = 0; i < tr_beam_adc.size(); i++)
      delete tr_beam_adc[i];
  }

  for (Int_t reac2 = kReacMin; reac2 <= kReacMax; reac2++)
    delete scatter[reac2];
  for (std::size_t i = 0; i < run_order.size(); i++)
    delete chain_by_run[run_order[i]];
}
