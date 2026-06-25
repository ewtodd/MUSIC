#include "StripSumScatter.hpp"

StripSumScatter::StripSumScatter() {
  for (Int_t i = 0; i < 64; i++) {
    m_yLo[i] = 0.0;
    m_yHi[i] = 0.0;
  }
}

StripSumScatter::~StripSumScatter() {
  std::map<Int_t, TH2F *>::iterator it;
  for (it = m_scatter.begin(); it != m_scatter.end(); ++it)
    delete it->second;
  m_scatter.clear();
}

Int_t StripSumScatter::ReacIndex(Int_t reac) {
  return reac - Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
}

Int_t StripSumScatter::YLoOf(Int_t reac) { return reac + 1; }

Int_t StripSumScatter::YHiOf(Int_t reac) { return TMath::Min(reac + 6, 17); }

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
  const Double_t kBeamFlatTol =
      Constants::STRIP_SUM_SCATTER_CONFIG.BEAM_FLAT_TOL;
  const Double_t kReacJumpMin =
      Constants::STRIP_SUM_SCATTER_CONFIG.REAC_JUMP_MIN;
  const Double_t kReacJumpMax =
      Constants::STRIP_SUM_SCATTER_CONFIG.REAC_JUMP_MAX;
  const Double_t kSmoothMaxStep =
      Constants::STRIP_SUM_SCATTER_CONFIG.REQUIRE_SMOOTHNESS_MAX_STEP;
  const Int_t kSmoothHiStrip =
      Constants::STRIP_SUM_SCATTER_CONFIG.REQUIRE_SMOOTHNESS_END_STRIP;
  const Double_t kStrip17Max = Constants::STRIP_SUM_SCATTER_CONFIG.STRIP_17_MAX;

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
  if (Constants::STRIP_SUM_SCATTER_CONFIG.REQUIRE_SMOOTHNESS)
    for (Int_t s = reac + 1; s <= kSmoothHiStrip; s++)
      if (TMath::Abs(ev.total[s] - ev.total[s - 1]) > kSmoothMaxStep)
        return kFALSE;
  return ev.total[17] < kStrip17Max;
}

Bool_t StripSumScatter::IsBeamFlat(const EnergyView &ev) {
  const Double_t kBeamFlatTol =
      Constants::STRIP_SUM_SCATTER_CONFIG.BEAM_FLAT_TOL;
  if (!AllStripsFired(ev))
    return kFALSE;
  for (Int_t s = 1; s <= 16; s++)
    if (TMath::Abs(ev.total[s] - 1.0) > kBeamFlatTol)
      return kFALSE;
  return kTRUE;
}

Bool_t StripSumScatter::IsPileup(const EnergyView &ev) {
  const Double_t kThresh = Constants::STRIP_SUM_SCATTER_CONFIG.PILEUP_THRESHOLD;
  Int_t n = 0;
  for (Int_t s = 1; s <= 16; s++)
    if (ev.total[s] >= kThresh && ++n >= 2)
      return kTRUE;
  return kFALSE;
}

Bool_t StripSumScatter::IsNoise(const EnergyView &ev) {
  const Double_t kThresh = Constants::STRIP_SUM_SCATTER_CONFIG.NOISE_THRESHOLD;
  Int_t n = 0;
  for (Int_t s = 1; s <= 16; s++)
    if (ev.total[s] <= kThresh && ++n >= 3)
      return kTRUE;
  return kFALSE;
}

Double_t StripSumScatter::SumRange(const Double_t *total, Int_t lo, Int_t hi) {
  Double_t sum = 0.0;
  for (Int_t s = lo; s <= hi; s++)
    sum += total[s];
  return sum;
}

std::vector<GateSpec> StripSumScatter::ActiveGates() {
  std::vector<GateSpec> gates;
  GateSpec g;
  g.sx = Constants::STRIP_SUM_SCATTER_CONFIG.GATE_STRIP_X;
  g.sy = Constants::STRIP_SUM_SCATTER_CONFIG.GATE_STRIP_Y;
  gates.push_back(g);
  if (Constants::STRIP_SUM_SCATTER_CONFIG.REQUIRE_GATE_S3_S4) {
    g.sx = 3;
    g.sy = 4;
    gates.push_back(g);
  }
  if (Constants::STRIP_SUM_SCATTER_CONFIG.REQUIRE_GATE_S5_S6) {
    g.sx = 5;
    g.sy = 6;
    gates.push_back(g);
  }
  return gates;
}

TString StripSumScatter::CacheName() {
  TString name = "StripSumScatter_cache";
  if (Constants::STRIP_SUM_SCATTER_CONFIG.REQUIRE_GATE_S3_S4)
    name += "_g34";
  if (Constants::STRIP_SUM_SCATTER_CONFIG.REQUIRE_GATE_S5_S6)
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
  const Double_t kGateNSigmaX =
      Constants::STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_X;
  const Double_t kGateNSigmaY =
      Constants::STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_Y;
  return BeamFitUtils::InEllipseXY(gate, g0, g1, kGateNSigmaX, kGateNSigmaY);
}

BeamFit2D
StripSumScatter::FindBeamGate(TChain *chain, Int_t sx, Int_t sy,
                              const std::vector<GateSpec> &prior_specs,
                              const std::vector<BeamFit2D> &prior_gates,
                              const TString &tag, const TString &subdir) {
  const Int_t kGateBins = Constants::STRIP_SUM_SCATTER_CONFIG.GATE_BINS;
  const Double_t kGateMin = Constants::STRIP_SUM_SCATTER_CONFIG.GATE_MIN;
  const Double_t kGateMax = Constants::STRIP_SUM_SCATTER_CONFIG.GATE_MAX;
  const Int_t kSeedHalfBins =
      Constants::STRIP_SUM_SCATTER_CONFIG.SEED_HALF_BINS;
  const Double_t kSeedFrac = 0.3;
  const Long64_t kSampleMaxPoints =
      Constants::STRIP_SUM_SCATTER_CONFIG.SAMPLE_MAX_POINTS;

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

  const Double_t kGateNSigmaX =
      Constants::STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_X;
  const Double_t kGateNSigmaY =
      Constants::STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_Y;

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

void StripSumScatter::DrawRegionMeanTraces(const TString &save_name,
                                           const TString &subdir,
                                           const std::vector<TGraph *> &beam,
                                           const std::vector<TGraph *> &aa,
                                           const std::vector<TGraph *> &an,
                                           Double_t y_min, Double_t y_max,
                                           const char *y_title) {
  std::lock_guard<std::mutex> lock(g_plot_mutex);
  TH2F *frame = new TH2F("h_region_mean_frame", Form(";Strip;%s", y_title), 18,
                         -0.5, 17.5, 100, y_min, y_max);
  frame->SetStats(0);
  frame->SetDirectory(nullptr);
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  frame->Draw();

  const std::vector<TGraph *> *regions[3] = {&beam, &aa, &an};
  Int_t colors[3] = {kGray + 2, kAzure + 2, kRed + 1};
  const char *labels[3] = {"Beam", "(#alpha,#alpha')", "(#alpha,n)"};
  std::vector<TGraphErrors *> means;
  TLegend *leg = PlottingUtils::AddLegend(0.725, 0.875, 0.70, 0.86);

  for (Int_t r = 0; r < 3; r++) {
    const std::vector<TGraph *> &tr = *regions[r];
    if (tr.empty())
      continue;
    Int_t npts = tr[0]->GetN();
    std::vector<Double_t> mean(npts, 0.0), m2(npts, 0.0);
    for (std::size_t t = 0; t < tr.size(); t++) {
      Double_t *yv = tr[t]->GetY();
      for (Int_t p = 0; p < npts; p++) {
        mean[p] += yv[p];
        m2[p] += yv[p] * yv[p];
      }
    }
    Double_t nt = Double_t(tr.size());
    TGraphErrors *ge = new TGraphErrors(npts);
    Double_t *xv = tr[0]->GetX();
    for (Int_t p = 0; p < npts; p++) {
      mean[p] /= nt;
      Double_t var = m2[p] / nt - mean[p] * mean[p];
      ge->SetPoint(p, xv[p], mean[p]);
      ge->SetPointError(p, 0.0, var > 0.0 ? TMath::Sqrt(var) : 0.0);
    }
    ge->SetLineColor(colors[r]);
    ge->SetLineWidth(3);
    ge->SetFillColorAlpha(colors[r], 0.15);
    ge->Draw("3 SAME"); // +-1 RMS band (all bands first, behind the lines)
    means.push_back(ge);
    leg->AddEntry(ge, labels[r], "l");
  }
  // Mean lines on top of every band.
  for (std::size_t i = 0; i < means.size(); i++)
    means[i]->Draw("LX SAME"); // mean line, no end caps
  leg->Draw();

  PlottingUtils::SaveFigure(c, save_name, subdir, PlotSaveOptions::kLINEAR);
  for (std::size_t i = 0; i < means.size(); i++)
    delete means[i];
  delete leg;
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

void StripSumScatter::SmoothTrace(const Double_t *in, Double_t *out,
                                  Int_t width) {
  Int_t half = width / 2;
  for (Int_t s = 0; s < 18; s++) {
    Int_t lo = TMath::Max(0, s - half);
    Int_t hi = TMath::Min(17, s + half);
    Double_t sum = 0.0;
    for (Int_t t = lo; t <= hi; t++)
      sum += in[t];
    out[s] = sum / Double_t(hi - lo + 1);
  }
}

// Savitzky-Golay smoothing: 3rd-degree polynomial, half-window of 2
// (5-point convolution). Uses standard SG coefficients that sum to 1.
// For a 5-point window with 3rd-degree polynomial, the smoothed value at
// the center uses coefficients: [-3, 12, 17, 12, -3] / 35.
// At edges, the window shrinks and coefficients are renormalised.
void StripSumScatter::SavitzkyGolay(const Double_t *in, Double_t *out) {
  static const Int_t K = 2; // half-width (5-point window)

  // Standard SG coefficients for 5-point, 3rd-degree polynomial (smoothed
  // value): These are translation-invariant - same for all center positions.
  static const Double_t sg_coeff[2 * K + 1] = {
      -3.0 / 35.0, // coefficient for t = s - 2
      12.0 / 35.0, // coefficient for t = s - 1
      17.0 / 35.0, // coefficient for t = s (center)
      12.0 / 35.0, // coefficient for t = s + 1
      -3.0 / 35.0  // coefficient for t = s + 2
  };

  for (Int_t s = 0; s < 18; s++) {
    Int_t lo = TMath::Max(0, s - K);
    Int_t hi = TMath::Min(17, s + K);
    Double_t val = 0.0;
    Double_t wsum = 0.0;

    // Apply SG coefficients for the positions within the clipped window
    for (Int_t t = lo; t <= hi; t++) {
      Int_t offset = t - s + K; // 0..4, position within 5-point window
      val += sg_coeff[offset] * in[t];
      wsum += sg_coeff[offset];
    }

    // Renormalise at edges where window shrinks
    out[s] = (wsum != 0.0) ? val / wsum : in[s];
  }
}

// CFD-style trigger finder: scan left-to-right for the first strip whose
// beam-subtracted signal exceeds both a fraction of the trace peak and a
// multiple of the beam sigma. Returns the strip index, or -1 if none fires.

Int_t StripSumScatter::FindTrigger(const Double_t *td, const Double_t *base,
                                   Double_t beam_sigma) {
  const Int_t s_lo = 2;
  const Int_t s_hi = 16;

  const Double_t frac = Constants::STRIP_SUM_SCATTER_CONFIG.TRIGGER_CFD_FRAC;

  const Double_t nsigma =
      Constants::STRIP_SUM_SCATTER_CONFIG.TRIGGER_NSIGMA * beam_sigma;

  Double_t peak_signal = -1.0e30;

  for (Int_t s = s_lo; s <= s_hi; s++) {
    Double_t signal = td[s] - base[s];
    if (signal > peak_signal)
      peak_signal = signal;
  }

  if (peak_signal < nsigma)
    return -1;

  Double_t thresh = frac * peak_signal;

  for (Int_t s = s_lo; s <= s_hi; s++) {
    Double_t signal = td[s] - base[s];

    if (signal >= thresh && signal >= nsigma)
      return s;
  }

  return -1;
}

// Build a TGraph trace from a Savitzky-Golay-smoothed set of per-strip totals.
TGraph *StripSumScatter::SmoothedTraceFromTotal(const Float_t *total) {
  Double_t td[18], sgd[18];
  for (Int_t s = 0; s < 18; s++)
    td[s] = Double_t(total[s]);
  SavitzkyGolay(td, sgd);
  return TraceCreator::BuildTraceFromTotals(sgd);
}

void StripSumScatter::ClusterVarHists(Int_t reac, TCutG *cut_aa, TCutG *cut_an,
                                      const TString &subdir) {
  const Int_t NV = 9;
  const Int_t NC = 3;
  const char *vkey[NV] = {"energy",      "peak3",      "plateau",
                          "tail",        "reacstrip",  "mult",
                          "trigtaildev", "reacslope3", "beamdev"};
  const char *vtitle[NV] = {
      "#Sigma_{all strips}(#DeltaE#minus1) [a.u.]",
      "#Sigma_{trig#pm1}#DeltaE (0 if no trigger) [a.u.]",
      "Plateau Excess  #Sigma_{trig+1..trig+POST}(#DeltaE#minus1) [a.u.]",
      "#DeltaE(s17) [a.u.]",
      "Trigger Strip",
      "Both-side Multiplicity (strips 1-16)",
      "|#DeltaE#minusbeam| at trigger + at s17 [a.u.]",
      "#DeltaE(reac+3) #minus #DeltaE(reac#minus3) [a.u.]",
      "RMS_{8-17}(#DeltaE#minusbeam) [a.u.]"};
  const char *clabel[NC] = {"beam", "(a,a')", "(a,n)"};

  const Int_t kXLo = Constants::STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::STRIP_SUM_SCATTER_CONFIG.X_HI;
  const Int_t kClusterSmoothWindow =
      Constants::STRIP_SUM_SCATTER_CONFIG.CLUSTER_SMOOTH_WINDOW;

  // Per (class, variable) value lists for raw traces.
  std::vector<Double_t> vals_raw[NC][NV];
  // Same structure, but cluster variables computed on Savitzky-Golay-smoothed
  // traces (SG kernel removes the L_odd/R_even sawtooth so peak/onset features
  // land on the real physical peak rather than an odd strip).
  std::vector<Double_t> vals_sg[NC][NV];
  UInt_t bit = (1u << ReacIndex(reac));

  // Per-strip beam baseline = mean over the beam-flat reservoir events. It
  // carries the L_odd/R_even sawtooth, so subtracting it removes that
  // systematic exactly before the reaction-onset jump search (better than
  // blurring it with smoothing).
  Double_t base[18];
  for (Int_t s = 0; s < 18; s++)
    base[s] = 0.0;
  Long64_t nbeam = 0;
  for (std::size_t k = 0; k < m_reservoir.size(); k++)
    if (m_reservoir[k].beam_flat) {
      for (Int_t s = 0; s < 18; s++)
        base[s] += Double_t(m_reservoir[k].total[s]);
      nbeam++;
    }
  for (Int_t s = 0; s < 18; s++)
    base[s] = (nbeam > 0) ? base[s] / Double_t(nbeam) : 1.0;

  // Pooled beam-noise RMS = sqrt(mean over beam-flat events and all strips of
  // (total - base)^2). The reaction-onset threshold is this many sigma (an
  // N-sigma discriminator), matching the Python pipeline -- so on flat beam
  // the excess does NOT cross on average and the event gets NO trigger.
  Double_t beam_sumsq = 0.0;
  Long64_t beam_npt = 0;
  for (std::size_t k = 0; k < m_reservoir.size(); k++)
    if (m_reservoir[k].beam_flat)
      for (Int_t s = 0; s < 18; s++) {
        Double_t d = Double_t(m_reservoir[k].total[s]) - base[s];
        beam_sumsq += d * d;
        beam_npt++;
      }
  Double_t beam_sigma =
      (beam_npt > 0) ? TMath::Sqrt(beam_sumsq / beam_npt) : 0.0;

  for (std::size_t k = 0; k < m_reservoir.size(); k++) {
    const TraceEvt &e = m_reservoir[k];
    Double_t td[18];
    for (Int_t s = 0; s < 18; s++)
      td[s] = Double_t(e.total[s]);
    Int_t cls = -1;
    if (e.beam_flat)
      cls = 0;
    else if (e.reac_mask & bit) {
      Double_t x = SumRange(td, kXLo, kXHi);
      Double_t y = SumRange(td, YLoOf(reac), YHiOf(reac));
      if (cut_aa && cut_aa->IsInside(x, y))
        cls = 1;
      else if (cut_an && cut_an->IsInside(x, y))
        cls = 2;
    }
    if (cls < 0)
      continue;
    // The same five variables the blind clustering uses (normed total, beam
    // at 1 per strip; guards 0/17 included).
    Double_t energy = 0.0;
    for (Int_t s = 0; s < 18; s++) {
      energy += td[s] - 1.0;
    }

    Int_t trigger_strip = FindTrigger(td, base, beam_sigma);
    Bool_t has_trig = (trigger_strip >= 0);
    // Plateau excess: beam-subtracted sum over a sliding window that tracks
    // the trigger (reacstrip+1 .. reacstrip+PLATEAU_POST). 0 when no trigger.
    // Out-of-range strips are dropped from the sum.
    const Int_t kPlateauPost = Constants::STRIP_SUM_SCATTER_CONFIG.PLATEAU_POST;
    Double_t plateau = 0.0;
    if (has_trig)
      for (Int_t d = 1; d <= kPlateauPost; d++) {
        Int_t s = trigger_strip + d;
        if (s >= 0 && s < 18)
          plateau += td[s] - 1.0;
      }
    // peak-3 sum centered on the TRIGGER strip (not the argmax). 0 when there
    // is no trigger, which is the discriminator: flat beam scores 0.
    Double_t peak3 = 0.0;
    if (has_trig)
      for (Int_t s = trigger_strip - 1; s <= trigger_strip + 1; s++)
        if (s >= 0 && s < 18)
          peak3 += td[s];
    // |deviation from beam| at the TRIGGER strip plus at the end strip s17 --
    // the (a,n) signature is a rise at the trigger AND a collapse at s17.
    Double_t trigtaildev =
        has_trig ? TMath::Abs(td[trigger_strip] - base[trigger_strip]) +
                       TMath::Abs(td[17] - base[17])
                 : 0.0;
    // Slopes ACROSS the trigger: dE(reac+n) - dE(reac-n). reac+-n share parity
    // so the L_odd/R_even sawtooth cancels. Only valid (and filled) when there
    // is a trigger AND both endpoints are in range (symmetric, no clamping).
    Bool_t ok3 =
        has_trig && (trigger_strip - 3 >= 0) && (trigger_strip + 3 <= 17);
    Double_t reacslope3 =
        ok3 ? td[trigger_strip + 3] - td[trigger_strip - 3] : 0.0;

    // How beam-like the BACK HALF (strips 8-17) is: RMS deviation of the trace
    // from the beam baseline over those strips. Subtracting base[] removes the
    // L/R sawtooth, so this is clean; LOW = beam-like (flat at beam level),
    // high for a reaction's plateau/collapse or the elevation of pileup.
    // Amplitude-aware (NOT max-normalized like the template-prune residual):
    // the beam has a fixed level, so a flat-but-elevated trace must not read as
    // beam. No trigger needed (fixed window), so always filled.
    Double_t beamdev = 0.0;
    Int_t n_bl = 0;
    for (Int_t s = 8; s <= 17; s++) {
      Double_t d = td[s] - base[s];
      beamdev += d * d;
      n_bl++;
    }
    beamdev = TMath::Sqrt(beamdev / Double_t(n_bl));
    Double_t v[NV] = {energy,
                      peak3,
                      plateau,
                      td[17], // raw end strip (was td[17] - 1.0)
                      Double_t(trigger_strip),
                      Double_t(e.both_mult),
                      trigtaildev,
                      reacslope3,
                      beamdev};
    // peak3 is filled even with no trigger (it scores 0 -- the discriminator).
    // The other trigger-centered vars are skipped when there is no trigger
    // (reacstrip) or the symmetric window runs off an edge (slopes /
    // jaggedness).
    Bool_t vok[NV] = {kTRUE, kTRUE,    kTRUE, kTRUE, kTRUE,
                      kTRUE, has_trig, ok3,   kTRUE};
    for (Int_t iv = 0; iv < NV; iv++)
      if (vok[iv])
        vals_raw[cls][iv].push_back(v[iv]);

    // Savitzky-Golay-smoothed trace: recompute trigger and cluster variables on
    // the SG-filtered copy. This removes the L_odd/R_even sawtooth so onset /
    // peak / slope features land on the real physical structure.
    Double_t sgd[18];
    SavitzkyGolay(td, sgd);
    Double_t energy_sg = 0.0;
    for (Int_t s = 0; s < 18; s++)
      energy_sg += sgd[s] - 1.0;
    Double_t ex_sg[18], sm_ex_sg[18];
    for (Int_t s = 0; s < 18; s++)
      ex_sg[s] = sgd[s] - base[s];
    SmoothTrace(ex_sg, sm_ex_sg, kClusterSmoothWindow);
    Int_t reacstrip_sg = FindTrigger(sgd, base, beam_sigma);
    Bool_t has_trig_sg = (reacstrip_sg >= 0);
    // Plateau excess on SG trace: sliding window over reacstrip_sg+1 ..
    // reacstrip_sg+PLATEAU_POST. 0 when no trigger; out-of-range strips
    // dropped.
    Double_t plateau_sg = 0.0;
    if (has_trig_sg)
      for (Int_t d = 1; d <= kPlateauPost; d++) {
        Int_t s = reacstrip_sg + d;
        if (s >= 0 && s < 18)
          plateau_sg += sgd[s] - 1.0;
      }
    Double_t peak3_sg = 0.0;
    if (has_trig_sg)
      for (Int_t s = reacstrip_sg - 1; s <= reacstrip_sg + 1; s++)
        if (s >= 0 && s < 18)
          peak3_sg += sgd[s];
    Double_t trigtaildev_sg =
        has_trig_sg ? TMath::Abs(sgd[reacstrip_sg] - base[reacstrip_sg]) +
                          TMath::Abs(sgd[17] - base[17])
                    : 0.0;
    Bool_t ok3_sg =
        has_trig_sg && (reacstrip_sg - 3 >= 0) && (reacstrip_sg + 3 <= 17);
    Double_t reacslope3_sg =
        ok3_sg ? sgd[reacstrip_sg + 3] - sgd[reacstrip_sg - 3] : 0.0;
    Double_t beamdev_sg = 0.0;
    Int_t n_bl_sg = 0;
    for (Int_t s = 8; s <= 17; s++) {
      Double_t d = sgd[s] - base[s];
      beamdev_sg += d * d;
      n_bl_sg++;
    }
    beamdev_sg = TMath::Sqrt(beamdev_sg / Double_t(n_bl_sg));
    Double_t v_sg[NV] = {energy_sg,
                         peak3_sg,
                         plateau_sg,
                         sgd[17],
                         Double_t(reacstrip_sg),
                         Double_t(e.both_mult),
                         trigtaildev_sg,
                         reacslope3_sg,
                         beamdev_sg};
    Bool_t vok_sg[NV] = {kTRUE, kTRUE,       kTRUE,  kTRUE, kTRUE,
                         kTRUE, has_trig_sg, ok3_sg, kTRUE};
    for (Int_t iv = 0; iv < NV; iv++)
      if (vok_sg[iv])
        vals_sg[cls][iv].push_back(v_sg[iv]);
  }

  std::cout << "cluster-var hists (reac " << reac
            << "): beam=" << vals_raw[0][0].size()
            << " (a,a')=" << vals_raw[1][0].size()
            << " (a,n)=" << vals_raw[2][0].size() << std::endl;

  std::vector<Int_t> colors = PlottingUtils::GetDefaultColors();
  // Two smoothing passes: raw traces, then Savitzky-Golay smoothed.
  const Int_t kNP = 2;
  const char *pass_label[kNP] = {"raw", "sg"};

  for (Int_t ip = 0; ip < kNP; ip++) {
    std::vector<Double_t>(*vals)[NC][NV] = (ip == 0) ? &vals_raw : &vals_sg;
    for (Int_t iv = 0; iv < NV; iv++) {
      Double_t lo = 1.0e30, hi = -1.0e30;
      Double_t mean[NC] = {0.0, 0.0, 0.0};
      for (Int_t ic = 0; ic < NC; ic++) {
        for (std::size_t m = 0; m < (*vals)[ic][iv].size(); m++) {
          lo = TMath::Min(lo, (*vals)[ic][iv][m]);
          hi = TMath::Max(hi, (*vals)[ic][iv][m]);
          mean[ic] += (*vals)[ic][iv][m];
        }
        if (!(*vals)[ic][iv].empty())
          mean[ic] /= Double_t((*vals)[ic][iv].size());
      }
      if (ip == 0) {
        std::cout << "  [" << pass_label[ip] << "] " << vkey[iv]
                  << ": mean beam=" << mean[0] << " (a,a')=" << mean[1]
                  << " (a,n)=" << mean[2] << " [range " << lo << ".." << hi
                  << "]" << std::endl;
      }
      Int_t nbins = 80;
      if (iv == 4) { // reaction strip: integer bins 0..17
        lo = -1.5;
        hi = 17.5;
        nbins = 19;
      } else if (iv == 5) { // both-channel multiplicity: integer bins 0..16
        lo = -0.5;
        hi = 16.5;
        nbins = 17;
      } else {
        if (!(hi > lo)) { // constant -> give it a drawable range
          lo -= 0.5;
          hi += 0.5;
        }
        Double_t pad = 0.05 * (hi - lo);
        lo -= pad;
        hi += pad;
      }

      TCanvas *c = PlottingUtils::GetConfiguredCanvas(kTRUE);
      TString axis = Form(";%s;Counts", vtitle[iv]);
      std::vector<TH1F *> hs;
      Double_t ymax = 0.0;
      for (Int_t ic = 0; ic < NC; ic++) {
        TH1F *h = new TH1F(
            Form("h_cv_%s_%s_c%d_r%d", vkey[iv], pass_label[ip], ic, reac),
            axis, nbins, lo, hi);
        h->SetDirectory(nullptr);
        for (std::size_t m = 0; m < (*vals)[ic][iv].size(); m++)
          h->Fill((*vals)[ic][iv][m]);

        PlottingUtils::ConfigureHistogram(h, colors[ic % Int_t(colors.size())],
                                          axis);
        h->SetStats(0);
        ymax = TMath::Max(ymax, h->GetMaximum());
        hs.push_back(h);
      }
      if (ymax <= 0.0)
        ymax = 1.0;
      TLegend *leg = PlottingUtils::AddLegend(0.775, 0.875, 0.70, 0.86);

      for (Int_t ic = 0; ic < NC; ic++) {
        if (ic == 0) {
          hs[0]->SetMaximum(3.0 * ymax);
          hs[0]->SetMinimum(1.0e-1);
          hs[0]->Draw("HIST");
        } else {
          hs[ic]->Draw("HIST SAME");
        }
        leg->AddEntry(hs[ic], clabel[ic], "l");
      }
      leg->Draw();

      TString sub_subdir = subdir + "/clusters_" + pass_label[ip];

      PlottingUtils::SaveFigure(
          c, Form("cluster_var_%s_%s_reac%d", vkey[iv], pass_label[ip], reac),
          sub_subdir, PlotSaveOptions::kLOG);
      for (std::size_t m = 0; m < hs.size(); m++)
        delete hs[m];
      delete c;
    }
  }
}

TString StripSumScatter::BuildFingerprint(const std::vector<Int_t> &run_order,
                                          std::map<Int_t, TChain *> &chains) {
  const Int_t kReacMin = Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax = Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Double_t kBeamFlatTol =
      Constants::STRIP_SUM_SCATTER_CONFIG.BEAM_FLAT_TOL;
  const Double_t kReacJumpMin =
      Constants::STRIP_SUM_SCATTER_CONFIG.REAC_JUMP_MIN;
  const Double_t kReacJumpMax =
      Constants::STRIP_SUM_SCATTER_CONFIG.REAC_JUMP_MAX;
  const Int_t kSmoothHiStrip =
      Constants::STRIP_SUM_SCATTER_CONFIG.REQUIRE_SMOOTHNESS_END_STRIP;
  const Double_t kSmoothMaxStep =
      Constants::STRIP_SUM_SCATTER_CONFIG.REQUIRE_SMOOTHNESS_MAX_STEP;
  const Double_t kStrip17Max = Constants::STRIP_SUM_SCATTER_CONFIG.STRIP_17_MAX;

  const Int_t kGateStripX = Constants::STRIP_SUM_SCATTER_CONFIG.GATE_STRIP_X;
  const Int_t kGateStripY = Constants::STRIP_SUM_SCATTER_CONFIG.GATE_STRIP_Y;
  const Double_t kGateNSigmaX =
      Constants::STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_X;
  const Double_t kGateNSigmaY =
      Constants::STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_Y;
  const Int_t kGateBins = Constants::STRIP_SUM_SCATTER_CONFIG.GATE_BINS;
  const Double_t kGateMin = Constants::STRIP_SUM_SCATTER_CONFIG.GATE_MIN;
  const Double_t kGateMax = Constants::STRIP_SUM_SCATTER_CONFIG.GATE_MAX;

  const Double_t kXMin = Constants::STRIP_SUM_SCATTER_CONFIG.XMIN;
  const Double_t kXMax = Constants::STRIP_SUM_SCATTER_CONFIG.XMAX;
  const Int_t kXBins = Constants::STRIP_SUM_SCATTER_CONFIG.XBINS;
  const Int_t kYBins = Constants::STRIP_SUM_SCATTER_CONFIG.YBINS;

  TString s = Form(
      "v8 reac[%d,%d] flatTol=%.3f jump[%.3f,%.3f] smooth=%d,%d "
      "step=%.3f s17=%.3f gate[s%d,s%d,%.2f,%.2f,%d,%.3f,%.3f] x[%.3f,%.3f,%d] "
      "ybins=%d",
      kReacMin, kReacMax, kBeamFlatTol, kReacJumpMin, kReacJumpMax,
      Int_t(Constants::STRIP_SUM_SCATTER_CONFIG.REQUIRE_SMOOTHNESS),
      kSmoothHiStrip, kSmoothMaxStep, kStrip17Max, kGateStripX, kGateStripY,
      kGateNSigmaX, kGateNSigmaY, kGateBins, kGateMin, kGateMax, kXMin, kXMax,
      kXBins, kYBins);
  // Active beam gates (also keyed by cache filename, but folded in here too so
  // a mismatch never silently reuses a stale same-named cache).
  std::vector<GateSpec> gates = ActiveGates();
  for (std::size_t i = 0; i < gates.size(); i++)
    s += Form(" g[s%d,s%d]", gates[i].sx, gates[i].sy);

  Double_t y_lo[64], y_hi[64];
  YBounds(y_lo, y_hi);
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++)
    s += Form(
        " y%d[%.3f,%.3f]", reac,
        y_lo[reac - Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN],
        y_hi[reac - Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN]);
  for (std::size_t i = 0; i < run_order.size(); i++) {
    Int_t run = run_order[i];
    s += Form(" r%d:%lld", run, chains[run]->GetEntries());
  }
  return s;
}

// Per-reaction-strip y-axis bounds straight from
// Constants::STRIP_SUM_SCATTER_CONFIG.Y_RANGE (tunable per dataset,
// per strip); strips absent from the map fall back to YMIN/YMAX. x stays
// fixed (strip-independent).
void StripSumScatter::YBounds(Double_t *y_lo, Double_t *y_hi) {
  const Int_t kReacMin = Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax = Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    Int_t ri = reac - kReacMin;
    std::map<Int_t, std::pair<Double_t, Double_t>>::const_iterator it =
        Constants::STRIP_SUM_SCATTER_CONFIG.Y_RANGE.find(reac);
    if (it != Constants::STRIP_SUM_SCATTER_CONFIG.Y_RANGE.end()) {
      y_lo[ri] = it->second.first;
      y_hi[ri] = it->second.second;
    } else {
      y_lo[ri] = Constants::STRIP_SUM_SCATTER_CONFIG.YMIN;
      y_hi[ri] = Constants::STRIP_SUM_SCATTER_CONFIG.YMAX;
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
  const Long64_t kSampleMaxPoints =
      Constants::STRIP_SUM_SCATTER_CONFIG.SAMPLE_MAX_POINTS;

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
  const Int_t kXLo = Constants::STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::STRIP_SUM_SCATTER_CONFIG.X_HI;

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
  Int_t y_lo = reac + 1;
  Int_t y_hi = TMath::Min(reac + 6, 17);
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
    Double_t y = SumRange(total, y_lo, y_hi);
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

// Per reaction strip, overlay TRACES_PER_CLASS sampled per-strip traces of each
// sim population in the experimental DrawRegionTraces style (beam grey, (a,a')
// azure, (a,n) red). The beam reference is the same for every strip. Sampled
// fresh each run (40 traces/file is trivial).
void StripSumScatter::SimTraceOverlay() {
  const Int_t kReacMin = Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax = Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Int_t kTracesPerRegion =
      Constants::STRIP_SUM_SCATTER_CONFIG.TRACES_PER_CLASS;

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
  const Int_t kReacMin = Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax = Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Int_t kXLo = Constants::STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::STRIP_SUM_SCATTER_CONFIG.X_HI;

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
  TString full = IO::GetRootFilesBaseDir() + TString("/") +
                 "StripSumScatter_simcache.root";
  if (gSystem->AccessPathName(full))
    return kFALSE;
  TFile *f = IO::OpenForReading("StripSumScatter_simcache.root");
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
  TFile *out = IO::OpenForWriting("StripSumScatter_simcache.root", "RECREATE");
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
void StripSumScatter::SimOverlay() {
  const Int_t kReacMin = Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax = Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Int_t kXLo = Constants::STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::STRIP_SUM_SCATTER_CONFIG.X_HI;

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
    std::map<Int_t, TH2F *>::const_iterator sit = m_scatter.find(r);
    if (sit == m_scatter.end())
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

Bool_t StripSumScatter::TryLoadCache(const TString &cacheName,
                                     const TString &fingerprint) {
  const Int_t kReacMin = Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax = Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;

  TString cache_full = IO::GetRootFilesBaseDir() + TString("/") + cacheName;
  if (gSystem->AccessPathName(cache_full)) {
    std::cout << "strip-sum-scatter: no cache file found; will rebuild."
              << std::endl;
    return kFALSE;
  }

  TFile *cf = IO::OpenForReading(cacheName);
  if (!cf || cf->IsZombie()) {
    if (cf)
      delete cf;
    std::cout << "strip-sum-scatter: cache file unreadable; rebuilding."
              << std::endl;
    return kFALSE;
  }

  TNamed *fp = dynamic_cast<TNamed *>(cf->Get("fingerprint"));
  if (!fp || fingerprint != fp->GetTitle()) {
    cf->Close();
    delete cf;
    std::cout << "strip-sum-scatter: cache present but stale; rebuilding."
              << std::endl;
    return kFALSE;
  }

  Bool_t ok = kTRUE;
  for (Int_t reac = kReacMin; reac <= kReacMax && ok; reac++) {
    TH2F *h = dynamic_cast<TH2F *>(cf->Get(Form("scatter_r%d", reac)));
    if (!h) {
      ok = kFALSE;
      break;
    }
    TH2F *hc = static_cast<TH2F *>(h->Clone());
    hc->SetDirectory(nullptr);
    m_scatter[reac] = hc;
  }

  TTree *tt = dynamic_cast<TTree *>(cf->Get("traces"));
  if (ok && tt) {
    TraceEvt e;
    tt->SetBranchAddress("total", e.total);
    tt->SetBranchAddress("total_adc", e.total_adc);
    tt->SetBranchAddress("reac_mask", &e.reac_mask);
    tt->SetBranchAddress("beam_flat", &e.beam_flat);
    tt->SetBranchAddress("both_mult", &e.both_mult);
    Long64_t nt = tt->GetEntries();
    m_reservoir.reserve(nt);
    for (Long64_t j = 0; j < nt; j++) {
      tt->GetEntry(j);
      m_reservoir.push_back(e);
    }
  }

  cf->Close();
  delete cf;

  if (ok)
    std::cout << "strip-sum-scatter: loaded cached scatters + "
              << m_reservoir.size() << " reservoir events (fingerprint match)."
              << std::endl;
  else
    std::cout << "strip-sum-scatter: cache partially corrupt; rebuilding."
              << std::endl;

  return ok;
}

void StripSumScatter::WriteCache(const TString &cacheName,
                                 const TString &fingerprint) {
  const Int_t kReacMin = Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax = Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;

  TFile *out = IO::OpenForWriting(cacheName, "RECREATE");
  if (!out || out->IsZombie()) {
    if (out)
      delete out;
    return;
  }
  out->cd();
  TNamed fp("fingerprint", fingerprint.Data());
  fp.Write();
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++)
    m_scatter[reac]->Write(Form("scatter_r%d", reac));

  TTree *tt = new TTree("traces", "strip-sum trace reservoir");
  TraceEvt e;
  tt->Branch("total", e.total, "total[18]/F");
  tt->Branch("total_adc", e.total_adc, "total_adc[18]/F");
  tt->Branch("reac_mask", &e.reac_mask, "reac_mask/i");
  tt->Branch("beam_flat", &e.beam_flat, "beam_flat/O");
  tt->Branch("both_mult", &e.both_mult, "both_mult/I");
  for (std::size_t k = 0; k < m_reservoir.size(); k++) {
    e = m_reservoir[k];
    tt->Fill();
  }
  tt->Write();
  out->Close();
  delete out;
  std::cout << "strip-sum-scatter: wrote cache " << cacheName << std::endl;
}

void StripSumScatter::FillScatters(const std::vector<Int_t> &runOrder,
                                   std::map<Int_t, TChain *> &chains) {
  const Int_t kReacMin = Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax = Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Int_t kXLo = Constants::STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::STRIP_SUM_SCATTER_CONFIG.X_HI;
  const Int_t kXBins = Constants::STRIP_SUM_SCATTER_CONFIG.XBINS;
  const Double_t kXMin = Constants::STRIP_SUM_SCATTER_CONFIG.XMIN;
  const Double_t kXMax = Constants::STRIP_SUM_SCATTER_CONFIG.XMAX;
  const Int_t kYBins = Constants::STRIP_SUM_SCATTER_CONFIG.YBINS;
  const Int_t kBeamReservoirCap =
      Constants::STRIP_SUM_SCATTER_CONFIG.TRACES_PER_CLASS * 10;

  // Compute y-bounds once.
  YBounds(m_yLo, m_yHi);

  // Allocate scatter histograms.
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    Int_t ri = ReacIndex(reac);
    TH2F *h = new TH2F(
        Form("scatter_r%d", reac),
        Form(";norm. #DeltaE strips %d#rightarrow%d [a.u.];norm. #DeltaE "
             "strips %d#rightarrow%d [a.u.]",
             kXLo, kXHi, YLoOf(reac), YHiOf(reac)),
        kXBins, kXMin, kXMax, kYBins, m_yLo[ri], m_yHi[ri]);
    h->SetDirectory(nullptr);
    h->SetStats(0);
    m_scatter[reac] = h;
  }

  std::vector<GateSpec> activeGates = ActiveGates();

  std::map<Int_t, std::vector<BeamFit2D>> gates;
  std::map<Int_t, Bool_t> gatesOk;
  for (std::size_t i = 0; i < runOrder.size(); i++) {
    Int_t run = runOrder[i];
    gatesOk[run] = kFALSE;
    if (!chains[run] || chains[run]->GetEntries() == 0)
      continue;
    std::vector<BeamFit2D> runGates;
    std::vector<GateSpec> priorSpecs;
    Bool_t allOk = kTRUE;
    for (std::size_t gi = 0; gi < activeGates.size(); gi++) {
      BeamFit2D g = FindBeamGate(
          chains[run], activeGates[gi].sx, activeGates[gi].sy, priorSpecs,
          runGates, Form("run%d", run), Form("strip_sum_scatter/run%d", run));
      if (g.ok)
        std::cout << "  run " << run << " beam gate s" << activeGates[gi].sx
                  << "/s" << activeGates[gi].sy << ": mu=(" << g.mu_x << ","
                  << g.mu_y << ")" << std::endl;
      else {
        std::cerr << "  run " << run << " beam gate s" << activeGates[gi].sx
                  << "/s" << activeGates[gi].sy << " failed; skipping run"
                  << std::endl;
        allOk = kFALSE;
      }
      runGates.push_back(g);
      priorSpecs.push_back(activeGates[gi]);
    }
    gates[run] = runGates;
    gatesOk[run] = allOk;
  }

  Long64_t totalGated = 0, totalSeen = 0;
  Int_t nBeamKept = 0;

  for (std::size_t i = 0; i < runOrder.size(); i++) {
    Int_t run = runOrder[i];
    TChain *chain = chains[run];
    if (!chain || !gatesOk[run])
      continue;
    const std::vector<BeamFit2D> &runGates = gates[run];

    EnergyView ev;
    ev.Attach(chain);
    EnableEventBranches(chain);
    Long64_t n = chain->GetEntries();
    Int_t nReac = kReacMax - kReacMin + 1;
    std::cout << "Run " << run << ": filling " << nReac
              << " reaction-strip scatters over " << n << " events..."
              << std::endl;

    for (Long64_t j = 0; j < n; j++) {
      chain->GetEntry(j);
      ev.Decode();
      totalSeen++;

      Bool_t passesAll = kTRUE;
      for (std::size_t gi = 0; gi < activeGates.size(); gi++)
        if (!PassesGate(runGates[gi], ev, activeGates[gi].sx,
                        activeGates[gi].sy)) {
          passesAll = kFALSE;
          break;
        }
      if (!passesAll)
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
        m_scatter[reac]->Fill(x, SumRange(ev.total, YLoOf(reac), YHiOf(reac)));
      }

      // Keep every reaction-passing event for traces; cap clean flat-beam
      // events (only ~TRACES_PER_CLASS are ever drawn). The two are mutually
      // exclusive -- a flat-beam event has no reaction jump.
      Bool_t beam = (mask == 0) && IsBeamFlat(ev);
      if (mask == 0 && !(beam && nBeamKept < kBeamReservoirCap))
        continue;
      if (beam)
        nBeamKept++;

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
      // Both-channel multiplicity: split strips (1-16) where both ends
      // FIRED. Read off the RAW ADC, not the calibrated ends -- the
      // short-end gains are 0 (uncalibrated, no sim anchor), so the
      // calibrated short ends are always zero; the raw ADC still carries
      // whether the channel fired.
      Int_t both = 0;
      for (Int_t s = 1; s <= 16; s++)
        if (ev.left_0_17_adc[s] > 0.0 && ev.rightdE_adc[s] > 0.0)
          both++;
      e.both_mult = both;
      e.reac_mask = mask;
      e.beam_flat = beam;
      m_reservoir.push_back(e);
      if (mask != 0)
        totalGated++;
    }
  }

  std::cout << "Built scatters: " << totalGated
            << " reaction events across strips " << kReacMin << "-" << kReacMax
            << " (" << totalSeen << " seen), reservoir " << m_reservoir.size()
            << " events." << std::endl;
}

void StripSumScatter::PlotScatters() {
  const Int_t kXLo = Constants::STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::STRIP_SUM_SCATTER_CONFIG.X_HI;
  const Int_t kReacMin = Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax = Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;

  std::lock_guard<std::mutex> lock(g_plot_mutex);
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(m_scatter[reac], c);
    m_scatter[reac]->GetYaxis()->SetTitleOffset(1.3);
    c->SetLeftMargin(0.18);
    PlottingUtils::SaveFigure(c,
                              Form("normsumE_reac%d_s%d_%d_vs_s%d_%d", reac,
                                   YLoOf(reac), YHiOf(reac), kXLo, kXHi),
                              "strip_sum_scatter", PlotSaveOptions::kLINEAR);
    delete c;
  }
}

void StripSumScatter::InteractiveOverlay(Int_t reac) {
  const Int_t kReacMin = Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax = Constants::STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Int_t kXLo = Constants::STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::STRIP_SUM_SCATTER_CONFIG.X_HI;
  const Int_t kTracesPerRegion =
      Constants::STRIP_SUM_SCATTER_CONFIG.TRACES_PER_CLASS;

  if (reac < kReacMin || reac > kReacMax) {
    std::cerr << "strip-sum-scatter: candidate reaction strip " << reac
              << " outside [" << kReacMin << "," << kReacMax
              << "]; skipping interactive overlay." << std::endl;
    return;
  }
  if (!gSystem->Getenv("DISPLAY")) {
    std::cerr << "strip-sum-scatter: no DISPLAY; skipping interactive "
                 "region-trace overlay (scatters already saved)."
              << std::endl;
    return;
  }

  Int_t app_argc = 1;
  char app_arg0[] = "strip-sum-scatter";
  char *app_argv[] = {app_arg0};
  TApplication app("strip-sum-scatter", &app_argc, app_argv);
  gROOT->SetBatch(kFALSE);

  TCanvas *cutCanvas = new TCanvas("c_strip_sum_regions",
                                   "Draw (a,n) then (a,a') regions", 900, 700);
  cutCanvas->SetLogz(kTRUE); // match the saved scatter's z-scale
  m_scatter[reac]->Draw("COLZ");
  cutCanvas->Update();
  TCutG *cutAn = PromptCut(cutCanvas, "region_an", "(a,n)");
  TCutG *cutAa = PromptCut(cutCanvas, "region_aa", "(a,a')");

  cutCanvas->GetListOfPrimitives()->Remove(cutAn);
  cutCanvas->GetListOfPrimitives()->Remove(cutAa);
  gROOT->SetEditorMode();
  delete cutCanvas;
  gSystem->ProcessEvents();
  gROOT->SetBatch(kTRUE);

  ClusterVarHists(reac, cutAa, cutAn, "strip_sum_scatter");

  std::vector<TGraph *> tr_an, tr_aa, tr_beam;
  // Same selected events, raw (un-normalized) ADC -- one entry per normed
  // trace, kept in lock-step so the two overlays show the identical events.
  std::vector<TGraph *> tr_an_adc, tr_aa_adc, tr_beam_adc;
  // Same selected events again, Savitzky-Golay smoothed (normed a.u. space),
  // for the with-smoothing overlay -- also in lock-step with the raw normed
  // traces, so the two a.u. overlays show the identical events.
  std::vector<TGraph *> tr_an_sg, tr_aa_sg, tr_beam_sg;
  UInt_t bit = (1u << ReacIndex(reac));

  for (std::size_t k = 0; k < m_reservoir.size(); k++) {
    if (tr_an.size() >= std::size_t(kTracesPerRegion) &&
        tr_aa.size() >= std::size_t(kTracesPerRegion) &&
        tr_beam.size() >= std::size_t(kTracesPerRegion))
      break;
    const TraceEvt &e = m_reservoir[k];
    if (e.beam_flat && tr_beam.size() < std::size_t(kTracesPerRegion)) {
      tr_beam.push_back(TraceFromTotal(e.total));
      tr_beam_adc.push_back(TraceFromTotal(e.total_adc));
      tr_beam_sg.push_back(SmoothedTraceFromTotal(e.total));
      continue;
    }
    if (!(e.reac_mask & bit))
      continue;
    Double_t td[18];
    for (Int_t s = 0; s < 18; s++)
      td[s] = Double_t(e.total[s]);
    Double_t x = SumRange(td, kXLo, kXHi);
    Double_t y = SumRange(td, YLoOf(reac), YHiOf(reac));
    if (cutAn && tr_an.size() < std::size_t(kTracesPerRegion) &&
        cutAn->IsInside(x, y)) {
      tr_an.push_back(TraceFromTotal(e.total));
      tr_an_adc.push_back(TraceFromTotal(e.total_adc));
      tr_an_sg.push_back(SmoothedTraceFromTotal(e.total));
    } else if (cutAa && tr_aa.size() < std::size_t(kTracesPerRegion) &&
               cutAa->IsInside(x, y)) {
      tr_aa.push_back(TraceFromTotal(e.total));
      tr_aa_adc.push_back(TraceFromTotal(e.total_adc));
      tr_aa_sg.push_back(SmoothedTraceFromTotal(e.total));
    }
  }

  std::cout << "Sampled traces: beam=" << tr_beam.size()
            << " (a,a')=" << tr_aa.size() << " (a,n)=" << tr_an.size()
            << std::endl;
  DrawRegionTraces(Form("region_traces_reac%d", reac), "strip_sum_scatter",
                   tr_beam, tr_aa, tr_an);
  DrawRegionMeanTraces(Form("region_mean_traces_reac%d", reac),
                       "strip_sum_scatter", tr_beam, tr_aa, tr_an);
  Double_t adc_y_lo = 0.0, adc_y_hi = 0.0;
  TraceYRange(tr_beam_adc, tr_aa_adc, tr_an_adc, adc_y_lo, adc_y_hi);
  DrawRegionTraces(Form("region_traces_reac%d_adc", reac), "strip_sum_scatter",
                   tr_beam_adc, tr_aa_adc, tr_an_adc, adc_y_lo, adc_y_hi,
                   "#DeltaE [ADC]");
  DrawRegionMeanTraces(Form("region_mean_traces_reac%d_adc", reac),
                       "strip_sum_scatter", tr_beam_adc, tr_aa_adc, tr_an_adc,
                       adc_y_lo, adc_y_hi, "#DeltaE [ADC]");
  DrawRegionTraces(Form("region_traces_reac%d_sg", reac), "strip_sum_scatter",
                   tr_beam_sg, tr_aa_sg, tr_an_sg);
  DrawRegionMeanTraces(Form("region_mean_traces_reac%d_sg", reac),
                       "strip_sum_scatter", tr_beam_sg, tr_aa_sg, tr_an_sg);

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
  for (std::size_t i = 0; i < tr_an_sg.size(); i++)
    delete tr_an_sg[i];
  for (std::size_t i = 0; i < tr_aa_sg.size(); i++)
    delete tr_aa_sg[i];
  for (std::size_t i = 0; i < tr_beam_sg.size(); i++)
    delete tr_beam_sg[i];

  delete cutAn;
  delete cutAa;
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

  // Build fingerprint and try cache.
  TString fingerprint = BuildFingerprint(run_order, chain_by_run);
  TString cache_name = CacheName();

  Bool_t loaded = TryLoadCache(cache_name, fingerprint);

  if (!loaded) {
    FillScatters(run_order, chain_by_run);
    WriteCache(cache_name, fingerprint);
  }

  // Batch plotting (always done).
  PlotScatters();

  // Optional sim overlays.
  if (Constants::STRIP_SUM_SCATTER_CONFIG.RERUN_SIM) {
    SimOverlay();
    SimTraceOverlay();
  }

  // Interactive region-trace overlay (requires DISPLAY).
  Int_t reac = Constants::STRIP_SUM_SCATTER_CONFIG.CANDIDATE_REAC_STRIP;
  InteractiveOverlay(reac);

  // Cleanup chains.
  for (std::size_t i = 0; i < run_order.size(); i++)
    delete chain_by_run[run_order[i]];
}
