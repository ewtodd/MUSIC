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
  return reac - Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
}

Int_t StripSumScatter::YLoOf(Int_t reac) { return reac + 1; }

Int_t StripSumScatter::YHiOf(Int_t reac) {
  const Int_t kPost =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.POST_TRIGGER_SUM_STRIPS;
  return TMath::Min(reac + kPost, 17);
}

void StripSumScatter::EnableEventBranches(TChain *chain) {
  chain->SetBranchStatus("*", 0);
  chain->SetBranchStatus("Left_0_17_dE", 1);
  chain->SetBranchStatus("RightdE", 1);
  chain->SetBranchStatus("Cathode", 1);
}

Bool_t StripSumScatter::AllStripsFired(const EnergyView &ev) {
  if (!Constants::cfg.IGNORE_STRIP_0 && !(ev.total[0] > 0.0))
    return kFALSE;
  if (!Constants::cfg.IGNORE_STRIP_17 && !(ev.total[17] > 0.0))
    return kFALSE;
  for (Int_t s = 1; s <= 16; s++)
    if (!(ev.total[s] > 0.0))
      return kFALSE;
  return kTRUE;
}

Bool_t StripSumScatter::PassesReaction(const EnergyView &ev, Int_t reac) {
  const Double_t kReacJumpMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REAC_JUMP_MIN;
  const Double_t kReacJumpMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REAC_JUMP_MAX;
  const Double_t kSmoothMaxStep =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REQUIRE_SMOOTHNESS_MAX_STEP;
  const Int_t kSmoothHiStrip =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REQUIRE_SMOOTHNESS_END_STRIP;
  const Double_t kEndStripMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.END_STRIP_MAX;

  if (!AllStripsFired(ev))
    return kFALSE;
  Double_t reac_jump = ev.total[reac] - ev.total[reac - 1];
  if (!(reac_jump > kReacJumpMin && reac_jump < kReacJumpMax))
    return kFALSE;
  if (!(ev.total[reac] > 1.0 + kReacJumpMin &&
        ev.total[reac] < 1.0 + kReacJumpMax))
    return kFALSE;
  if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REQUIRE_SMOOTHNESS)
    for (Int_t s = reac + 1; s <= kSmoothHiStrip; s++)
      if (TMath::Abs(ev.total[s] - ev.total[s - 1]) > kSmoothMaxStep)
        return kFALSE;
  if (Constants::cfg.IGNORE_STRIP_17)
    return ev.total[16] < kEndStripMax;
  Int_t end_strip =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REQUIRE_STRIP_16_BELOW_BEAM ? 16
                                                                          : 17;
  return ev.total[end_strip] < kEndStripMax;
}

Bool_t StripSumScatter::IsPureBeam(const EnergyView &ev,
                                   const BeamEllipses &be) {
  if (!be.ok)
    return kFALSE;
  if (!AllStripsFired(ev))
    return kFALSE;
  // Must pass BOTH the entrance AND exit ellipses.
  if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PURE_BEAM_GATE ==
      StripSumScatterConfig::PURE_BEAM_GATE_S1_S2) {
    if (!PassesGate(be.s1_s2, ev, 1, 2))
      return kFALSE;
  } else {
    if (!PassesGate(be.s0_s1, ev, 0, 1))
      return kFALSE;
  }
  if (be.use_s15_s16) {
    if (!PassesGate(be.s15_s16, ev, 15, 16))
      return kFALSE;
  } else {
    if (!PassesGate(be.s16_s17, ev, 16, 17))
      return kFALSE;
  }
  return kTRUE;
}

Bool_t StripSumScatter::IsPileup(const EnergyView &ev) {
  if (!Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REJECT_PILEUP)
    return kFALSE;
  const Double_t kThresh =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PILEUP_THRESHOLD;
  const Int_t kMinStrips =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PILEUP_MIN_STRIPS;
  Int_t n = 0;
  for (Int_t s = 1; s <= 16; s++)
    if (ev.total[s] >= kThresh && ++n >= kMinStrips)
      return kTRUE;
  return kFALSE;
}

Bool_t StripSumScatter::IsNoise(const EnergyView &ev) {
  if (!Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REJECT_NOISE)
    return kFALSE;
  const Double_t kThresh =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.NOISE_THRESHOLD;
  const Int_t kMinStrips =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.NOISE_MIN_STRIPS;
  Int_t n = 0;
  for (Int_t s = 1; s <= 16; s++)
    if (ev.total[s] <= kThresh && ++n >= kMinStrips)
      return kTRUE;
  return kFALSE;
}

Bool_t StripSumScatter::IsHighStrip(const EnergyView &ev) {
  if (!Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REJECT_HIGH_STRIP)
    return kFALSE;
  const Double_t kCap =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.HIGH_STRIP_THRESHOLD;
  for (Int_t s = 1; s <= 16; s++)
    if (ev.total[s] > kCap)
      return kTRUE;
  return kFALSE;
}

Bool_t StripSumScatter::IsOffbeam(const EnergyView &ev) {
  const Double_t kDist = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.OFFBEAM_DIST;
  const Int_t kMinStrips =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.OFFBEAM_MIN_STRIPS;
  // Use flat normed beam level (1.0) as reference, matching Python which
  // computes beam reference after these pre-filters.
  const Double_t kBeamRef = 1.0;
  Int_t n = 0;
  for (Int_t s = 1; s <= 16; s++)
    if (TMath::Abs(ev.total[s] - kBeamRef) >= kDist && ++n >= kMinStrips)
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
  g.sx = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_STRIP_X;
  g.sy = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_STRIP_Y;
  gates.push_back(g);
  if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REQUIRE_GATE_S3_S4) {
    g.sx = 3;
    g.sy = 4;
    gates.push_back(g);
  }
  if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REQUIRE_GATE_S5_S6) {
    g.sx = 5;
    g.sy = 6;
    gates.push_back(g);
  }
  return gates;
}

TString StripSumScatter::CacheName() {
  TString name = "StripSumScatter_cache";
  if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REQUIRE_GATE_S3_S4)
    name += "_g34";
  if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REQUIRE_GATE_S5_S6)
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
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_X;
  const Double_t kGateNSigmaY =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_Y;
  return BeamFitUtils::InEllipseXY(gate, g0, g1, kGateNSigmaX, kGateNSigmaY);
}

BeamFit2D
StripSumScatter::FindBeamGate(TChain *chain, Int_t sx, Int_t sy,
                              const std::vector<GateSpec> &prior_specs,
                              const std::vector<BeamFit2D> &prior_gates,
                              const TString &tag, const TString &subdir) {
  const Int_t kGateBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_BINS;
  const Double_t kGateMin = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_MIN;
  const Double_t kGateMax = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_MAX;
  const Int_t kSeedHalfBins =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.SEED_HALF_BINS;
  const Double_t kSeedFrac = 0.3;
  const Long64_t kSampleMaxPoints =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.SAMPLE_MAX_POINTS;

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
    for (Int_t gi = 0; gi < prior_specs.size(); gi++)
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
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_X;
  const Double_t kGateNSigmaY =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_Y;

  std::lock_guard<std::mutex> lock(g_plot_mutex);
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureAndDraw2DHistogram(h, c);
  // Correlated 2D Gaussian ellipse matching the InEllipseXY gate.
  Double_t sxx = out.sigma_x * out.sigma_x;
  Double_t syy = out.sigma_y * out.sigma_y;
  Double_t sxy = out.rho * out.sigma_x * out.sigma_y;
  Double_t sum = sxx + syy;
  Double_t diff = sxx - syy;
  Double_t det = TMath::Sqrt(diff * diff + 4.0 * sxy * sxy);
  Double_t lambda1 = 0.5 * (sum + det);
  Double_t lambda2 = 0.5 * (sum - det);
  Double_t theta = 0.5 * TMath::ATan2(2.0 * sxy, diff) * 180.0 / TMath::Pi();
  TEllipse *e =
      new TEllipse(out.mu_x, out.mu_y, kGateNSigmaX * TMath::Sqrt(lambda1),
                   kGateNSigmaX * TMath::Sqrt(lambda2), 0, 360, theta);
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

TString StripSumScatter::BuildFingerprint(const std::vector<Int_t> &run_order,
                                          std::map<Int_t, TChain *> &chains) {
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Double_t kReacJumpMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REAC_JUMP_MIN;
  const Double_t kReacJumpMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REAC_JUMP_MAX;
  const Int_t kSmoothHiStrip =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REQUIRE_SMOOTHNESS_END_STRIP;
  const Double_t kSmoothMaxStep =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REQUIRE_SMOOTHNESS_MAX_STEP;
  const Double_t kEndStripMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.END_STRIP_MAX;

  const Int_t kGateStripX =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_STRIP_X;
  const Int_t kGateStripY =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_STRIP_Y;
  const Double_t kGateNSigmaX =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_X;
  const Double_t kGateNSigmaY =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_Y;
  const Int_t kGateBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_BINS;
  const Double_t kGateMin = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_MIN;
  const Double_t kGateMax = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_MAX;

  const Double_t kXMin = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.XMIN;
  const Double_t kXMax = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.XMAX;
  const Int_t kXBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.XBINS;
  const Int_t kYBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.YBINS;

  // Pre-fill event filters: any change here must invalidate the cached
  // scatters/reservoir (they were filled under the old filter settings).
  const Bool_t kRejNoise = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REJECT_NOISE;
  const Double_t kNoiseThresh =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.NOISE_THRESHOLD;
  const Int_t kNoiseMinStrips =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.NOISE_MIN_STRIPS;
  const Bool_t kRejPileup =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REJECT_PILEUP;
  const Double_t kPileupThresh =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PILEUP_THRESHOLD;
  const Int_t kPileupMinStrips =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PILEUP_MIN_STRIPS;
  const Bool_t kRejOffbeam =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REJECT_OFFBEAM;
  const Double_t kOffbeamDist =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.OFFBEAM_DIST;
  const Int_t kOffbeamMinStrips =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.OFFBEAM_MIN_STRIPS;
  const Bool_t kRejHighStrip =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REJECT_HIGH_STRIP;
  const Double_t kHighStripCap =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.HIGH_STRIP_THRESHOLD;

  // Remaining knobs that change the cached scatter/reservoir contents:
  // x-axis strip span (X_LO..X_HI), beam-classification settings (entrance
  // gate choice, exit strip17 handling), the reservoir cap, and the
  // reaction/beam event predicates' top-level toggles.
  const Int_t kXLo = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_HI;
  const Int_t kPureBeamGate =
      Int_t(Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PURE_BEAM_GATE);
  const Int_t kTracesPerClass =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TRACES_PER_CLASS;
  const Bool_t kReqS16BelowBeam =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REQUIRE_STRIP_16_BELOW_BEAM;
  const Bool_t kIgnStrip0 = Constants::cfg.IGNORE_STRIP_0;
  const Bool_t kIgnStrip17 = Constants::cfg.IGNORE_STRIP_17;
  const Bool_t kIgnShort = Constants::cfg.IGNORE_SHORT_STRIPS;

  TString s = Form(
      "v12 reac[%d,%d] jump[%.3f,%.3f] smooth=%d,%d "
      "step=%.3f s17=%.3f gate[s%d,s%d,%.2f,%.2f,%d,%.3f,%.3f] x[%.3f,%.3f,%d] "
      "ybins=%d filt[noise:%d,%.3f,%d pileup:%d,%.3f,%d offbeam:%d,%.3f,%d "
      "high:%d,%.3f] "
      "xsum[%d,%d] beamgate=%d traces=%d s16below=%d ign[s0:%d,s17:%d,short:%d]",
      kReacMin, kReacMax, kReacJumpMin, kReacJumpMax,
      Int_t(Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REQUIRE_SMOOTHNESS),
      kSmoothHiStrip, kSmoothMaxStep, kEndStripMax, kGateStripX, kGateStripY,
      kGateNSigmaX, kGateNSigmaY, kGateBins, kGateMin, kGateMax, kXMin, kXMax,
      kXBins, kYBins, kRejNoise, kNoiseThresh, kNoiseMinStrips, kRejPileup,
      kPileupThresh, kPileupMinStrips, kRejOffbeam, kOffbeamDist,
      kOffbeamMinStrips, kRejHighStrip, kHighStripCap, kXLo, kXHi,
      kPureBeamGate, kTracesPerClass, kReqS16BelowBeam, kIgnStrip0, kIgnStrip17,
      kIgnShort);
  // Active beam gates (also keyed by cache filename, but folded in here too so
  // a mismatch never silently reuses a stale same-named cache).
  std::vector<GateSpec> gates = ActiveGates();
  for (Int_t i = 0; i < gates.size(); i++)
    s += Form(" g[s%d,s%d]", gates[i].sx, gates[i].sy);

  Double_t y_lo[64], y_hi[64];
  YBounds(y_lo, y_hi);
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++)
    s += Form(
        " y%d[%.3f,%.3f]", reac,
        y_lo[reac - Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN],
        y_hi[reac -
             Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN]);
  for (Int_t i = 0; i < run_order.size(); i++) {
    Int_t run = run_order[i];
    s += Form(" r%d:%lld", run, chains[run]->GetEntries());
  }
  return s;
}

// Per-reaction-strip y-axis bounds straight from
// Constants::cfg.STRIP_SUM_SCATTER_CONFIG.Y_RANGE (tunable per dataset,
// per strip); strips absent from the map fall back to YMIN/YMAX. x stays
// fixed (strip-independent).
void StripSumScatter::YBounds(Double_t *y_lo, Double_t *y_hi) {
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    Int_t ri = reac - kReacMin;
    std::map<Int_t, std::pair<Double_t, Double_t>>::const_iterator it =
        Constants::cfg.STRIP_SUM_SCATTER_CONFIG.Y_RANGE.find(reac);
    if (it != Constants::cfg.STRIP_SUM_SCATTER_CONFIG.Y_RANGE.end()) {
      y_lo[ri] = it->second.first;
      y_hi[ri] = it->second.second;
    } else {
      y_lo[ri] = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.YMIN;
      y_hi[ri] = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.YMAX;
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

Bool_t StripSumScatter::TryLoadCache(const TString &cacheName,
                                     const TString &fingerprint) {
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;

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

  TNamed *fp = static_cast<TNamed *>(cf->Get("fingerprint"));
  if (!fp || fingerprint != fp->GetTitle()) {
    cf->Close();
    delete cf;
    std::cout << "strip-sum-scatter: cache present but stale; rebuilding."
              << std::endl;
    return kFALSE;
  }

  Bool_t ok = kTRUE;
  for (Int_t reac = kReacMin; reac <= kReacMax && ok; reac++) {
    TH2F *h = static_cast<TH2F *>(cf->Get(Form("scatter_r%d", reac)));
    if (!h) {
      ok = kFALSE;
      break;
    }
    TH2F *hc = static_cast<TH2F *>(h->Clone());
    hc->SetDirectory(nullptr);
    m_scatter[reac] = hc;
  }

  TTree *tt = static_cast<TTree *>(cf->Get("traces"));
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
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;

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
  for (Int_t k = 0; k < m_reservoir.size(); k++) {
    e = m_reservoir[k];
    tt->Fill();
  }
  tt->Write();
  out->Close();
  delete out;
  std::cout << "strip-sum-scatter: wrote cache " << cacheName << std::endl;
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
  if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.RERUN_SIM) {
    SimOverlay();
    SimTraceOverlay();
  }

  // Interactive region-trace overlay (requires DISPLAY).
  Int_t reac = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.CANDIDATE_REAC_STRIP;
  InteractiveOverlay(reac);

  // Cleanup chains.
  for (Int_t i = 0; i < run_order.size(); i++)
    delete chain_by_run[run_order[i]];
}
