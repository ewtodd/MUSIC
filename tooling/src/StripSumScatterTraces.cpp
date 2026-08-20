#include "StripSumScatter.hpp"

void StripSumScatter::DrawTraceSet(const std::vector<TGraph *> &traces,
                                   Int_t color) {
  for (Int_t i = 0; i < Int_t(traces.size()); i++) {
    traces[i]->SetLineColor(color);
    traces[i]->SetLineWidth(1);
    traces[i]->Draw("L SAME");
  }
}

TGraph *StripSumScatter::TraceFromTotal(const Float_t *total) {
  Double_t td[18];
  for (Int_t s = 0; s < 18; s++)
    td[s] = Double_t(total[s]);
  return EventsSummary::BuildTraceFromTotals(td);
}

void StripSumScatter::DrawRegionTraces(const TString &save_name,
                                       const TString &subdir,
                                       const std::vector<TGraph *> &beam,
                                       const std::vector<TGraph *> &aa,
                                       const std::vector<TGraph *> &an,
                                       Double_t y_min, Double_t y_max,
                                       const char *y_title) {
  std::lock_guard<std::mutex> lock(g_plot_mutex);
  Int_t s_lo = Constants::cfg.IGNORE_STRIP_0 ? 1 : 0;
  Int_t s_hi = Constants::cfg.IGNORE_STRIP_17 ? 16 : 17;
  TH2F *frame =
      new TH2F("h_region_trace_frame", Form(";Strip;%s", y_title),
               s_hi - s_lo + 1, s_lo - 0.5, s_hi + 0.5, 100, y_min, y_max);
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
    for (Int_t t = 0; t < Int_t(tr.size()); t++) {
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
  for (Int_t i = 0; i < Int_t(means.size()); i++)
    means[i]->Draw("LX SAME"); // mean line, no end caps
  leg->Draw();

  PlottingUtils::SaveFigure(c, save_name, subdir, PlotSaveOptions::kLINEAR);
  for (Int_t i = 0; i < Int_t(means.size()); i++)
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
    for (Int_t i = 0; i < Int_t(v.size()); i++) {
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
  TCutG *cut = static_cast<TCutG *>(c->WaitPrimitive("CUTG", "CutG"));
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

// Savitzky-Golay: 3rd-degree, 5-point window, coeff [-3,12,17,12,-3]/35.
// At edges the window shrinks and the coefficients are renormalised.
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

// CFD-style trigger: first strip (left-to-right) whose beam-subtracted
// signal exceeds a peak fraction and an N-sigma gate; -1 if none fires.

Int_t StripSumScatter::FindTrigger(const Double_t *td, const Double_t *base,
                                   Double_t beam_sigma) {
  const Int_t s_lo = 2;
  const Int_t s_hi = 16;

  const Double_t frac =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TRIGGER_CFD_FRAC;

  const Double_t nsigma =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TRIGGER_NSIGMA * beam_sigma;

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
  return EventsSummary::BuildTraceFromTotals(sgd);
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

  const Int_t kXLo = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_HI;
  const Int_t kClusterSmoothWindow =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.CLUSTER_SMOOTH_WINDOW;
  const Bool_t kScatterSavgol =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.SCATTER_SAVGOL;

  // Per (class, variable) value lists for raw traces.
  std::vector<Double_t> vals_raw[NC][NV];
  // Same, but cluster variables on SG-smoothed traces (the kernel removes the
  // L_odd/R_even sawtooth so features land on the real physical peak).
  std::vector<Double_t> vals_sg[NC][NV];
  UInt_t bit = (1u << ReacIndex(reac));

  // Per-strip beam baseline = mean of the beam-flat reservoir events. It
  // carries the L_odd/R_even sawtooth, so subtraction removes it exactly before
  // the onset search (better than blurring).
  Double_t base[18];
  for (Int_t s = 0; s < 18; s++)
    base[s] = 0.0;
  Long64_t nbeam = 0;
  for (Int_t k = 0; k < Int_t(m_reservoir.size()); k++)
    if (m_reservoir[k].beam_flat) {
      for (Int_t s = 0; s < 18; s++)
        base[s] += Double_t(m_reservoir[k].total[s]);
      nbeam++;
    }
  for (Int_t s = 0; s < 18; s++)
    base[s] = (nbeam > 0) ? base[s] / Double_t(nbeam) : 1.0;

  // Pooled beam-noise RMS = sqrt(mean (total-base)^2 over beam-flat events and
  // all strips); the onset gate is N sigma of this, so flat beam does not
  // trigger on average (matches Python).
  Double_t beam_sumsq = 0.0;
  Long64_t beam_npt = 0;
  for (Int_t k = 0; k < Int_t(m_reservoir.size()); k++)
    if (m_reservoir[k].beam_flat)
      for (Int_t s = 0; s < 18; s++) {
        Double_t d = Double_t(m_reservoir[k].total[s]) - base[s];
        beam_sumsq += d * d;
        beam_npt++;
      }
  Double_t beam_sigma =
      (beam_npt > 0) ? TMath::Sqrt(beam_sumsq / beam_npt) : 0.0;

  // Count triggers over ALL reservoir events (not just classified subset) so
  // the numbers are directly comparable to the Python pipeline.
  Long64_t triggered = 0, no_trigger = 0;
  Double_t td_all[18];
  for (Int_t k = 0; k < Int_t(m_reservoir.size()); k++) {
    for (Int_t s = 0; s < 18; s++)
      td_all[s] = Double_t(m_reservoir[k].total[s]);
    if (FindTrigger(td_all, base, beam_sigma) >= 0)
      triggered++;
    else
      no_trigger++;
  }
  const Double_t reac_onset_gate =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TRIGGER_NSIGMA * beam_sigma;
  const Double_t cf_frac =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TRIGGER_CFD_FRAC;
  std::cout << "  beam reference: mean+RMS of " << nbeam
            << " pure-beam events (fitted s0,s1 & s16,s17 ellipses); "
            << Form("noise sigma=%.4f", beam_sigma) << std::endl;
  std::cout << Form("  reaction onset: gate %g-sigma = %.4f, CF fraction %g; "
                    "triggered %lld of %lld (no trigger: %lld)",
                    Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TRIGGER_NSIGMA,
                    reac_onset_gate, cf_frac, (long long)triggered,
                    (long long)m_reservoir.size(), (long long)no_trigger)
            << std::endl;

  for (Int_t k = 0; k < Int_t(m_reservoir.size()); k++) {
    const TraceEvt &e = m_reservoir[k];
    Double_t td[18];
    for (Int_t s = 0; s < 18; s++)
      td[s] = Double_t(e.total[s]);
    Int_t cls = -1;
    if (e.beam_flat)
      cls = 0;
    else if (e.reac_mask & bit) {
      // Cuts are drawn on the scatter, which uses SG-smoothed sums when
      // SCATTER_SAVGOL; region classification must use that same space (cluster
      // vars below still use raw).
      Double_t td_sg[18];
      const Double_t *coords = td;
      if (kScatterSavgol) {
        SavitzkyGolay(td, td_sg);
        coords = td_sg;
      }
      Double_t x = SumRange(coords, kXLo, kXHi);
      Double_t y = SumRange(coords, YLoOf(reac), YHiOf(reac));
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
    // Post-trigger excess: beam-subtracted sum over trigger+1 .. trigger+
    // POST_TRIGGER_SUM_STRIPS; 0 when no trigger, out-of-range strips dropped.
    const Int_t kPostTrig =
        Constants::cfg.STRIP_SUM_SCATTER_CONFIG.POST_TRIGGER_SUM_STRIPS;
    Double_t plateau = 0.0;
    if (has_trig)
      for (Int_t d = 1; d <= kPostTrig; d++) {
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
    // Slope across the trigger: dE(reac+3) - dE(reac-3); the endpoints share
    // parity so the sawtooth cancels. Filled only with an in-range trigger (no
    // clamping).
    Bool_t ok3 =
        has_trig && (trigger_strip - 3 >= 0) && (trigger_strip + 3 <= 17);
    Double_t reacslope3 =
        ok3 ? td[trigger_strip + 3] - td[trigger_strip - 3] : 0.0;

    // Beam-likeness of the back half (strips 8-17): RMS of (trace-base) there
    // (base subtraction kills the sawtooth). LOW = beam-flat, high = reaction
    // or pileup. Amplitude-aware: a flat-but-elevated trace must not read as
    // beam; fixed window, so always filled.
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
    // peak3 fills even without a trigger (it scores 0 -- the discriminator);
    // other trigger-centered vars skip without a trigger or window off-edge.
    Bool_t vok[NV] = {kTRUE, kTRUE,    kTRUE, kTRUE, kTRUE,
                      kTRUE, has_trig, ok3,   kTRUE};
    for (Int_t iv = 0; iv < NV; iv++)
      if (vok[iv])
        vals_raw[cls][iv].push_back(v[iv]);

    // SG-smoothed pass: recompute trigger + cluster variables on the smoothed
    // copy (the kernel removes the sawtooth so features land on real
    // structure).
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
    // Post-trigger excess on SG trace (same sliding window as above); 0 when
    // no trigger, out-of-range strips dropped.
    Double_t plateau_sg = 0.0;
    if (has_trig_sg)
      for (Int_t d = 1; d <= kPostTrig; d++) {
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
    if (ip == 1 && Constants::cfg.STRIP_SUM_SCATTER_CONFIG.SKIP_SAVGOL_PLOTS)
      continue;
    std::vector<Double_t>(*vals)[NC][NV] = (ip == 0) ? &vals_raw : &vals_sg;
    for (Int_t iv = 0; iv < NV; iv++) {
      Double_t lo = 1.0e30, hi = -1.0e30;
      Double_t mean[NC] = {0.0, 0.0, 0.0};
      for (Int_t ic = 0; ic < NC; ic++) {
        for (Int_t m = 0; m < Int_t((*vals)[ic][iv].size()); m++) {
          lo = TMath::Min(lo, (*vals)[ic][iv][m]);
          hi = TMath::Max(hi, (*vals)[ic][iv][m]);
          mean[ic] += (*vals)[ic][iv][m];
        }
        if (!(*vals)[ic][iv].empty())
          mean[ic] /= Double_t((*vals)[ic][iv].size());
      }
      // No class has values for this variable (e.g. no triggers anywhere): skip
      // it instead of building a TH1F over the un-pulled [1e30, -1e30] range.
      Bool_t all_empty = kTRUE;
      for (Int_t ic = 0; ic < NC; ic++)
        if (!(*vals)[ic][iv].empty())
          all_empty = kFALSE;
      if (all_empty) {
        std::cout << "  [" << pass_label[ip] << "] " << vkey[iv]
                  << ": no events in any class; skipping" << std::endl;
        continue;
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
        for (Int_t m = 0; m < Int_t((*vals)[ic][iv].size()); m++)
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
      for (Int_t m = 0; m < Int_t(hs.size()); m++)
        delete hs[m];
      delete c;
    }
  }
}
