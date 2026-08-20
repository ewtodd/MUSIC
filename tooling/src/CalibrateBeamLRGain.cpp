#include "CalibrateBeamInternal.hpp"

const Long64_t kSampleCap = 20000;
// (L, R) pair cap — larger than kSampleCap because the shoulder slice hits
// only events where the long side reads low.
const Long64_t kPairCap = 100000;

// Shoulder slice collects (short, long) pairs where LONG side reads below
// absolute cap (superset of relative cut: long < 0.35 × long beam peak).
const Double_t kGmSliceLongCap = 1000.0; // ADC, absolute long-side cap
// Relative cut and band bounds (fractions of long beam peak); upper edge
// below short-side pileup (2× full charge) so pileup never wins the mode.
const Double_t kGmShoulderLongFrac = 0.35;
const Double_t kGmShoulderBandLoFrac = 0.40;
const Double_t kGmShoulderBandHiFrac = 1.25;

void CollectAnchorSamplesOneSubfile(const FileSpec &spec,
                                    const std::vector<ChannelCal> &chans,
                                    const BeamFit2D &beam,
                                    const BeamFit2D &beam0vGrid,
                                    std::vector<std::vector<Float_t>> &samples,
                                    StripPairSamples pairs[18]) {
  Int_t n_chans = Int_t(chans.size());
  samples.assign(n_chans, std::vector<Float_t>());
  if (!beam.ok)
    return;

  TString sub = FileSet::EventsName(spec) + ".root";
  TFile *sf = IO::OpenForReading(sub);
  if (!sf || sf->IsZombie()) {
    if (sf)
      sf->Close();
    return;
  }
  TTree *tree = static_cast<TTree *>(sf->Get("events"));
  if (!tree) {
    sf->Close();
    delete sf;
    return;
  }
  // Guard strips (S) and left ends (L) → Left_0_17_dE; right ends → RightdE.
  UShort_t left_0_17_adc[18], rightdE_adc[18];
  Short_t cathode_adc = 0, grid_adc = 0;
  tree->SetBranchAddress("Left_0_17_dE", left_0_17_adc);
  tree->SetBranchAddress("RightdE", rightdE_adc);
  tree->SetBranchAddress("Cathode", &cathode_adc);
  tree->SetBranchAddress("Grid", &grid_adc);

  Long64_t n = tree->GetEntries();
  for (Long64_t j = 0; j < n; j++) {
    tree->GetEntry(j);
    // UNGATED L/R pairs for check_LR notebook; shoulder slice needs reaction
    // events that a beam gate would remove (charge mostly on short end).
    if (pairs) {
      for (Int_t s = 1; s <= 16; s++) {
        Int_t lv = Int_t(left_0_17_adc[s]);
        Int_t rv = Int_t(rightdE_adc[s]);
        if (lv > 0 && rv > 0 && Long64_t(pairs[s].l.size()) < kPairCap) {
          pairs[s].l.push_back(Float_t(lv));
          pairs[s].r.push_back(Float_t(rv));
        }
        // (short, long) pairs where LONG reads low; collected independently so
        // rare slice events aren't crowded out.
        Bool_t l_is_long = (LongSide(s) == 'L');
        Int_t long_v = l_is_long ? lv : rv;
        Int_t short_v = l_is_long ? rv : lv;
        if (short_v > 0 && Double_t(long_v) < kGmSliceLongCap &&
            Long64_t(pairs[s].slice_short.size()) < kPairCap) {
          pairs[s].slice_short.push_back(Float_t(short_v));
          pairs[s].slice_long.push_back(Float_t(long_v));
        }
      }
    }
    Double_t x = Double_t(left_0_17_adc[1]) + Double_t(rightdE_adc[1]);
    Double_t y = Double_t(left_0_17_adc[2]) + Double_t(rightdE_adc[2]);
    if (x <= 0 || y <= 0)
      continue;
    if (!BeamFitUtils::InEllipseXY(beam, x, y, kEllipseNSigmaX,
                                   kEllipseNSigmaY))
      continue;
    if (beam0vGrid.ok) {
      Double_t g0 = Double_t(grid_adc);
      Double_t s0 = Double_t(left_0_17_adc[0]);
      if (g0 <= 0 || s0 <= 0)
        continue;
      if (!BeamFitUtils::InEllipseXY(beam0vGrid, g0, s0, kEllipseNSigmaX,
                                     kEllipseNSigmaY))
        continue;
    }
    for (Int_t i = 0; i < n_chans; i++) {
      if (Long64_t(samples[i].size()) >= kSampleCap)
        continue;
      const ChannelCal &c = chans[i];
      Int_t v = 0;
      if (c.side == 'S' || c.side == 'L')
        v = Int_t(left_0_17_adc[c.strip]);
      else if (c.side == 'R')
        v = Int_t(rightdE_adc[c.strip]);
      else if (c.side == 'C')
        v = Int_t(cathode_adc);
      if (v > 0)
        samples[i].push_back(Float_t(v));
    }
  }
  sf->Close();
  delete sf;
}

// L/R gain match: PASS 1 — LONG = beam peak, SHORT = full-charge shoulder mode.
// Gains = 1/anchor/target; PASS 2 — eSum fine-tunes short-side gain to median.
const Int_t kGmBins = 512;
const Long64_t kGmMinSlice = 100; // min entries in the shoulder slice
const Double_t kGmEsumLo = 0.8;   // a.u. eSum peak search window (pass 2)
const Double_t kGmEsumHi = 2.5;
const Double_t kGmCorrThresh = 0.05; // a.u. (notebook: 30/1000 = 3% of TARGET)

// Histogram-mode peak finder: bin v over [lo, hi] with kGmBins bins, skip first
// skip_frac bins (avoid threshold pile), return max-bin centre; 0 if empty.
Double_t GmFindPeak(const std::vector<Float_t> &v, Double_t lo, Double_t hi,
                    Double_t skip_frac) {
  if (v.empty() || hi <= lo)
    return 0.0;
  std::vector<Long64_t> h(kGmBins, 0);
  Double_t bw = (hi - lo) / kGmBins;
  for (Int_t j = 0; j < Int_t(v.size()); j++) {
    Double_t x = Double_t(v[j]);
    if (x < lo || x >= hi)
      continue;
    Int_t b = Int_t((x - lo) / bw);
    if (b >= 0 && b < kGmBins)
      h[b]++;
  }
  Int_t skip = Int_t(kGmBins * skip_frac);
  Int_t bmax = -1;
  Long64_t vmax = 0;
  for (Int_t b = skip; b < kGmBins; b++) {
    if (h[b] > vmax) {
      vmax = h[b];
      bmax = b;
    }
  }
  if (bmax < 0 || vmax <= 0)
    return 0.0;
  return lo + (bmax + 0.5) * bw;
}

void ComputeLRGainMatch(std::vector<ChannelCal> &chans,
                        const StripPairSamples pairs[18]) {
  Int_t idx_l[18], idx_r[18];
  for (Int_t s = 0; s < 18; s++) {
    idx_l[s] = -1;
    idx_r[s] = -1;
  }
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    if (chans[i].strip >= 1 && chans[i].strip <= 16) {
      if (chans[i].side == 'L')
        idx_l[chans[i].strip] = i;
      else if (chans[i].side == 'R')
        idx_r[chans[i].strip] = i;
    }
  }

  // Pass 1: LONG anchor = beam peak; SHORT = full-charge shoulder mode (LONG
  // reads low), searched in [0.40, 1.25] × long_peak.
  Bool_t matched[18] = {kFALSE};
  Double_t short_anchor_adc[18] = {0};
  std::vector<Double_t> anchors_odd, anchors_even;
  for (Int_t s = 1; s <= 16; s++) {
    if (idx_l[s] < 0 || idx_r[s] < 0)
      continue;
    Bool_t l_is_long = (LongSide(s) == 'L');
    ChannelCal &c_long = chans[l_is_long ? idx_l[s] : idx_r[s]];
    if (!IsCalibrated(c_long)) {
      std::cerr << "  strip " << s
                << ": long side uncalibrated; skipping L/R gain match"
                << std::endl;
      continue;
    }

    const StripPairSamples &p = pairs[s];
    Double_t long_peak = c_long.fit_adc;
    std::vector<Float_t> rel;
    for (Int_t j = 0; j < Int_t(p.slice_long.size()); j++) {
      if (p.slice_long[j] > 0 &&
          Double_t(p.slice_long[j]) < kGmShoulderLongFrac * long_peak)
        rel.push_back(p.slice_short[j]);
    }

    Double_t peak_short = 0.0;
    if (Long64_t(rel.size()) >= kGmMinSlice) {
      Double_t lo = kGmShoulderBandLoFrac * long_peak;
      Double_t hi = kGmShoulderBandHiFrac * long_peak;
      peak_short = GmFindPeak(rel, lo, hi, 0.0);
      // Reject anchors hugging band edge: edge-hugging mode indicates
      // partial-charge events — treat as unmeasured and fall back to median.
      if (peak_short < 1.05 * lo || peak_short > 0.95 * hi)
        peak_short = 0.0;
    }

    if (peak_short > 0) {
      short_anchor_adc[s] = peak_short;
      if (l_is_long)
        anchors_odd.push_back(peak_short);
      else
        anchors_even.push_back(peak_short);
      std::cout << "  strip " << s
                << " short_anchor=" << Form("%.1f", peak_short)
                << " ADC (shoulder, n=" << rel.size() << ")" << std::endl;
    }
  }

  Double_t median_anchor_odd = 0.0, median_anchor_even = 0.0;
  if (!anchors_odd.empty()) {
    std::sort(anchors_odd.begin(), anchors_odd.end());
    Int_t ms = Int_t(anchors_odd.size());
    median_anchor_odd =
        (ms % 2 == 1) ? anchors_odd[ms / 2]
                      : 0.5 * (anchors_odd[ms / 2 - 1] + anchors_odd[ms / 2]);
  }
  if (!anchors_even.empty()) {
    std::sort(anchors_even.begin(), anchors_even.end());
    Int_t ms = Int_t(anchors_even.size());
    median_anchor_even =
        (ms % 2 == 1) ? anchors_even[ms / 2]
                      : 0.5 * (anchors_even[ms / 2 - 1] + anchors_even[ms / 2]);
  }

  for (Int_t s = 1; s <= 16; s++) {
    if (idx_l[s] < 0 || idx_r[s] < 0)
      continue;
    Bool_t l_is_long = (LongSide(s) == 'L');
    ChannelCal &c_long = chans[l_is_long ? idx_l[s] : idx_r[s]];
    ChannelCal &c_short = chans[l_is_long ? idx_r[s] : idx_l[s]];
    if (!IsCalibrated(c_long))
      continue;
    if (short_anchor_adc[s] <= 0) {
      Double_t fallback = l_is_long ? median_anchor_odd : median_anchor_even;
      // Cross-parity last resort: if this parity measured no anchors at all,
      // the other parity's scale is still the same electronics family.
      if (fallback <= 0)
        fallback = l_is_long ? median_anchor_even : median_anchor_odd;
      if (fallback <= 0) {
        std::cerr << "  strip " << s
                  << ": no anchor and no median fallback; keeping "
                     "independent gains"
                  << std::endl;
        continue;
      }
      short_anchor_adc[s] = fallback;
      std::cout << "  strip " << s << " short_anchor=median fallback "
                << Form("%.1f", fallback) << " ADC" << std::endl;
    }
    c_long.gain = 1.0 / c_long.fit_adc;
    c_short.gain = 1.0 / short_anchor_adc[s];
    matched[s] = kTRUE;
    std::cout << "  strip " << s
              << " L/R match: long peak=" << Form("%.1f", c_long.fit_adc)
              << " ADC  short anchor=" << Form("%.1f", short_anchor_adc[s])
              << " ADC" << std::endl;
  }

  // ── Pass 2: per-strip eSum peak, median reference, short-side correction ──
  Double_t esum_peak[18] = {0};
  std::vector<Double_t> peaks_for_median;
  for (Int_t s = 1; s <= 16; s++) {
    if (!matched[s])
      continue;
    const StripPairSamples &p = pairs[s];
    Double_t gl = Gain(chans[idx_l[s]]);
    Double_t gr = Gain(chans[idx_r[s]]);
    std::vector<Float_t> esum;
    esum.reserve(p.l.size());
    for (Int_t j = 0; j < Int_t(p.l.size()); j++)
      esum.push_back(Float_t(gl * Double_t(p.l[j]) + gr * Double_t(p.r[j])));
    esum_peak[s] = GmFindPeak(esum, kGmEsumLo, kGmEsumHi, 0.0);
    if (esum_peak[s] > 0)
      peaks_for_median.push_back(esum_peak[s]);
  }
  if (peaks_for_median.size() < 4) {
    std::cerr << "  L/R gain match pass 2: too few strips with eSum peaks ("
              << peaks_for_median.size() << "); skipping short-side correction"
              << std::endl;
  } else {
    std::sort(peaks_for_median.begin(), peaks_for_median.end());
    Int_t m = Int_t(peaks_for_median.size());
    Double_t esum_ref =
        (m % 2 == 1)
            ? peaks_for_median[m / 2]
            : 0.5 * (peaks_for_median[m / 2 - 1] + peaks_for_median[m / 2]);
    std::cout << "  L/R gain match: eSum reference (median) = "
              << Form("%.4f", esum_ref) << " a.u." << std::endl;

    for (Int_t s = 1; s <= 16; s++) {
      if (!matched[s] || esum_peak[s] <= 0)
        continue;
      Double_t delta = esum_peak[s] - esum_ref;
      Double_t short_now = esum_peak[s] - 1.0;
      Double_t short_want = esum_ref - 1.0;
      Bool_t l_is_long = (LongSide(s) == 'L');
      ChannelCal &c_short = chans[l_is_long ? idx_r[s] : idx_l[s]];
      if (TMath::Abs(delta) > kGmCorrThresh && short_now > 0 &&
          short_want > 0) {
        Double_t corr = short_want / short_now;
        c_short.gain *= corr;
        std::cout << "  strip " << s
                  << " eSum peak=" << Form("%.4f", esum_peak[s])
                  << " delta=" << Form("%+.4f", delta) << " -> short gain x"
                  << Form("%.4f", corr) << std::endl;
      } else {
        std::cout << "  strip " << s
                  << " eSum peak=" << Form("%.4f", esum_peak[s])
                  << " delta=" << Form("%+.4f", delta) << " OK" << std::endl;
      }
    }
  }
}

// Cathode: median + IQR (asymmetric tail). All other channels: robust
// mode-seeded Gaussian fit over peak core; fall back to robust mode on fit
// failure.
void ReduceToAnchors(std::vector<ChannelCal> &chans,
                     std::vector<std::vector<Float_t>> &samples,
                     std::vector<TF1 *> &fits_out, const TString &run_label,
                     const StripPairSamples pairs[18]) {
  Int_t n_chans = Int_t(chans.size());
  fits_out.assign(n_chans, nullptr);

  for (Int_t i = 0; i < n_chans; i++) {
    ChannelCal &c = chans[i];
    std::vector<Float_t> &v = samples[i];
    c.n_samples = Long64_t(v.size());

    if (Long64_t(v.size()) < kMinSamples) {
      c.fit_adc = 0;
      c.fit_sigma_adc = 0;
      continue;
    }
    if (c.side == 'C') {
      // Cathode: median + IQR (asymmetric tail, no clean peak).
      c.fit_adc = Median(v);
      c.fit_sigma_adc = InterquartileRange(v) / 1.349;
    } else {
      // Mode-fallback: robust mode, never tail-biased sample mean.
      Double_t peak = 0, sig = 0;
      TF1 *fit = nullptr;
      TString fname =
          Form("f_peak_gaus_%s_%s", c.name.Data(), run_label.Data());
      if (FitBeamPeakGaussian(v, fname, peak, sig, fit)) {
        c.fit_adc = peak;
        c.fit_sigma_adc = sig;
        fits_out[i] = fit;
      } else {
        // Fit failed; anchor on robust mode (still "calibrated"),
        // flagged with long/short tag — long-side fallback is a miscalibration
        // risk.
        Double_t mode = 0.0, rsigma = 0.0;
        RobustPeakSeed(v, mode, rsigma);
        if (rsigma > 0.5 * mode) {
          // Sanity gate: no discernible peak (rsigma > 0.5*mode) → gain 0,
          // don't anchor on a smeared noise feature.
          c.fit_adc = 0;
          c.fit_sigma_adc = 0;
          std::cerr << "  [uncalibrated no-peak] " << c.name
                    << ": no peak-like feature (mode=" << Form("%.1f", mode)
                    << ", rsigma=" << Form("%.1f", rsigma)
                    << "); gain 0 (n=" << c.n_samples << ")" << std::endl;
          continue;
        }
        c.fit_adc = mode;
        c.fit_sigma_adc = rsigma;
        TString kind = (c.side == 'S')                 ? "guard"
                       : (c.side == LongSide(c.strip)) ? "long"
                                                       : "short";
        std::cerr << "  [fit-fallback " << kind << "] " << c.name
                  << ": peak fit failed; using mode anchor "
                  << Form("%.1f", c.fit_adc) << " ADC (n=" << c.n_samples << ")"
                  << std::endl;
      }
    }
    std::cout << "  " << c.name << " anchor[ADC]=" << c.fit_adc
              << " sig=" << c.fit_sigma_adc << " (n=" << c.n_samples << ")"
              << std::endl;
  }

  // After all channels are fitted, run LR gain matching: long-beam peak +
  // short-shoulder anchors, then eSum median alignment on the short side.
  if (pairs)
    ComputeLRGainMatch(chans, pairs);
}
