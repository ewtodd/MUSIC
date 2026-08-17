#include "CalibrateBeamInternal.hpp"

const Long64_t kSampleCap = 20000;
// (L, R) pair cap for the L/R gain-matching passes. Larger than kSampleCap
// because the short-side "shoulder" anchor is found in a narrow slice of the
// pairs (long side reads low) that only a small fraction of events populate.
const Long64_t kPairCap = 100000;

// Shoulder slice constants for the short-side full-charge anchor. The slice
// collects (short, long) pairs where the LONG side reads below an absolute
// cap (a superset of the relative cut applied later: long < 0.35 x long beam
// peak, which never exceeds ~700 ADC for these beam peaks). Defined before
// CollectAnchorSamplesOneSubfile because the slice is selected during sample
// collection.
const Double_t kGmSliceLongCap = 1000.0; // ADC, absolute long-side cap
// Relative cut and search band (fractions of the long side's beam peak) used
// in ComputeLRGainMatch. The band's upper edge sits below the short side's
// pileup reading (2x full charge) so pileup never wins the mode.
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
  // Raw ADC, pre-calibration. Guard strips (S) and left ends (L) live in
  // Left_0_17_dE; right ends in RightdE. Strip totals are L+R; the gate uses
  // the strip1 total (L1+R1) and the strip2 total (L2+R2).
  UShort_t left_0_17_adc[18], rightdE_adc[18];
  Short_t cathode_adc = 0, grid_adc = 0;
  tree->SetBranchAddress("Left_0_17_dE", left_0_17_adc);
  tree->SetBranchAddress("RightdE", rightdE_adc);
  tree->SetBranchAddress("Cathode", &cathode_adc);
  tree->SetBranchAddress("Grid", &grid_adc);

  Long64_t n = tree->GetEntries();
  for (Long64_t j = 0; j < n; j++) {
    tree->GetEntry(j);
    // (L, R) pairs for strips 1–16 where both ends fired — UNGATED, matching
    // the check_LR notebook which runs on all events. The short-side anchor
    // ("shoulder") is found in a slice where the LONG side reads low
    // (~300-350 ADC), i.e. events where the charge went mostly to the short
    // end. Those are reaction/off-position events that a beam gate would
    // remove, so the pairs must not be beam-gated.
    if (pairs) {
      for (Int_t s = 1; s <= 16; s++) {
        Int_t lv = Int_t(left_0_17_adc[s]);
        Int_t rv = Int_t(rightdE_adc[s]);
        if (lv > 0 && rv > 0 && Long64_t(pairs[s].l.size()) < kPairCap) {
          pairs[s].l.push_back(Float_t(lv));
          pairs[s].r.push_back(Float_t(rv));
        }
        // Shoulder slice: (short, long) pairs where the LONG side reads low,
        // i.e. events where the charge went mostly to the short end. Collected
        // independently of the capped pairs above so rare slice events are
        // never crowded out by beam statistics. The relative cut (long < 0.35
        // x beam peak) is applied in ComputeLRGainMatch once the beam peak is
        // known; kGmSliceLongCap is an absolute superset of that cut.
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

// L/R gain matching, following the check_LR notebook's two-pass recipe
// (37Cl_an_check_LR.ipynb, "Save 2-pass calibration" cell) in spirit, but
// with the short-side anchor fixed to the FULL-CHARGE response:
//
// Pass 1 — per-side anchors:
//   * LONG side anchor  = beam peak of the long side (histogram mode). The
//     per-channel anchor from ReduceToAnchors is exactly this (the beam
//     dominates the spectrum), so it is reused: gain_long = 1/long_peak.
//     Anchoring the long side at its beam peak keeps long-only events at
//     1.0 a.u., matching every downstream threshold that assumes a flat beam.
//   * SHORT side anchor = the short side's full-charge response: the mode of
//     the short side among events where the LONG side reads low
//     (< 0.35 x long beam peak). Those events carry most of the charge on
//     the short end, so the short side reads the same charge the long side
//     reads at its own full-charge scale. Anchoring there makes the strip
//     total gL*L + gR*R partition-independent: an event that routes a
//     different fraction of its charge to the short end reads the same
//     total, instead of being scaled by the anchor error. (The notebook's
//     2D-correlation ridge was tried first; on this data it locks onto the
//     beam blob's top edge / pileup instead of the full-charge point, so the
//     shoulder slice is the primary method for both parities.)
//   * gain = 1/anchor per side (TARGET = 1.0 a.u.).
//
// Pass 2 — per-strip eSum alignment via the SHORT side only:
//   * eSum = gain_L·L + gain_R·R over the pairs; per-strip beam peak found
//     in (0.8, 2.5) a.u.
//   * reference = median of the per-strip peaks.
//   * strips with |peak − ref| > 0.05 get gain_short *= (ref−1)/(peak−1),
//     which moves the strip's eSum peak onto the reference without touching
//     the long side. With full-charge short anchors this correction is small
//     (the anchors already share one electronics scale); it only fine-tunes
//     residual per-strip differences.
//
// Pairs are collected UNGATED (all events with both ends firing) because the
// shoulder slice needs reaction/off-position events that the beam gate
// removes. Strips whose slice is still too small (upstream strips see few
// reaction products, so slice statistics grow with strip number) fall back to
// the MEDIAN anchor of the same parity (short side = R on odd strips, L on
// even strips — each parity's preamps share one electronics scale), and pass
// 2 then fine-tunes each strip individually using its own eSum peak (full
// beam statistics).
const Int_t kGmBins = 512;
const Long64_t kGmMinSlice = 100; // min entries in the shoulder slice
const Double_t kGmEsumLo = 0.8;   // a.u. eSum peak search window (pass 2)
const Double_t kGmEsumHi = 2.5;
const Double_t kGmCorrThresh = 0.05; // a.u. (notebook: 30/1000 = 3% of TARGET)

// Histogram-mode peak finder, mirroring the notebook's find_peak(): histogram
// `v` over [lo, hi] with kGmBins bins, skip the first skip_frac of bins (to
// avoid the threshold pile), return the max-bin centre. Returns 0 when empty.
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

  // ── Pass 1: per-side peak anchors ──
  // LONG side anchor = the beam peak of the long side (from ReduceToAnchors).
  // SHORT side anchor = the short side's FULL-CHARGE response from the
  // shoulder slice: the mode of the short side among events where the LONG
  // side reads low (< 0.35 x long beam peak). Those events carry most of the
  // charge on the short end, so the short side reads the same charge the
  // long side reads at its full-charge scale. The mode is searched in
  // [0.40, 1.25] x long_peak: the lower edge sits above the short side's
  // beam-mode pile, the upper edge below its pileup reading (2x full charge).
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
      // Reject anchors hugging the band edge: a trustworthy full-charge mode
      // sits comfortably inside the band, so an edge-hugging mode means the
      // slice was dominated by partial-charge events (or the band is wrong
      // for this strip) — treat it as unmeasured and fall back to the median.
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

// Cathode uses median + IQR/1.349 (asymmetric tail not as clean and the user
// prefers to keep cathode on the existing approach). All other channels
// (S guard strips + L/R long anodes) use a robust mode-seeded Gaussian fit
// over the peak core (mode ± 2σ), anchored to the fitted centroid. Fall back
// to the robust mode itself on fit failure.
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
      // Primary: robust mode-seeded Gaussian fit of the peak core. Last
      // resort: the robust mode itself (peak-like), never the tail-biased
      // sample mean.
      Double_t peak = 0, sig = 0;
      TF1 *fit = nullptr;
      TString fname =
          Form("f_peak_gaus_%s_%s", c.name.Data(), run_label.Data());
      if (FitBeamPeakGaussian(v, fname, peak, sig, fit)) {
        c.fit_adc = peak;
        c.fit_sigma_adc = sig;
        fits_out[i] = fit;
      } else {
        // Fit failed; anchor on the robust mode. Still "calibrated", but
        // no fit curve is drawn -- flag it, tagged long/short, since a
        // long-side fallback is a real miscalibration risk.
        Double_t mode = 0.0, rsigma = 0.0;
        RobustPeakSeed(v, mode, rsigma);
        if (rsigma > 0.5 * mode) {
          // Same peak-like sanity gate as the fitter: no discernible peak,
          // leave the channel uncalibrated (gain 0 -> reads 0 a.u.) instead
          // of anchoring on a smear's mode.
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

  // After all per-channel peaks are fitted, run the check_LR-style two-pass
  // L/R gain matching for strips 1–16: long-side beam peak + short-side
  // shoulder anchors, then a per-strip eSum median alignment applied to the
  // short side only. Sets the ChannelCal::gain overrides; strips where the
  // shoulder cannot be found keep the independent 1/fit_adc gains.
  if (pairs)
    ComputeLRGainMatch(chans, pairs);
}
