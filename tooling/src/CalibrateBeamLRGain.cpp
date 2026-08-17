#include "CalibrateBeamInternal.hpp"

const Long64_t kSampleCap = 20000;
// (L, R) pair cap for the L/R gain-matching passes. Larger than kSampleCap
// because the short-side "shoulder" anchor is found in a narrow slice of the
// pairs (long side ≈ 300-350 ADC) that only a small fraction of events
// populate.
const Long64_t kPairCap = 100000;

// check_LR gain-match slice constants (37Cl_an_check_LR.ipynb, cell "Save
// 2-pass calibration"). Defined before CollectAnchorSamplesOneSubfile because
// the shoulder slice is selected during sample collection.
const Double_t kGmOffsetLongL = 350.0; // ADC slice centre when L is long
const Double_t kGmOffsetLongR = 300.0; // ADC slice centre when R is long
const Double_t kGmOffsetWin = 80.0;    // ADC slice half-width
const Double_t kGmPeakWin = 150.0;     // ADC near-peak gate half-width

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
        // Shoulder slice: short-side value when the long side reads ≈ the
        // offset. Collected independently of the capped pairs above so rare
        // slice events are never crowded out by beam statistics.
        Bool_t l_is_long = (LongSide(s) == 'L');
        Int_t long_v = l_is_long ? lv : rv;
        Int_t short_v = l_is_long ? rv : lv;
        Double_t offset = l_is_long ? kGmOffsetLongL : kGmOffsetLongR;
        if (short_v > 0 &&
            TMath::Abs(Double_t(long_v) - offset) < kGmOffsetWin &&
            Long64_t(pairs[s].shoulder.size()) < kPairCap)
          pairs[s].shoulder.push_back(Float_t(short_v));
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

// L/R gain matching, replicating the check_LR notebook's two-pass recipe
// (37Cl_an_check_LR.ipynb, "Save 2-pass calibration" cell) exactly, minus the
// eta/position correction:
//
// Pass 1 — per-side anchors:
//   * LONG side anchor  = beam peak of the long side (histogram mode). The
//     per-channel anchor from ReduceToAnchors is exactly this (the beam
//     dominates the spectrum), so it is reused.
//   * SHORT side anchor = "shoulder": the mode of the short side among events
//     where the LONG side reads low (≈350 ADC for L-long strips, ≈300 for
//     R-long, ±80). Those are events where the charge went mostly to the
//     short end, so the short side is measured at the same charge scale as
//     the long side's beam peak. Search window (500, 2000) ADC.
//   * gain = TARGET / anchor per side, with TARGET = 1.0 a.u. (the notebook
//     uses 1000 ADC; only the overall scale differs).
//
// Pass 2 — per-strip eSum alignment via the SHORT side only:
//   * eSum = gain_L·L + gain_R·R over the pairs; per-strip beam peak found
//     in (0.8, 2.5) a.u. (notebook: TARGET=1000, eSum ≈1200 →
//     (1100/1000, 2500/1000) = (1.1, 2.5); widened slightly here).
//   * reference = median of the per-strip peaks.
//   * strips with |peak − ref| > 0.05 get gain_short *= (ref−1)/(peak−1),
//     which moves the strip's eSum peak onto the reference without touching
//     the long side.
//
// Pairs are collected UNGATED (all events with both ends firing) because the
// shoulder slice needs reaction/off-position events that the beam gate
// removes. Strips whose slice is still too small (upstream strips see few
// reaction products, so slice statistics grow with strip number) fall back to
// the MEDIAN shoulder of the strips that measured one — the short-side
// full-charge response is the same electronics scale on every strip, and
// pass 2 then fine-tunes each strip individually using its own eSum peak
// (full beam statistics).
const Double_t kGmShortLo = 500.0; // ADC short-side shoulder search window
const Double_t kGmShortHi = 2000.0;
const Int_t kGmBins = 512;
const Double_t kGmSkipFrac = 0.05;
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

  // ── Pass 1: per-side peak anchors, matching the notebook cell 7 ──
  // 2D correlation ridge on every strip: build a 2D histogram of
  //   (short, long), find the most probable short value given the long
  //   beam peak, then extract the short-side anchor. The short side is
  //   L on even strips and R on odd strips. Strips whose correlation has
  //   too few counts fall back to the shoulder slice, then to the median.
  const Int_t kNCorrBins = 256;
  const Double_t kCorrRange = 2000.0;
  const Int_t kCorrMinCounts = 20;

  Bool_t matched[18] = {kFALSE};
  Double_t short_anchor_adc[18] = {0};
  std::vector<Double_t> anchors_found;
  for (Int_t s = 1; s <= 16; s++) {
    if (idx_l[s] < 0 || idx_r[s] < 0)
      continue;
    Bool_t l_is_long = (LongSide(s) == 'L');
    ChannelCal &c_long = chans[l_is_long ? idx_l[s] : idx_r[s]];
    ChannelCal &c_short = chans[l_is_long ? idx_r[s] : idx_l[s]];
    if (!IsCalibrated(c_long)) {
      std::cerr << "  strip " << s
                << ": long side uncalibrated; skipping L/R gain match"
                << std::endl;
      continue;
    }
    Double_t peak_short = 0.0;

    // 2D correlation ridge for the short-side anchor, used on every strip
    // (even and odd alike): build a 2D histogram of (short, long), scan from
    // high short downward for the ridge column, then peak-find the short side
    // near the long beam peak. Anchoring the short side to the long side's
    // beam peak is more robust against pileup than the shoulder slice, so it
    // is the primary method for both parities. Falls back to the shoulder
    // slice when the correlation has too few counts to be meaningful.
    const StripPairSamples &p = pairs[s];
    std::vector<Float_t> short_v, long_v;
    if (l_is_long) {
      short_v = p.r;
      long_v = p.l;
    } else {
      short_v = p.l;
      long_v = p.r;
    }

    if (Long64_t(short_v.size()) > 100 && Long64_t(long_v.size()) > 100) {
      TH2F *h2d =
          new TH2F(Form("h2d_lr_corr_s%d", s), ";short [ADC];long [ADC]",
                   kNCorrBins, 0.0, kCorrRange, kNCorrBins, 0.0, kCorrRange);
      h2d->SetDirectory(nullptr);
      for (Int_t j = 0; j < Int_t(long_v.size()); j++)
        h2d->Fill(Double_t(short_v[j]), Double_t(long_v[j]));

      // Notebook: scan from high short downward, find first short bin whose
      // column (projection onto long) has any bin >= 20 counts.
      Double_t rough_short = 0.0;
      for (Int_t ix = kNCorrBins; ix >= 1; ix--) {
        Double_t col_max = 0.0;
        for (Int_t iy = 1; iy <= kNCorrBins; iy++) {
          Double_t v = h2d->GetBinContent(ix, iy);
          if (v > col_max)
            col_max = v;
        }
        if (col_max >= kCorrMinCounts) {
          rough_short = h2d->GetXaxis()->GetBinCenter(ix);
          break;
        }
      }

      // Find offset_short: short peak in events where long ≈ long beam peak
      Double_t long_peak = c_long.fit_adc;
      std::vector<Float_t> short_near_long;
      for (Int_t j = 0; j < Int_t(long_v.size()); j++) {
        if (TMath::Abs(Double_t(long_v[j]) - long_peak) < kGmPeakWin)
          short_near_long.push_back(short_v[j]);
      }
      Double_t offset_short = 0.0;
      if (short_near_long.size() > 50)
        offset_short = GmFindPeak(short_near_long, 50.0, 600.0, 0.05);

      if (rough_short > 0) {
        peak_short =
            GmFindPeak(short_v, rough_short * 0.85, rough_short * 1.15, 0.0);
      } else {
        Double_t lo = (offset_short > 0) ? offset_short + 100.0 : 500.0;
        peak_short = GmFindPeak(short_v, lo, 2000.0, 0.05);
      }
      delete h2d;
    }

    if (peak_short <= 0) {
      // Fall back to shoulder if 2D correlation fails
      const std::vector<Float_t> &slice = pairs[s].shoulder;
      if (Long64_t(slice.size()) >= kGmMinSlice)
        peak_short = GmFindPeak(slice, kGmShortLo, kGmShortHi, kGmSkipFrac);
    }

    if (peak_short <= 0)
      continue;
    short_anchor_adc[s] = peak_short;
    anchors_found.push_back(peak_short);
    std::cout << "  strip " << s << " short_anchor=" << Form("%.1f", peak_short)
              << " ADC (method=" << (l_is_long ? "2Dcorr" : "shoulder") << ")"
              << std::endl;
  }

  Double_t median_anchor = 0.0;
  if (!anchors_found.empty()) {
    std::sort(anchors_found.begin(), anchors_found.end());
    Int_t ms = Int_t(anchors_found.size());
    median_anchor =
        (ms % 2 == 1)
            ? anchors_found[ms / 2]
            : 0.5 * (anchors_found[ms / 2 - 1] + anchors_found[ms / 2]);
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
      if (median_anchor <= 0) {
        std::cerr << "  strip " << s
                  << ": no anchor and no median fallback; keeping "
                     "independent gains"
                  << std::endl;
        continue;
      }
      short_anchor_adc[s] = median_anchor;
      std::cout << "  strip " << s << " short_anchor=median fallback "
                << Form("%.1f", median_anchor) << " ADC" << std::endl;
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
