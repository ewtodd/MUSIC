#include "CalibrateBeamInternal.hpp"

const Long64_t kSampleCap = 20000;

void CollectAnchorSamplesOneSubfile(
    const FileSpec &spec, const std::vector<ChannelCal> &chans,
    const BeamFit2D &beam, const BeamFit2D &beam0vGrid,
    std::vector<std::vector<Float_t>> &samples) {
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

// Fill `out[s]` (s = 1..16) with the beam peak of one anode family on every
// strip. `base` holds the direct measurements on the strips where that family
// is the long side (odd strips for L, even strips for R); the in-between
// strips get the linear interpolation of the two nearest measured neighbours,
// with linear extrapolation at the detector ends. Zero (unmeasured) base
// entries are skipped by the neighbour search.
void InterpolateParityPeak(const Double_t base[18], Double_t out[18],
                           Bool_t measured_on_odd) {
  Int_t parity = measured_on_odd ? 1 : 0;
  for (Int_t s = 1; s <= 16; s++)
    out[s] = base[s];

  for (Int_t s = 1; s <= 16; s++) {
    if (s % 2 == parity)
      continue;
    Int_t lo = -1, hi = -1;
    for (Int_t t = s - 1; t >= 1; t--)
      if (t % 2 == parity && base[t] > 0) {
        lo = t;
        break;
      }
    for (Int_t t = s + 1; t <= 16; t++)
      if (t % 2 == parity && base[t] > 0) {
        hi = t;
        break;
      }

    if (lo >= 0 && hi >= 0) {
      Double_t frac = Double_t(s - lo) / Double_t(hi - lo);
      out[s] = base[lo] + frac * (base[hi] - base[lo]);
    } else if (lo >= 0) {
      // Past the last measured strip: linear extrapolation from the last two.
      Int_t lo2 = -1;
      for (Int_t t = lo - 2; t >= 1; t--)
        if (t % 2 == parity && base[t] > 0) {
          lo2 = t;
          break;
        }
      out[s] = (lo2 >= 0)
                   ? base[lo] + (base[lo] - base[lo2]) *
                                    (Double_t(s - lo) / Double_t(lo - lo2))
                   : base[lo];
    } else if (hi >= 0) {
      // Before the first measured strip: linear extrapolation from the first
      // two.
      Int_t hi2 = -1;
      for (Int_t t = hi + 2; t <= 16; t++)
        if (t % 2 == parity && base[t] > 0) {
          hi2 = t;
          break;
        }
      out[s] = (hi2 >= 0)
                   ? base[hi] - (base[hi2] - base[hi]) *
                                    (Double_t(hi - s) / Double_t(hi2 - hi))
                   : base[hi];
    }
  }
}

// L/R gain matching: put the two ends of every strip on one charge scale by
// anchoring each end at its own beam peak. The long side's beam peak is
// measured directly (it carries most of the beam deposit). The short side's
// beam peak comes from the SAME physical end measured where that end is the
// long side (the L preamp on odd strips, the R preamp on even strips) and is
// interpolated onto the short strips; the short-side channels sit below
// threshold at the beam, so their own beam sample is not usable. Anchoring
// both ends at their beam peaks makes g_L/g_R = R_peak/L_peak, so
// total = g_L*L + g_R*R is partition-independent.
void ComputeLRGainMatch(std::vector<ChannelCal> &chans) {
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

  // Long-side beam peaks, one per anode family.
  Double_t L_peak[18] = {0};
  Double_t R_peak[18] = {0};
  for (Int_t s = 1; s <= 16; s++) {
    if (s % 2 == 1 && idx_l[s] >= 0 && IsCalibrated(chans[idx_l[s]]))
      L_peak[s] = chans[idx_l[s]].fit_adc;
    if (s % 2 == 0 && idx_r[s] >= 0 && IsCalibrated(chans[idx_r[s]]))
      R_peak[s] = chans[idx_r[s]].fit_adc;
  }

  Double_t L_peak_all[18] = {0};
  Double_t R_peak_all[18] = {0};
  InterpolateParityPeak(L_peak, L_peak_all, kTRUE);
  InterpolateParityPeak(R_peak, R_peak_all, kFALSE);

  for (Int_t s = 1; s <= 16; s++) {
    if (idx_l[s] < 0 || idx_r[s] < 0)
      continue;
    Double_t l_peak = L_peak_all[s];
    Double_t r_peak = R_peak_all[s];
    if (l_peak <= 0 || r_peak <= 0) {
      std::cerr << "  strip " << s
                << ": missing beam-peak anchor; keeping independent gains"
                << std::endl;
      continue;
    }
    chans[idx_l[s]].gain = 1.0 / l_peak;
    chans[idx_r[s]].gain = 1.0 / r_peak;
    std::cout << "  strip " << s
              << " L/R match: L peak=" << Form("%.1f", l_peak)
              << " ADC  R peak=" << Form("%.1f", r_peak)
              << " ADC  g_L/g_R=" << Form("%.4f", r_peak / l_peak) << std::endl;
  }
}

// Cathode: median + IQR (asymmetric tail). All other channels: robust
// mode-seeded Gaussian fit over peak core; fall back to robust mode on fit
// failure.
void ReduceToAnchors(std::vector<ChannelCal> &chans,
                     std::vector<std::vector<Float_t>> &samples,
                     std::vector<TF1 *> &fits_out, const TString &run_label) {
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

  // After all per-channel peaks are fitted, match the L/R gains: each end is
  // anchored at its own beam peak (short side interpolated from the opposite
  // parity where that end is the long side).
  ComputeLRGainMatch(chans);
}
