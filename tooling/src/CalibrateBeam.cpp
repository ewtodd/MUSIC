#include "CalibrateBeam.hpp"
#include <Rtypes.h>
#include <TH2D.h>
#include <utility>

const Int_t kMaxChannels = 36;
// Raw total of one strip. Strips 0 and 17 are unsegmented and read once;
// strips 1-16 are read at two ends, held in arrays of 16 indexed by strip - 1,
// and their total is the sum of both.
inline Double_t StripTotalAdc(const UShort_t *l, const UShort_t *r, UShort_t s0,
                              UShort_t s17, Int_t s) {
  if (s <= 0)
    return Double_t(s0);
  if (s >= 17)
    return Double_t(s17);
  return Double_t(l[s - 1]) + Double_t(r[s - 1]);
}

/// The strip whose total is paired with strip s to gate it: the one before it,
/// except strips 0 and 1 which share the (0, 1) pair because there is nothing
/// before strip 0.
inline Int_t GatePartner(Int_t s) { return s <= 1 ? 0 : s - 1; }
const Long64_t kMinSamples = 200;
const Long64_t kSampleCap = 20000;
/// (L, R) pair cap for the L/R gain-matching passes. Larger than kSampleCap
/// because the short-side "shoulder" anchor is found in a narrow slice of the
/// pairs (long side ≈ 300-350 ADC) that only a small fraction of events
/// populate.
const Long64_t kPairCap = 100000;

// Defined below; forward-declared for the anchor reduction.
void RobustPeakSeed(const std::vector<Float_t> &v, Double_t &mode,
                    Double_t &sigma);

std::vector<ChannelCal> CalibrateBeam::BuildChannels() {
  std::vector<ChannelCal> chans;
  ChannelCal c{};
  c.name = "Strip0";
  c.side = 'S';
  c.strip = 0;
  chans.push_back(c);
  c.name = "Strip17";
  c.side = 'S';
  c.strip = 17;
  chans.push_back(c);
  for (Int_t s = 1; s <= 16; s++) {
    c.name = Form("L%d", s);
    c.side = 'L';
    c.strip = s;
    chans.push_back(c);
    c.name = Form("R%d", s);
    c.side = 'R';
    c.strip = s;
    chans.push_back(c);
  }
  c.name = "Cathode";
  c.side = 'C';
  c.strip = -1;
  chans.push_back(c);
  // The grid: anchored on its modal beam peak like a strip, so the analysis
  // can gate on it in beam units (GATE_AXIS_GRID).
  c.name = "Grid";
  c.side = 'G';
  c.strip = -1;
  chans.push_back(c);
  return chans;
}

Char_t LongSide(Int_t strip) { return (strip % 2 == 0) ? 'R' : 'L'; }

Bool_t IsCalibrated(const ChannelCal &c) { return c.fit_adc > 0; }

Double_t Gain(const ChannelCal &c) {
  return c.gain >= 0.0 ? c.gain : 1.0 / c.fit_adc;
}

/// The beam gate for one strip: a 2D Gaussian on the (strip sx, strip sy) raw
/// totals. Gating a strip on itself and its neighbour selects single-particle
/// events in THAT strip; the gate is on totals, so it says nothing about how
/// the charge divides between the two ends, which is what the anchor fit reads.
BeamFit2D FindBeamGateStrips(const FileSpec &spec, Int_t sx, Int_t sy,
                             const TString &run_label,
                             const TString &plot_subdir,
                             Bool_t save_plot = kTRUE) {
  BeamFit2D out;

  TString sub = FileSet::EventsName(spec) + ".root";
  TFile *sf = IO::OpenForReading(sub);
  if (!sf || sf->IsZombie()) {
    if (sf)
      delete sf;
    return out;
  }
  TTree *tree = static_cast<TTree *>(sf->Get("events"));
  if (!tree) {
    sf->Close();
    delete sf;
    return out;
  }
  // Raw ADC, read before any calibration exists (this fit produces it), so read
  // the branches directly. strip1 total = L1 + R1; strip2 total = L2 + R2.
  UShort_t leftdE_adc[16], rightdE_adc[16], strip0_adc = 0, strip17_adc = 0;
  tree->SetBranchAddress("LeftdE", leftdE_adc);
  tree->SetBranchAddress("RightdE", rightdE_adc);
  tree->SetBranchAddress("Strip0dE", &strip0_adc);
  tree->SetBranchAddress("Strip17dE", &strip17_adc);

  /// Use 1024 bins (previously 256) so each bin is narrower (16 ADC for 37Cl).
  /// The sigma floor in ComputeMoments is 2 * bin_width; at the old 64 ADC/bin
  /// it clamped sigma to 128 ADC, hiding the true beam width and washing out
  /// correlation (rho ≈ 0.09 even though strip1/strip2 track the same beam).
  const Int_t kBeamGateNBins = 1024;
  TH2F *h =
      new TH2F(Form("h2_gate_s%d_s%d_%s", sx, sy, run_label.Data()),
               Form(";Strip%d #DeltaE [ADC];Strip%d #DeltaE [ADC]", sx, sy),
               kBeamGateNBins, 0.0, Constants::ActiveStripEMaxAdc(),
               kBeamGateNBins, 0.0, Constants::ActiveStripEMaxAdc());
  h->SetDirectory(nullptr);
  Long64_t n = tree->GetEntries();
  // The same events unbinned (up to a cap), for the clipped moments below.
  const Long64_t kMaxPoints = 2000000;
  std::vector<std::pair<Float_t, Float_t>> pts;
  pts.reserve(std::min(n, kMaxPoints));
  for (Long64_t j = 0; j < n; j++) {
    tree->GetEntry(j);
    const Double_t tx =
        StripTotalAdc(leftdE_adc, rightdE_adc, strip0_adc, strip17_adc, sx);
    const Double_t ty =
        StripTotalAdc(leftdE_adc, rightdE_adc, strip0_adc, strip17_adc, sy);
    if (tx > 0.0 && ty > 0.0) {
      h->Fill(tx, ty);
      if (Long64_t(pts.size()) < kMaxPoints)
        pts.push_back(std::make_pair(Float_t(tx), Float_t(ty)));
    }
  }
  sf->Close();
  delete sf;

  if (h->GetEntries() < 100) {
    Constants::DetailErr() << "  " << run_label
                           << ": too few events for the strip " << sy
                           << " beam gate (strips " << sx << " vs " << sy << ")"
                           << std::endl;
    delete h;
    return out;
  }

  /// With 1024 bins instead of 256, scale from 10→40 to keep the ADC seed
  /// window roughly 640 ADC (10 bins × 16384/256 = 640; 40 bins × 16384/1024 =
  /// 640).
  const Int_t kSeedHalfBins = 40;
  Moments2D m = BeamFitUtils::SeedSpotMoments(h, kSeedHalfBins);
  if (m.weight <= 0) {
    Constants::DetailErr() << "  " << run_label
                           << ": no bins above beam seed threshold"
                           << std::endl;
    delete h;
    return out;
  }
  // Re-center iteratively (moments inside ±2.5σ of the centroid); at high
  // rate (run84: 45 kHz) the ~2x-beam pileup blob passes the seed cut.
  const Int_t kMomentRefineIters = 4;
  const Double_t kMomentRefineNSigma = 2.5;
  const Double_t bw_x = h->GetXaxis()->GetBinWidth(1);
  const Double_t bw_y = h->GetYaxis()->GetBinWidth(1);
  Int_t bx, by, bz;
  h->GetMaximumBin(bx, by, bz);
  const Double_t peak_val = h->GetBinContent(bx, by);
  for (Int_t iter = 0; iter < kMomentRefineIters; iter++) {
    Int_t wlo_bx = std::max(
        1, h->GetXaxis()->FindBin(m.mu_x - kMomentRefineNSigma * m.sigma_x));
    Int_t whi_bx = std::min(
        h->GetNbinsX(),
        h->GetXaxis()->FindBin(m.mu_x + kMomentRefineNSigma * m.sigma_x));
    Int_t wlo_by = std::max(
        1, h->GetYaxis()->FindBin(m.mu_y - kMomentRefineNSigma * m.sigma_y));
    Int_t whi_by = std::min(
        h->GetNbinsY(),
        h->GetYaxis()->FindBin(m.mu_y + kMomentRefineNSigma * m.sigma_y));
    Moments2D m_ref = BeamFitUtils::ComputeMoments(
        h, wlo_bx, whi_bx, wlo_by, whi_by, 0.30 * peak_val, bw_x, bw_y);
    if (m_ref.weight <= 0)
      break;
    m = m_ref;
  }
  // The width BEAM_GATE_NSIGMA is in is a fitted sigma: the spot events,
  // finely binned around the seed, fitted with a correlated Gaussian on a
  // pedestal. Fine binning matters — over the full ADC range the spot is a
  // few bins wide and the fit would be fitting the binning.
  BeamFitUtils::SpotFit spot = BeamFitUtils::FitSpotFromPoints(h, pts, m);
  out = spot.fit;
  Constants::Detail() << "  beam gate strip " << sy << " (strips " << sx
                      << " vs " << sy << "): "
                      << (spot.fit_used ? "fit"
                                        : "CLIPPED MOMENTS (fit failed)")
                      << " mu=(" << out.mu_x << "," << out.mu_y << ") sigma=("
                      << out.sigma_x << "," << out.sigma_y
                      << ") rho=" << out.rho << " chi2/ndf=" << spot.chi2_ndf
                      << std::endl;

  if (save_plot) {
    TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(h, cv);
    BeamFitUtils::DrawEllipse(out, Constants::cfg.BEAM_GATE_NSIGMA,
                              kViolet + 2);
    if (Constants::SavePlots())
      PlottingUtils::SaveFigure(cv, Form("beam_gate_s%02d", sy),
                                plot_subdir + "/beam_gate",
                                PlotSaveOptions::kLINEAR);
    delete cv;
  }
  delete h;
  return out;
}

// Median by partial sort; v need not be sorted, but is only partially
// ordered afterwards, so a sort is still needed for other order statistics.
template <class T> inline Double_t Median(std::vector<T> &v) {
  if (v.empty())
    return 0.0;
  Int_t n = Int_t(v.size());
  std::nth_element(v.begin(), v.begin() + n / 2, v.end());
  Double_t med = Double_t(v[n / 2]);
  if (n % 2 == 0) {
    std::nth_element(v.begin(), v.begin() + n / 2 - 1, v.end());
    med = 0.5 * (med + Double_t(v[n / 2 - 1]));
  }
  return med;
}

/// Robust peak/width estimate used to seed the beam-peak Gaussian fit, and as
/// the fallback anchor when the fit fails. Histograms only the
/// 5th-95th-percentile core so outlier ADC values can't stretch the binning,
/// takes the modal bin centre as the peak and IQR/1.349 as the width.
/// Peak-like: unlike the raw sample mean it is not pulled up by the
/// straggling beam-dE tail.
void RobustPeakSeed(const std::vector<Float_t> &v, Double_t &mode,
                    Double_t &sigma) {
  mode = 0.0;
  sigma = 0.0;
  Int_t n = Int_t(v.size());
  if (n < 4)
    return;
  std::vector<Float_t> s(v);
  std::sort(s.begin(), s.end());
  Int_t i_lo = Int_t(0.05 * n);
  Int_t i_hi = Int_t(0.95 * n);
  Int_t i_q1 = Int_t(0.25 * n);
  Int_t i_q3 = Int_t(0.75 * n);
  if (i_hi >= n)
    i_hi = n - 1;
  Double_t p_lo = Double_t(s[i_lo]);
  Double_t p_hi = Double_t(s[i_hi]);
  Double_t med = Double_t(s[n / 2]);
  sigma = (Double_t(s[i_q3]) - Double_t(s[i_q1])) / 1.349;
  if (!(sigma > 0.0))
    sigma = 0.05 * (med > 0.0 ? med : 1.0);
  if (!(p_hi > p_lo)) {
    mode = med;
    return;
  }
  const Int_t nbins = 100;
  TH1F h("robust_peak_seed_h", "", nbins, p_lo, p_hi);
  h.SetDirectory(nullptr);
  for (Int_t j = 0; j < n; j++)
    if (Double_t(s[j]) >= p_lo && Double_t(s[j]) <= p_hi)
      h.Fill(Double_t(s[j]));
  mode = h.GetBinCenter(h.GetMaximumBin());
  if (!(mode > 0.0))
    mode = med;
}

/// Paired (L, R) raw-ADC samples for one strip, collected UNGATED from events
/// where both ends fire, plus the uncapped "shoulder" slice: short-side values
/// from events where the LONG side reads low (charge went mostly to the short
/// end). The slice is collected separately because those events are rare — the
/// capped pair vectors fill up with beam events long before enough slice
/// events arrive.
struct StripPairSamples {
  std::vector<Float_t> l;
  std::vector<Float_t> r;
  // Beam-gated, long-triggered pairs: keep events with silent short end
  // (recorded 0) -- the population the decode sums; `l`/`r` need both ends.
  std::vector<Float_t> gated_long;
  std::vector<Float_t> gated_short;
};

// check_LR gain-match slice constants (37Cl_an_check_LR.ipynb); defined
// before CollectAnchorSamplesOneSubfile, which picks the shoulder slice.

/// Every strip carries its own beam gate, on its own total and its neighbour's
/// (GatePartner). A strip is therefore anchored on events that were beam-like
/// THERE, rather than on events that were beam-like at strips 1 and 2 and
/// whatever they happened to be doing further down the chamber.
void CollectAnchorSamplesOneSubfile(const FileSpec &spec,
                                    const std::vector<ChannelCal> &chans,
                                    const BeamFit2D gate[18],
                                    std::vector<std::vector<Float_t>> &samples,
                                    StripPairSamples pairs[18]) {
  Int_t n_chans = Int_t(chans.size());
  samples.assign(n_chans, std::vector<Float_t>());
  Bool_t any = kFALSE;
  for (Int_t s = 0; s <= 17; s++)
    any = any || gate[s].ok;
  if (!any)
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
  // Raw ADC, pre-calibration: the two ends of strips 1-16 in LeftdE/RightdE,
  // the unsegmented strips in Strip0dE/Strip17dE. Strip totals are L+R; the
  // gate uses the strip1/strip2 totals.
  UShort_t leftdE_adc[16], rightdE_adc[16], strip0_adc = 0, strip17_adc = 0;
  Short_t cathode_adc = 0, grid_adc = 0;
  tree->SetBranchAddress("LeftdE", leftdE_adc);
  tree->SetBranchAddress("RightdE", rightdE_adc);
  tree->SetBranchAddress("Strip0dE", &strip0_adc);
  tree->SetBranchAddress("Strip17dE", &strip17_adc);
  tree->SetBranchAddress("Cathode", &cathode_adc);
  tree->SetBranchAddress("Grid", &grid_adc);

  Long64_t n = tree->GetEntries();
  for (Long64_t j = 0; j < n; j++) {
    tree->GetEntry(j);
    // (L, R) pairs for strips 1–16, UNGATED (check_LR runs on all events):
    // a beam gate would remove the shoulder events (long side ~300-350 ADC).
    if (pairs) {
      for (Int_t s = 1; s <= 16; s++) {
        Int_t lv = Int_t(leftdE_adc[s - 1]);
        Int_t rv = Int_t(rightdE_adc[s - 1]);
        if (lv > 0 && rv > 0 && Long64_t(pairs[s].l.size()) < kPairCap) {
          pairs[s].l.push_back(Float_t(lv));
          pairs[s].r.push_back(Float_t(rv));
        }
      }
    }
    // Which strips this event is beam-like in, each judged by its own gate.
    Bool_t pass[18] = {kFALSE};
    for (Int_t s = 1; s <= 17; s++) {
      if (!gate[s].ok)
        continue;
      const Double_t x = StripTotalAdc(leftdE_adc, rightdE_adc, strip0_adc,
                                       strip17_adc, GatePartner(s));
      const Double_t y =
          StripTotalAdc(leftdE_adc, rightdE_adc, strip0_adc, strip17_adc, s);
      if (x <= 0.0 || y <= 0.0)
        continue;
      pass[s] = BeamFitUtils::InEllipseXY(gate[s], x, y,
                                          Constants::cfg.BEAM_GATE_NSIGMA);
    }
    // Strip 0 shares strip 1's gate, so share its verdict: GatePartner(0) is
    // 0, and a self-test lands every event outside, zeroing its gain.
    pass[0] = pass[1];
    if (pairs) {
      for (Int_t s = 1; s <= 16; s++) {
        if (!pass[s] || Long64_t(pairs[s].gated_long.size()) >= kPairCap)
          continue;
        Bool_t l_is_long = (LongSide(s) == 'L');
        Int_t long_v =
            l_is_long ? Int_t(leftdE_adc[s - 1]) : Int_t(rightdE_adc[s - 1]);
        Int_t short_v =
            l_is_long ? Int_t(rightdE_adc[s - 1]) : Int_t(leftdE_adc[s - 1]);
        if (long_v > 0) {
          pairs[s].gated_long.push_back(Float_t(long_v));
          pairs[s].gated_short.push_back(Float_t(short_v > 0 ? short_v : 0));
        }
      }
    }
    for (Int_t i = 0; i < n_chans; i++) {
      if (Long64_t(samples[i].size()) >= kSampleCap)
        continue;
      const ChannelCal &c = chans[i];
      // The cathode and the grid have no strip of their own, so they ride on
      // strip 1's gate.
      const Int_t s = (c.side == 'C' || c.side == 'G') ? 1 : c.strip;
      if (s < 0 || s > 17 || !pass[s])
        continue;
      Int_t v = 0;
      if (c.side == 'S')
        v = Int_t(c.strip == 0 ? strip0_adc : strip17_adc);
      else if (c.side == 'L')
        v = Int_t(leftdE_adc[c.strip - 1]);
      else if (c.side == 'R')
        v = Int_t(rightdE_adc[c.strip - 1]);
      else if (c.side == 'C')
        v = Int_t(cathode_adc);
      else if (c.side == 'G')
        v = Int_t(grid_adc);
      if (v > 0)
        samples[i].push_back(Float_t(v));
    }
  }
  sf->Close();
  delete sf;
}

/// L/R gain matching, after the check_LR notebook's recipe
/// (37Cl_an_check_LR.ipynb, "Save 2-pass calibration" cell) minus the
/// eta/position correction, with a short-end offset read off the ridge line:
///
/// Pass 1 — per-side anchors and the short-end offset:
///   * a fixed deposit divided between the two ends traces the line
///         long/C_long + (short - S0)/C_short = 1
///     in the (short, long) plane. S0 is an additive offset the short end's
///     energy carries whenever that end fired: on 37Cl the R-chain short ends
///     read 100-230 ADC more than the charge they collected (R5 231, R3 137,
///     R13 127 ADC in run 97; R9 about 0), and the L-chain short ends a few
///     tens of ADC, which is what a trigger threshold looks like from here.
///     An end that did not fire reads 0 and carries no offset. Without S0
///     the strip total is bimodal: events without a short hit sum to C_long,
///     events with one to C_long + S0 x C_long/C_short, and the second peak
///     is the high-side tail the odd strips showed.
///   * LONG anchor C_long = the long end's own modal beam peak
///     (ReduceToAnchors): most beam events leave the short end silent, so
///     that peak is the full deposit on the long end to within the charge a
///     silent short end can hide, which is below its trigger threshold.
///   * the robust line through the beam-gated pairs (RidgeShortAnchor) has
///     slope -k = -C_long/C_short and intercept a = C_long + k S0, so
///         C_short = C_long / k,   S0 = (a - C_long) / k.
///     Both come from the same line; S0 is the sideways shift of the ridge
///     relative to the long-only peak.
///   * where the ridge is not measurable the short anchor is the parity's
///     median C_short/C_long (a preamp gain ratio, so it transfers across
///     strips) times C_long, or the global median, and S0 is the parity's
///     median offset in units of C_long, or 0.
///   * gain = TARGET / anchor per side, with TARGET = 1.0 a.u. (the notebook
///     uses 1000 ADC; only the overall scale differs); the decode applies
///     gain x (ADC - S0) to an end that fired (EnergyView::Decode).
///
/// Pass 2 — check, do not correct: with the offset removed the summed beam
/// peak of the pairs on the ridge (short inside the fit window) equals the
/// long-only peak, 1.0 a.u., by construction; a departure exposes a bad ridge
/// fit or a poor fallback. Censored events (short end 0) sum to the long end
/// alone and are left out so they cannot mask a bad fit.
///   * the old pass 2 rescaled both gains by 1/peak, which left the L/R ratio
///     untouched, so a wrong short anchor misallocated charge within a strip
///     but no longer shifts the strip total.
///   * when IGNORE_SHORT_STRIPS is set the decode keeps only the long end, so
///     the long end alone is normalised to 1.0 and pass 1 already provides it.
///
/// Pairs are collected UNGATED (all events with both ends firing) because the
/// ridge needs the off-centre crossings that a beam gate removes. Strips whose
/// ridge cannot be fitted fall back to the MEDIAN anchor of the strips that
/// could; that fallback assumes a common electronics scale across strips, which
/// the measured ridge slopes can be used to check rather than assume.
const Int_t kGmBins = 512;
/// Ridge-slope short anchor. In the long-vs-short plane a fixed deposit divided
/// between the two ends of a strip traces
///     ADC_long/C_long + ADC_short/C_short = 1,
/// a line of slope -C_long/C_short. The slope is the same for every deposit
/// energy, so the ridge direction alone gives the anchor ratio and C_short need
/// never be observed directly -- which matters because the beam is collimated
/// onto the long end and rarely deposits its full charge on the short one.
/// The plane also contains pile-up bands at 2x, 3x the single-particle deposit;
/// the 2-particle band falls into the single-particle window once
/// short > ~0.55*C_short, so the fit is capped well below that.
///
/// The fit window also starts well above the short end's threshold. Just above
/// threshold the short reading carries no position information: on the ends
/// that are censored for most beam events (the L end of the even strips on
/// 37Cl) the pairs that survive pile up in a flat blob at 0.04-0.15 C_long
/// with the long end at its full deposit, and on the other ends the slice
/// medians bend upward there too. A slice fit that sees that blob comes out
/// nearly flat (slopes of -0.1 to -0.5 on strips 2 and 4, ratios of 2-10) and
/// puts the normed band far below -1. Only the falling ridge beyond the blob
/// is fitted; it is sparse on the censored ends, so the slice count adapts to
/// the events available.
const Double_t kRidgeShortMinFrac =
    0.15; // start of the short window, as a fraction of C_long
const Double_t kRidgeShortMaxFrac =
    0.5;                            // cap on short, as a fraction of C_long
const Double_t kRidgeBandLo = 0.30; // single-particle band, x C_long
const Double_t kRidgeBandHi = 1.45;
const Int_t kRidgeSlices = 60;   // at most this many slices of the window
const Int_t kRidgeMinSlices = 6; // at least this many, for sparse ridges
const Long64_t kRidgeMinPerSlice = 200;
// Absolute noise floor per slice when the per-slice minimum is relaxed for
// thinly-spread ridges: below this a slice median is dominated by shot noise.
const Long64_t kRidgeNoiseFloor = 20;
const Int_t kRidgeMinPts = 6;
// C_short/C_long is a preamp gain ratio, so it is order unity (0.9-1.3 on
// 37Cl); outside this range the ridge fit has failed in a way the intercept
// test cannot see, typically a flat fit through the threshold blob.
const Double_t kRidgeRatioLo = 0.50;
const Double_t kRidgeRatioHi = 2.00;
const Double_t kGmEsumLo = 0.8; // a.u. eSum peak search window (pass 2)
const Double_t kGmEsumHi = 2.5;

/// Everything SaveRidgeFitPlots needs to redraw a fit: the slice medians it was
/// fitted to, the window they were taken from, and the resulting line. Filled
/// even when the fit is later rejected, so a bad fit can be looked at.
struct RidgeFit {
  std::vector<Double_t> x, y, ey; // slice centres, medians, median errors
  Double_t lo = 0.0, hi = 0.0;    // short-axis fit window
  Double_t slope = 0.0, intercept = 0.0;
  Double_t c_long = 0.0, c_short = 0.0;
  Bool_t fitted = kFALSE; // a line was fitted (it may still be rejected)
  // Diagnostic: how far a strip got before the fit, so a "not measurable"
  // log line can say which check rejected it instead of only slope/intercept.
  Long64_t n_gated = 0;        // pairs handed to RidgeShortAnchor
  Long64_t n_short_window = 0; // stayed inside the short window
  Long64_t n_long_band = 0;    // also inside the long band
  Int_t n_slices = 0;          // slices with >= kRidgeMinPerSlice
  Long64_t min_per_slice = 0;  // adaptive per-slice minimum used this strip
  const char *fail = "";       // the check that returned 0, if any
};

/// Theil-Sen robust line: slope = median of pairwise slopes, intercept = median
/// residual. Unlike OLS ("pol1" via a chi-square fit) a single outlying slice
/// median cannot bend the slope: it corrects the intercept, complementing
/// the slice median's own robustness. Returns kFALSE when the slice set
/// is too small or degenerate (all-same x).
Bool_t TheilSenLine(const std::vector<Double_t> &x,
                    const std::vector<Double_t> &y, Double_t &slope,
                    Double_t &intercept) {
  const Int_t n = Int_t(x.size());
  if (n < 2)
    return kFALSE;
  std::vector<Double_t> slopes;
  slopes.reserve(n * (n - 1) / 2);
  for (Int_t i = 0; i < n; i++)
    for (Int_t j = i + 1; j < n; j++) {
      Double_t dx = x[j] - x[i];
      if (dx == 0.0)
        continue;
      slopes.push_back((y[j] - y[i]) / dx);
    }
  if (slopes.empty())
    return kFALSE;
  slope = Median(slopes);
  std::vector<Double_t> resid;
  resid.reserve(n);
  for (Int_t i = 0; i < n; i++)
    resid.push_back(y[i] - slope * x[i]);
  intercept = Median(resid);
  return kTRUE;
}

/// Robust line through the single-particle charge-sharing ridge, returning
/// the fitted line's short-axis crossing, -intercept/slope (nonzero means the
/// ridge was measurable), and the line itself through
/// `slope_out`/`intercept_out`. The caller turns slope and intercept into the
/// short anchor and the short-end offset against the long end's modal peak
/// (see ComputeLRGainMatch). `c_long` is only the seed that places the short
/// window and the long band; the line's intercept must land near it. Returns
/// 0 when the ridge is not measurable.
Double_t RidgeShortAnchor(const std::vector<Float_t> &v_short,
                          const std::vector<Float_t> &v_long, Double_t c_long,
                          Double_t &slope_out, Double_t &intercept_out,
                          RidgeFit &dbg) {
  slope_out = 0.0;
  intercept_out = 0.0;
  dbg.c_long = c_long;
  dbg.n_gated = Long64_t(v_short.size());
  if (c_long <= 0 || v_short.size() != v_long.size() || v_short.size() < 500) {
    dbg.fail = c_long <= 0                       ? "c_long<=0"
               : v_short.size() != v_long.size() ? "size mismatch"
                                                 : "fewer than 500 gated pairs";
    return 0.0;
  }
  const Double_t hi = kRidgeShortMaxFrac * c_long;
  const Double_t lo = kRidgeShortMinFrac * c_long;
  // Pairs inside the short window and the long band, kept so the slice count
  // can be chosen from how many there are before they are binned.
  std::vector<std::pair<Double_t, Double_t>> in_window;
  for (Int_t j = 0; j < Int_t(v_short.size()); j++) {
    Double_t sh = Double_t(v_short[j]), lg = Double_t(v_long[j]);
    if (sh < lo || sh >= hi)
      continue;
    dbg.n_short_window++;
    if (lg <= kRidgeBandLo * c_long || lg >= kRidgeBandHi * c_long)
      continue;
    dbg.n_long_band++;
    in_window.push_back(std::make_pair(sh, lg));
  }
  // Slice count from the population: kRidgeMinPerSlice per slice on average,
  // between kRidgeMinSlices (sparse censored ends) and kRidgeSlices.
  Int_t n_slices = Int_t(Long64_t(in_window.size()) / kRidgeMinPerSlice);
  n_slices = TMath::Max(kRidgeMinSlices, TMath::Min(kRidgeSlices, n_slices));
  std::vector<std::vector<Double_t>> slice(n_slices);
  for (Int_t j = 0; j < Int_t(in_window.size()); j++) {
    Int_t b = Int_t((in_window[j].first - lo) / (hi - lo) * n_slices);
    if (b >= 0 && b < n_slices)
      slice[b].push_back(in_window[j].second);
  }
  std::vector<Double_t> x, y, ey;
  {
    // Per-slice min: strict floor first; if kRidgeMinPts slices don't fill,
    // halve & retry down to the noise floor (the censored ends are sparse).
    Long64_t min_per_slice = kRidgeMinPerSlice;
    for (;;) {
      x.clear();
      y.clear();
      ey.clear();
      for (Int_t b = 0; b < n_slices; b++) {
        if (Long64_t(slice[b].size()) < min_per_slice)
          continue;
        std::sort(slice[b].begin(), slice[b].end());
        const Int_t nb = Int_t(slice[b].size());
        // The IQR reads the sorted slice, so it comes before the median,
        // which only partially orders it.
        const Double_t iqr = slice[b][nb * 3 / 4] - slice[b][nb / 4];
        const Double_t med = Median(slice[b]);
        x.push_back(lo + (b + 0.5) * (hi - lo) / n_slices);
        y.push_back(med);
        ey.push_back(1.253 * (iqr / 1.349) / TMath::Sqrt(Double_t(nb)));
      }
      if (Int_t(x.size()) >= kRidgeMinPts || min_per_slice <= kRidgeNoiseFloor)
        break;
      // Halve the floor and retry.
      min_per_slice = TMath::Max(kRidgeNoiseFloor, min_per_slice / 2);
    }
    dbg.min_per_slice = min_per_slice;
  }
  dbg.x = x;
  dbg.y = y;
  dbg.ey = ey;
  dbg.lo = lo;
  dbg.hi = hi;
  dbg.n_slices = Int_t(x.size());
  if (Int_t(x.size()) < kRidgeMinPts) {
    dbg.fail = "fewer than kRidgeMinPts filled slices";
    return 0.0;
  }
  // TheilSen: median of pairwise slopes; outlying slice medians (2/3-particle
  // bands, low-long background) can't bend the slope as a chi-square pol1.
  Double_t slope = 0.0, inter = 0.0;
  if (!TheilSenLine(x, y, slope, inter)) {
    dbg.fail = "degenerate slice set (TheilSenLine)";
    return 0.0;
  }
  slope_out = slope;
  intercept_out = inter;
  dbg.slope = slope;
  dbg.intercept = inter;
  dbg.fitted = kTRUE;
  if (slope < 0)
    dbg.c_short = -inter / slope;
  if (slope >= 0) {
    dbg.fail = "slope >= 0";
    return 0.0;
  }
  // Line must cross the long axis near the modal seed: above it by the
  // short-end offset times the slope (up to ~20% on the R-chain short ends,
  // see ComputeLRGainMatch); a large departure means the band selection
  // missed the single-particle ridge; window loose.
  if (inter < 0.80 * c_long || inter > 1.35 * c_long) {
    dbg.fail = "intercept outside [0.80, 1.35]*C_long";
    return 0.0;
  }
  // The line's short-axis crossing, as the "measurable" flag; the caller
  // works from slope and intercept directly.
  return -inter / slope;
}

/// One plot per strip under `<plot_subdir>/ridge`, named `ridge_s<NN>`: the
/// beam-gated long-vs-short plane, the slice medians the fit was actually given
/// (black), and the fitted line (violet). Drawn for every strip, including the
/// ones whose fit was rejected, so a bad ridge can be seen rather than inferred
/// from the slope in the log.
void SaveRidgeFitPlots(const StripPairSamples pairs[18], const RidgeFit dbg[18],
                       const TString &plot_subdir) {
  TString subdir = plot_subdir + "/ridge";
  for (Int_t s = 1; s <= 16; s++) {
    const RidgeFit &d = dbg[s];
    const StripPairSamples &p = pairs[s];
    if (d.c_long <= 0 || p.gated_short.empty())
      continue;
    Bool_t l_is_long = (LongSide(s) == 'L');
    // A little past the fit window, so the cap the fit stops at is visible.
    Double_t xhi = (d.hi > 0 ? d.hi : kRidgeShortMaxFrac * d.c_long) * 1.35;
    TH2F *h = new TH2F(Form("h_ridge_s%d", s),
                       Form(";Strip %d Short (%c) #DeltaE [ADC];Strip %d Long "
                            "(%c) #DeltaE [ADC]",
                            s, l_is_long ? 'R' : 'L', s, l_is_long ? 'L' : 'R'),
                       200, 0.0, xhi, 200, 0.0, 1.60 * d.c_long);
    for (Int_t j = 0; j < Int_t(p.gated_short.size()); j++)
      h->Fill(Double_t(p.gated_short[j]), Double_t(p.gated_long[j]));
    TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::Configure2DHistogram(h, cv);
    h->Draw("COLZ");
    TGraphErrors *g = nullptr;
    if (!d.x.empty()) {
      g = new TGraphErrors(Int_t(d.x.size()), &d.x[0], &d.y[0], nullptr,
                           &d.ey[0]);
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.8);
      g->SetMarkerColor(kBlack);
      g->SetLineColor(kBlack);
      g->Draw("P SAME");
    }
    // Extended to short = 0 on purpose: where the line lands above the
    // long-only peak is the short-end offset times the slope.
    TLine *lf = nullptr;
    if (d.fitted) {
      lf = new TLine(0.0, d.intercept, xhi, d.intercept + d.slope * xhi);
      lf->SetLineColor(kViolet + 2);
      lf->SetLineWidth(2);
      lf->Draw();
    }
    // The short window the slices were taken from, so the threshold blob left
    // of it can be seen to be excluded.
    TLine *lw[2] = {nullptr, nullptr};
    if (d.hi > 0) {
      const Double_t ymax = 1.60 * d.c_long;
      lw[0] = new TLine(d.lo, 0.0, d.lo, ymax);
      lw[1] = new TLine(d.hi, 0.0, d.hi, ymax);
      for (Int_t k = 0; k < 2; k++) {
        lw[k]->SetLineColor(kGray + 2);
        lw[k]->SetLineStyle(2);
        lw[k]->Draw();
      }
    }
    if (Constants::SavePlots())
      PlottingUtils::SaveFigure(cv, Form("ridge_s%02d", s), subdir,
                                PlotSaveOptions::kLINEAR);
    delete cv;
    delete g;
    delete lf;
    delete lw[0];
    delete lw[1];
    delete h;
  }
}

/// Histogram-mode peak finder, mirroring the notebook's find_peak(): histogram
/// `v` over [lo, hi] with kGmBins bins, skip the first skip_frac of bins (to
/// avoid the threshold pile), return the max-bin centre. Returns 0 when empty.
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
                        const StripPairSamples pairs[18],
                        const TString &plot_subdir) {
  RidgeFit ridge_dbg[18];
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

  // Pass 1: LONG anchor = modal peak (ReduceToAnchors); from the ridge line
  // SHORT anchor = C_long/k and the short-end offset S0 = (a - C_long)/k.
  // Unmeasurable: short = parity median ratio x long, S0 = parity median
  // offset x long.

  Bool_t matched[18] = {kFALSE};
  Double_t long_anchor_adc[18] = {0};
  Double_t short_anchor_adc[18] = {0};
  Double_t short_offset_adc[18] = {0};
  Bool_t ridge_ok[18] = {kFALSE};
  std::vector<Double_t> ratios_found[2], offsets_found[2];
  for (Int_t s = 1; s <= 16; s++) {
    if (idx_l[s] < 0 || idx_r[s] < 0)
      continue;
    Bool_t l_is_long = (LongSide(s) == 'L');
    ChannelCal &c_long = chans[l_is_long ? idx_l[s] : idx_r[s]];
    ChannelCal &c_short = chans[l_is_long ? idx_r[s] : idx_l[s]];
    if (!IsCalibrated(c_long)) {
      Constants::DetailErr()
          << "  strip " << s
          << ": long side uncalibrated; skipping L/R gain match" << std::endl;
      continue;
    }
    Double_t peak_short = 0.0, peak_long = 0.0, offset = 0.0;

    {
      const StripPairSamples &p = pairs[s];
      // Beam-gated ridge fit: the gate removes pile-up/junk while keeping
      // off-centre crossings; ungated's 2/3-particle bands steepen the slope.
      const std::vector<Float_t> &v_short = p.gated_short;
      const std::vector<Float_t> &v_long = p.gated_long;
      Double_t slope = 0.0, inter = 0.0;
      const Double_t crossing = RidgeShortAnchor(
          v_short, v_long, c_long.fit_adc, slope, inter, &ridge_dbg[s]);
      if (crossing > 0) {
        // The line: long = a - k short. C_short = C_long/k; the offset is
        // the short reading at which the line reaches the long-only peak.
        const Double_t k = -slope;
        peak_long = c_long.fit_adc;
        peak_short = peak_long / k;
        offset = (inter - peak_long) / k;
        Constants::Detail()
            << "  strip " << s << " ridge slope=" << Form("%.3f", slope)
            << " intercept=" << Form("%.1f", inter)
            << " ADC; long anchor=" << Form("%.1f", peak_long)
            << " (modal)  short anchor=" << Form("%.1f", peak_short)
            << " ADC  short offset=" << Form("%.1f", offset) << " ADC ("
            << Form("%.3f", offset / peak_long) << " x C_long)" << std::endl;
      } else
        Constants::DetailErr()
            << "  strip " << s << ": ridge not measurable "
            << "(slope=" << Form("%.3f", slope)
            << ", intercept=" << Form("%.1f", inter) << ")"
            << (ridge_dbg[s].fail[0] ? Form(" [%s]", ridge_dbg[s].fail) : "")
            << " gated=" << Form("%lld", ridge_dbg[s].n_gated)
            << " short_win=" << Form("%lld", ridge_dbg[s].n_short_window)
            << " long_band=" << Form("%lld", ridge_dbg[s].n_long_band)
            << " slices=" << ridge_dbg[s].n_slices
            << " min_per_slice=" << ridge_dbg[s].min_per_slice << std::endl;
    }

    if (peak_short <= 0 || peak_long <= 0)
      continue;
    // C_short/C_long = -1/slope is a preamp-gain ratio, so order unity.
    // Catches ridges too flat for the intercept test: ratio absurd.
    Double_t ratio = peak_short / peak_long;
    if (ratio < kRidgeRatioLo || ratio > kRidgeRatioHi) {
      Constants::DetailErr()
          << "  strip " << s << ": ridge ratio " << Form("%.2f", ratio)
          << " outside [" << kRidgeRatioLo << ", " << kRidgeRatioHi
          << "]; rejecting anchor " << Form("%.1f", peak_short) << " ADC"
          << std::endl;
      continue;
    }
    long_anchor_adc[s] = peak_long;
    short_anchor_adc[s] = peak_short;
    short_offset_adc[s] = offset;
    ridge_ok[s] = kTRUE;
    ratios_found[s % 2].push_back(ratio);
    offsets_found[s % 2].push_back(offset / peak_long);
    c_short.ridge_ratio = ratio;
    c_short.ridge_offset = offset / peak_long;
  }

  SaveRidgeFitPlots(pairs, ridge_dbg, plot_subdir);

  // Fall back on the median RATIO (preamp property, transfers across strips;
  // ADC anchors don't), per parity: L/R preamps differ (odd ~1.2, even ~0.9).
  Double_t median_ratio[2] = {0.0, 0.0};
  std::vector<Double_t> all_ratios;
  for (Int_t par = 0; par < 2; par++) {
    std::vector<Double_t> &v = ratios_found[par];
    all_ratios.insert(all_ratios.end(), v.begin(), v.end());
    if (v.size() < 2)
      continue;
    median_ratio[par] = Median(v);
  }
  Double_t global_ratio = 0.0;
  if (!all_ratios.empty())
    global_ratio = Median(all_ratios);
  Constants::Detail() << "  ridge ratio medians: odd="
                      << Form("%.3f", median_ratio[1])
                      << " even=" << Form("%.3f", median_ratio[0])
                      << " all=" << Form("%.3f", global_ratio) << std::endl;
  // The offset fallback, likewise per parity (it is a chain property), in
  // units of C_long; 0 where the parity measured none.
  Double_t median_offset[2] = {0.0, 0.0};
  for (Int_t par = 0; par < 2; par++) {
    std::vector<Double_t> &v = offsets_found[par];
    if (v.size() < 2)
      continue;
    median_offset[par] = Median(v);
  }
  Constants::Detail() << "  ridge offset medians (x C_long): odd="
                      << Form("%.3f", median_offset[1])
                      << " even=" << Form("%.3f", median_offset[0])
                      << std::endl;

  for (Int_t s = 1; s <= 16; s++) {
    if (idx_l[s] < 0 || idx_r[s] < 0)
      continue;
    Bool_t l_is_long = (LongSide(s) == 'L');
    ChannelCal &c_long = chans[l_is_long ? idx_l[s] : idx_r[s]];
    ChannelCal &c_short = chans[l_is_long ? idx_r[s] : idx_l[s]];
    if (!IsCalibrated(c_long))
      continue;
    if (short_anchor_adc[s] <= 0) {
      // No ridge: the modal peak is the best long anchor available.
      long_anchor_adc[s] = c_long.fit_adc;
      Double_t r = median_ratio[s % 2] > 0 ? median_ratio[s % 2] : global_ratio;
      if (r <= 0) {
        Constants::DetailErr() << "  strip " << s
                               << ": no ridge and no ratio fallback; keeping "
                                  "independent gains"
                               << std::endl;
        continue;
      }
      short_anchor_adc[s] = r * long_anchor_adc[s];
      short_offset_adc[s] = median_offset[s % 2] * long_anchor_adc[s];
      Constants::Detail() << "  strip " << s << " short_anchor=ratio fallback "
                          << Form("%.3f", r)
                          << " x C_long = " << Form("%.1f", short_anchor_adc[s])
                          << " ADC"
                          << (median_ratio[s % 2] > 0
                                  ? ""
                                  : " (global, parity had none)")
                          << "  short offset=parity median "
                          << Form("%.3f", median_offset[s % 2])
                          << " x C_long = " << Form("%.1f", short_offset_adc[s])
                          << " ADC" << std::endl;
    }
    c_long.gain = 1.0 / long_anchor_adc[s];
    c_short.gain = 1.0 / short_anchor_adc[s];
    c_short.offset_adc = short_offset_adc[s];
    matched[s] = kTRUE;
    Constants::Detail() << "  strip " << s << " L/R match: long anchor="
                        << Form("%.1f", long_anchor_adc[s])
                        << " ADC  short anchor="
                        << Form("%.1f", short_anchor_adc[s])
                        << " ADC  short offset="
                        << Form("%.1f", short_offset_adc[s]) << " ADC"
                        << (ridge_ok[s] ? "" : " (fallback)") << std::endl;
  }

  // Pass 2: check, do not correct. On the pairs inside the ridge fit window
  // (short between kRidgeShortMinFrac and kRidgeShortMaxFrac of the long
  // anchor, so on the ridge and clear of the threshold blob) eSum, with the
  // short offset removed, peaks at the long-only peak, 1.0 a.u., by
  // construction; departures expose bad ridge fits or a poor fallback —
  // rescaling would hide that failure.
  Double_t esum_peak[18] = {0};
  for (Int_t s = 1; s <= 16; s++) {
    if (!matched[s])
      continue;
    const StripPairSamples &p = pairs[s];
    Bool_t l_is_long = (LongSide(s) == 'L');
    Double_t g_long = Gain(chans[l_is_long ? idx_l[s] : idx_r[s]]);
    Double_t g_short = Gain(chans[l_is_long ? idx_r[s] : idx_l[s]]);
    std::vector<Float_t> esum;
    esum.reserve(p.gated_long.size());
    const Double_t sh_lo = kRidgeShortMinFrac * long_anchor_adc[s];
    const Double_t sh_hi = kRidgeShortMaxFrac * long_anchor_adc[s];
    for (Int_t j = 0; j < Int_t(p.gated_long.size()); j++) {
      if (p.gated_short[j] < sh_lo || p.gated_short[j] >= sh_hi)
        continue;
      esum.push_back(Float_t(
          g_long * Double_t(p.gated_long[j]) +
          g_short * (Double_t(p.gated_short[j]) - short_offset_adc[s])));
    }
    esum_peak[s] = GmFindPeak(esum, kGmEsumLo, kGmEsumHi, 0.0);
  }
  for (Int_t s = 1; s <= 16; s++) {
    if (!matched[s])
      continue;
    if (esum_peak[s] <= 0) {
      Constants::DetailErr()
          << "  strip " << s << ": no summed beam peak in (" << kGmEsumLo
          << ", " << kGmEsumHi << ") a.u." << std::endl;
      continue;
    }
    Double_t dev = esum_peak[s] - 1.0;
    Constants::Detail() << "  strip " << s
                        << " summed beam peak=" << Form("%.4f", esum_peak[s])
                        << " a.u. (" << Form("%+.1f%%", 100.0 * dev) << ")"
                        << (TMath::Abs(dev) > 0.05 ? "  <-- check ridge fit"
                                                   : "")
                        << std::endl;
  }
}

/// Cathode uses median + IQR/1.349 (asymmetric tail, no clean peak). Every
/// other channel is anchored on the modal bin of its beam-peak histogram
/// (RobustPeakSeed: 100 bins over the 5th-95th percentile of the samples), with
/// no fit on top and no width: nothing downstream uses one.
void ReduceToAnchors(std::vector<ChannelCal> &chans,
                     std::vector<std::vector<Float_t>> &samples,
                     const StripPairSamples pairs[18],
                     const TString &plot_subdir) {
  Int_t n_chans = Int_t(chans.size());

  for (Int_t i = 0; i < n_chans; i++) {
    ChannelCal &c = chans[i];
    std::vector<Float_t> &v = samples[i];
    c.n_samples = Long64_t(v.size());

    if (Long64_t(v.size()) < kMinSamples) {
      c.fit_adc = 0;
      continue;
    }
    if (c.side == 'C') {
      // Cathode: the median (asymmetric tail, no clean peak).
      c.fit_adc = Median(v);
    } else {
      // The modal bin, never the tail-biased mean.
      Double_t mode = 0.0, rsigma = 0.0;
      RobustPeakSeed(v, mode, rsigma);
      c.fit_adc = mode;
    }
    Constants::Detail() << "  " << c.name << " anchor[ADC]=" << c.fit_adc
                        << " (n=" << c.n_samples << ")" << std::endl;
  }

  // L/R match (strips 1-16): both anchors from the charge-sharing ridge, the
  // modal peak above being the seed and the fallback.
  if (pairs)
    ComputeLRGainMatch(chans, pairs, plot_subdir);
}

/// Writes a one-row `calibration` tree into the open file `dst` (typically a
/// per-subfile .cal.root). `align` carries the per-strip alignment factors;
/// pass nullptr when they have not been computed yet (identity is written).
void WriteCalibrationTree(TFile *dst, const std::vector<ChannelCal> &chans,
                          const StripAlignmentResult *align) {
  dst->cd();
  if (TObject *old = dst->Get("calibration"))
    old->Delete();
  TTree *cal = new TTree("calibration", "Per-channel normMUSIC calibration");
  Float_t gain[kMaxChannels] = {0};
  Float_t fit_adc[kMaxChannels] = {0};
  Long64_t fit_n[kMaxChannels] = {0};
  Bool_t ok[kMaxChannels] = {0};
  // GainLeft[k]/GainRight[k] scale LeftdE[k]/RightdE[k], the two ends of
  // strip k+1; GainStrip0/GainStrip17 the unsegmented strips. EnergyView
  // reads these to calibrate on the fly.
  Float_t gain_left[16] = {0}, gain_right[16] = {0};
  Float_t gain_strip0 = 0.0f, gain_strip17 = 0.0f;
  Float_t gain_cathode = 0.0f;
  Float_t gain_grid = 0.0f;
  // OffsetLeft[k]/OffsetRight[k]: ADC subtracted from LeftdE[k]/RightdE[k]
  // before the gain, only when that end fired (see ComputeLRGainMatch: the
  // short-end offset). Zero on the long ends.
  Float_t offset_left[16] = {0}, offset_right[16] = {0};
  // Per-strip ridge ratio and offset (x C_long) measured in THIS subfile, 0
  // where the ridge was not measurable. AggregateRidgeRatiosForRun medians
  // these across a run.
  Float_t ridge_ratio[18] = {0};
  Float_t ridge_offset[18] = {0};
  Float_t long_anchor[18] = {0};
  Int_t n_actual = TMath::Min(Int_t(chans.size()), kMaxChannels);
  for (Int_t k = 0; k < n_actual; k++) {
    const ChannelCal &c = chans[k];
    ok[k] = IsCalibrated(c) || c.gain > 0;
    gain[k] = ok[k] ? Float_t(Gain(c)) : 0.0f;
    fit_adc[k] = Float_t(c.fit_adc);
    fit_n[k] = c.n_samples;
    if (c.strip >= 1 && c.strip <= 16 && (c.side == 'L' || c.side == 'R')) {
      if (c.ridge_ratio > 0) {
        ridge_ratio[c.strip] = Float_t(c.ridge_ratio);
        ridge_offset[c.strip] = Float_t(c.ridge_offset);
      }
      // The long anchor in use (the modal peak): AggregateRidgeRatiosForRun
      // rebuilds the short gain and offset from it.
      if (c.side == LongSide(c.strip) && ok[k] && gain[k] > 0)
        long_anchor[c.strip] = Float_t(1.0 / gain[k]);
    }
    if (c.side == 'S' && c.strip == 0)
      gain_strip0 = gain[k];
    else if (c.side == 'S' && c.strip == 17)
      gain_strip17 = gain[k];
    else if (c.side == 'L' && c.strip >= 1 && c.strip <= 16) {
      gain_left[c.strip - 1] = gain[k];
      offset_left[c.strip - 1] = Float_t(c.offset_adc);
    } else if (c.side == 'R' && c.strip >= 1 && c.strip <= 16) {
      gain_right[c.strip - 1] = gain[k];
      offset_right[c.strip - 1] = Float_t(c.offset_adc);
    } else if (c.side == 'C')
      gain_cathode = gain[k];
    else if (c.side == 'G')
      gain_grid = gain[k];
  }
  cal->Branch("Gain", gain, Form("Gain[%d]/F", kMaxChannels));
  cal->Branch("Ok", ok, Form("Ok[%d]/O", kMaxChannels));
  cal->Branch("FitADC", fit_adc, Form("FitADC[%d]/F", kMaxChannels));
  cal->Branch("FitN", fit_n, Form("FitN[%d]/L", kMaxChannels));
  cal->Branch("GainLeft", gain_left, "GainLeft[16]/F");
  cal->Branch("GainRight", gain_right, "GainRight[16]/F");
  cal->Branch("GainStrip0", &gain_strip0, "GainStrip0/F");
  cal->Branch("GainStrip17", &gain_strip17, "GainStrip17/F");
  cal->Branch("GainCathode", &gain_cathode, "GainCathode/F");
  cal->Branch("GainGrid", &gain_grid, "GainGrid/F");
  cal->Branch("OffsetLeft", offset_left, "OffsetLeft[16]/F");
  cal->Branch("OffsetRight", offset_right, "OffsetRight[16]/F");
  cal->Branch("RidgeRatio", ridge_ratio, "RidgeRatio[18]/F");
  cal->Branch("RidgeOffset", ridge_offset, "RidgeOffset[18]/F");
  cal->Branch("LongAnchor", long_anchor, "LongAnchor[18]/F");

  // Per-strip two-point alignment: EnergyView multiplies every end of a strip
  // by StripFactor after the gain and adds StripOffset to the long end (or the
  // unsegmented value), so beam = 1 and pile-up = 2 on the strip total.
  // Defaults 1 and 0 (identity).
  Float_t strip_factor[18], strip_offset[18];
  for (Int_t s = 0; s < 18; s++) {
    strip_factor[s] = 1.0f;
    strip_offset[s] = 0.0f;
  }
  if (align && align->ok)
    for (Int_t s = 0; s < 18; s++) {
      strip_factor[s] = Float_t(align->factors[s]);
      strip_offset[s] = Float_t(align->offsets[s]);
    }
  cal->Branch("StripFactor", strip_factor, "StripFactor[18]/F");
  cal->Branch("StripOffset", strip_offset, "StripOffset[18]/F");
  cal->Fill();
  cal->Write("calibration", TObject::kOverwrite);
}

/// Writes the per-channel gain table (tree "calibration") into the subfile's
/// own events file. No per-event calibrated tree is produced: downstream
/// readers recover a.u. on the fly via gain x raw ADC (EnergyView), so the
/// raw events tree plus this one-row gain table fully determine every
/// calibrated value.
void WriteCalibrationToEvents(const FileSpec &spec,
                              const std::vector<ChannelCal> &chans,
                              const StripAlignmentResult *align) {
  TString events_subpath = FileSet::EventsName(spec) + ".root";
  TFile *f = IO::OpenForWriting(events_subpath, "UPDATE");
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open " << events_subpath << " to write calibration"
              << std::endl;
    if (f)
      delete f;
    return;
  }
  WriteCalibrationTree(f, chans, align);
  Constants::Detail() << "  wrote calibration into " << events_subpath
                      << std::endl;
  f->Close();
  delete f;
}

// Per-channel calibrated overlay for one subfile, via AttachCalSidecar (the
// same path downstream macros use). The sidecar must already be on disk.
void SaveDynamicRangeOverlay(const FileSpec &spec, const TString &plot_subdir,
                             const TString &file_label) {
  const Int_t kNStrips = 18;
  const Double_t emin = Constants::cfg.STRIP_DE_MIN_NORMED;
  const Double_t emax = Constants::cfg.STRIP_DE_MAX_NORMED;
  // 0.025 a.u. bins: a few per beam-peak width, so eighteen overlaid curves
  // read as curves rather than as noise.
  const Int_t nbins = TMath::Max(20, Int_t((emax - emin) / 0.025 + 0.5));
  TH1D *h[kNStrips];
  for (Int_t s = 0; s < kNStrips; s++) {
    h[s] = new TH1D(Form("h_dynrange_%s_S%d", file_label.Data(), s),
                    ";#DeltaE [a.u.];Counts", nbins, emin, emax);
    h[s]->SetDirectory(nullptr);
  }
  TString sub = FileSet::EventsName(spec) + ".root";
  TFile *sf = IO::OpenForReading(sub);
  if (!sf || sf->IsZombie()) {
    if (sf)
      sf->Close();
    for (Int_t s = 0; s < kNStrips; s++)
      delete h[s];
    return;
  }
  TTree *tree = static_cast<TTree *>(sf->Get("events"));
  if (!tree) {
    sf->Close();
    delete sf;
    for (Int_t s = 0; s < kNStrips; s++)
      delete h[s];
    return;
  }
  EnergyView ev;
  ev.Attach(tree);
  if (!ev.is_normed) {
    sf->Close();
    delete sf;
    for (Int_t s = 0; s < kNStrips; s++)
      delete h[s];
    return;
  }
  // The grid alone on the same axis, so its beam classes read against the
  // strips' range.
  TH1D *hg = new TH1D(Form("h_dynrange_%s_grid", file_label.Data()),
                      ";Grid #DeltaE [a.u.];Counts", nbins, emin, emax);
  hg->SetDirectory(nullptr);
  Long64_t n = tree->GetEntries();
  for (Long64_t j = 0; j < n; j++) {
    tree->GetEntry(j);
    ev.Decode();
    for (Int_t s = 0; s < kNStrips; s++) {
      Double_t v = ev.Total(s);
      if (v > 0)
        h[s]->Fill(v);
    }
    if (ev.grid > 0.0)
      hg->Fill(ev.grid);
  }
  sf->Close();
  delete sf;
  std::vector<Int_t> colors = PlottingUtils::GetDefaultColors();
  Double_t y_top = 0;
  for (Int_t s = 0; s < kNStrips; s++) {
    Double_t m = h[s]->GetMaximum();
    if (m > y_top)
      y_top = m;
  }
  TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
  cv->SetRightMargin(0.20);
  Bool_t first = kTRUE;
  for (Int_t s = 0; s < kNStrips; s++) {
    Int_t color = colors[s % Int_t(colors.size())];
    h[s]->SetLineColor(color);
    h[s]->SetLineWidth(2);
    h[s]->SetMaximum(1.15 * y_top);
    h[s]->Draw(first ? "HIST" : "HIST SAME");
    first = kFALSE;
  }
  TLegend *leg = PlottingUtils::AddLegend(0.81, 0.99, 0.10, 0.95);
  for (Int_t s = 0; s < kNStrips; s++)
    leg->AddEntry(h[s], Form("S%d", s), "l");
  leg->Draw();
  if (Constants::SavePlots())
    PlottingUtils::SaveFigure(cv, "dynamic_range_check", plot_subdir,
                              PlotSaveOptions::kLOG);
  delete cv;
  delete leg;
  TCanvas *cg = PlottingUtils::GetConfiguredCanvas(kTRUE);
  PlottingUtils::ConfigureAndDrawHistogram(hg, kBlue + 1);
  if (Constants::SavePlots())
    PlottingUtils::SaveFigure(cg, "grid_spectrum", plot_subdir,
                              PlotSaveOptions::kLOG);
  delete cg;
  delete hg;
  for (Int_t s = 0; s < kNStrips; s++)
    delete h[s];
}

/// Peak of a smoothed projection in [lo, hi]: the highest bin of the window,
/// refined to sub-bin precision by the parabola through it and its two
/// neighbours (kept within half a bin). No fit. Returns 0 when the window is
/// empty or its highest bin sits on the window edge, which means the peak
/// itself lies outside the window. The window is what selects the feature:
/// picking the local maximum nearest a target instead let a ripple of the
/// smoothed tail at 2.1-2.2 a.u. win over the real pile-up peak at 1.8 on
/// the strips whose beam sits low, which threw their alignment by 20%.
Double_t SmoothedPeakIn(TH1D *proj, Double_t lo, Double_t hi) {
  Int_t b_lo = proj->FindBin(lo);
  Int_t b_hi = proj->FindBin(hi);
  Int_t b_max = -1;
  Double_t val_max = 0;
  for (Int_t b = b_lo; b <= b_hi; b++) {
    Double_t v = proj->GetBinContent(b);
    if (v > val_max) {
      val_max = v;
      b_max = b;
    }
  }
  if (b_max < 0 || b_max == b_lo || b_max == b_hi)
    return 0.0;
  Double_t peak = proj->GetBinCenter(b_max);
  if (b_max > 1 && b_max < proj->GetNbinsX()) {
    const Double_t ym = proj->GetBinContent(b_max - 1);
    const Double_t y0 = proj->GetBinContent(b_max);
    const Double_t yp = proj->GetBinContent(b_max + 1);
    const Double_t denom = ym - 2.0 * y0 + yp;
    if (denom < 0.0) {
      Double_t shift = 0.5 * (ym - yp) / denom;
      shift = TMath::Max(-0.5, TMath::Min(0.5, shift));
      peak += shift * proj->GetBinWidth(b_max);
    }
  }
  return peak;
}

/// Per-strip two-point alignment of the strip total. Decodes every event with
/// the per-channel gains already on disk, projects each strip's total, and
/// locates two peaks in the smoothed projection: the beam (B, the highest
/// bin in [0.6, 1.5]) and the two-particle pile-up (P, the highest bin in
/// [1.7 B, 2.3 B], a window placed from the beam found so a strip whose beam
/// sits at 0.9 is searched around 1.8 and the 3-particle shoulder at 3 B
/// stays out).
/// Two beam particles in one event deposit exactly twice what one does, so
/// the total is put on the scale that has B at 1.0 and P at 2.0:
///
///     total' = factor * total + offset,
///     factor = 1 / (P - B),  offset = 1 - factor * B.
///
/// A beam-only anchor cannot separate a gain from a pedestal; the pile-up
/// peak is the second point that fixes both, and the summed spectra then show
/// the beam at 1 and the pile-up at 2 on every strip instead of pile-up peaks
/// spread over 1.8-2.2 by per-channel baseline offsets. Where no pile-up peak
/// is found the strip keeps a factor-only alignment, 1 / B, with no offset.
/// The offset is a per-strip quantity; EnergyView carries it on the long end
/// (the unsegmented value on strips 0 and 17), so a strip's total is still
/// the plain sum of its ends.
///
/// Uses ALL events (not beam-gated): the beam dominates every strip's
/// histogram by a wide margin, and the pile-up peak is the next feature. The
/// peaks come from 0.005 a.u. bins with parabolic sub-bin refinement
/// (SmoothedPeakIn); a percent-level alignment needs peaks located well
/// inside a percent, which the earlier 0.05 a.u. bins could not do.
StripAlignmentResult FindStripCentroidAlignment(const FileSpec &spec,
                                                const TString &plot_subdir,
                                                const TString &file_label) {
  const Int_t kNStrips = 18;
  const Int_t kNHistBins = 2000; // 0.005 a.u. bins
  const Double_t kHistMin = 0.0;
  const Double_t kHistMax = 10.0;
  const Int_t kSmoothTimes = 5;
  const Double_t kBeamLo = 0.6, kBeamHi = 1.5;   // beam peak window
  const Double_t kPileLoX = 1.7, kPileHiX = 2.3; // 2-particle window, x beam
  StripAlignmentResult result;
  for (Int_t s = 0; s < kNStrips; s++) {
    result.factors[s] = 1.0;
    result.offsets[s] = 0.0;
    result.centroids[s] = 0.0;
    result.pileups[s] = 0.0;
  }

  TString sub = FileSet::EventsName(spec) + ".root";
  TFile *sf = IO::OpenForReading(sub);
  if (!sf || sf->IsZombie()) {
    if (sf)
      sf->Close();
    delete sf;
    return result;
  }
  TTree *tree = static_cast<TTree *>(sf->Get("events"));
  if (!tree) {
    sf->Close();
    delete sf;
    return result;
  }
  EnergyView ev;
  ev.Attach(tree);
  if (!ev.is_normed) {
    std::cerr << "  " << file_label
              << ": calibration tree not found -- skipping alignment"
              << std::endl;
    sf->Close();
    delete sf;
    return result;
  }

  TH2D *h2 = new TH2D(Form("h2_strip_align_%s", file_label.Data()),
                      ";Strip number;#DeltaE [a.u.]", kNStrips, -0.5,
                      kNStrips - 0.5, kNHistBins, kHistMin, kHistMax);
  h2->SetDirectory(nullptr);

  Long64_t n = tree->GetEntries();
  Long64_t n_used = 0;
  for (Long64_t j = 0; j < n; j++) {
    tree->GetEntry(j);
    ev.Decode();
    Bool_t any = kFALSE;
    for (Int_t s = 0; s < kNStrips; s++) {
      Double_t v = ev.Total(s);
      if (v <= 0)
        continue;
      h2->Fill(Double_t(s), v);
      any = kTRUE;
    }
    if (any)
      n_used++;
  }
  sf->Close();
  delete sf;

  Constants::Detail() << "  strip alignment: " << n_used << " events decoded"
                      << std::endl;

  Bool_t beam_ok[kNStrips] = {kFALSE};
  Bool_t pile_ok[kNStrips] = {kFALSE};
  Int_t valid_strips = 0;

  for (Int_t s = 0; s < kNStrips; s++) {
    Int_t bin_ix = s + 1;
    TH1D *proj = h2->ProjectionY(
        Form("hproj_align_%s_s%d", file_label.Data(), s), bin_ix, bin_ix);
    proj->SetDirectory(nullptr);
    Long64_t n_entries = Long64_t(proj->GetEntries());
    if (n_entries < kMinSamples) {
      Constants::DetailErr()
          << "  strip " << s << ": too few entries for alignment (" << n_entries
          << ")" << std::endl;
      delete proj;
      continue;
    }
    proj->Smooth(kSmoothTimes);

    Double_t beam = SmoothedPeakIn(proj, kBeamLo, kBeamHi);
    if (beam <= 0) {
      Constants::DetailErr()
          << "  strip " << s << ": no beam peak in [" << kBeamLo << ", "
          << kBeamHi << "] a.u.; left unaligned" << std::endl;
      delete proj;
      continue;
    }
    // The 2-particle peak sits near twice the beam; a window placed from the
    // beam keeps the 3-particle shoulder and the inter-peak tail out, and a
    // peak on the window edge is rejected by SmoothedPeakIn.
    Double_t pile = SmoothedPeakIn(proj, kPileLoX * beam, kPileHiX * beam);
    delete proj;

    result.centroids[s] = beam;
    beam_ok[s] = kTRUE;
    if (pile > 0) {
      result.pileups[s] = pile;
      pile_ok[s] = kTRUE;
      result.factors[s] = 1.0 / (pile - beam);
      result.offsets[s] = 1.0 - result.factors[s] * beam;
    } else {
      result.factors[s] = 1.0 / beam;
      result.offsets[s] = 0.0;
    }
    valid_strips++;
    Constants::Detail() << "  strip " << s << " beam=" << Form("%.4f", beam)
                        << " pileup="
                        << (pile > 0 ? Form("%.4f", pile)
                                     : "none (factor only)")
                        << " a.u.  factor=" << Form("%.4f", result.factors[s])
                        << " offset=" << Form("%+.4f", result.offsets[s])
                        << " (n=" << n_entries << ")" << std::endl;
  }
  result.ok = (valid_strips >= 4) ? kTRUE : kFALSE;

  // Diagnostic plot: the pre-alignment totals with the two peaks each strip
  // was aligned on (beam orange, pile-up violet).
  {
    TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(h2, cv);
    TGraph *g_beam_plot = new TGraph(kNStrips);
    TGraph *g_pile_plot = new TGraph(kNStrips);
    Int_t nb = 0, npk = 0;
    for (Int_t s = 0; s < kNStrips; s++) {
      if (beam_ok[s])
        g_beam_plot->SetPoint(nb++, Double_t(s), result.centroids[s]);
      if (pile_ok[s])
        g_pile_plot->SetPoint(npk++, Double_t(s), result.pileups[s]);
    }
    g_beam_plot->Set(nb);
    g_pile_plot->Set(npk);
    if (nb > 0) {
      g_beam_plot->SetMarkerStyle(20);
      g_beam_plot->SetMarkerColor(kOrange);
      g_beam_plot->Draw("P SAME");
    }
    if (npk > 0) {
      g_pile_plot->SetMarkerStyle(21);
      g_pile_plot->SetMarkerColor(kViolet + 2);
      g_pile_plot->Draw("P SAME");
    }
    if (Constants::SavePlots())
      PlottingUtils::SaveFigure(cv, "strip_alignment_check", plot_subdir,
                                PlotSaveOptions::kLINEAR);
    delete cv;
    delete g_beam_plot;
    delete g_pile_plot;
  }

  delete h2;

  return result;
}

void CalibrateBeam::CalibrateBeamOneSubfile(
    const FileSpec &spec, const std::vector<ChannelCal> &chans_template) {
  TString file_label = FileSet::FileLabel(spec);
  TString plot_subdir = "beam_calibration/" + file_label;
  Constants::Detail() << "Beam calibration: " << file_label << std::endl;

  // One gate per strip: strip total vs its neighbour's. A strip whose gate
  // cannot be fitted is excluded from calibration, not the whole subfile.
  BeamFit2D gate[18];
  Int_t n_gates = 0;
  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    for (Int_t s = 1; s <= 17; s++) {
      gate[s] =
          FindBeamGateStrips(spec, GatePartner(s), s, file_label, plot_subdir);
      if (gate[s].ok)
        n_gates++;
      else
        Constants::DetailErr()
            << "  " << file_label << ": strip " << s
            << " beam gate failed; that strip is not calibrated here"
            << std::endl;
    }
    // Strip 0 shares strip 1's gate: there is no strip before it, so the
    // (0, 1) pair is the only one available to either.
    gate[0] = gate[1];
  }
  if (n_gates == 0) {
    std::cerr << "  " << file_label << ": every per-strip beam gate failed"
              << std::endl;
    return;
  }

  std::vector<ChannelCal> chans = chans_template;
  std::vector<std::vector<Float_t>> samples;
  StripPairSamples pairs[18];
  CollectAnchorSamplesOneSubfile(spec, chans, gate, samples, pairs);
  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    ReduceToAnchors(chans, samples, pairs, plot_subdir);
  }

  // Uncalibrated channels silently get gain 0 (drop out of strip totals), so
  // log which/why; grep "[uncalibrated long]" — long side should never starve.
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    const ChannelCal &c = chans[i];
    if (IsCalibrated(c))
      continue;
    TString kind;
    if (c.side == 'C')
      kind = "cathode";
    else if (c.side == 'G')
      kind = "grid";
    else if (c.side == 'S')
      kind = "unsegmented";
    else
      kind = (c.side == LongSide(c.strip)) ? "long" : "short";
    TString why;
    if (c.n_samples < kMinSamples)
      why =
          Form("too few beam samples (%lld < %lld)", c.n_samples, kMinSamples);
    else
      why = Form("bad exp anchor (fit_adc=%.1f)", c.fit_adc);
    Constants::DetailErr() << "  [uncalibrated " << kind << "] " << c.name
                           << " -> gain 0; " << why
                           << " (beam n=" << c.n_samples << ")" << std::endl;
  }

  for (Int_t i = 0; i < Int_t(chans.size()); i++)
    if (IsCalibrated(chans[i]))
      Constants::Detail() << "  " << chans[i].name << " gain=" << Gain(chans[i])
                          << " a.u./ADC" << std::endl;

  StripAlignmentResult align;

  // Write initial calibration tree: per-channel gains. The alignment step
  // needs this on disk so EnergyView can decode events in a.u.
  WriteCalibrationToEvents(spec, chans, &align);

  // Per-strip two-point alignment on the strip total: beam peak to 1.0,
  // two-particle pile-up peak to 2.0 (factor on every end, offset on the
  // long end).
  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    align = FindStripCentroidAlignment(spec, plot_subdir, file_label);
  }
  if (align.ok) {
    WriteCalibrationToEvents(spec, chans, &align);
  }

  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    SaveDynamicRangeOverlay(spec, plot_subdir, file_label);
  }

  // One line per subfile whatever the sample: what was measured and what fell
  // back. The per-strip and per-channel lines above are detail.
  Int_t n_cal = 0, n_uncal = 0, n_ridge = 0, n_aligned = 0;
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    if (IsCalibrated(chans[i]))
      n_cal++;
    else
      n_uncal++;
    if (chans[i].ridge_ratio > 0.0)
      n_ridge++;
  }
  for (Int_t s = 1; s <= 16; s++)
    if (align.ok && align.factors[s] != 1.0)
      n_aligned++;
  // One string, one write, so concurrent workers' lines do not interleave.
  const TString line =
      Form("[calibration] %s: gates %d/17, channels %d calibrated (%d not), "
           "ridge measured on %d/16 strips, alignment on %d/16",
           file_label.Data(), n_gates, n_cal, n_uncal, n_ridge, n_aligned);
  std::cout << line << std::endl;
}

/// Replace each subfile's short-side gain with one built from the run-level
/// median ridge ratio.
///
/// C_short/C_long is a ratio of preamp gains, so it is fixed per channel and
/// does not vary subfile to subfile -- the measured scatter is ~1%, well below
/// the ~7% spread between strips. But a single subfile often cannot fit the
/// ridge on the low-occupancy short ends (on 87Rb the even strips fit in only
/// ~45% of subfiles, and 37Cl's upstream strips almost never do), so
/// per-subfile fitting leaves a large fraction of strips on a fallback.
///
/// Aggregating fixes that without re-reading any event data: take the median
/// ratio and the median short-end offset (both in units of C_long) per strip
/// over the subfiles that did measure the ridge, then rebuild every
/// subfile's short gain as 1/(ratio * C_long) and its short offset as
/// offset * C_long, using that subfile's own C_long so per-subfile gain drift
/// is preserved. Strips that no subfile could fit keep whatever the
/// per-subfile fallback gave them.
void CalibrateBeam::AggregateRidgeRatiosForRun(
    Int_t run, const std::vector<FileSpec> &specs) {
  std::vector<Double_t> per_strip[18], per_strip_off[18];
  for (Int_t k = 0; k < Int_t(specs.size()); k++) {
    TString sub = FileSet::EventsName(specs[k]) + ".root";
    TFile *cf = IO::OpenForReading(sub);
    if (!cf || cf->IsZombie()) {
      if (cf)
        delete cf;
      continue;
    }
    TTree *t = static_cast<TTree *>(cf->Get("calibration"));
    if (!t || t->GetEntries() < 1 || !t->GetBranch("RidgeRatio") ||
        !t->GetBranch("RidgeOffset")) {
      cf->Close();
      delete cf;
      continue;
    }
    Float_t rr[18] = {0}, ro[18] = {0};
    t->SetBranchAddress("RidgeRatio", rr);
    t->SetBranchAddress("RidgeOffset", ro);
    t->GetEntry(0);
    for (Int_t s = 1; s <= 16; s++)
      if (rr[s] > 0) {
        per_strip[s].push_back(Double_t(rr[s]));
        per_strip_off[s].push_back(Double_t(ro[s]));
      }
    cf->Close();
    delete cf;
  }

  // The ratio and the offset (both in units of C_long) travel together: a
  // subfile that measured the ridge measured both.
  Double_t med[18] = {0}, med_off[18] = {0};
  std::cout << "Run " << run << " ridge-ratio aggregation:" << std::endl;
  for (Int_t s = 1; s <= 16; s++) {
    std::vector<Double_t> &v = per_strip[s];
    std::vector<Double_t> &vo = per_strip_off[s];
    if (v.size() < 3) {
      std::cerr << "  strip " << s << ": only " << v.size()
                << " subfiles measured the ridge; leaving per-subfile gains"
                << std::endl;
      continue;
    }
    std::sort(v.begin(), v.end());
    std::sort(vo.begin(), vo.end());
    Int_t m = Int_t(v.size());
    // The IQR reads the sorted vectors, so it comes before the medians,
    // which only partially order them.
    Double_t lo = v[m / 4], hi = v[(3 * m) / 4];
    Double_t lo_off = vo[m / 4], hi_off = vo[(3 * m) / 4];
    med[s] = Median(v);
    med_off[s] = Median(vo);
    std::cout << "  strip " << s << " ratio=" << Form("%.4f", med[s])
              << "  IQR " << Form("%.4f", lo) << "-" << Form("%.4f", hi)
              << "  offset=" << Form("%.4f", med_off[s]) << " x C_long  IQR "
              << Form("%.4f", lo_off) << "-" << Form("%.4f", hi_off)
              << "  from " << m << " subfiles" << std::endl;
  }

  Int_t n_rewritten = 0;
  for (Int_t k = 0; k < Int_t(specs.size()); k++) {
    TString sub = FileSet::EventsName(specs[k]) + ".root";
    TFile *cf = IO::OpenForWriting(sub, "UPDATE");
    if (!cf || cf->IsZombie()) {
      if (cf)
        delete cf;
      continue;
    }
    TTree *t = static_cast<TTree *>(cf->Get("calibration"));
    if (!t || t->GetEntries() < 1 || !t->GetBranch("LongAnchor") ||
        !t->GetBranch("OffsetLeft")) {
      cf->Close();
      delete cf;
      continue;
    }
    Float_t gl[16] = {0}, gr[16] = {0}, ol[16] = {0}, orr[16] = {0},
            la[18] = {0}, rr[18] = {0}, ro[18] = {0};
    t->SetBranchAddress("GainLeft", gl);
    t->SetBranchAddress("GainRight", gr);
    t->SetBranchAddress("OffsetLeft", ol);
    t->SetBranchAddress("OffsetRight", orr);
    t->SetBranchAddress("LongAnchor", la);
    t->SetBranchAddress("RidgeRatio", rr);
    t->SetBranchAddress("RidgeOffset", ro);
    t->GetEntry(0);
    Bool_t changed = kFALSE;
    for (Int_t s = 1; s <= 16; s++) {
      if (med[s] <= 0 || la[s] <= 0)
        continue;
      // Short anchor and offset from the run medians and this subfile's own
      // long anchor, so per-subfile gain drift is preserved.
      Double_t anchor = med[s] * Double_t(la[s]);
      if (anchor <= 0)
        continue;
      const Float_t g = Float_t(1.0 / anchor);
      const Float_t off = Float_t(med_off[s] * Double_t(la[s]));
      if (LongSide(s) == 'L') {
        gr[s - 1] = g;
        orr[s - 1] = off;
      } else {
        gl[s - 1] = g;
        ol[s - 1] = off;
      }
      rr[s] = Float_t(med[s]);
      ro[s] = Float_t(med_off[s]);
      changed = kTRUE;
    }
    if (changed) {
      TTree *nt = t->CloneTree(0);
      nt->Fill();
      cf->cd();
      nt->Write("calibration", TObject::kOverwrite);
      n_rewritten++;
    }
    cf->Close();
    delete cf;
  }
  std::cout << "  rewrote short gains and offsets in " << n_rewritten
            << " subfiles" << std::endl;
}

void CalibrateBeam::Run(const TString &file_label) {
  const TString project_root = Paths::DatasetDir();
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                Paths::ResultsDir() + "/plots",
                                Paths::ResultsDir() + "/root_files");
  gROOT->SetBatch(kTRUE);

  std::vector<FileSpec> specs;
  if (file_label.IsNull()) {
    specs = FileSet::BuildProcessedFileSpecs();
    if (specs.empty()) {
      std::cerr << "No file specs from FileSet::BuildProcessedFileSpecs()"
                << std::endl;
      return;
    }
  } else {
    FileSpec s = FileSet::ResolveFileSpec(file_label);
    if (s.run < 0) {
      std::cerr << "Could not resolve file label '" << file_label << "'"
                << std::endl;
      return;
    }
    specs.push_back(s);
  }

  std::vector<ChannelCal> chans = CalibrateBeam::BuildChannels();

  Int_t n_specs = Int_t(specs.size());

  std::set<Int_t> runs;
  for (Int_t k = 0; k < n_specs; k++)
    runs.insert(specs[k].run);

  Int_t n_workers =
      TMath::Min(Int_t(std::thread::hardware_concurrency()), n_specs);
  n_workers = TMath::Min(n_workers, Constants::cfg.MAX_FUSED_WORKERS);
  if (n_workers < 1)
    n_workers = 1;
  std::cout << "calibrate-beam: " << n_specs << " subfiles on " << n_workers
            << " workers" << std::endl;

  std::queue<Int_t> work;
  for (Int_t k = 0; k < n_specs; k++)
    work.push(k);
  std::mutex work_mutex;

  std::vector<std::thread> workers;
  for (Int_t w = 0; w < n_workers; w++) {
    workers.emplace_back([&]() {
      while (true) {
        Int_t k;
        {
          std::lock_guard<std::mutex> lk(work_mutex);
          if (work.empty())
            return;
          k = work.front();
          work.pop();
        }
        // The first files of the epoch draw (a single named file always).
        Constants::SetPlotsThisFile(Constants::InPlotSample(k));
        CalibrateBeamOneSubfile(specs[k], chans);
      }
    });
  }
  for (Int_t w = 0; w < Int_t(workers.size()); w++)
    workers[w].join();
}
