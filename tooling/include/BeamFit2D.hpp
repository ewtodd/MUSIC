#ifndef BEAM_FIT_2D_HPP
#define BEAM_FIT_2D_HPP

#include <Rtypes.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

/**
 * @brief A fitted 2-D Gaussian beam spot.
 *
 * Check #ok before reading anything else; a failed fit leaves the parameters at
 * zero, which is indistinguishable from a genuine result otherwise.
 */
struct BeamFit2D {
  Double_t amp = 0;                  ///< Peak amplitude.
  Double_t mu_x = 0, mu_y = 0;       ///< Centroid.
  Double_t sigma_x = 0, sigma_y = 0; ///< Widths along each axis.
  Double_t rho = 0;                  ///< Correlation coefficient, in [-1, 1].
  Bool_t ok = kFALSE;                ///< Whether the fit converged.
};

/**
 * @brief Second moments of a 2-D distribution over a bin range.
 *
 * The cheap alternative to a fit: computed directly from the histogram
 * contents, and used to seed BeamFit2D or when a fit is not warranted.
 */
struct Moments2D {
  Double_t mu_x = 0, mu_y = 0;       ///< Weighted centroid.
  Double_t sigma_x = 0, sigma_y = 0; ///< Weighted RMS widths.
  Double_t rho = 0;                  ///< Correlation coefficient, in [-1, 1].
  Double_t weight = 0;               ///< Total weight included; zero means the
                                     ///< range held nothing above threshold.
};

/// @brief Helpers for the 2-D beam-spot description.
class BeamFitUtils {
public:
  /**
   * @brief Whether a point lies inside the fitted ellipse at @p n sigma: the
   *        correlation-aware Mahalanobis chi2 below n^2, one level, since a
   *        density contour has one. A Gaussian keeps 1 - exp(-n^2/2) inside.
   * @param b Fitted spot; `ok` must be true.
   * @param n Level, in sigma.
   */
  static Bool_t InEllipseXY(const BeamFit2D &b, Double_t x, Double_t y,
                            Double_t n);

  /**
   * @brief Second moments of a histogram region above a threshold.
   * @param h      Histogram to measure. Must not be null.
   * @param lo_bx  First x bin, inclusive.
   * @param hi_bx  Last x bin, inclusive.
   * @param lo_by  First y bin, inclusive.
   * @param hi_by  Last y bin, inclusive.
   * @param thresh Bin contents at or below this are ignored, which keeps
   *               background out of the widths.
   * @param bw_x   x bin width, for converting bin indices to physical units.
   * @param bw_y   y bin width.
   * @return The moments; `weight == 0` when nothing cleared @p thresh.
   */
  static Moments2D ComputeMoments(TH2F *h, Int_t lo_bx, Int_t hi_bx,
                                  Int_t lo_by, Int_t hi_by, Double_t thresh,
                                  Double_t bw_x, Double_t bw_y);

  /**
   * @brief Sigma-clipped centroid and covariance of a point sample, the 2-D
   *        counterpart of the strip-width clip: each pass keeps the points
   *        within @p clip (Mahalanobis) of the current estimate and takes
   *        their mean and covariance, scaled up by the variance a 2-D Gaussian
   *        loses outside @p clip. The widths are then a Gaussian's sigma for
   *        a Gaussian spot, and the clipped RMS of anything else; no binning
   *        or peak-fraction threshold enters.
   * @param pts   The (x, y) sample.
   * @param seed  Starting estimate, e.g. ComputeMoments over the core. Its
   *              widths may be too narrow; the passes grow them.
   * @param clip  Mahalanobis radius kept per pass.
   * @return The converged moments; the seed when fewer than 10 points are
   *         kept.
   */
  static Moments2D
  ClippedMoments(const std::vector<std::pair<Float_t, Float_t>> &pts,
                 const Moments2D &seed, Double_t clip = 3.0);

  /**
   * @brief Fit a correlated 2-D Gaussian on a flat pedestal to @p h: ROOT's
   *        bigaus written out, Poisson likelihood, Minuit limits from @p seed.
   * @param h             Finely binned spot reaching past the fit window; a
   *                      couple of bins per sigma would fit the binning.
   * @param seed          Sets the window and the limits, e.g. ClippedMoments.
   * @param window_nsigma Half-width fitted, in the seed's sigma; narrow keeps a
   *                      shoulder out of the width (+4 percent at 2, +10 at 4).
   * @param chi2_ndf      Reduced chi2 over the whole of @p h, so what the
   *                      window excluded is what it notices; -1 if no fit ran.
   * @return `ok == kFALSE` when the fit failed or ran away: keep the seed.
   */
  static BeamFit2D FitSpot(TH2F *h, const Moments2D &seed,
                           Double_t window_nsigma = 2.0,
                           Double_t *chi2_ndf = nullptr);

  /**
   * @brief The 1-D counterpart: a Gaussian on a flat pedestal, fitted over
   *        @p window_nsigma of the seed and judged over the whole of @p h.
   * @param fix_mu Hold the centre at @p seed_mu, for a distribution with more
   *               than one peak.
   * @return `kFALSE` when the fit failed or ran away; @p mu and @p sigma are
   *         then untouched.
   */
  static Bool_t FitPeak(TH1F *h, Double_t seed_mu, Double_t seed_sigma,
                        Double_t window_nsigma, Bool_t fix_mu, Double_t &mu,
                        Double_t &sigma, Double_t *chi2_ndf = nullptr);
};

#endif
