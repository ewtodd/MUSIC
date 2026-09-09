#ifndef BEAM_FIT_2D_HPP
#define BEAM_FIT_2D_HPP

#include <Rtypes.h>
#include <TH2F.h>
#include <TString.h>
#include <algorithm>
#include <cmath>

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
   * @brief Whether a point lies inside the fitted beam ellipse.
   * @param b  Fitted spot. Must have `ok == kTRUE` to be meaningful.
   * @param x  Point x.
   * @param y  Point y.
   * @param nx Ellipse half-width in units of `sigma_x`.
   * @param ny Ellipse half-height in units of `sigma_y`.
   * @return `kTRUE` if inside the correlation-aware ellipse.
   */
  static Bool_t InEllipseXY(const BeamFit2D &b, Double_t x, Double_t y,
                            Double_t nx, Double_t ny);

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
};

#endif
