#ifndef REGION_CUTS_HPP
#define REGION_CUTS_HPP

#include <TCutG.h>
#include <TH2F.h>
#include <TString.h>

/**
 * @brief On-disk home of the per-strip (a,n) and (a,a') region cuts.
 *
 * Shared by `strip-sum-scatter`, which draws or loads them, and
 * `compute-regions`, which fits them — so the two tools can never disagree
 * about where a cut lives.
 *
 * One cut per file, written `RECREATE`. Saving one strip's region therefore
 * cannot touch another's, and never touches the scatter cache: cuts are applied
 * at read time, so a new cut is picked up on the next read with no cache
 * rebuild.
 */
namespace RegionCutStore {
/// @brief Directory holding the region cut files.
TString Dir();
/// @brief Path to one cut's file.
/// @param name Cut name, e.g. the (a,n) or (a,a') label.
/// @param reac Reaction strip index.
TString Path(const char *name, Int_t reac);
/**
 * @brief Write a strip's two region cuts.
 * @param reac   Reaction strip index.
 * @param cut_an The (a,n) region.
 * @param cut_aa The (a,a') region.
 * @param n_an_assigned Events the fit attributed to the reaction component.
 *        When >= 0 this is stored beside the (a,n) cut, for the cross section
 *        to prefer over a geometric count of the region. Hand-drawn cuts have
 *        none, hence the `-1` default.
 */
void Save(Int_t reac, TCutG *cut_an, TCutG *cut_aa,
          Double_t n_an_assigned = -1.0);
/**
 * @brief Load a cut, from either storage generation.
 *
 * Tries the per-cut file first, then the pre-split combined `RegionCuts.root`,
 * so cuts drawn before the split keep working until they are next replaced.
 *
 * @param name Cut name.
 * @param reac Reaction strip index.
 * @return The cut, or null if neither file holds one. The caller owns it.
 */
TCutG *Load(const char *name, Int_t reac);
/**
 * @brief Load only a hand-drawn cut, never a fitted one.
 *
 * Sources are interactive `strip-sum-scatter` and the pre-split
 * `RegionCuts.root`. Used where a measurement must rest on a human's judgement
 * rather than on the mixture fit — the tag efficiency bootstrap, for instance.
 *
 * @param name Cut name.
 * @param reac Reaction strip index.
 * @return The cut, or null. The caller owns it.
 */
TCutG *LoadDrawn(const char *name, Int_t reac);
/// @brief The attributed count stored with a cut.
/// @param name Cut name.
/// @param reac Reaction strip index.
/// @return The count, or `-1` when the cut carries none.
Double_t LoadAssigned(const char *name, Int_t reac);
} // namespace RegionCutStore

/**
 * @brief One bivariate Gaussian component of a strip's scatter.
 *
 * Fields are in ROOT's `bigaus` parameter order, so the struct can be handed
 * straight to a `TF2` built from that formula.
 */
struct Gauss2D {
  Double_t amp = 0.0; ///< Amplitude.
  Double_t mx = 0.0;  ///< Mean in x.
  Double_t sx = 0.0;  ///< Width in x.
  Double_t my = 0.0;  ///< Mean in y.
  Double_t sy = 0.0;  ///< Width in y.
  Double_t rho = 0.0; ///< x-y correlation, in [-1, 1].
};

/**
 * @brief Result of fitting the two-component mixture to one strip's scatter.
 *
 * Check #ok before reading anything; #why says what went wrong when it is
 * false.
 */
struct RegionFit {
  Gauss2D beam; ///< The beam-like population, pinned from its own core.
  Gauss2D reac; ///< The reaction population above the ridge.
  Double_t x_lo = 0.0, x_hi = 0.0, y_lo = 0.0, y_hi = 0.0; ///< Fit window.
  /// Events the mixture attributes to each component.
  ///
  /// Where the island sits on the beam's tail, this is the count a person
  /// drawing a tight region gets, and what a geometric cut cannot give: the
  /// region's ellipse also holds whatever tail lies under it, and the beam's
  /// real tail is far heavier than its Gaussian, so no model subtraction
  /// recovers it.
  Double_t n_beam = 0.0, n_reac = 0.0;
  /// `kFALSE` for a beam-only fit, as used by the ridge-band region mode.
  /// #reac is then a copy of #beam with zero amplitude, and nothing draws or
  /// counts it.
  Bool_t has_reac = kTRUE;
  Bool_t ok = kFALSE; ///< Whether the fit succeeded. Check this first.
  TString why;        ///< Why the fit failed; set only when #ok is false.
};

namespace RegionCutStore {
/**
 * @brief Store the fitted components beside the (a,n) cut.
 *
 * Kept so a later step can weight events by the mixture's own posterior rather
 * than by a geometric cut.
 *
 * @param reac Reaction strip index.
 * @param fit  Fit to store.
 */
void SaveFit(Int_t reac, const RegionFit &fit);

/**
 * @brief Read back a stored fit.
 * @param reac Reaction strip index.
 * @param[out] fit Set on success; untouched otherwise.
 * @return `kFALSE` when the store holds no fit for that strip.
 */
Bool_t LoadFit(Int_t reac, RegionFit &fit);
} // namespace RegionCutStore

/**
 * @brief Finding regions by fitting a two-component mixture to the scatter.
 *
 * This is ApJ 983:142's "two Gaussian peaks" carried into two dimensions. The
 * two populations separate along the beam ridge's diagonal, which 1-D
 * projections wash out, and the reaction island is a few times 1e-4 of the
 * scatter.
 *
 * The fit is staged accordingly: the beam-like component is pinned from its own
 * core, the reaction component is seeded at the first local maximum above the
 * ridge, and the mixture is fitted **only above the core**, so that the island
 * rather than the 1e5-count peak decides the second component.
 */
namespace RegionCutFinder {
/**
 * @brief Fit the two-component mixture to one strip's scatter.
 * @param scatter Scatter histogram for the strip. Must not be null.
 * @param reac    Reaction strip index.
 * @param x_lo    Fit window, lower x.
 * @param x_hi    Fit window, upper x.
 * @param y_lo    Fit window, lower y.
 * @param y_hi    Fit window, upper y.
 * @return The fit; check RegionFit::ok, and RegionFit::why when it is false.
 */
RegionFit FitMixture(TH2F *scatter, Int_t reac, Double_t x_lo, Double_t x_hi,
                     Double_t y_lo, Double_t y_hi);

/**
 * @brief Closed polygon of a component's Mahalanobis contour.
 * @param name   Name for the returned cut.
 * @param g      Component to contour.
 * @param nsigma Contour level, in Mahalanobis sigma.
 * @param npts   Polygon vertices; more gives a smoother ellipse.
 * @return The cut. The caller owns it.
 */
TCutG *EllipseCut(const char *name, const Gauss2D &g, Double_t nsigma,
                  Int_t npts = 64);

/**
 * @brief Fit the beam component alone, with no reaction component.
 *
 * Pinned from its core exactly as FitMixture() does, but returns a fit with
 * `has_reac == kFALSE`. For the ridge-band region mode, where the (a,n)
 * population is not a compact island and so cannot be modelled as a second
 * Gaussian.
 *
 * @param scatter Scatter histogram. Must not be null.
 * @param reac    Reaction strip index.
 * @param x_lo    Fit window, lower x.
 * @param x_hi    Fit window, upper x.
 * @param y_lo    Fit window, lower y.
 * @param y_hi    Fit window, upper y.
 * @return The beam-only fit.
 */
RegionFit FitBeam(TH2F *scatter, Int_t reac, Double_t x_lo, Double_t x_hi,
                  Double_t y_lo, Double_t y_hi);

/**
 * @brief How far a point sits above the beam ridge, in conditional sigma.
 *
 * The residual of y about the ridge line, in units of the beam's conditional
 * width `sigma_y * sqrt(1 - rho^2)`.
 *
 * @param beam Fitted beam component.
 * @param x    Point x.
 * @param y    Point y.
 * @return The height in sigma; negative below the ridge.
 */
Double_t AboveRidge(const Gauss2D &beam, Double_t x, Double_t y);

/**
 * @brief Closed polygon of a band above the beam ridge.
 *
 * The (a,n) region in `AN_REGION_RIDGE_BAND` mode. Spans the window in x and is
 * clipped to the window in y.
 *
 * @param name    Name for the returned cut.
 * @param beam    Fitted beam component defining the ridge.
 * @param nsig_lo Lower band edge, in conditional sigma above the ridge.
 * @param nsig_hi Upper band edge.
 * @param x_lo    Window, lower x.
 * @param x_hi    Window, upper x.
 * @param y_lo    Window, lower y.
 * @param y_hi    Window, upper y.
 * @return The cut. The caller owns it.
 */
TCutG *RidgeBandCut(const char *name, const Gauss2D &beam, Double_t nsig_lo,
                    Double_t nsig_hi, Double_t x_lo, Double_t x_hi,
                    Double_t y_lo, Double_t y_hi);

/**
 * @brief Count scatter events inside a cut.
 * @param scatter Scatter histogram.
 * @param cut     Region to count within.
 * @param fit     Fit supplying the window; only bins inside it are considered.
 * @return Summed contents of the bins whose centre lies inside @p cut.
 */
Double_t CountInside(TH2F *scatter, TCutG *cut, const RegionFit &fit);

/**
 * @brief Save the diagnostic figures for a strip's regions.
 *
 * The scatter with both regions drawn, then each projection with the fitted
 * mixture's projection overlaid — so a region that landed somewhere silly is
 * visible rather than merely wrong.
 *
 * @param scatter Scatter histogram.
 * @param reac    Reaction strip index.
 * @param fit     Fit to overlay.
 * @param an      The (a,n) region.
 * @param aa      The (a,a') region.
 * @param subdir  Plot subdirectory under the plots base.
 */
void SaveFigures(TH2F *scatter, Int_t reac, const RegionFit &fit, TCutG *an,
                 TCutG *aa, const TString &subdir);
} // namespace RegionCutFinder

#endif
