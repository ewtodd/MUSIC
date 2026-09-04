#ifndef REGION_CUTS_HPP
#define REGION_CUTS_HPP

#include <TCutG.h>
#include <TH2F.h>
#include <TString.h>

// Where each reaction strip's (a,n) and (a,a') region cuts live, and how they
// are named. Shared by strip-sum-scatter (which draws or loads them) and
// compute-regions (which fits them), so the two tools can never disagree about
// the file. One cut per file, written RECREATE, so saving one strip's region
// can never touch another's and never touches the scatter cache: cuts are
// applied at read time, so a new cut is picked up on the next read.
namespace RegionCutStore {
TString Dir();
TString Path(const char *name, Int_t reac);
// n_an_assigned, when >= 0, is stored beside the (a,n) cut: the events the
// fit attributed to the reaction component, for cross-section to prefer over
// a geometric count of the region. Drawn cuts have none.
void Save(Int_t reac, TCutG *cut_an, TCutG *cut_aa,
          Double_t n_an_assigned = -1.0);
// The per-cut file first, then the pre-split combined RegionCuts.root, so cuts
// drawn before the split keep working until they are next replaced.
TCutG *Load(const char *name, Int_t reac);
// Only a cut drawn by hand (interactive strip-sum-scatter or the pre-split
// RegionCuts.root), never one compute-regions fitted.
TCutG *LoadDrawn(const char *name, Int_t reac);
// The attributed count saved with a cut, or -1 if there is none.
Double_t LoadAssigned(const char *name, Int_t reac);
} // namespace RegionCutStore

// One bivariate Gaussian component of a strip's scatter: amplitude, means,
// widths and the x-y correlation. ROOT's "bigaus" parameter order.
struct Gauss2D {
  Double_t amp = 0.0, mx = 0.0, sx = 0.0, my = 0.0, sy = 0.0, rho = 0.0;
};

struct RegionFit {
  Gauss2D beam; // the beam-like population, pinned from its own core
  Gauss2D reac; // the reaction population above the ridge
  Double_t x_lo = 0.0, x_hi = 0.0, y_lo = 0.0, y_hi = 0.0; // fit window
  // Events the mixture attributes to each component. Where the island sits
  // on the beam's tail this is the count a person drawing a tight region
  // gets, and what a geometric cut cannot: the region's ellipse also holds
  // whatever tail lies under it, and the beam's real tail is far heavier than
  // its Gaussian, so no model subtraction recovers it.
  Double_t n_beam = 0.0, n_reac = 0.0;
  Bool_t ok = kFALSE;
  TString why; // set when !ok
};

namespace RegionCutStore {
// The fitted components, stored beside the (a,n) cut so a later step can
// weight events by the mixture's own posterior instead of a geometric cut.
void SaveFit(Int_t reac, const RegionFit &fit);
Bool_t LoadFit(Int_t reac, RegionFit &fit);
} // namespace RegionCutStore

// Regions by a two-component bivariate Gaussian mixture fitted to the scatter
// itself. This is ApJ 983:142's "two Gaussian peaks" carried into 2D: the two
// populations are separated along the beam ridge's diagonal, which 1D
// projections wash out, and the reaction island is a few 1e-4 of the scatter.
// The fit is staged accordingly: the beam-like component is pinned from its
// own core, the reaction component is seeded at the first local maximum above
// the ridge, and the mixture is fitted only above the core so the island, not
// the 1e5-count peak, decides the second component.
namespace RegionCutFinder {
RegionFit FitMixture(TH2F *scatter, Int_t reac, Double_t x_lo, Double_t x_hi,
                     Double_t y_lo, Double_t y_hi);

// Closed polygon of the component's nsigma Mahalanobis contour.
TCutG *EllipseCut(const char *name, const Gauss2D &g, Double_t nsigma,
                  Int_t npts = 64);

// Scatter events whose bin centre lies inside the cut (window only).
Double_t CountInside(TH2F *scatter, TCutG *cut, const RegionFit &fit);

// One figure per plot, through PlottingUtils like every other figure: the
// scatter with both regions drawn, then each projection with the fitted
// mixture's projection overlaid so a region that landed somewhere silly is
// visible.
void SaveFigures(TH2F *scatter, Int_t reac, const RegionFit &fit, TCutG *an,
                 TCutG *aa, const TString &subdir);
} // namespace RegionCutFinder

#endif
