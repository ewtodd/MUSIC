// compute-regions: fit each reaction strip's (a,n) and (a,a') regions from the
// scatter cache and save them as cuts, exactly where hand-drawn ones go.
//
// A read-only pass over StripSumScatter_cache.root: it never fills, never
// touches the cache, and needs no DISPLAY. strip-sum-scatter then loads the
// saved cuts like any drawn ones. Re-run it after any refill.
#include "Constants.hpp"
#include "InitUtils.hpp"
#include "Paths.hpp"
#include "RegionCuts.hpp"
#include "StripSumScatter.hpp"
#include <TFile.h>
#include <TH2F.h>
#include <TNamed.h>
#include <TROOT.h>
#include <iostream>

int main() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                Paths::ResultsDir() + "/plots",
                                Paths::ResultsDir() + "/root_files");
  gROOT->SetBatch(kTRUE);
  const StripSumScatterConfig &C = Constants::cfg.STRIP_SUM_SCATTER_CONFIG;

  // Same resolution as StripSumScatter::TryLoadCache: the bare name lives
  // under the dataset's root_files directory.
  TString cache =
      IO::GetRootFilesBaseDir() + "/" + StripSumScatter::CacheName();
  TFile f(cache, "READ");
  if (f.IsZombie()) {
    std::cerr << "compute-regions: no scatter cache at " << cache
              << "; run strip-sum-scatter first" << std::endl;
    return 1;
  }
  if (TNamed *fp = static_cast<TNamed *>(f.Get("fingerprint"))) {
    std::cout << "compute-regions: cache " << cache << std::endl;
    std::cout << "  fingerprint " << fp->GetTitle() << std::endl;
  }

  const Bool_t band =
      C.AN_REGION_MODE == StripSumScatterConfig::AN_REGION_RIDGE_BAND;
  if (band)
    std::cout << "compute-regions: beam component per strip; (a,n) region = "
                 "the band "
              << C.AN_RIDGE_NSIGMA_LO << ".." << C.AN_RIDGE_NSIGMA_HI
              << " conditional sigma above the ridge, (a,a') the beam's "
              << C.AA_REGION_NSIGMA << " sigma ellipse; reac "
              << C.REACTION_STRIP_MIN << ".." << C.REACTION_STRIP_MAX
              << std::endl;
  else
    std::cout << "compute-regions: bivariate Gaussian mixture per strip; "
                 "(a,n) region at "
              << C.AN_REGION_NSIGMA << " sigma, (a,a') at "
              << C.AA_REGION_NSIGMA << " sigma; reac " << C.REACTION_STRIP_MIN
              << ".." << C.REACTION_STRIP_MAX << std::endl;

  Int_t nOk = 0, nTried = 0;
  for (Int_t reac = C.REACTION_STRIP_MIN; reac <= C.REACTION_STRIP_MAX;
       reac++) {
    TH2F *h = static_cast<TH2F *>(f.Get(Form("scatter_r%d", reac)));
    if (!h)
      continue;
    nTried++;
    // Same window the scatter is displayed in.
    Double_t y_lo = C.Y_DISPLAY_MIN, y_hi = C.Y_DISPLAY_MAX;
    std::map<Int_t, std::pair<Double_t, Double_t>>::const_iterator it =
        C.Y_DISPLAY_RANGE.find(reac);
    if (it != C.Y_DISPLAY_RANGE.end()) {
      y_lo = it->second.first;
      y_hi = it->second.second;
    }
    RegionFit fit =
        band ? RegionCutFinder::FitBeam(h, reac, C.X_DISPLAY_MIN,
                                        C.X_DISPLAY_MAX, y_lo, y_hi)
             : RegionCutFinder::FitMixture(h, reac, C.X_DISPLAY_MIN,
                                           C.X_DISPLAY_MAX, y_lo, y_hi);
    if (!fit.ok) {
      std::cout << "  [region] reac " << reac << ": no region -- " << fit.why
                << std::endl;
      continue;
    }
    TCutG *an = band ? RegionCutFinder::RidgeBandCut(
                           "region_an", fit.beam, C.AN_RIDGE_NSIGMA_LO,
                           C.AN_RIDGE_NSIGMA_HI, C.X_DISPLAY_MIN,
                           C.X_DISPLAY_MAX, y_lo, y_hi)
                     : RegionCutFinder::EllipseCut("region_an", fit.reac,
                                                   C.AN_REGION_NSIGMA);
    TCutG *aa =
        RegionCutFinder::EllipseCut("region_aa", fit.beam, C.AA_REGION_NSIGMA);
    const Double_t nAn = RegionCutFinder::CountInside(h, an, fit);
    const Double_t nAa = RegionCutFinder::CountInside(h, aa, fit);
    if (band) {
      // The band's count is geometric by construction; no attributed count
      // is stored, so cross-section takes what lies inside the region.
      fit.n_reac = nAn;
      std::cout << Form("  [region] reac %2d: beam (%.3f, %.3f) s(%.3f, %.3f) "
                        "rho %+.2f | (a,n) band %.1f..%.1f sigma above the "
                        "ridge: %.0f in region; (a,a') in region %.0f",
                        reac, fit.beam.mx, fit.beam.my, fit.beam.sx,
                        fit.beam.sy, fit.beam.rho, C.AN_RIDGE_NSIGMA_LO,
                        C.AN_RIDGE_NSIGMA_HI, nAn, nAa)
                << std::endl;
    } else {
      // Two counts: what the fit attributes to the reaction, and what lies
      // inside the drawn region. They agree where the island is well off the
      // ridge and diverge where the region also covers beam tail.
      std::cout
          << Form("  [region] reac %2d: beam (%.3f, %.3f) s(%.3f, %.3f) "
                  "rho %+.2f | (a,n) (%.3f, %.3f) s(%.3f, %.3f) rho %+.2f "
                  "| attributed %.0f, in region %.0f; (a,a') in region %.0f",
                  reac, fit.beam.mx, fit.beam.my, fit.beam.sx, fit.beam.sy,
                  fit.beam.rho, fit.reac.mx, fit.reac.my, fit.reac.sx,
                  fit.reac.sy, fit.reac.rho, fit.n_reac, nAn, nAa)
          << std::endl;
    }
    RegionCutStore::Save(reac, an, aa, band ? -1.0 : fit.n_reac);
    RegionCutStore::SaveFit(reac, fit);
    RegionCutFinder::SaveFigures(h, reac, fit, an, aa, "compute_regions");
    delete an;
    delete aa;
    nOk++;
  }
  f.Close();
  std::cout << "compute-regions: regions for " << nOk << " of " << nTried
            << " reaction strips" << std::endl;
  return nOk > 0 ? 0 : 1;
}
