/// compute-regions: fit each reaction strip's (a,n) and (a,a') regions from the
/// scatter cache and save them as cuts, exactly where hand-drawn ones go.
///
/// In the mixture mode this is a read-only pass over
/// StripSumScatter_cache.root: it never fills, never touches the cache, and
/// needs no DISPLAY; strip-sum-scatter then loads the saved cuts like any
/// drawn ones, so re-run it after any refill. In the all-tagged mode there is
/// nothing to draw or fit, so this is the only binary needed: it measures the
/// beam and loads or builds the cache itself (StripSumScatter::Prepare),
/// saves the whole build window as each strip's region, and draws every
/// strip's tagged traces against the beam under plots/compute_regions.
#include "Constants.hpp"
#include "InitUtils.hpp"
#include "Paths.hpp"
#include "RegionCuts.hpp"
#include "StripSumScatter.hpp"
#include <TFile.h>
#include <TH2F.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>

int main() {
  Constants::ActivateAnalysisEpoch();
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                Paths::ResultsDir() + "/plots",
                                Paths::ResultsDir() + "/root_files");
  gROOT->SetBatch(kTRUE);
  const StripSumScatterConfig &C = Constants::cfg.STRIP_SUM_SCATTER_CONFIG;
  const Bool_t all_tagged =
      C.AN_REGION_MODE == StripSumScatterConfig::AN_REGION_ALL_TAGGED;

  // All-tagged: build or load the cache here, so strip-sum-scatter need not
  // run; the scatter object also keeps the trace reservoir for the overlays.
  StripSumScatter scatter;
  if (all_tagged && !scatter.Prepare()) {
    std::cerr << "compute-regions: could not build the scatter cache"
              << std::endl;
    return 1;
  }

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

  if (C.AN_REGION_MODE == StripSumScatterConfig::AN_REGION_ALL_TAGGED)
    std::cout << "compute-regions: every tagged event is (a,n), no fit; reac "
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
    const StripSumScatterConfig::AnRegionMode mode = C.AN_REGION_MODE;
    if (mode == StripSumScatterConfig::AN_REGION_ALL_TAGGED) {
      // Every tagged event is the reaction: region = whole build window,
      // attributed count = scatter's content; no (a,a') -- unlinks stale.
      TCutG *an = new TCutG("region_an", 5);
      an->SetPoint(0, ScatterBuildRange::kXMin, ScatterBuildRange::kYMin);
      an->SetPoint(1, ScatterBuildRange::kXMax, ScatterBuildRange::kYMin);
      an->SetPoint(2, ScatterBuildRange::kXMax, ScatterBuildRange::kYMax);
      an->SetPoint(3, ScatterBuildRange::kXMin, ScatterBuildRange::kYMax);
      an->SetPoint(4, ScatterBuildRange::kXMin, ScatterBuildRange::kYMin);
      const Double_t nAn = h->Integral();
      std::cout << Form("  [region] reac %2d: all tagged events are (a,n): "
                        "%.0f",
                        reac, nAn)
                << std::endl;
      RegionCutStore::Save(reac, an, nullptr, nAn);
      gSystem->Unlink(RegionCutStore::Path("region_aa", reac));
      delete an;
      nOk++;
      continue;
    }
    // Same window the scatter is displayed in.
    Double_t y_lo = C.Y_DISPLAY_MIN, y_hi = C.Y_DISPLAY_MAX;
    std::map<Int_t, std::pair<Double_t, Double_t>>::const_iterator it =
        C.Y_DISPLAY_RANGE.find(reac);
    if (it != C.Y_DISPLAY_RANGE.end()) {
      y_lo = it->second.first;
      y_hi = it->second.second;
    }
    RegionFit fit = RegionCutFinder::FitMixture(h, reac, C.X_DISPLAY_MIN,
                                                C.X_DISPLAY_MAX, y_lo, y_hi);
    if (!fit.ok) {
      std::cout << "  [region] reac " << reac << ": no region -- " << fit.why
                << std::endl;
      continue;
    }
    TCutG *an =
        RegionCutFinder::EllipseCut("region_an", fit.reac, C.AN_REGION_NSIGMA);
    TCutG *aa =
        RegionCutFinder::EllipseCut("region_aa", fit.beam, C.AA_REGION_NSIGMA);
    const Double_t nAn = RegionCutFinder::CountInside(h, an, fit);
    const Double_t nAa = RegionCutFinder::CountInside(h, aa, fit);
    {
      // Two counts: fit-attributed vs in-region; agree where the island
      // is off the ridge; diverge where the region also covers beam tail.
      std::cout
          << Form("  [region] reac %2d: beam (%.3f, %.3f) s(%.3f, %.3f) "
                  "rho %+.2f | (a,n) (%.3f, %.3f) s(%.3f, %.3f) rho %+.2f "
                  "| attributed %.0f, in region %.0f; (a,a') in region %.0f",
                  reac, fit.beam.mx, fit.beam.my, fit.beam.sx, fit.beam.sy,
                  fit.beam.rho, fit.reac.mx, fit.reac.my, fit.reac.sx,
                  fit.reac.sy, fit.reac.rho, fit.n_reac, nAn, nAa)
          << std::endl;
    }
    RegionCutStore::Save(reac, an, aa, fit.n_reac);
    RegionCutStore::SaveFit(reac, fit);
    RegionCutFinder::SaveFigures(h, reac, fit, an, aa, "compute_regions");
    delete an;
    delete aa;
    nOk++;
  }
  f.Close();
  std::cout << "compute-regions: regions for " << nOk << " of " << nTried
            << " reaction strips" << std::endl;
  if (all_tagged) {
    std::cout << "compute-regions: tagged traces per reaction strip against "
                 "the beam"
              << std::endl;
    scatter.DrawAllTaggedTraces("compute_regions");
  }
  return nOk > 0 ? 0 : 1;
}
