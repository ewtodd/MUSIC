#ifndef STRIP_SUM_SCATTER_HPP
#define STRIP_SUM_SCATTER_HPP

#include "BeamFit2D.hpp"
#include "Constants.hpp"
#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Normalization.hpp"
#include "PlottingUtils.hpp"
#include "RemixSim.hpp"
#include "TraceCreator.hpp"
#include <Rtypes.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TCutG.h>
#include <TEllipse.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

class StripSumScatter {
public:
  static void Run();

private:
  static constexpr Int_t kSmoothHiStrip =
      12; // smoothness checked through this (absolute) strip (upstream k<13)
  static constexpr Double_t kBeamFlatTol =
      1.0; // |E - NORM| tol for the pre-reaction strips (0..reac-1)
  static constexpr Double_t kReacJumpMin =
      0.1; // reaction dE jump / excess over NORM, lower bound
  static constexpr Double_t kReacJumpMax =
      2.0; // reaction dE jump / excess over NORM, upper bound
  static constexpr Double_t kSmoothMaxStep =
      1.2; // max |dE| step between adjacent post-reaction strips
  static constexpr Double_t kStrip17Max =
      1; // Strip17 upper bound (upstream S17R[1])
  static constexpr Double_t kLongStripMinValue = 0.85;

  // Candidate reaction strips: every scatter is computed in one pass over the
  // data, one TH2F per reaction strip in [MIN, MAX].
  static constexpr Int_t kReacMin = Constants::REACTION_STRIP_MIN;
  static constexpr Int_t kReacMax = Constants::REACTION_STRIP_MAX;
  static constexpr Int_t kNReac = kReacMax - kReacMin + 1;

  // Max sampled per-strip traces drawn per region (beam / (a,a') / (a,n)).
  static constexpr Int_t kTracesPerRegion = 40;
  // Cap on beam-flat events kept in the trace reservoir (only ~40 are drawn).
  static constexpr Int_t kBeamReservoirCap = 400;

  // Sum windows: x = all long anodes (1-16); y = the downstream reaction window
  // reac+1 .. reac+6, clamped to the last strip (17).
  static constexpr Int_t kXLo = 1;
  static constexpr Int_t kXHi = 16;

  static constexpr Int_t kGateStripX = 1;
  static constexpr Int_t kGateStripY = 2;
  static constexpr Double_t kGateNSigmaX = 2.0;
  static constexpr Double_t kGateNSigmaY = 2.0;

  static constexpr Double_t kGateMin = 0.0;
  static constexpr Double_t kGateMax = 3.0;
  static constexpr Int_t kGateBins = 240;
  static constexpr Int_t kSeedHalfBins = 40;
  static constexpr Double_t kSeedFrac = 0.30;

  inline static const Double_t kXMin = Constants::STRIP_SUM_XMIN;
  inline static const Double_t kXMax = Constants::STRIP_SUM_XMAX;
  static constexpr Int_t kXBins = Constants::STRIP_SUM_XBINS;
  static constexpr Int_t kYBins = Constants::STRIP_SUM_YBINS;

  // Sampled-pass budget for the beam-gate fit and the y-bound auto-range. Both
  // only need to characterize a population, not visit every event.
  static constexpr Long64_t kSampleMaxPoints = 2000000;

  // Cache file (under root_files): all per-strip scatters + the trace
  // reservoir, stamped with a fingerprint of the cut constants and input entry
  // counts. Its name is gate-config dependent (see CacheName()).
  static constexpr const char *kSimCacheName = "StripSumScatter_simcache.root";

  // One beam gate: an ellipse fit on strip sx vs strip sy. The strip1-vs-strip2
  // gate (kGateStripX/Y) is always active; strips 3/4 and 5/6 are opt-in via
  // Constants::STRIP_SUM_GATE_S3_S4 / _S5_S6.
  struct GateSpec {
    Int_t sx;
    Int_t sy;
  };

  struct TraceEvt {
    Float_t total[18];
    Float_t total_adc[18]; // raw (un-normalized) ADC sum per strip, same event
    UInt_t reac_mask;
    Bool_t beam_flat;
  };

  struct SimPop {
    TString file;
    TString label;
  };

  static Int_t ReacIndex(Int_t reac) { return reac - kReacMin; }
  static Int_t YLoOf(Int_t reac) { return reac + 1; }
  static Int_t YHiOf(Int_t reac) { return TMath::Min(reac + 6, 17); }

  static void EnableEventBranches(TChain *chain);
  static Bool_t AllStripsFired(const EnergyView &ev);
  static Bool_t PassesReaction(const EnergyView &ev, Int_t reac);
  static Bool_t IsBeamFlat(const EnergyView &ev);
  static Bool_t IsPileup(const EnergyView &ev);
  static Bool_t IsNoise(const EnergyView &ev);
  static Double_t SumRange(const Double_t *total, Int_t lo, Int_t hi);
  static std::vector<GateSpec> ActiveGates();

  // Per-gate-config cache file: StripSumScatter_cache[_g34][_g56].root, so each
  // on/off combination keeps its own cache and toggling never clobbers another.
  static TString CacheName();
  static Bool_t PassesGate(const BeamFit2D &gate, const EnergyView &ev,
                           Int_t sx, Int_t sy);

  static BeamFit2D FindBeamGate(TChain *chain, Int_t sx, Int_t sy,
                                const std::vector<GateSpec> &prior_specs,
                                const std::vector<BeamFit2D> &prior_gates,
                                const TString &tag, const TString &subdir);

  static void DrawTraceSet(const std::vector<TGraph *> &traces, Int_t color);
  static TGraph *TraceFromTotal(const Float_t *total);
  static void DrawRegionTraces(const TString &save_name, const TString &subdir,
                               const std::vector<TGraph *> &beam,
                               const std::vector<TGraph *> &aa,
                               const std::vector<TGraph *> &an,
                               Double_t y_min = Constants::STRIP_DE_MIN_NORMED,
                               Double_t y_max = Constants::STRIP_DE_MAX_NORMED,
                               const char *y_title = "#DeltaE [a.u.]");

  static void TraceYRange(const std::vector<TGraph *> &beam,
                          const std::vector<TGraph *> &aa,
                          const std::vector<TGraph *> &an, Double_t &y_min,
                          Double_t &y_max);
  static TCutG *PromptCut(TCanvas *c, const char *name, const char *label);

  static TString
  SimFingerprint(const std::vector<RemixSim::SimFileSpec> &specs);
  static TString BuildFingerprint(const std::vector<Int_t> &run_order,
                                  std::map<Int_t, TChain *> &chains);
  static void YBounds(Double_t *y_lo, Double_t *y_hi);
  static TString PrettyLabel(const TString &tag);
  static Bool_t SimBeamGains(Double_t *gain);
  static void SimTotal(const Float_t *left, const Float_t *right,
                       const Double_t *gain, Double_t *total);
  static TGraph *SimPopScatter(const TString &file, Int_t reac,
                               const Double_t *gain, Long64_t max_points);
  static std::vector<TGraph *>
  SimPopTraces(const TString &file, const Double_t *gain, Long64_t max_traces);
  static void SimTraceOverlay();
  static Bool_t LoadSimCache(const TString &fp,
                             std::map<Int_t, std::vector<TGraph *>> &by_strip);
  static void
  WriteSimCache(const TString &fp,
                const std::map<Int_t, std::vector<TGraph *>> &by_strip);
  static void SimOverlay(const std::map<Int_t, TH2F *> &scatter);
};

#endif
