#ifndef STRIP_SUM_SCATTER_HPP
#define STRIP_SUM_SCATTER_HPP

#include "BeamFit2D.hpp"
#include "Constants.hpp"
#include "EventsSummary.hpp"
#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Normalization.hpp"
#include "PlottingUtils.hpp"
#include "RemixSim.hpp"
#include <Rtypes.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TCutG.h>
#include <TEllipse.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TKey.h>
#include <TLegend.h>
#include <TMath.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

struct GateSpec {
  Int_t sx;
  Int_t sy;
};

// Beam classification ellipses: entrance (s0,s1 or s1,s2, per config's
// PURE_BEAM_GATE) AND exit (s16,s17 or s15,s16). An event is "pure beam"
// only if it passes both ellipses.
namespace ScatterBuildRange {
// Every strip-sum scatter is built over this fixed window, so the display
// windows can be retuned without refilling. The axes get different ceilings
// because they sum different numbers of strips: x runs over X_LO..X_HI (16
// strips, so it reaches ~20), while y covers only POST_TRIGGER_SUM_STRIPS
// after the trigger (6, so it never approaches 40). Giving y its own, lower
// ceiling doubles its resolution at the same bin count -- which matters,
// because y is the axis the reaction populations separate along and its
// display windows are much narrower than x's.
const Double_t kXMin = 0.0;
const Double_t kXMax = 40.0;
const Double_t kYMin = 0.0;
const Double_t kYMax = 20.0;
} // namespace ScatterBuildRange

struct BeamEllipses {
  BeamFit2D s0_s1;
  BeamFit2D s1_s2;
  BeamFit2D s16_s17;
  BeamFit2D s15_s16;
  Bool_t ok;
  Bool_t use_s15_s16;
};

struct TraceEvt {
  Float_t total[18];
  Float_t total_adc[18]; // raw (un-normalized) ADC sum per strip
  // Calibrated a.u. of each half of a split strip, kept separately and
  // independently of IGNORE_SHORT_STRIPS (which is a decode-time switch and
  // zeroes one half of EnergyView::left/right). Reconstructed from the raw
  // ADC and the per-channel gains, so one calibration can be rendered under
  // either decode.
  Float_t long_au[18];
  Float_t short_au[18];
  UInt_t reac_mask;
  Bool_t beam_flat;
  Int_t both_mult; // # split strips (1-16) with BOTH ends above threshold
  // Timestamp of the Grid hit that seeded the event, from the events tree.
  // Unique per event and independent of how runs are split into files, so a
  // cached event can be joined back to its source record. 0 when the source
  // predates the SeedTs branch.
  ULong64_t seed_ts;
};

// Per-run results, so both fill phases can run one worker per run and merge
// afterwards in run order -- which keeps the threaded result identical to the
// sequential one.
struct SingleRunFitResult {
  BeamEllipses pure_beam;
  std::vector<BeamFit2D> series_gates;
  Bool_t ok;
  SingleRunFitResult() : ok(kFALSE) {}
};

struct SingleRunFillResult {
  // Private clones, one per reaction strip; owned by the caller after merge.
  std::vector<TH2F *> scatters;
  std::vector<TraceEvt> reservoir;
  Long64_t gated;
  Long64_t seen;
  // Events surviving every cut applied before reaction tagging: the beam
  // gates, pileup/noise/offbeam rejection and the both-ends multiplicity.
  // Every tagged event passed exactly this selection, so it is the
  // denominator in which those efficiencies cancel.
  Long64_t normed;
  // Per-strip denominator, indexed by ReacIndex: beam particles that reached
  // that strip under exactly the conditions a reaction there would have had to
  // satisfy -- every strip fired, and strips 1..reac-1 beam-like. Those
  // conditions are strip-dependent (strip 5 must have four strips upstream
  // look like beam where strip 3 needs two), so a single denominator would
  // leave that efficiency uncancelled and bias the excitation function.
  std::vector<Long64_t> normed_at;
  // Events tagged at each reaction strip, indexed by ReacIndex.
  std::vector<Long64_t> tagged;
  SingleRunFillResult() : gated(0), seen(0), normed(0) {}
};

struct SimPop {
  TString file;
  TString label;
};

class StripSumScatter {
public:
  StripSumScatter();
  ~StripSumScatter();

  // Main entry point called from main_strip_sum_scatter.cpp
  void Run();

  // The scatter cache this configuration reads and writes; compute-regions
  // reads the same file, so the name lives here and nowhere else.
  static TString CacheName();

  // The scatter plane's coordinates of one event's normed strip totals, for
  // reaction strip `reac`: x the sum over X_LO..X_HI, y the sum over the
  // post-trigger window. Every consumer of a region cut goes through this, so
  // the fill, the overlay and the cross section can never disagree about
  // where an event sits.
  static void PlaneXY(const Double_t *total, Int_t reac, Double_t &x,
                      Double_t &y);
  // The post-trigger window summed onto y for reaction strip `reac`:
  // YLoOf..YHiOf inclusive, per the POST_TRIGGER_SUM_STRIPS /
  // POST_WINDOW_LAST_STRIP / POST_WINDOW_STRIPS rule.
  static Int_t YLoOf(Int_t reac);
  static Int_t YHiOf(Int_t reac);
  // The beam's noise, measured by MeasureBeamNoise and stamped in the cache:
  // JumpSigma is the sigma of total[s] - total[s-1] (jump_sigma_s<N>),
  // StripSigma the sigma of total[s] itself (strip_sigma_s<N>). Set once
  // before any tagging and read-only after, so the worker threads share them
  // safely.
  static Double_t JumpSigma(Int_t strip);
  static Double_t StripSigma(Int_t strip);
  static void SetJumpSigma(const Double_t *sigma);
  static void SetStripSigma(const Double_t *sigma);
  // Minimum jump for a tag at `reac`: REAC_JUMP_NSIGMA times JumpSigma(reac).
  static Double_t JumpMin(Int_t reac);

private:
  static Double_t s_jumpSigma[18];
  static Double_t s_stripSigma[18];
  std::map<Int_t, TH2F *> m_scatter;
  std::vector<TraceEvt> m_reservoir;
  // Normalization counts, merged over every run and persisted in the cache so
  // a cross section can be taken from it without a second pass over the data.
  Long64_t m_nSeen;
  Long64_t m_nNormed;
  std::vector<Long64_t> m_normedAt;
  std::vector<Long64_t> m_tagged;
  Double_t m_yLo[64];
  Double_t m_yHi[64];

  const TString kSimCacheName = "StripSumScatter_simcache.root";

  static Int_t ReacIndex(Int_t reac);

  Bool_t TryLoadCache(const TString &cacheName, const TString &fingerprint);
  void WriteCache(const TString &cacheName, const TString &fingerprint);

  // Fresh, empty scatters for the configured plane and build range.
  void AllocateScatters();
  // Fill the scatters from the reservoir alone. The reservoir keeps every
  // tagged event, so a plane or axes change needs only this and not a
  // 25-minute pass over the events files.
  void ReprojectFromReservoir();
  void FillScatters(const std::vector<Int_t> &runOrder,
                    std::map<Int_t, TChain *> &chains);

  void PlotScatters();

  void InteractiveOverlay(Int_t reac);

  static void EnableEventBranches(TChain *chain);
  static Bool_t AllStripsFired(const EnergyView &ev);
  static Bool_t IsPureBeam(const EnergyView &ev, const BeamEllipses &be);
  // Sigma-clipped width of each strip-to-strip difference and of each strip's
  // deposit over a capped sample of `chain`, after the cheap pre-tag cuts.
  // False when too few events survive to measure them.
  static Bool_t MeasureBeamNoise(TChain *chain, Double_t *jump_sigma,
                                 Double_t *strip_sigma);
  // Strips 1..reac-1 within BEAM_UPSTREAM_NSIGMA of the beam, or the
  // requirement is off. Shared by the tag and its per-strip denominator.
  static Bool_t BeamUpstreamOf(const EnergyView &ev, Int_t reac);
  static Bool_t IsPileup(const EnergyView &ev);
  static Bool_t IsNoise(const EnergyView &ev);
  static Bool_t IsOffbeam(const EnergyView &ev);
  static Bool_t IsParityAsymmetric(const EnergyView &ev);
  static Double_t SumRange(const Double_t *total, Int_t lo, Int_t hi);
  static std::vector<GateSpec> ActiveGates();

  static Bool_t PassesGate(const BeamFit2D &gate, const EnergyView &ev,
                           Int_t sx, Int_t sy);

  static BeamFit2D FindBeamGate(TChain *chain, Int_t sx, Int_t sy,
                                const std::vector<GateSpec> &prior_specs,
                                const std::vector<BeamFit2D> &prior_gates,
                                const TString &tag, const TString &subdir);

  static void DrawTraceSet(const std::vector<TGraph *> &traces, Int_t color);
  void DrawAltDecodeRegionTraces(Int_t reac, TCutG *cutAn, TCutG *cutAa);
  static TGraph *TraceFromTotal(const Float_t *total);
  static void DrawRegionTraces(const TString &save_name, const TString &subdir,
                               const std::vector<TGraph *> &beam,
                               const std::vector<TGraph *> &aa,
                               const std::vector<TGraph *> &an, Double_t y_min,
                               Double_t y_max, const char *y_title);

  static void DrawRegionMeanTraces(const TString &save_name,
                                   const TString &subdir,
                                   const std::vector<TGraph *> &beam,
                                   const std::vector<TGraph *> &aa,
                                   const std::vector<TGraph *> &an,
                                   Double_t y_min, Double_t y_max,
                                   const char *y_title);

  static void TraceYRange(const std::vector<TGraph *> &beam,
                          const std::vector<TGraph *> &aa,
                          const std::vector<TGraph *> &an, Double_t &y_min,
                          Double_t &y_max);
  static SingleRunFitResult
  FitRunGates(Int_t run, TChain *chain,
              const std::vector<GateSpec> &activeGates);
  static SingleRunFillResult FillRunScatters(
      Int_t run, TChain *chain, const std::vector<GateSpec> &activeGates,
      const std::vector<BeamFit2D> &runGates, const BeamEllipses &runBeam);
  static TCutG *PromptCut(TCanvas *c, const char *name, const char *label);
  static void SaveRegionCuts(Int_t reac, TCutG *cut_an, TCutG *cut_aa);
  static TCutG *LoadRegionCut(const char *name, Int_t reac);

  static void SmoothTrace(const Double_t *in, Double_t *out, Int_t width);

  // Savitzky-Golay smoothing: 3rd-degree polynomial, half-window of 2
  // (5-point convolution). Uses standard SG coefficients [-3,12,17,12,-3]/35.
  // At edges, the window shrinks and coefficients are renormalised.
  static void SavitzkyGolay(const Double_t *in, Double_t *out);

  // CFD-style trigger finder: locate the first strip whose beam-subtracted
  // signal (td[s]-1) exceeds both a fraction of the trace peak and a multiple
  // of the beam sigma. Returns the strip index, or -1 if no trigger fires.
  static Int_t FindTrigger(const Double_t *td, const Double_t *base,
                           Double_t beam_sigma);

  // Build a TGraph from Savitzky-Golay-smoothed per-strip totals. Input is
  // the raw normed array; smoothing is applied internally before graph build.
  static TGraph *SmoothedTraceFromTotal(const Float_t *total);

  void ClusterVarHists(Int_t reac, TCutG *cut_aa, TCutG *cut_an,
                       const TString &subdir);

  static TString
  SimFingerprint(const std::vector<RemixSim::SimFileSpec> &specs);
  static TString BuildFingerprint(const std::vector<Int_t> &run_order,
                                  std::map<Int_t, TChain *> &chains);
  static void YBounds(Double_t *y_lo, Double_t *y_hi);
  static TString PrettyLabel(const TString &tag);

public:
  // Shared with tag-efficiency, which pushes bootstrapped traces through the
  // same tag as the scatter.
  static Bool_t PassesReaction(const EnergyView &ev, Int_t reac);

private:
  static Bool_t SimBeamGains(Double_t *gain);
  static void SimTotal(const Float_t *left, const Float_t *right,
                       const Double_t *gain, Double_t *total);
  static TGraph *SimPopScatter(const TString &file, Int_t reac,
                               const Double_t *gain, Long64_t max_points);
  static std::vector<TGraph *>
  SimPopTraces(const TString &file, const Double_t *gain, Long64_t max_traces);
  void SimTraceOverlay();
  static Bool_t LoadSimCache(const TString &fp,
                             std::map<Int_t, std::vector<TGraph *>> &by_strip);
  static void
  WriteSimCache(const TString &fp,
                const std::map<Int_t, std::vector<TGraph *>> &by_strip);
  void SimOverlay();
};

#endif
