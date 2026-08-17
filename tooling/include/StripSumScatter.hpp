#ifndef STRIP_SUM_SCATTER_HPP
#define STRIP_SUM_SCATTER_HPP

#include "BeamFit2D.hpp"
#include "Constants.hpp"
#include "EventsSummary.hpp"
#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "InteractiveEditorX11Guard.hpp"
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
#include <TVirtualX.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

struct GateSpec {
  Int_t sx;
  Int_t sy;
};

// Beam classification ellipses: entrance (s0,s1 or s1,s2, per config's
// PURE_BEAM_GATE) AND exit (s16,s17 or s15,s16). An event is "pure beam"
// only if it passes both ellipses.
struct BeamEllipses {
  BeamFit2D s0_s1;
  BeamFit2D s1_s2;
  BeamFit2D s16_s17;
  BeamFit2D s15_s16;
  Bool_t ok = kFALSE;
  Bool_t use_s15_s16;
};

struct TraceEvt {
  Float_t total[18];
  Float_t total_adc[18]; // raw (un-normalized) ADC sum per strip
  UInt_t reac_mask;
  Bool_t beam_flat;
  Int_t both_mult; // # split strips (1-16) with BOTH ends above threshold
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

private:
  std::map<Int_t, TH2F *> m_scatter;
  std::vector<TraceEvt> m_reservoir;
  Double_t m_yLo[64];
  Double_t m_yHi[64];

  static Int_t ReacIndex(Int_t reac);
  static Int_t YLoOf(Int_t reac);
  static Int_t YHiOf(Int_t reac);

  Bool_t TryLoadCache(const TString &cacheName, const TString &fingerprint);
  void WriteCache(const TString &cacheName, const TString &fingerprint);

  void FillScatters(const std::vector<Int_t> &runOrder,
                    std::map<Int_t, TChain *> &chains);

  // Per-run result of the beam-ellipse + series-gate fits (FillScatters
  // phase 1, one entry per run, computed on worker threads).
  struct RunGateFit {
    BeamEllipses beam;
    std::vector<BeamFit2D> gates;
    Bool_t ok = kFALSE;
  };

  // Per-run result of the event fill (FillScatters phase 2, one entry per
  // run): the run's private scatter histograms (merged afterwards) and its
  // trace reservoir slice.
  struct FillRunResult {
    std::vector<TH2F *> scatters;
    std::vector<TraceEvt> reservoir;
    Long64_t seen = 0;       // events read from the chain
    Long64_t gated = 0;      // events passing all filters with a reaction mask
    Long64_t n_rej_gate = 0; // rejected by the series beam gates
    Long64_t n_rej_pileup = 0;    // rejected by IsPileup
    Long64_t n_rej_noise = 0;     // rejected by IsNoise
    Long64_t n_rej_highstrip = 0; // rejected by IsHighStrip
    Long64_t n_rej_offbeam = 0;   // rejected by IsOffbeam
  };

  static RunGateFit FitRunGates(Int_t run, TChain *chain);
  static FillRunResult FillRunScatters(
      Int_t run, TChain *chain, const std::vector<GateSpec> &active_gates,
      const std::vector<BeamFit2D> &run_gates, const BeamEllipses &run_beam,
      const Double_t *y_lo, const Double_t *y_hi);

  void PlotScatters();

  void InteractiveOverlay(Int_t reac);

  static void EnableEventBranches(TChain *chain);
  static Bool_t AllStripsFired(const EnergyView &ev);
  static Bool_t PassesReaction(const EnergyView &ev, Int_t reac);
  static Bool_t IsPureBeam(const EnergyView &ev, const BeamEllipses &be);
  static Bool_t IsPileup(const EnergyView &ev);
  static Bool_t IsNoise(const EnergyView &ev);
  static Bool_t IsHighStrip(const EnergyView &ev);
  static Bool_t IsOffbeam(const EnergyView &ev);
  static Double_t SumRange(const Double_t *total, Int_t lo, Int_t hi);
  static std::vector<GateSpec> ActiveGates();

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
  static TCutG *PromptCut(TCanvas *c, const char *name, const char *label);

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
