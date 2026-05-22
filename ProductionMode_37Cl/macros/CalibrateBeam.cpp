#include "Constants.hpp"
#include "DrawCut.hpp"
#include "FittingUtils.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Normalization.hpp"
#include "Pipeline.hpp"
#include "PlottingUtils.hpp"
#include "RooFitUtils.hpp"
#include <TApplication.h>
#include <TCanvas.h>
#include <TCutG.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TKey.h>
#include <TLegend.h>
#include <TList.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVectorD.h>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <queue>
#include <thread>
#include <toml++/toml.hpp>
#include <vector>

namespace CalibrateBeamNS {

const Int_t kSimBins = 300;
const Int_t kNPeaks = 4; // beam + 3 pile-up multiplicities
const Double_t kFitHalfWindowSigmas = 5.0;
const Int_t kFitMinHalfWindowBins = 15;
const Long64_t kMinEntriesForFit = 200;
const Int_t kMaxChannels = 35;
const Long64_t kMaxEventsPerJob = 200000;

struct PeakFit {
  Double_t mean = 0;
  Double_t mean_err = 0;
  Double_t sigma = 0;
  Double_t sigma_err = 0;
  Double_t amp = 0;
  Double_t reduced_chi2 = -1;
  Bool_t ok = kFALSE;
};

struct ChannelCal {
  TString name; // "Strip0", "Strip17", "L1".."R16", "Cathode"
  Char_t side;  // 'S' (unsplit guard) / 'L' / 'R' / 'C'
  Int_t strip;  // 0..17 for strip channels; -1 for cathode

  TH1D *sim_hist_mev = nullptr;

  PeakFit sim_fit;
  PeakFit peak_fits[kNPeaks];
  Long64_t peak_n[kNPeaks] = {0, 0, 0, 0};

  Double_t quad_a = 0, quad_b = 0, quad_c = 0;
  Double_t quad_chi2_ndf = -1;
  Bool_t quad_ok = kFALSE;

  // Multiplier on peak_sigma_ADC so this channel's k=1 MeV sigma equals the
  // long-side k=1 MeV mean (via the local pol2 derivative). 1.0 for long
  // sides by construction; != 1.0 for Strip0/Strip17/Cathode.
  Double_t sigma_match_scale = 1.0;
};

inline TString SimBeamPath(const TString &project_root) {
  TString p = project_root + "/../Remix_37Cl/root_files/traces_37Cl_beam.root";
  return TString(gSystem->ExpandPathName(p.Data()));
}

inline TString CalSummarySubpath(Int_t run) {
  return Form("Calibration_Run%d.root", run);
}

inline Char_t LongSide(Int_t strip) { return (strip % 2 == 1) ? 'L' : 'R'; }
inline Char_t ShortSide(Int_t strip) { return (strip % 2 == 1) ? 'R' : 'L'; }

inline Bool_t IsShortSide(const ChannelCal &c) {
  if (c.side != 'L' && c.side != 'R')
    return kFALSE;
  if (c.strip < 1 || c.strip > 16)
    return kFALSE;
  return c.side == ShortSide(c.strip);
}

inline TString ChannelHistName(const ChannelCal &c) {
  if (c.side == 'C')
    return "h2_totalE_vs_cathode";
  if (IsShortSide(c))
    return Form("h2_totalE_vs_short_s%d", c.strip);
  return Form("h2_totalE_vs_stripE_s%d", c.strip);
}

inline TString PeakCutName(const ChannelCal &c, Int_t k) {
  return Form("peak_%s_k%d", c.name.Data(), k);
}

std::vector<ChannelCal> BuildChannels() {
  std::vector<ChannelCal> chans;
  ChannelCal c{};
  if (!Constants::IGNORE_STRIP_0) {
    c.name = "Strip0";
    c.side = 'S';
    c.strip = 0;
    chans.push_back(c);
  }
  if (!Constants::IGNORE_STRIP_17) {
    c.name = "Strip17";
    c.side = 'S';
    c.strip = 17;
    chans.push_back(c);
  }
  for (Int_t s = 1; s <= 16; s++) {
    Char_t long_side = LongSide(s);
    c.name = Form("%c%d", long_side, s);
    c.side = long_side;
    c.strip = s;
    chans.push_back(c);
    if (!Constants::IGNORE_SHORT_STRIPS) {
      Char_t short_side = ShortSide(s);
      c.name = Form("%c%d", short_side, s);
      c.side = short_side;
      c.strip = s;
      chans.push_back(c);
    }
  }
  c.name = "Cathode";
  c.side = 'C';
  c.strip = -1;
  chans.push_back(c);
  return chans;
}

// Pull a TH2F written by EventBuilder. The builder calls
// `c->Write(h2->GetName(), ...)` so file->Get(name) returns either the
// histogram (older files) or a canvas wrapping it.
TH2F *ExtractH2FromCanvasInFile(TFile *file, const TString &name) {
  TObject *obj = file->Get(name);
  if (!obj) {
    std::cerr << "  not found in " << file->GetName() << ": " << name
              << std::endl;
    return nullptr;
  }
  if (obj->InheritsFrom(TH2F::Class()))
    return static_cast<TH2F *>(obj);
  if (!obj->InheritsFrom(TCanvas::Class())) {
    std::cerr << "  object " << name << " is neither TH2F nor TCanvas"
              << std::endl;
    return nullptr;
  }
  TCanvas *canv = static_cast<TCanvas *>(obj);
  TList *prims = canv->GetListOfPrimitives();
  for (Int_t i = 0; i < prims->GetSize(); i++) {
    TObject *p = prims->At(i);
    if (p && p->InheritsFrom(TH2F::Class()))
      return static_cast<TH2F *>(p);
  }
  std::cerr << "  no TH2F primitive in canvas " << name << std::endl;
  return nullptr;
}

Bool_t LoadSimHistograms(const TString &sim_path,
                         std::vector<ChannelCal> &chans,
                         const TString &run_label) {
  TFile *f = TFile::Open(sim_path);
  if (!f || f->IsZombie()) {
    std::cerr << "Could not open sim beam file: " << sim_path << std::endl;
    if (f)
      f->Close();
    return kFALSE;
  }
  TTree *t = static_cast<TTree *>(f->Get("events_MeV"));
  if (!t) {
    std::cerr << "Sim file has no 'events_MeV' tree: " << sim_path << std::endl;
    f->Close();
    return kFALSE;
  }

  // InitUtils sets TH1::AddDirectory(kFALSE) for MT safety, but TTree::Draw
  // needs the new histogram registered in gDirectory so we can fetch it by
  // name. Flip it back on for this load and restore after.
  Bool_t saved_add_dir = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kTRUE);

  for (Int_t k = 0; k < Int_t(chans.size()); k++) {
    ChannelCal &c = chans[k];
    TString expr;
    TString cut;
    if (c.side == 'S') {
      expr = Form("TotaldE[%d]", c.strip);
      cut = Form("TotaldE[%d]>0", c.strip);
    } else if (c.side == 'L') {
      expr = Form("LeftdE[%d]", c.strip);
      cut = Form("LeftdE[%d]>0", c.strip);
    } else if (c.side == 'R') {
      expr = Form("RightdE[%d]", c.strip);
      cut = Form("RightdE[%d]>0", c.strip);
    } else if (c.side == 'C') {
      expr = "Cathode";
      cut = "Cathode>0";
    } else {
      ::Fatal("CalibrateBeamNS::LoadSimHistograms",
              "Unknown channel side '%c' for '%s'", c.side, c.name.Data());
    }

    TString probe_name = "h_simprobe_" + PlottingUtils::GetRandomName();
    t->Draw(expr + ">>" + probe_name + "(200,0,0)", cut, "goff");
    TH1F *probe = static_cast<TH1F *>(gDirectory->Get(probe_name));
    if (!probe || probe->GetEntries() == 0) {
      std::cerr << "  sim Draw produced no histogram for " << c.name
                << std::endl;
      delete probe;
      continue;
    }
    Double_t hi = probe->GetXaxis()->GetXmax() * 1.5;
    delete probe;

    TString hname = "sim_" + run_label + "_" + c.name;
    TH1D *h = new TH1D(hname,
                       Form("Sim %s (%s);#DeltaE [MeV];Counts", c.name.Data(),
                            run_label.Data()),
                       kSimBins, 0.0, hi);
    t->Project(hname, expr, cut);
    h->SetDirectory(nullptr);
    c.sim_hist_mev = h;
  }
  TH1::AddDirectory(saved_add_dir);
  f->Close();
  return kTRUE;
}

Bool_t EstimateFitRange(TH1 *h, Double_t &lo, Double_t &hi) {
  if (h->GetEntries() < kMinEntriesForFit)
    return kFALSE;
  Int_t bmax = h->GetMaximumBin();
  Double_t xmax = h->GetBinCenter(bmax);
  Double_t ymax = h->GetBinContent(bmax);
  Double_t bin_w = h->GetBinWidth(bmax);
  Double_t half = ymax * 0.5;
  Int_t bL = bmax;
  while (bL > 1 && h->GetBinContent(bL) > half)
    bL--;
  Int_t bR = bmax;
  while (bR < h->GetNbinsX() && h->GetBinContent(bR) > half)
    bR++;
  Double_t fwhm = h->GetBinCenter(bR) - h->GetBinCenter(bL);
  Double_t sigma = std::max(fwhm / 2.355, bin_w * 2.0);
  Double_t halfwin =
      std::max(kFitHalfWindowSigmas * sigma, kFitMinHalfWindowBins * bin_w);
  lo = std::max(h->GetXaxis()->GetXmin(), xmax - halfwin);
  hi = std::min(h->GetXaxis()->GetXmax(), xmax + halfwin);
  return kTRUE;
}

PeakFit PackFitResult(const FitResult &r) {
  PeakFit p;
  if (r.valid && !r.peaks.empty()) {
    p.mean = r.peaks[0].mu;
    p.mean_err = r.peaks[0].mu_error;
    p.sigma = r.peaks[0].sigma;
    p.sigma_err = r.peaks[0].sigma_error;
    p.amp = r.peaks[0].gaus_amplitude;
    p.reduced_chi2 = r.reduced_chi2;
    p.ok = kTRUE;
  }
  return p;
}

PeakFit FitGaussian(TH1 *h, Bool_t flat_bg, const TString &peak_name,
                    const TString &input_name) {
  PeakFit p;
  Double_t lo = 0, hi = 0;
  if (!EstimateFitRange(h, lo, hi))
    return p;
  FittingUtils fitter(h, lo, hi,
                      /*use_flat_background=*/flat_bg,
                      /*use_step=*/kFALSE,
                      /*use_low_exp_tail=*/kFALSE,
                      /*use_low_lin_tail=*/kFALSE,
                      /*use_high_exp_tail=*/kFALSE);
  return PackFitResult(fitter.FitSinglePeak(input_name, peak_name));
}

// Per-cut fit cache. Avoids re-scanning every subfile for cuts whose
// Gaussian fit hasn't changed. Keyed by cut name, stored as a TVectorD
// inside `peak_fits/` of the cal-summary file. Invalidated when a cut is
// (re)drawn in DrawAllPeakCuts.
void SaveCachedPeakFit(TFile *cal_summary, const TString &cut_name,
                       const PeakFit &pf, Long64_t n_events) {
  TDirectory *cwd = gDirectory;
  TDirectory *d = cal_summary->GetDirectory("peak_fits");
  if (!d)
    d = cal_summary->mkdir("peak_fits");
  d->cd();
  TVectorD v(8);
  v[0] = pf.mean;
  v[1] = pf.mean_err;
  v[2] = pf.sigma;
  v[3] = pf.sigma_err;
  v[4] = pf.amp;
  v[5] = pf.reduced_chi2;
  v[6] = pf.ok ? 1.0 : 0.0;
  v[7] = Double_t(n_events);
  v.Write(cut_name, TObject::kOverwrite);
  // Flush the dir's key index so an interrupted run can still read fits back.
  d->SaveSelf(kTRUE);
  if (cwd)
    cwd->cd();
}

Bool_t LoadCachedPeakFit(TFile *cal_summary, const TString &cut_name,
                         PeakFit &pf, Long64_t &n_events) {
  TDirectory *d = cal_summary->GetDirectory("peak_fits");
  if (!d)
    return kFALSE;
  TVectorD *v = static_cast<TVectorD *>(d->Get(cut_name));
  if (!v || v->GetNoElements() < 8) {
    delete v;
    return kFALSE;
  }
  pf.mean = (*v)[0];
  pf.mean_err = (*v)[1];
  pf.sigma = (*v)[2];
  pf.sigma_err = (*v)[3];
  pf.amp = (*v)[4];
  pf.reduced_chi2 = (*v)[5];
  pf.ok = ((*v)[6] > 0.5);
  n_events = Long64_t((*v)[7]);
  delete v;
  return kTRUE;
}

void InvalidateCachedPeakFit(TFile *cal_summary, const TString &cut_name) {
  TDirectory *d = cal_summary->GetDirectory("peak_fits");
  if (!d)
    return;
  d->Delete(cut_name + ";*");
  d->SaveSelf(kTRUE);
}

// Aggregate per-channel TH2s from each subfile's embedded canvases into
// Calibration_Run{N}.root and CLOSE cleanly. Returns kTRUE on success.
// The file is closed before returning so every subsequent phase opens
// it fresh; a Ctrl-C between phases (or between cut draws) then always
// lands on a fully-closed, fully-recoverable file.
Bool_t BuildOrOpenRunSummaryFile(Int_t run, const std::vector<FileSpec> &specs,
                                 const std::vector<ChannelCal> &chans) {
  TString cal_path = CalSummarySubpath(run);
  TFile *cal_summary = IO::OpenForWriting(cal_path, "UPDATE");
  if (!cal_summary || cal_summary->IsZombie()) {
    std::cerr << "Cannot open " << cal_path << std::endl;
    if (cal_summary)
      delete cal_summary;
    return kFALSE;
  }

  Bool_t all_present = kTRUE;
  for (Int_t k = 0; k < Int_t(chans.size()); k++) {
    TString name = ChannelHistName(chans[k]);
    if (!cal_summary->FindKey(name)) {
      all_present = kFALSE;
      break;
    }
  }
  if (all_present) {
    std::cout << "  cal summary already has aggregated h2s; skipping rebuild"
              << std::endl;
    if (!cal_summary->GetDirectory("cuts"))
      cal_summary->mkdir("cuts");
    cal_summary->Close();
    delete cal_summary;
    return kTRUE;
  }

  std::vector<TH2F *> agg(chans.size(), nullptr);
  Bool_t initialised = kFALSE;
  for (Int_t s = 0; s < Int_t(specs.size()); s++) {
    TString events_subpath = EventsName(specs[s]) + ".root";
    TFile *sf = IO::OpenForReading(events_subpath);
    if (!sf || sf->IsZombie()) {
      std::cerr << "  cannot open subfile " << events_subpath << std::endl;
      if (sf)
        sf->Close();
      continue;
    }
    for (Int_t k = 0; k < Int_t(chans.size()); k++) {
      TString name = ChannelHistName(chans[k]);
      TH2F *h = ExtractH2FromCanvasInFile(sf, name);
      if (!h)
        continue;
      if (!agg[k]) {
        agg[k] = static_cast<TH2F *>(h->Clone(name));
        agg[k]->SetDirectory(nullptr);
      } else {
        agg[k]->Add(h);
      }
    }
    initialised = kTRUE;
    sf->Close();
    delete sf;
  }
  if (!initialised) {
    std::cerr << "Run " << run << ": no subfile histograms could be read"
              << std::endl;
    cal_summary->Close();
    delete cal_summary;
    return kFALSE;
  }

  cal_summary->cd();
  for (Int_t k = 0; k < Int_t(chans.size()); k++) {
    if (agg[k])
      agg[k]->Write(ChannelHistName(chans[k]), TObject::kOverwrite);
    delete agg[k];
  }
  if (!cal_summary->GetDirectory("cuts"))
    cal_summary->mkdir("cuts");
  cal_summary->Close();
  delete cal_summary;
  return kTRUE;
}

// L1 (the first long-L odd strip) seeds all other long-L cuts; R2 seeds the
// other long-R cuts. Returns kTRUE iff `c` should be derived from such a
// template, with the source cut name in `out_template_name`.
Bool_t CutTemplateFor(const ChannelCal &c, Int_t k,
                      TString &out_template_name) {
  if (IsShortSide(c) || c.side == 'S' || c.side == 'C')
    return kFALSE;
  if (c.side == 'L' && c.strip > 1) {
    out_template_name = Form("peak_L1_k%d", k);
    return kTRUE;
  }
  if (c.side == 'R' && c.strip > 2) {
    out_template_name = Form("peak_R2_k%d", k);
    return kTRUE;
  }
  return kFALSE;
}

// Open/close the cal-summary around every cut so an interrupt at any
// moment lands on a cleanly-closed file. ROOT's `SaveSelf(kTRUE)` flushes
// a TDirectory's keys but not the file-level FREE/header records, so
// keeping the file open across many cuts (the old model) left the
// /cuts/ key index unrecoverable after a Ctrl-C.
void DrawAllPeakCuts(const TString &cal_path,
                     const std::vector<ChannelCal> &chans, Int_t run) {
  TString label = Form("run%d", run);
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    Int_t n_k = IsShortSide(chans[i]) ? 1 : kNPeaks;
    for (Int_t k = 1; k <= n_k; k++) {
      TString cut_name = PeakCutName(chans[i], k);
      TFile *cal_summary = IO::OpenForWriting(cal_path, "UPDATE");
      if (!cal_summary || cal_summary->IsZombie()) {
        std::cerr << "  cannot reopen " << cal_path << " for cut '" << cut_name
                  << "'" << std::endl;
        if (cal_summary)
          delete cal_summary;
        continue;
      }
      TDirectory *cuts_dir = cal_summary->GetDirectory("cuts");

      PeakFit cached_fit;
      Long64_t cached_n = 0;
      Bool_t has_cached_fit =
          LoadCachedPeakFit(cal_summary, cut_name, cached_fit, cached_n);
      Bool_t bad_cached_fit =
          has_cached_fit && (!cached_fit.ok || cached_fit.mean <= 0);
      if (!Constants::BEAM_CAL_REDRAW_CUTS && cuts_dir &&
          cuts_dir->Get(cut_name) && !bad_cached_fit) {
        std::cout << "  reusing existing cut '" << cut_name << "'" << std::endl;
        cal_summary->Close();
        delete cal_summary;
        continue;
      }
      if (bad_cached_fit) {
        std::cout << "  cached fit for '" << cut_name
                  << "' failed (mu=" << cached_fit.mean
                  << ", ok=" << cached_fit.ok
                  << "); reopening cut interactively" << std::endl;
        if (cuts_dir && cuts_dir->Get(cut_name)) {
          cuts_dir->Delete(cut_name + ";*");
          cuts_dir->SaveSelf(kTRUE);
        }
      }
      InvalidateCachedPeakFit(cal_summary, cut_name);
      TString tmpl_name;
      Bool_t used_template = kFALSE;
      if (CutTemplateFor(chans[i], k, tmpl_name)) {
        TCutG *tmpl =
            cuts_dir ? static_cast<TCutG *>(cuts_dir->Get(tmpl_name)) : nullptr;
        if (tmpl) {
          DrawCutNS::DrawCutOnFileFromTemplate(
              cal_summary, ChannelHistName(chans[i]), cut_name, tmpl, label);
          used_template = kTRUE;
        } else {
          std::cerr << "  template '" << tmpl_name << "' missing; manual draw"
                    << std::endl;
        }
      }
      if (!used_template)
        DrawCutNS::DrawCutOnFile(cal_summary, ChannelHistName(chans[i]),
                                 cut_name, label);
      cal_summary->Close();
      delete cal_summary;
    }
  }
}

// Copy every `peak_<chan>_kN` cut from the cal-summary's /cuts/ into each
// subfile's Events_*.root:/cuts/. This is what makes the cuts visible to
// DrawCut and RegionTraces downstream.
void MirrorCutsToEventsFiles(TFile *cal_summary,
                             const std::vector<FileSpec> &specs,
                             const std::vector<ChannelCal> &chans) {
  TDirectory *src_cuts = cal_summary->GetDirectory("cuts");
  if (!src_cuts) {
    std::cerr << "  no /cuts/ in cal summary; nothing to mirror" << std::endl;
    return;
  }
  for (Int_t s = 0; s < Int_t(specs.size()); s++) {
    TString events_subpath = EventsName(specs[s]) + ".root";
    TFile *ef = IO::OpenForWriting(events_subpath, "UPDATE");
    if (!ef || ef->IsZombie()) {
      std::cerr << "  cannot open events file for cut mirror: "
                << events_subpath << std::endl;
      if (ef)
        delete ef;
      continue;
    }
    TDirectory *dst_cuts = ef->GetDirectory("cuts");
    if (!dst_cuts)
      dst_cuts = ef->mkdir("cuts");
    for (Int_t i = 0; i < Int_t(chans.size()); i++) {
      Int_t n_k = IsShortSide(chans[i]) ? 1 : kNPeaks;
      for (Int_t k = 1; k <= n_k; k++) {
        TString cut_name = PeakCutName(chans[i], k);
        TCutG *src = static_cast<TCutG *>(src_cuts->Get(cut_name));
        if (!src)
          continue;
        dst_cuts->cd();
        TCutG *copy = static_cast<TCutG *>(src->Clone(cut_name));
        copy->Write(cut_name, TObject::kOverwrite);
        delete copy;
      }
    }
    dst_cuts->SaveSelf(kTRUE);
    ef->Close();
    delete ef;
  }
}

// Per-event values used inside cuts/fits. Pulled out so callers can use
// EnergyView's raw ADC arrays without going through Decode() (Decode flattens
// total[] under IGNORE_SHORT_STRIPS, which would mismatch the cut Y axis
// from EventBuilder, see EventBuilder.cpp:106).
inline Double_t CutXFromEvent(const ChannelCal &c, const EnergyView &ev) {
  if (c.side == 'C')
    return Double_t(ev.cathode_adc);
  if (IsShortSide(c)) {
    return (c.side == 'L') ? Double_t(ev.leftdE_adc[c.strip])
                           : Double_t(ev.rightdE_adc[c.strip]);
  }
  return Double_t(ev.totaldE_adc[c.strip]);
}

inline Double_t FitValueFromEvent(const ChannelCal &c, const EnergyView &ev) {
  if (c.side == 'S')
    return Double_t(ev.totaldE_adc[c.strip]);
  if (c.side == 'L')
    return Double_t(ev.leftdE_adc[c.strip]);
  if (c.side == 'R')
    return Double_t(ev.rightdE_adc[c.strip]);
  if (c.side == 'C')
    return Double_t(ev.cathode_adc);
  return 0.0;
}

inline Double_t RawEventTotalFromAdc(const EnergyView &ev) {
  Double_t t = 0.0;
  for (Int_t s = 0; s < 18; s++)
    t += Double_t(ev.totaldE_adc[s]);
  return t;
}

PeakFit FitGaussianFromEvents(const std::vector<Double_t> &events,
                              Double_t adc_max, const TString &peak_name,
                              const TString &input_name) {
  PeakFit p;
  if (Long64_t(events.size()) < kMinEntriesForFit)
    return p;
  TH1D h_range(TString("h_range_") + peak_name, "", 400, 0.0, adc_max);
  h_range.SetDirectory(nullptr);
  for (Int_t i = 0; i < Int_t(events.size()); i++)
    h_range.Fill(events[i]);
  Double_t lo = 0, hi = 0;
  if (!EstimateFitRange(&h_range, lo, hi))
    return p;
  Float_t bin_w = Float_t((hi - lo) / 80.0);
  RooFitUtils fitter(events, Float_t(lo), Float_t(hi), bin_w,
                     /*use_flat_background=*/kTRUE,
                     /*use_step=*/kFALSE,
                     /*use_low_exp_tail=*/kFALSE,
                     /*use_low_lin_tail=*/kFALSE,
                     /*use_high_exp_tail=*/kFALSE);
  return PackFitResult(fitter.FitSinglePeak(input_name, peak_name));
}

struct CutJob {
  Int_t channel_idx;
  Int_t k;
  TCutG *cut;
  TString cut_name;
  Double_t adc_max;
  std::vector<Double_t> events;
};

void CollectInCutEventsWorker(
    Int_t worker_id, const std::vector<FileSpec> *specs,
    const std::vector<ChannelCal> *chans, const std::vector<CutJob> *jobs,
    std::vector<std::vector<std::vector<Double_t>>> *per_worker,
    std::atomic<Long64_t> *counts, std::queue<Int_t> *work,
    std::mutex *work_mutex) {
  Int_t n_jobs = Int_t(jobs->size());
  while (kTRUE) {
    Int_t spec_idx;
    {
      std::lock_guard<std::mutex> lk(*work_mutex);
      if (work->empty())
        return;
      spec_idx = work->front();
      work->pop();
    }
    TString events_subpath = EventsName((*specs)[spec_idx]) + ".root";
    TFile *sf = IO::OpenForReading(events_subpath);
    if (!sf || sf->IsZombie()) {
      if (sf)
        sf->Close();
      continue;
    }
    TTree *tree = static_cast<TTree *>(sf->Get("events"));
    if (!tree) {
      sf->Close();
      delete sf;
      continue;
    }
    EnergyView ev;
    ev.Attach(tree);
    Long64_t n = tree->GetEntries();
    std::vector<std::vector<Double_t>> &local = (*per_worker)[worker_id];
    for (Long64_t j = 0; j < n; j++) {
      tree->GetEntry(j);
      Double_t cy = RawEventTotalFromAdc(ev);
      Bool_t cathode_zero = (ev.cathode_adc <= 0);
      for (Int_t ji = 0; ji < n_jobs; ji++) {
        const CutJob &job = (*jobs)[ji];
        const ChannelCal &c = (*chans)[job.channel_idx];
        if (c.side == 'C' && cathode_zero)
          continue;
        if (counts[ji].load(std::memory_order_relaxed) >= kMaxEventsPerJob)
          continue;
        Double_t cx = CutXFromEvent(c, ev);
        if (!job.cut->IsInside(cx, cy))
          continue;
        Double_t v = FitValueFromEvent(c, ev);
        if (v <= 0)
          continue;
        Long64_t prev = counts[ji].fetch_add(1, std::memory_order_relaxed);
        if (prev < kMaxEventsPerJob)
          local[ji].push_back(v);
      }
    }
    sf->Close();
    delete sf;
    Bool_t all_full = kTRUE;
    for (Int_t ji = 0; ji < n_jobs; ji++) {
      if (counts[ji].load(std::memory_order_relaxed) < kMaxEventsPerJob) {
        all_full = kFALSE;
        break;
      }
    }
    if (all_full) {
      std::lock_guard<std::mutex> lk(*work_mutex);
      while (!work->empty())
        work->pop();
      return;
    }
  }
}

void CollectInCutEvents(const std::vector<FileSpec> &specs,
                        const std::vector<ChannelCal> &chans,
                        std::vector<CutJob> &jobs) {
  if (jobs.empty())
    return;

  ROOT::EnableThreadSafety();

  Int_t n_specs = Int_t(specs.size());
  Int_t n_workers =
      TMath::Min(Int_t(std::thread::hardware_concurrency()), n_specs);
  n_workers = TMath::Min(n_workers, Constants::MAX_FUSED_WORKERS);
  if (n_workers < 1)
    n_workers = 1;

  std::cout << "  in-cut collect: " << n_specs << " subfiles, " << jobs.size()
            << " jobs, " << n_workers << " workers, cap " << kMaxEventsPerJob
            << " events/job" << std::endl;

  std::queue<Int_t> work;
  for (Int_t k = 0; k < n_specs; k++)
    work.push(k);
  std::mutex work_mutex;

  std::vector<std::vector<std::vector<Double_t>>> per_worker(n_workers);
  for (Int_t w = 0; w < n_workers; w++)
    per_worker[w].resize(jobs.size());

  Int_t n_jobs = Int_t(jobs.size());
  std::atomic<Long64_t> *counts = new std::atomic<Long64_t>[n_jobs];
  for (Int_t ji = 0; ji < n_jobs; ji++)
    counts[ji].store(0);

  std::vector<std::thread> workers;
  for (Int_t w = 0; w < n_workers; w++)
    workers.emplace_back(CollectInCutEventsWorker, w, &specs, &chans, &jobs,
                         &per_worker, counts, &work, &work_mutex);
  for (Int_t w = 0; w < Int_t(workers.size()); w++)
    workers[w].join();

  for (Int_t ji = 0; ji < n_jobs; ji++) {
    Long64_t total = 0;
    for (Int_t w = 0; w < n_workers; w++)
      total += Long64_t(per_worker[w][ji].size());
    jobs[ji].events.reserve(jobs[ji].events.size() + total);
    for (Int_t w = 0; w < n_workers; w++) {
      std::vector<Double_t> &src = per_worker[w][ji];
      for (Int_t i = 0; i < Int_t(src.size()); i++)
        jobs[ji].events.push_back(src[i]);
    }
  }
  delete[] counts;
}

// MeV-space sigma matching for Strip0/Strip17/Cathode. Uses the local pol2
// derivative at peak_mu_ADC to convert each long-side k=1 sigma_ADC into
// sigma_MeV, takes the mean, then rescales each S/C channel's peak_sigma_ADC
// so its k=1 sigma_MeV matches that reference.
void ScaleSigmasToLongSideMeanMeV(std::vector<ChannelCal> &chans) {
  Double_t sum_mev = 0.0;
  Int_t n = 0;
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    const ChannelCal &c = chans[i];
    if (c.side != 'L' && c.side != 'R')
      continue;
    if (c.strip < 1 || c.strip > 16)
      continue;
    if (IsShortSide(c))
      continue;
    if (!c.quad_ok)
      continue;
    const PeakFit &k1 = c.peak_fits[0];
    if (!k1.ok || k1.sigma <= 0 || k1.mean <= 0)
      continue;
    Double_t slope = c.quad_b + 2.0 * c.quad_c * k1.mean;
    Double_t sigma_mev = k1.sigma * slope;
    if (sigma_mev <= 0)
      continue;
    sum_mev += sigma_mev;
    n++;
  }
  if (n == 0) {
    std::cerr << "  ScaleSigmasToLongSideMeanMeV: no long-side k=1 fits; "
                 "skipping rescale"
              << std::endl;
    return;
  }
  Double_t ref_mev = sum_mev / Double_t(n);
  std::cout << "  long-side k=1 mean peak_sigma_MeV = " << ref_mev
            << " (n=" << n << ")" << std::endl;
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    ChannelCal &c = chans[i];
    if (c.side != 'S' && c.side != 'C') {
      c.sigma_match_scale = 1.0;
      continue;
    }
    if (!c.quad_ok) {
      c.sigma_match_scale = 1.0;
      continue;
    }
    const PeakFit &k1 = c.peak_fits[0];
    if (!k1.ok || k1.sigma <= 0 || k1.mean <= 0) {
      c.sigma_match_scale = 1.0;
      continue;
    }
    Double_t slope = c.quad_b + 2.0 * c.quad_c * k1.mean;
    Double_t current_sigma_mev = k1.sigma * slope;
    if (current_sigma_mev <= 0) {
      c.sigma_match_scale = 1.0;
      continue;
    }
    Double_t factor = ref_mev / current_sigma_mev;
    c.sigma_match_scale = factor;
    for (Int_t k = 0; k < kNPeaks; k++) {
      if (c.peak_fits[k].ok && c.peak_fits[k].sigma > 0)
        c.peak_fits[k].sigma *= factor;
    }
    std::cout << "  " << c.name << " MeV sigma-rescale factor=" << factor
              << " (orig sigma_MeV=" << current_sigma_mev
              << ", target=" << ref_mev << ")" << std::endl;
  }
}

// Long-side / unsplit / cathode: pol2 fit on up to 4 (peak_mu_k, k*sim_mu)
// data points plus an anchor (0, 0). Short-side: linear from k=1 only.
void FitQuadratic(ChannelCal &c) {
  Double_t sim_mu_mev = c.sim_fit.mean;
  if (sim_mu_mev <= 0) {
    c.quad_ok = kFALSE;
    return;
  }

  if (IsShortSide(c)) {
    const PeakFit &pf = c.peak_fits[0];
    if (!pf.ok || pf.mean <= 0) {
      c.quad_ok = kFALSE;
      std::cerr << "  " << c.name << ": short-side k=1 fit failed (ok=" << pf.ok
                << ", mean=" << pf.mean << "); no linear cal" << std::endl;
      return;
    }
    c.quad_a = 0.0;
    c.quad_b = sim_mu_mev / pf.mean;
    c.quad_c = 0.0;
    c.quad_chi2_ndf = -1;
    c.quad_ok = kTRUE;
    std::cout << "  " << c.name << " linear: peak_mu=" << pf.mean
              << " ADC, slope=" << c.quad_b << std::endl;
    return;
  }

  TGraph g;
  Int_t np = 0;
  // Zero-point anchor: forces the pol2 toward MeV(0 ADC) = 0 so the cal
  // doesn't dip negative at low ADC.
  g.SetPoint(np++, 0.0, 0.0);
  std::cout << "  " << c.name << " quad points: (0, 0)";
  for (Int_t k = 1; k <= kNPeaks; k++) {
    const PeakFit &pf = c.peak_fits[k - 1];
    if (!pf.ok || pf.mean <= 0)
      continue;
    g.SetPoint(np, pf.mean, Double_t(k) * sim_mu_mev);
    std::cout << " (" << pf.mean << ", " << Double_t(k) * sim_mu_mev << ")";
    np++;
  }
  std::cout << std::endl;
  if (np < 4) {
    c.quad_ok = kFALSE;
    std::cerr << "  " << c.name << ": only " << (np - 1)
              << " valid peak fits (need 3+ besides the zero anchor)"
              << std::endl;
    return;
  }

  TF1 f("f_quad", "pol2", 0.0, 16384.0);
  // "W" forces unit weights; without errors the default chi2 fit is
  // ill-defined for exact-point sets.
  TFitResultPtr r = g.Fit(&f, "QSW");
  if (!r.Get() || !r->IsValid()) {
    c.quad_ok = kFALSE;
    std::cerr << "  " << c.name << ": pol2 fit invalid" << std::endl;
    return;
  }
  c.quad_a = f.GetParameter(0);
  c.quad_b = f.GetParameter(1);
  c.quad_c = f.GetParameter(2);
  c.quad_chi2_ndf =
      (f.GetNDF() > 0) ? f.GetChisquare() / Double_t(f.GetNDF()) : -1;
  c.quad_ok = kTRUE;
}

inline Double_t ApplyQuad(const ChannelCal &c, Double_t adc) {
  return c.quad_a + c.quad_b * adc + c.quad_c * adc * adc;
}

// Emit a [detector.eres] TOML snippet for the simulator. Keys are Cathode,
// S0, S17, L1..L16, R1..R16; values are per-channel k=1 sigma in MeV after
// the long-side-mean rescale. Channels with no valid fit get -1 (sim
// interprets that as no noise on that channel).
void WriteEresTomlForRun(Int_t run, const std::vector<ChannelCal> &chans) {
  Double_t eres_vals[35];
  for (Int_t i = 0; i < 35; i++)
    eres_vals[i] = -1.0;

  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    const ChannelCal &c = chans[i];
    if (!c.quad_ok)
      continue;
    const PeakFit &k1 = c.peak_fits[0];
    if (!k1.ok || k1.mean <= 0 || k1.sigma <= 0)
      continue;
    Double_t slope = c.quad_b + 2.0 * c.quad_c * k1.mean;
    Double_t sigma_mev = k1.sigma * slope;
    if (sigma_mev <= 0)
      continue;
    Int_t idx = -1;
    if (c.side == 'C')
      idx = 0;
    else if (c.side == 'S' && c.strip == 0)
      idx = 1;
    else if (c.side == 'S' && c.strip == 17)
      idx = 2;
    else if (c.side == 'L' && c.strip >= 1 && c.strip <= 16)
      idx = 3 + (c.strip - 1);
    else if (c.side == 'R' && c.strip >= 1 && c.strip <= 16)
      idx = 19 + (c.strip - 1);
    if (idx >= 0)
      eres_vals[idx] = sigma_mev;
  }

  toml::table eres_tbl;
  eres_tbl.insert("Cathode", eres_vals[0]);
  eres_tbl.insert("S0", eres_vals[1]);
  eres_tbl.insert("S17", eres_vals[2]);
  for (Int_t s = 1; s <= 16; s++) {
    std::string key = "L" + std::to_string(s);
    eres_tbl.insert(key, eres_vals[3 + (s - 1)]);
  }
  for (Int_t s = 1; s <= 16; s++) {
    std::string key = "R" + std::to_string(s);
    eres_tbl.insert(key, eres_vals[19 + (s - 1)]);
  }
  toml::table detector_tbl;
  detector_tbl.insert("eres", eres_tbl);
  toml::table root_tbl;
  root_tbl.insert("detector", detector_tbl);

  TString out_subpath = Form("Calibration_Run%d_eres.toml", run);
  TString out_full = IO::GetRootFilesBaseDir() + "/" + out_subpath;
  std::ofstream f(out_full.Data());
  if (!f) {
    std::cerr << "Cannot write eres TOML: " << out_full << std::endl;
    return;
  }
  f << root_tbl << std::endl;
  std::cout << "  wrote eres TOML: " << out_full << std::endl;
}

void WriteCalibrationTreeToSummary(TFile *cal_summary,
                                   const std::vector<ChannelCal> &chans) {
  cal_summary->cd();
  if (TObject *old = cal_summary->Get("calibration"))
    old->Delete();
  TTree *cal = new TTree("calibration", "Per-channel pol2 calibration");
  // Fixed-size branches so the schema is stable regardless of how
  // IGNORE_SHORT_STRIPS shrinks chans. Unused entries stay zero.
  Float_t quad_a[kMaxChannels] = {0}, quad_b[kMaxChannels] = {0},
          quad_c[kMaxChannels] = {0}, quad_chi2[kMaxChannels] = {0};
  Float_t sim_mu[kMaxChannels] = {0}, sim_sigma[kMaxChannels] = {0};
  Float_t sigma_match[kMaxChannels] = {0};
  Float_t peak_mu[kNPeaks][kMaxChannels] = {{0}},
          peak_sigma[kNPeaks][kMaxChannels] = {{0}};
  Long64_t peak_n[kNPeaks][kMaxChannels] = {{0}};
  Bool_t quad_ok[kMaxChannels] = {0};
  Int_t n_actual = TMath::Min(Int_t(chans.size()), kMaxChannels);
  for (Int_t k = 0; k < n_actual; k++) {
    const ChannelCal &c = chans[k];
    quad_a[k] = Float_t(c.quad_a);
    quad_b[k] = Float_t(c.quad_b);
    quad_c[k] = Float_t(c.quad_c);
    quad_chi2[k] = Float_t(c.quad_chi2_ndf);
    sim_mu[k] = Float_t(c.sim_fit.mean);
    sim_sigma[k] = Float_t(c.sim_fit.sigma);
    sigma_match[k] = Float_t(c.sigma_match_scale);
    for (Int_t p = 0; p < kNPeaks; p++) {
      peak_mu[p][k] = Float_t(c.peak_fits[p].mean);
      peak_sigma[p][k] = Float_t(c.peak_fits[p].sigma);
      peak_n[p][k] = c.peak_n[p];
    }
    quad_ok[k] = c.quad_ok;
  }
  cal->Branch("QuadA", quad_a, Form("QuadA[%d]/F", kMaxChannels));
  cal->Branch("QuadB", quad_b, Form("QuadB[%d]/F", kMaxChannels));
  cal->Branch("QuadC", quad_c, Form("QuadC[%d]/F", kMaxChannels));
  cal->Branch("QuadChi2NDF", quad_chi2,
              Form("QuadChi2NDF[%d]/F", kMaxChannels));
  cal->Branch("QuadOK", quad_ok, Form("QuadOK[%d]/O", kMaxChannels));
  cal->Branch("SimMu_MeV", sim_mu, Form("SimMu_MeV[%d]/F", kMaxChannels));
  cal->Branch("SimSigma_MeV", sim_sigma,
              Form("SimSigma_MeV[%d]/F", kMaxChannels));
  cal->Branch("SigmaMatchScale", sigma_match,
              Form("SigmaMatchScale[%d]/F", kMaxChannels));
  for (Int_t p = 0; p < kNPeaks; p++) {
    cal->Branch(Form("PeakMu_k%d_ADC", p + 1), peak_mu[p],
                Form("PeakMu_k%d_ADC[%d]/F", p + 1, kMaxChannels));
    cal->Branch(Form("PeakSigma_k%d_ADC", p + 1), peak_sigma[p],
                Form("PeakSigma_k%d_ADC[%d]/F", p + 1, kMaxChannels));
    cal->Branch(Form("PeakN_k%d", p + 1), peak_n[p],
                Form("PeakN_k%d[%d]/L", p + 1, kMaxChannels));
  }
  cal->Fill();
  cal->Write("calibration", TObject::kOverwrite);
}

// Write Events_RunXX[_M].cal.root with a single TTree `events_cal` holding
// LeftdEMeV/RightdEMeV/TotaldEMeV/CathodeMeV. Always RECREATE — no in-place
// rewrite of the events tree, so duplicate-branch corruption is impossible
// and an interrupted run leaves the events file untouched.
void WriteCalSidecar(const FileSpec &spec,
                     const std::vector<ChannelCal> &chans) {
  TString events_subpath = EventsName(spec) + ".root";
  TFile *src = IO::OpenForReading(events_subpath);
  if (!src || src->IsZombie()) {
    std::cerr << "Cannot open " << events_subpath << std::endl;
    if (src)
      delete src;
    return;
  }
  TTree *events = static_cast<TTree *>(src->Get("events"));
  if (!events) {
    std::cerr << "No events tree in " << events_subpath << std::endl;
    src->Close();
    delete src;
    return;
  }
  EnergyView ev;
  ev.Attach(events);
  if (ev.is_mev) {
    std::cerr << "WARNING: " << events_subpath
              << " already exposes MeV (stale in-place calibration?); cal "
                 "sidecar will still be written from raw ADC branches"
              << std::endl;
  }
  // Re-bind to the raw ADC branches even if Attach picked up a stale MeV
  // friend, so the sidecar is always built from authoritative ADC values.
  events->SetBranchAddress("LeftdE", ev.leftdE_adc);
  events->SetBranchAddress("RightdE", ev.rightdE_adc);
  events->SetBranchAddress("TotaldE", ev.totaldE_adc);
  events->SetBranchAddress("Cathode", &ev.cathode_adc);

  const ChannelCal *cal_S[18] = {nullptr};
  const ChannelCal *cal_L[18] = {nullptr};
  const ChannelCal *cal_R[18] = {nullptr};
  const ChannelCal *cal_C = nullptr;
  for (Int_t k = 0; k < Int_t(chans.size()); k++) {
    const ChannelCal &c = chans[k];
    if (c.side == 'S')
      cal_S[c.strip] = &c;
    else if (c.side == 'L' && c.strip >= 1 && c.strip <= 16)
      cal_L[c.strip] = &c;
    else if (c.side == 'R' && c.strip >= 1 && c.strip <= 16)
      cal_R[c.strip] = &c;
    else if (c.side == 'C')
      cal_C = &c;
  }

  TString sidecar_subpath = CalSidecarName(spec);
  TFile *dst = IO::OpenForWriting(sidecar_subpath, "RECREATE");
  if (!dst || dst->IsZombie()) {
    std::cerr << "Cannot create cal sidecar " << sidecar_subpath << std::endl;
    if (dst)
      delete dst;
    src->Close();
    delete src;
    return;
  }
  dst->cd();
  TTree *cal_tree = new TTree("events_cal", "Per-event calibrated energies");
  Float_t leftdE_MeV[18], rightdE_MeV[18], totaldE_MeV[18];
  Float_t cathode_MeV;
  cal_tree->Branch("LeftdEMeV", leftdE_MeV, "LeftdEMeV[18]/F");
  cal_tree->Branch("RightdEMeV", rightdE_MeV, "RightdEMeV[18]/F");
  cal_tree->Branch("TotaldEMeV", totaldE_MeV, "TotaldEMeV[18]/F");
  cal_tree->Branch("CathodeMeV", &cathode_MeV, "CathodeMeV/F");

  Long64_t n = events->GetEntries();
  std::cout << "  writing cal sidecar " << sidecar_subpath << " (" << n
            << " events)" << std::endl;
  for (Long64_t j = 0; j < n; j++) {
    events->GetEntry(j);
    for (Int_t s = 0; s < 18; s++) {
      leftdE_MeV[s] = 0.0f;
      rightdE_MeV[s] = 0.0f;
      totaldE_MeV[s] = 0.0f;
    }
    for (Int_t s = 0; s < 18; s++) {
      if (s == 0 || s == 17) {
        const ChannelCal *cS = cal_S[s];
        if (cS && cS->quad_ok && ev.totaldE_adc[s] > 0)
          totaldE_MeV[s] = Float_t(ApplyQuad(*cS, Double_t(ev.totaldE_adc[s])));
      } else {
        const ChannelCal *cL = cal_L[s];
        const ChannelCal *cR = cal_R[s];
        if (cL && cL->quad_ok && ev.leftdE_adc[s] > 0)
          leftdE_MeV[s] = Float_t(ApplyQuad(*cL, Double_t(ev.leftdE_adc[s])));
        if (cR && cR->quad_ok && ev.rightdE_adc[s] > 0)
          rightdE_MeV[s] = Float_t(ApplyQuad(*cR, Double_t(ev.rightdE_adc[s])));
        totaldE_MeV[s] = leftdE_MeV[s] + rightdE_MeV[s];
      }
    }
    cathode_MeV = (cal_C && cal_C->quad_ok && ev.cathode_adc > 0)
                      ? Float_t(ApplyQuad(*cal_C, Double_t(ev.cathode_adc)))
                      : 0.0f;
    cal_tree->Fill();
  }
  cal_tree->Write("events_cal", TObject::kOverwrite);
  dst->Close();
  delete dst;
  src->Close();
  delete src;
}

// Diagnostic for the MeV sigma rescale: for Strip0/Strip17/Cathode, dump
// the in-cut event distribution twice -- raw ADC (pre-rescale data) and
// MeV via the channel's pol2 (post-rescale data). No fits overlaid.
void SaveRescaledStripHistograms(const std::vector<FileSpec> &specs,
                                 const std::vector<ChannelCal> &chans,
                                 TFile *cal_summary, const TString &run_label,
                                 const TString &plot_subdir) {
  TDirectory *cuts_dir = cal_summary->GetDirectory("cuts");
  if (!cuts_dir) {
    std::cerr << "  SaveRescaledStripHistograms: no cuts/ dir" << std::endl;
    return;
  }
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    const ChannelCal &c = chans[i];
    if (c.side != 'S' && c.side != 'C')
      continue;
    if (!c.quad_ok)
      continue;
    for (Int_t k = 1; k <= kNPeaks; k++) {
      TString cut_name = PeakCutName(c, k);
      TCutG *cut = static_cast<TCutG *>(cuts_dir->Get(cut_name));
      if (!cut) {
        std::cerr << "  rescale-diag: missing cut " << cut_name << std::endl;
        continue;
      }
      Double_t adc_max = 16384.0;
      TString adc_name =
          Form("h_rescale_adc_%s_%s_k%d", run_label.Data(), c.name.Data(), k);
      TH1D *h_adc =
          new TH1D(adc_name,
                   Form("%s k=%d in-cut (ADC);%s #DeltaE [ADC];Counts",
                        c.name.Data(), k, c.name.Data()),
                   400, 0.0, adc_max);
      h_adc->SetDirectory(nullptr);

      Double_t mev_max =
          (c.sim_fit.mean > 0) ? Double_t(k) * c.sim_fit.mean * 2.0 : 60.0;
      TString mev_name =
          Form("h_rescale_mev_%s_%s_k%d", run_label.Data(), c.name.Data(), k);
      TH1D *h_mev =
          new TH1D(mev_name,
                   Form("%s k=%d in-cut (MeV);%s #DeltaE [MeV];Counts",
                        c.name.Data(), k, c.name.Data()),
                   400, 0.0, mev_max);
      h_mev->SetDirectory(nullptr);

      Long64_t n_in = 0;
      for (Int_t s = 0; s < Int_t(specs.size()); s++) {
        TString events_subpath = EventsName(specs[s]) + ".root";
        TFile *sf = IO::OpenForReading(events_subpath);
        if (!sf || sf->IsZombie()) {
          if (sf)
            sf->Close();
          continue;
        }
        TTree *tree = static_cast<TTree *>(sf->Get("events"));
        if (!tree) {
          sf->Close();
          delete sf;
          continue;
        }
        EnergyView ev;
        ev.Attach(tree);
        Long64_t n = tree->GetEntries();
        for (Long64_t j = 0; j < n; j++) {
          tree->GetEntry(j);
          if (c.side == 'C' && ev.cathode_adc <= 0)
            continue;
          Double_t cx = CutXFromEvent(c, ev);
          Double_t cy = RawEventTotalFromAdc(ev);
          if (!cut->IsInside(cx, cy))
            continue;
          Double_t v_adc = FitValueFromEvent(c, ev);
          if (v_adc <= 0)
            continue;
          h_adc->Fill(v_adc);
          h_mev->Fill(ApplyQuad(c, v_adc));
          n_in++;
        }
        sf->Close();
        delete sf;
      }
      std::cout << "  rescale-diag " << c.name << " k=" << k << ": " << n_in
                << " in-cut events" << std::endl;

      TCanvas *cv_adc = PlottingUtils::GetConfiguredCanvas(kFALSE);
      h_adc->Draw("HIST");
      if (Constants::SAVE_PLOTS)
        PlottingUtils::SaveFigure(cv_adc,
                                  Form("rescale_adc_%s_k%d", c.name.Data(), k),
                                  plot_subdir, PlotSaveOptions::kLINEAR);
      delete cv_adc;

      TCanvas *cv_mev = PlottingUtils::GetConfiguredCanvas(kFALSE);
      h_mev->Draw("HIST");
      if (Constants::SAVE_PLOTS)
        PlottingUtils::SaveFigure(cv_mev,
                                  Form("rescale_mev_%s_k%d", c.name.Data(), k),
                                  plot_subdir, PlotSaveOptions::kLINEAR);
      delete cv_mev;

      delete h_adc;
      delete h_mev;
    }
  }
}

void SaveQuadFitPlots(const std::vector<ChannelCal> &chans,
                      const TString &plot_subdir, const TString &run_label) {
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    const ChannelCal &c = chans[i];
    if (!c.quad_ok)
      continue;
    Int_t k_hi = IsShortSide(c) ? 1 : kNPeaks;
    TGraphErrors *g = new TGraphErrors();
    Int_t np = 0;
    for (Int_t k = 1; k <= k_hi; k++) {
      const PeakFit &pf = c.peak_fits[k - 1];
      if (!pf.ok || pf.mean <= 0)
        continue;
      g->SetPoint(np, pf.mean, Double_t(k) * c.sim_fit.mean);
      g->SetPointError(np, std::max(1.0, pf.mean_err), 0.0);
      np++;
    }
    Double_t xhi = 1.1 * c.peak_fits[k_hi - 1].mean;
    if (xhi <= 0)
      xhi = 16384.0;
    PlottingUtils::ConfigureGraph(
        g, kBlack, Form(";%s #DeltaE [ADC];#DeltaE [MeV]", c.name.Data()));
    g->SetMarkerStyle(20);
    TF1 *f = new TF1(Form("f_quad_%s", c.name.Data()), "pol2", 0.0, xhi);
    f->SetParameters(c.quad_a, c.quad_b, c.quad_c);
    f->SetLineColor(kRed + 1);
    f->SetLineWidth(2);
    TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
    g->Draw("AP");
    f->Draw("L SAME");
    if (Constants::SAVE_PLOTS)
      PlottingUtils::SaveFigure(cv, "quad_" + c.name, plot_subdir,
                                PlotSaveOptions::kLINEAR);
    delete cv;
    delete f;
    delete g;
  }
}

// Reads the just-written cal sidecars via the canonical AttachCalSidecar
// path, so we exercise the same code downstream macros use.
void SaveDynamicRangeOverlay(const std::vector<FileSpec> &specs,
                             const std::vector<ChannelCal> &chans,
                             const TString &plot_subdir,
                             const TString &run_label) {
  const Int_t n_chans = Int_t(chans.size());
  const Int_t nbins = 300;
  const Double_t emin = 0.0;
  const Double_t emax = 50;
  std::vector<TH1D *> h(n_chans, nullptr);
  for (Int_t i = 0; i < n_chans; i++) {
    TString hname =
        Form("h_dynrange_%s_%s", run_label.Data(), chans[i].name.Data());
    h[i] = new TH1D(hname, ";#DeltaE [MeV];Counts", nbins, emin, emax);
    h[i]->SetDirectory(nullptr);
  }

  for (Int_t s = 0; s < Int_t(specs.size()); s++) {
    TString events_subpath = EventsName(specs[s]) + ".root";
    TFile *sf = IO::OpenForReading(events_subpath);
    if (!sf || sf->IsZombie()) {
      if (sf)
        sf->Close();
      continue;
    }
    TTree *tree = static_cast<TTree *>(sf->Get("events"));
    if (!tree) {
      sf->Close();
      delete sf;
      continue;
    }
    TFile *cal_file = AttachCalSidecar(tree, specs[s]);
    EnergyView ev;
    ev.Attach(tree);
    if (!ev.is_mev) {
      sf->Close();
      delete sf;
      if (cal_file) {
        cal_file->Close();
        delete cal_file;
      }
      continue;
    }
    Long64_t n = tree->GetEntries();
    for (Long64_t j = 0; j < n; j++) {
      tree->GetEntry(j);
      ev.Decode();
      for (Int_t i = 0; i < n_chans; i++) {
        const ChannelCal &c = chans[i];
        if (IsShortSide(c))
          continue;
        Double_t v = 0;
        if (c.side == 'C')
          v = ev.cathode;
        else if (c.side == 'S')
          v = ev.total[c.strip];
        else if (c.side == 'L')
          v = ev.left[c.strip];
        else if (c.side == 'R')
          v = ev.right[c.strip];
        if (v > 0)
          h[i]->Fill(v);
      }
    }
    sf->Close();
    delete sf;
    if (cal_file) {
      cal_file->Close();
      delete cal_file;
    }
  }

  std::vector<Int_t> colors = PlottingUtils::GetDefaultColors();
  Double_t y_top = 0;
  for (Int_t i = 0; i < n_chans; i++) {
    if (IsShortSide(chans[i]))
      continue;
    Double_t m = h[i]->GetMaximum();
    if (m > y_top)
      y_top = m;
  }
  TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
  cv->SetRightMargin(0.20);
  Bool_t first = kTRUE;
  for (Int_t i = 0; i < n_chans; i++) {
    if (IsShortSide(chans[i]))
      continue;
    Int_t color = colors[i % Int_t(colors.size())];
    h[i]->SetLineColor(color);
    h[i]->SetLineWidth(2);
    h[i]->SetMaximum(1.15 * y_top);
    h[i]->Draw(first ? "HIST" : "HIST SAME");
    first = kFALSE;
  }
  TLegend *leg = PlottingUtils::AddLegend(0.81, 0.99, 0.10, 0.95);
  for (Int_t i = 0; i < n_chans; i++) {
    if (IsShortSide(chans[i]))
      continue;
    leg->AddEntry(h[i], chans[i].name.Data(), "l");
  }
  leg->Draw();
  if (Constants::SAVE_PLOTS)
    PlottingUtils::SaveFigure(cv, "dynamic_range_check", plot_subdir,
                              PlotSaveOptions::kLOG);
  delete cv;
  delete leg;
  for (Int_t i = 0; i < n_chans; i++)
    delete h[i];
}

void CalibrateBeamOneRun(Int_t run, const std::vector<FileSpec> &specs,
                         const TString &sim_beam_path) {
  TString run_label = Form("run%d", run);
  TString plot_subdir = "beam_calibration/" + run_label;

  std::cout << "Beam calibration: run " << run << " (" << specs.size()
            << " subfiles)" << std::endl;

  std::vector<ChannelCal> chans = BuildChannels();

  TString cal_path = CalSummarySubpath(run);

  if (!BuildOrOpenRunSummaryFile(run, specs, chans)) {
    std::cerr << "Run " << run << ": failed to build cal summary; skipping"
              << std::endl;
    return;
  }

  // DrawCutOnFile flips batch off for interactivity; force it back on after
  // so downstream fit/sidecar/plot work runs headlessly. DrawAllPeakCuts
  // opens/closes cal_summary around each cut so an interrupt is always safe.
  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    DrawAllPeakCuts(cal_path, chans, run);
    gROOT->SetBatch(kTRUE);
  }

  if (!LoadSimHistograms(sim_beam_path, chans, run_label)) {
    std::cerr << "Run " << run << ": sim histograms unavailable; skipping"
              << std::endl;
    return;
  }
  for (Int_t k = 0; k < Int_t(chans.size()); k++) {
    ChannelCal &c = chans[k];
    if (!c.sim_hist_mev) {
      std::cerr << "  " << c.name << " sim: NO HISTOGRAM" << std::endl;
      continue;
    }
    c.sim_fit = FitGaussian(c.sim_hist_mev, /*flat_bg=*/kTRUE, "sim_" + c.name,
                            run_label);
    std::cout << "  " << c.name << " sim: mu=" << c.sim_fit.mean
              << " MeV, sigma=" << c.sim_fit.sigma
              << " (entries=" << c.sim_hist_mev->GetEntries()
              << (c.sim_fit.ok ? "" : ", FIT FAILED") << ")" << std::endl;
  }

  // Fit-phase open: non-interactive, fast, single open/close.
  TFile *cal_summary = IO::OpenForWriting(cal_path, "UPDATE");
  if (!cal_summary || cal_summary->IsZombie()) {
    std::cerr << "Run " << run << ": cannot reopen cal summary for fit phase"
              << std::endl;
    if (cal_summary)
      delete cal_summary;
    for (Int_t k = 0; k < Int_t(chans.size()); k++)
      delete chans[k].sim_hist_mev;
    return;
  }
  TDirectory *cuts_dir = cal_summary->GetDirectory("cuts");
  if (!cuts_dir) {
    std::cerr << "Run " << run << ": no cuts/ dir in cal summary" << std::endl;
    cal_summary->Close();
    delete cal_summary;
    for (Int_t k = 0; k < Int_t(chans.size()); k++)
      delete chans[k].sim_hist_mev;
    return;
  }
  std::vector<CutJob> jobs;
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    ChannelCal &c = chans[i];
    Int_t n_k = IsShortSide(c) ? 1 : kNPeaks;
    for (Int_t k = 1; k <= n_k; k++) {
      TString cut_name = PeakCutName(c, k);
      Long64_t cached_n = 0;
      if (LoadCachedPeakFit(cal_summary, cut_name, c.peak_fits[k - 1],
                            cached_n)) {
        c.peak_n[k - 1] = cached_n;
        std::cout << "  " << c.name << " k=" << k
                  << ": cached mu=" << c.peak_fits[k - 1].mean
                  << " ADC, sigma=" << c.peak_fits[k - 1].sigma
                  << " (n=" << cached_n << ")" << std::endl;
        continue;
      }
      TCutG *cut = static_cast<TCutG *>(cuts_dir->Get(cut_name));
      if (!cut) {
        std::cerr << "  missing cut " << cut_name << "; skipping" << std::endl;
        continue;
      }
      CutJob job;
      job.channel_idx = i;
      job.k = k;
      job.cut = cut;
      job.cut_name = cut_name;
      job.adc_max = (c.side == 'L' || c.side == 'R') ? 5500.0 : 16384.0;
      jobs.push_back(job);
    }
  }

  CollectInCutEvents(specs, chans, jobs);

  for (Int_t ji = 0; ji < Int_t(jobs.size()); ji++) {
    CutJob &job = jobs[ji];
    ChannelCal &c = chans[job.channel_idx];
    Long64_t n_events = Long64_t(job.events.size());
    c.peak_fits[job.k - 1] = FitGaussianFromEvents(
        job.events, job.adc_max, c.name + "_" + job.cut_name, run_label);
    c.peak_n[job.k - 1] = n_events;
    SaveCachedPeakFit(cal_summary, job.cut_name, c.peak_fits[job.k - 1],
                      n_events);
    std::cout << "  " << c.name << " k=" << job.k
              << ": mu=" << c.peak_fits[job.k - 1].mean
              << " ADC, sigma=" << c.peak_fits[job.k - 1].sigma
              << " (n=" << n_events << ")"
              << (c.peak_fits[job.k - 1].ok ? "" : " [FIT FAILED]")
              << std::endl;
  }

  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    ChannelCal &c = chans[i];
    if (!c.sim_fit.ok) {
      c.quad_ok = kFALSE;
      std::cerr << "  " << c.name << ": no sim mu; cannot fit quad"
                << std::endl;
      continue;
    }
    FitQuadratic(c);
    if (c.quad_ok)
      std::cout << "  " << c.name << " quad: a=" << c.quad_a
                << ", b=" << c.quad_b << ", c=" << c.quad_c
                << ", chi2/ndf=" << c.quad_chi2_ndf << std::endl;
  }

  ScaleSigmasToLongSideMeanMeV(chans);

  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    SaveRescaledStripHistograms(specs, chans, cal_summary, run_label,
                                plot_subdir);
  }

  WriteCalibrationTreeToSummary(cal_summary, chans);
  WriteEresTomlForRun(run, chans);
  MirrorCutsToEventsFiles(cal_summary, specs, chans);
  cal_summary->Close();
  delete cal_summary;

  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    SaveQuadFitPlots(chans, plot_subdir, run_label);
  }

  for (Int_t s = 0; s < Int_t(specs.size()); s++)
    WriteCalSidecar(specs[s], chans);

  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    SaveDynamicRangeOverlay(specs, chans, plot_subdir, run_label);
  }

  for (Int_t k = 0; k < Int_t(chans.size()); k++)
    delete chans[k].sim_hist_mev;

  std::cout << "Run " << run << " calibration complete." << std::endl;
}

} // namespace CalibrateBeamNS

void CalibrateBeam(TString file_label = "") {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");
  // InitUtils sets batch ON; we need it OFF for interactive cut drawing.
  gROOT->SetBatch(kFALSE);

  // Standalone GUI needs a TApplication so canvases/event loop work.
  if (!gApplication) {
    static Int_t app_argc = 1;
    static char app_arg0[] = "calibrate-beam";
    static char *app_argv[] = {app_arg0};
    new TApplication("calibrate-beam", &app_argc, app_argv);
  }

  std::vector<FileSpec> specs;
  if (file_label.IsNull()) {
    specs = BuildFileSpecs();
    if (specs.empty()) {
      std::cerr << "No file specs from BuildFileSpecs()" << std::endl;
      return;
    }
  } else {
    FileSpec s = ResolveFileSpec(file_label);
    if (s.run < 0) {
      std::cerr << "Could not resolve file label '" << file_label << "'"
                << std::endl;
      return;
    }
    specs.push_back(s);
  }

  std::map<Int_t, std::vector<FileSpec>> run_to_specs;
  for (Int_t k = 0; k < Int_t(specs.size()); k++)
    run_to_specs[specs[k].run].push_back(specs[k]);

  TString sim_beam_path = CalibrateBeamNS::SimBeamPath(project_root);

  for (std::map<Int_t, std::vector<FileSpec>>::const_iterator it =
           run_to_specs.begin();
       it != run_to_specs.end(); ++it) {
    CalibrateBeamNS::CalibrateBeamOneRun(it->first, it->second, sim_beam_path);
  }
}
