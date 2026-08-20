#include "StripSumScatter.hpp"

// Per-strip gains = 1/mean beam total so sim beam lands at 1 a.u. like the
// per-channel data normalization. No-beam strips keep gain 0 (drop out).
Bool_t StripSumScatter::SimBeamGains(Double_t *gain) {
  const Long64_t kSampleMaxPoints =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.SAMPLE_MAX_POINTS;

  for (Int_t s = 0; s < 18; s++)
    gain[s] = 0.0;
  // Prefer the ERES beam file so the normalized eres beam lands at 1 a.u.
  // like the plotted eres populations; fall back to the standard file.
  TString file;
  Bool_t file_eres = kFALSE;
  std::vector<RemixSim::SimFileSpec> specs = RemixSim::BuildFileSpecs();
  for (Int_t i = 0; i < Int_t(specs.size()); i++) {
    TString base = RemixSim::TagWithoutStrip(specs[i].tag);
    base.ReplaceAll("_eres", "");
    if (base == "beam") {
      file = RemixSim::SimRootPath(specs[i]);
      file_eres = RemixSim::IsEresTag(specs[i].tag);
      break;
    }
  }
  if (file.Length() == 0)
    file =
        Paths::DatasetDir() + "/sim_root_files/" + Constants::cfg.SIM_BEAM_FILE;
  if (!file_eres)
    std::cout << "strip-sum-scatter: beam gains from non-eres file " << file
              << std::endl;
  TFile *f = IO::OpenForReading(file);
  if (!f || f->IsZombie()) {
    std::cerr << "strip-sum-scatter: cannot open sim beam file " << file
              << "; sim overlay stays in raw sim units." << std::endl;
    if (f)
      delete f;
    return kFALSE;
  }
  TTree *t = static_cast<TTree *>(f->Get("events_MeV"));
  if (!t) {
    std::cerr << "strip-sum-scatter: no events_MeV tree in sim beam file "
              << file << "; sim overlay stays in raw sim units." << std::endl;
    f->Close();
    delete f;
    return kFALSE;
  }
  Float_t left[18] = {0}, right[18] = {0};
  t->SetBranchAddress("Left_0_17_dE", left);
  t->SetBranchAddress("RightdE", right);
  Long64_t n = t->GetEntries();
  Long64_t stride = FileSet::SampleStride(n, kSampleMaxPoints);
  Double_t sum[18] = {0};
  Long64_t cnt[18] = {0};
  // Unit gains so SimTotal yields the raw per-strip beam total (IGNORE_SHORT
  // aware), which is exactly the quantity these gains will later normalize.
  Double_t unit[18];
  for (Int_t s = 0; s < 18; s++)
    unit[s] = 1.0;
  for (Long64_t j = 0; j < n; j += stride) {
    t->GetEntry(j);
    Double_t total[18];
    SimTotal(left, right, unit, total);
    for (Int_t s = 0; s < 18; s++)
      if (total[s] > 0.0) {
        sum[s] += total[s];
        cnt[s]++;
      }
  }
  f->Close();
  delete f;
  Int_t n_set = 0;
  for (Int_t s = 0; s < 18; s++) {
    if (cnt[s] > 0 && sum[s] > 0.0) {
      gain[s] = 1.0 / (sum[s] / Double_t(cnt[s]));
      n_set++;
    }
  }
  if (n_set == 0)
    return kFALSE;
  std::cout << "strip-sum-scatter: sim per-strip beam normalization to 1 a.u. ("
            << n_set << " strips)." << std::endl;
  return kTRUE;
}

void StripSumScatter::SimTotal(const Float_t *left, const Float_t *right,
                               const Double_t *gain, Double_t *total) {
  for (Int_t s = 0; s < 18; s++)
    total[s] = gain[s] * (Double_t(left[s]) + Double_t(right[s]));
  if (Constants::cfg.IGNORE_SHORT_STRIPS)
    for (Int_t s = 1; s <= 16; s++)
      total[s] =
          gain[s] * ((s % 2) != 0 ? Double_t(left[s]) : Double_t(right[s]));
}

TGraph *StripSumScatter::SimPopScatter(const TString &file, Int_t reac,
                                       const Double_t *gain,
                                       Long64_t max_points) {
  const Int_t kXLo = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_HI;

  TFile *f = IO::OpenForReading(file);
  if (!f || f->IsZombie()) {
    if (f)
      delete f;
    return nullptr;
  }
  TTree *t = static_cast<TTree *>(f->Get("events_MeV"));
  if (!t) {
    std::cerr << "  no events_MeV tree in " << file << std::endl;
    f->Close();
    delete f;
    return nullptr;
  }
  // Same y window as the data scatter (YLoOf/YHiOf), so sim points land in
  // the data's coordinate space instead of a hardcoded, stale one.
  Int_t y_lo = YLoOf(reac);
  Int_t y_hi = YHiOf(reac);
  Float_t left[18] = {0}, right[18] = {0};
  t->SetBranchAddress("Left_0_17_dE", left);
  t->SetBranchAddress("RightdE", right);
  Long64_t n = t->GetEntries();
  Long64_t stride = FileSet::SampleStride(n, max_points);
  TGraph *g = new TGraph();
  Long64_t k = 0;
  for (Long64_t j = 0; j < n; j += stride) {
    t->GetEntry(j);
    Double_t total[18];
    SimTotal(left, right, gain, total);
    Double_t x = SumRange(total, kXLo, kXHi);
    Double_t y = SumRange(total, y_lo, y_hi);
    if (x > 0.0)
      g->SetPoint(k++, x, y);
  }
  g->Set(k);
  f->Close();
  delete f;
  return g;
}

std::vector<TGraph *> StripSumScatter::SimPopTraces(const TString &file,
                                                    const Double_t *gain,
                                                    Long64_t max_traces) {
  std::vector<TGraph *> traces;
  TFile *f = IO::OpenForReading(file);
  if (!f || f->IsZombie()) {
    if (f)
      delete f;
    return traces;
  }
  TTree *t = static_cast<TTree *>(f->Get("events_MeV"));
  if (!t) {
    f->Close();
    delete f;
    return traces;
  }
  Float_t left[18] = {0}, right[18] = {0};
  t->SetBranchAddress("Left_0_17_dE", left);
  t->SetBranchAddress("RightdE", right);
  Long64_t n = t->GetEntries();
  Long64_t stride = FileSet::SampleStride(n, max_traces);
  for (Long64_t j = 0; j < n && Int_t(traces.size()) < max_traces;
       j += stride) {
    t->GetEntry(j);
    Double_t total[18];
    SimTotal(left, right, gain, total);
    traces.push_back(EventsSummary::BuildTraceFromTotals(total));
  }
  f->Close();
  delete f;
  return traces;
}

// Per-strip overlay of sampled sim traces in experimental style (beam grey,
// (a,a') azure, (a,n) red). Beam reference same for every strip; sampled each
// run.
void StripSumScatter::SimTraceOverlay() {
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Int_t kTracesPerRegion =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TRACES_PER_CLASS;

  std::vector<RemixSim::SimFileSpec> specs = RemixSim::BuildFileSpecs();
  if (specs.empty())
    return;
  std::map<Int_t, TString> aa_file, an_file; // reaction strip -> sim file
  std::vector<TString> beam_files;
  for (Int_t i = 0; i < Int_t(specs.size()); i++) {
    TString base = RemixSim::TagWithoutStrip(specs[i].tag);
    base.ReplaceAll("_eres", "");
    TString file = RemixSim::SimRootPath(specs[i]);
    Int_t strip = RemixSim::ReactionStripOf(specs[i].tag);
    if (base == "beam")
      beam_files.push_back(file);
    else if (strip >= kReacMin && strip <= kReacMax) {
      if (base == "aa")
        aa_file[strip] = file;
      else if (base == "an")
        an_file[strip] = file;
    }
  }

  Double_t gain[18];
  if (!SimBeamGains(gain))
    for (Int_t s = 0; s < 18; s++)
      gain[s] = 1.0;

  std::vector<TGraph *> beam_traces;
  for (Int_t i = 0; i < Int_t(beam_files.size()) &&
                    Int_t(beam_traces.size()) < kTracesPerRegion;
       i++) {
    std::vector<TGraph *> t = SimPopTraces(
        beam_files[i], gain, kTracesPerRegion - Int_t(beam_traces.size()));
    for (Int_t k = 0; k < Int_t(t.size()); k++)
      beam_traces.push_back(t[k]);
  }

  for (Int_t r = kReacMin; r <= kReacMax; r++) {
    std::vector<TGraph *> aa_traces, an_traces;
    if (aa_file.find(r) != aa_file.end())
      aa_traces = SimPopTraces(aa_file[r], gain, kTracesPerRegion);
    if (an_file.find(r) != an_file.end())
      an_traces = SimPopTraces(an_file[r], gain, kTracesPerRegion);
    if (aa_traces.empty() && an_traces.empty())
      continue;
    DrawRegionTraces(Form("sim_region_traces_reac%d", r), "sim_scatter",
                     beam_traces, aa_traces, an_traces, 0.6, 1.6,
                     "#DeltaE [a.u.]");
    for (Int_t i = 0; i < Int_t(aa_traces.size()); i++)
      delete aa_traces[i];
    for (Int_t i = 0; i < Int_t(an_traces.size()); i++)
      delete an_traces[i];
  }
  for (Int_t i = 0; i < Int_t(beam_traces.size()); i++)
    delete beam_traces[i];
}

// Fingerprint: per-eres file size+mtime (cheap, no open) + reaction-strip range
// + x window. Regenerating the sim (new mtimes) or changing windows invalidates
// the cache.
TString StripSumScatter::SimFingerprint(
    const std::vector<RemixSim::SimFileSpec> &specs) {
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Int_t kXLo = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_HI;

  TString s = Form("v3 reac[%d,%d] x[%d,%d]", kReacMin, kReacMax, kXLo, kXHi);
  for (Int_t i = 0; i < Int_t(specs.size()); i++) {
    TString f = RemixSim::SimRootPath(specs[i]);
    Long_t id = 0, flags = 0, mtime = 0;
    Long64_t size = -1;
    if (gSystem->GetPathInfo(f, &id, &size, &flags, &mtime) != 0) {
      size = -1;
      mtime = 0;
    }
    s += Form(" %s:%lld:%ld", specs[i].tag.Data(), (long long)size, mtime);
  }
  return s;
}

Bool_t StripSumScatter::LoadSimCache(
    const TString &fp, std::map<Int_t, std::vector<TGraph *>> &by_strip) {
  TString full = IO::GetRootFilesBaseDir() + TString("/") +
                 "StripSumScatter_simcache.root";
  if (gSystem->AccessPathName(full))
    return kFALSE;
  TFile *f = IO::OpenForReading("StripSumScatter_simcache.root");
  if (!f || f->IsZombie()) {
    if (f)
      delete f;
    return kFALSE;
  }
  TNamed *cfp = static_cast<TNamed *>(f->Get("sim_fingerprint"));
  if (!cfp || fp != cfp->GetTitle()) {
    f->Close();
    delete f;
    return kFALSE;
  }
  TIter next(f->GetListOfKeys());
  TKey *key;
  while ((key = static_cast<TKey *>(next()))) {
    TString name = key->GetName();
    if (!name.BeginsWith("simg_r"))
      continue;
    TString rest = name(6, name.Length() - 6); // after "simg_r": <strip>_p<idx>
    Int_t us = rest.Index("_p");
    if (us < 0)
      continue;
    Int_t r = TString(rest(0, us)).Atoi();
    TGraph *g = static_cast<TGraph *>(f->Get(name));
    if (!g)
      continue;
    by_strip[r].push_back(static_cast<TGraph *>(g->Clone()));
  }
  f->Close();
  delete f;
  return kTRUE;
}

void StripSumScatter::WriteSimCache(
    const TString &fp, const std::map<Int_t, std::vector<TGraph *>> &by_strip) {
  TFile *out = IO::OpenForWriting("StripSumScatter_simcache.root", "RECREATE");
  if (!out || out->IsZombie()) {
    if (out)
      delete out;
    return;
  }
  out->cd();
  TNamed cfp("sim_fingerprint", fp.Data());
  cfp.Write();
  std::map<Int_t, std::vector<TGraph *>>::const_iterator it;
  for (it = by_strip.begin(); it != by_strip.end(); ++it)
    for (Int_t i = 0; i < Int_t(it->second.size()); i++)
      it->second[i]->Write(Form("simg_r%d_p%d", it->first, Int_t(i)));
  out->Close();
  delete out;
}

// Sim overlay: each sim population as a labelled cloud on its reaction
// strip's data-scatter axes (beam overlays every strip); fingerprint-cached, so
// re-runs reload instead of rescanning sim files.
void StripSumScatter::SimOverlay() {
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Int_t kXLo = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_HI;

  std::vector<RemixSim::SimFileSpec> specs = RemixSim::BuildFileSpecs();
  if (specs.empty()) {
    std::cerr
        << "strip-sum-scatter: no sim control files; skipping sim overlay."
        << std::endl;
    return;
  }
  TString fp = SimFingerprint(specs);

  std::map<Int_t, std::vector<TGraph *>>
      by_strip; // strip -> graphs (title=label)
  Bool_t loaded = LoadSimCache(fp, by_strip);
  if (!loaded) {
    std::map<Int_t, std::vector<SimPop>> reacted;
    std::vector<SimPop> refs;
    for (Int_t i = 0; i < Int_t(specs.size()); i++) {
      SimPop p;
      p.file = RemixSim::SimRootPath(specs[i]);
      p.label = PrettyLabel(specs[i].tag);
      Int_t strip = RemixSim::ReactionStripOf(specs[i].tag);
      if (strip < 0)
        refs.push_back(p);
      else
        reacted[strip].push_back(p);
    }
    Double_t gain[18];
    if (!SimBeamGains(gain))
      for (Int_t s = 0; s < 18; s++)
        gain[s] = 1.0;
    const Long64_t kSimMaxPoints = 25000;
    for (Int_t r = kReacMin; r <= kReacMax; r++) {
      std::vector<SimPop> group = reacted[r];
      for (Int_t i = 0; i < Int_t(refs.size()); i++)
        group.push_back(refs[i]);
      for (Int_t i = 0; i < Int_t(group.size()); i++) {
        TGraph *g = SimPopScatter(group[i].file, r, gain, kSimMaxPoints);
        if (!g || g->GetN() == 0) {
          if (g)
            delete g;
          continue;
        }
        g->SetTitle(group[i].label);
        by_strip[r].push_back(g);
      }
    }
    Int_t n_graphs = 0;
    std::map<Int_t, std::vector<TGraph *>>::const_iterator cit;
    for (cit = by_strip.begin(); cit != by_strip.end(); ++cit)
      n_graphs += Int_t(cit->second.size());
    if (n_graphs == 0) {
      std::cerr << "strip-sum-scatter: no sim data found (regenerate "
                   "sim_root_files); skipping sim overlay."
                << std::endl;
      return;
    }
    WriteSimCache(fp, by_strip);
    std::cout << "strip-sum-scatter: built + cached sim overlay (" << n_graphs
              << " population graphs)." << std::endl;
  } else {
    std::cout
        << "strip-sum-scatter: loaded cached sim overlay (fingerprint match)."
        << std::endl;
  }

  std::map<Int_t, std::vector<TGraph *>>::iterator it;
  const Double_t kXMin = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.XMIN;
  const Double_t kXMax = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.XMAX;
  for (it = by_strip.begin(); it != by_strip.end(); ++it) {
    Int_t r = it->first;
    std::map<Int_t, TH2F *>::const_iterator sit = m_scatter.find(r);
    if (sit == m_scatter.end())
      continue;
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    TH2F *ref = sit->second;
    // Frame the sim overlay with the same display window as the data scatter
    // (x XMIN/XMAX, y Y_RANGE/YMIN/YMAX): the data hists use the wide
    // ScatterBuildRange, so their axis extents don't match the plotted window.
    Double_t y_lo[64], y_hi[64];
    YBounds(y_lo, y_hi);
    Int_t ri = ReacIndex(r);
    TH2F *frame = new TH2F(Form("sim_frame_r%d", r), "", 10, kXMin, kXMax, 10,
                           y_lo[ri], y_hi[ri]);
    frame->SetStats(0);
    frame->GetXaxis()->SetTitle(ref->GetXaxis()->GetTitle());
    frame->GetYaxis()->SetTitle(ref->GetYaxis()->GetTitle());
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    c->SetLeftMargin(0.18);
    frame->Draw();
    // Match the experimental region-traces legend placement (top-right).
    TLegend *leg = PlottingUtils::AddLegend(0.725, 0.875, 0.70, 0.86);
    for (Int_t i = 0; i < Int_t(it->second.size()); i++) {
      TGraph *g = it->second[i];
      // Match the experimental region-trace colours (DrawRegionTraces): beam
      // grey, (a,a') azure, (a,n) red -- keyed off the population label.
      TString lab = g->GetTitle();
      Int_t color = kBlack;
      if (lab == "Beam")
        color = kGray + 2;
      else if (lab == "(#alpha,#alpha')")
        color = kAzure + 2;
      else if (lab == "(#alpha,n)")
        color = kRed + 1;
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.3);
      g->SetMarkerColorAlpha(color, 0.35);
      g->SetLineColor(color);
      g->Draw("P SAME");
      leg->AddEntry(g, g->GetTitle(), "p");
    }
    leg->Draw();
    PlottingUtils::SaveFigure(c,
                              Form("sim_normsumE_reac%d_s%d_%d_vs_s%d_%d", r,
                                   YLoOf(r), YHiOf(r), kXLo, kXHi),
                              "sim_scatter", PlotSaveOptions::kLINEAR);
    delete leg;
    delete c;
    delete frame;
  }

  std::map<Int_t, std::vector<TGraph *>>::iterator dit;
  for (dit = by_strip.begin(); dit != by_strip.end(); ++dit)
    for (Int_t i = 0; i < Int_t(dit->second.size()); i++)
      delete dit->second[i];
}
