#include "TemplateMatch.hpp"
#include "Constants.hpp"
#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "Paths.hpp"
#include "PlottingUtils.hpp"
#include "RegionCuts.hpp"
#include "StripSumScatter.hpp"
#include "TagEfficiency.hpp"
#include <TCanvas.h>
#include <TCutG.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TTree.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <mutex>
#include <thread>

// ---------------------------------------------------------------------------
// The classifier and the data pass live on StripSumScatter, so they share the
// gates, the event-level cuts and the noise with the tag.
// ---------------------------------------------------------------------------

Int_t StripSumScatter::TemplateClassify(
    const Double_t *total, const std::vector<std::vector<Double_t>> &templates,
    Double_t delta_chi2, Double_t *chi2_beam_out, Double_t *chi2_best_out) {
  const Int_t s_lo = Constants::cfg.IGNORE_STRIP_0 ? 1 : 0;
  const Int_t s_hi = Constants::cfg.IGNORE_STRIP_17 ? 16 : 17;
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  Double_t chi2_beam = 0.0;
  for (Int_t s = s_lo; s <= s_hi; s++) {
    const Double_t sig = StripSigma(s);
    if (!(sig > 0.0))
      continue;
    const Double_t d = (total[s] - 1.0) / sig;
    chi2_beam += d * d;
  }
  Int_t best = -1;
  Double_t best_chi2 = 0.0;
  for (Int_t k = 0; k < Int_t(templates.size()); k++) {
    const std::vector<Double_t> &t = templates[k];
    if (t.empty())
      continue;
    Double_t chi2 = 0.0;
    for (Int_t s = s_lo; s <= s_hi; s++) {
      const Double_t sig = StripSigma(s);
      if (!(sig > 0.0))
        continue;
      const Double_t d = (total[s] - t[s]) / sig;
      chi2 += d * d;
    }
    if (best < 0 || chi2 < best_chi2) {
      best = k;
      best_chi2 = chi2;
    }
  }
  if (chi2_beam_out)
    *chi2_beam_out = chi2_beam;
  if (chi2_best_out)
    *chi2_best_out = best < 0 ? chi2_beam : best_chi2;
  if (best < 0 || chi2_beam - best_chi2 < delta_chi2)
    return -1;
  return kReacMin + best;
}

StripSumScatter::TemplatePassResult StripSumScatter::TemplateMatchPass(
    const std::vector<Int_t> &runOrder, std::map<Int_t, TChain *> &chains,
    const std::vector<std::vector<Double_t>> &templates, Double_t delta_chi2) {
  const StripSumScatterConfig &C = Constants::cfg.STRIP_SUM_SCATTER_CONFIG;
  const Int_t nReac = C.REACTION_STRIP_MAX - C.REACTION_STRIP_MIN + 1;
  TemplatePassResult total;
  total.matched.assign(nReac, 0);
  total.matched_and_tagged.assign(nReac, 0);

  std::vector<GateSpec> activeGates = ActiveGates();
  const Int_t nRuns = Int_t(runOrder.size());
  std::vector<TChain *> chainVec(nRuns);
  for (Int_t i = 0; i < nRuns; i++)
    chainVec[i] = chains[runOrder[i]];
  Int_t n_workers =
      TMath::Min(Int_t(std::thread::hardware_concurrency()), nRuns);
  n_workers = TMath::Min(n_workers, C.MAX_STRIP_SUM_WORKERS);
  if (n_workers < 1)
    n_workers = 1;

  // The same per-run gates the fill used, refitted here rather than read
  // back, so this pass never depends on what the cache happened to store.
  std::vector<SingleRunFitResult> fits(nRuns);
  RunIndexedParallel(nRuns, n_workers, [&](Int_t i) {
    fits[i] = FitRunGates(runOrder[i], chainVec[i], activeGates);
  });

  struct Task {
    Int_t run_idx;
    TString path;
  };
  std::vector<Task> tasks;
  {
    std::map<Int_t, Int_t> idx_of_run;
    for (Int_t i = 0; i < nRuns; i++)
      idx_of_run[runOrder[i]] = i;
    std::vector<FileSpec> specs = FileSet::BuildProcessedFileSpecs();
    for (Int_t k = 0; k < Int_t(specs.size()); k++) {
      std::map<Int_t, Int_t>::const_iterator it = idx_of_run.find(specs[k].run);
      if (it == idx_of_run.end() || !fits[it->second].ok)
        continue;
      TString full = IO::GetRootFilesBaseDir() + "/" +
                     FileSet::EventsName(specs[k]) + ".root";
      if (gSystem->AccessPathName(full))
        continue;
      Task t;
      t.run_idx = it->second;
      t.path = full;
      tasks.push_back(t);
    }
  }
  const Int_t nTasks = Int_t(tasks.size());
  Int_t workers =
      TMath::Min(Int_t(std::thread::hardware_concurrency()), nTasks);
  workers = TMath::Min(workers, C.MAX_STRIP_SUM_WORKERS);
  if (workers < 1)
    workers = 1;
  std::cout << "template-match: classifying " << nTasks << " files on "
            << workers << " workers, delta chi2 " << delta_chi2 << std::endl;

  std::vector<TemplatePassResult> parts(nTasks);
  std::mutex log_mutex;
  RunIndexedParallel(nTasks, workers, [&](Int_t t) {
    TemplatePassResult &out = parts[t];
    out.matched.assign(nReac, 0);
    out.matched_and_tagged.assign(nReac, 0);
    const Int_t i = tasks[t].run_idx;
    const std::vector<BeamFit2D> &runGates = fits[i].series_gates;
    TChain ch("events");
    ch.Add(tasks[t].path);
    if (ch.GetEntries() == 0)
      return;
    EnergyView ev;
    ev.Attach(&ch);
    EnableEventBranches(&ch);
    const Long64_t n = ch.GetEntries();
    for (Long64_t j = 0; j < n; j++) {
      ch.GetEntry(j);
      ev.Decode();
      out.seen++;
      Bool_t passesAll = kTRUE;
      for (Int_t gi = 0; gi < Int_t(activeGates.size()); gi++)
        if (!PassesGate(runGates[gi], ev, activeGates[gi].sx,
                        activeGates[gi].sy)) {
          passesAll = kFALSE;
          break;
        }
      if (!passesAll || IsPileup(ev) || IsNoise(ev))
        continue;
      if (C.REJECT_OFFBEAM && IsOffbeam(ev))
        continue;
      if (IsParityAsymmetric(ev))
        continue;
      if (C.BOTH_MULT_MAX >= 0) {
        const Int_t hi = TMath::Min(16, C.BOTH_MULT_COUNT_TO);
        Int_t nboth = 0;
        for (Int_t s = 1; s <= hi; s++)
          if (ev.left_0_17_adc[s] > 0.0 && ev.rightdE_adc[s] > 0.0)
            nboth++;
        if (nboth > C.BOTH_MULT_MAX)
          continue;
      }
      if (!AllStripsFired(ev))
        continue;
      out.considered++;
      const Int_t r = TemplateClassify(ev.total, templates, delta_chi2);
      if (r < 0)
        continue;
      out.matched[ReacIndex(r)]++;
      if (RejectReason(ev, r) == kTagPass)
        out.matched_and_tagged[ReacIndex(r)]++;
    }
    std::lock_guard<std::mutex> lk(log_mutex);
    std::cout << "  " << gSystem->BaseName(tasks[t].path) << ": " << out.seen
              << " events, " << out.considered << " considered" << std::endl;
  });
  for (Int_t t = 0; t < nTasks; t++) {
    total.seen += parts[t].seen;
    total.considered += parts[t].considered;
    for (Int_t k = 0; k < nReac; k++) {
      total.matched[k] += parts[t].matched[k];
      total.matched_and_tagged[k] += parts[t].matched_and_tagged[k];
    }
  }
  return total;
}

// ---------------------------------------------------------------------------
// The driver
// ---------------------------------------------------------------------------

namespace {

Bool_t LoadSigmas(TFile &cache) {
  Double_t jump[18] = {0}, strip[18] = {0};
  Int_t found = 0;
  for (Int_t s = 0; s < 18; s++) {
    if (TParameter<Double_t> *p = static_cast<TParameter<Double_t> *>(
            cache.Get(Form("jump_sigma_s%d", s)))) {
      jump[s] = p->GetVal();
      found++;
    }
    if (TParameter<Double_t> *p = static_cast<TParameter<Double_t> *>(
            cache.Get(Form("strip_sigma_s%d", s))))
      strip[s] = p->GetVal();
  }
  StripSumScatter::SetJumpSigma(jump);
  StripSumScatter::SetStripSigma(strip);
  return found > 0;
}

Long64_t ReadCount(TFile &f, const TString &key) {
  TParameter<Long64_t> *p = dynamic_cast<TParameter<Long64_t> *>(f.Get(key));
  return p ? p->GetVal() : -1;
}

} // namespace

void TemplateMatch::DrawTemplate(Int_t reac, const std::vector<Double_t> &tpl,
                                 Double_t n_events) {
  const Int_t s_lo = Constants::cfg.IGNORE_STRIP_0 ? 1 : 0;
  const Int_t s_hi = Constants::cfg.IGNORE_STRIP_17 ? 16 : 17;
  TH2F *frame =
      new TH2F(Form("h_tpl_frame_%d", reac), ";Strip;#DeltaE [a.u.]",
               s_hi - s_lo + 1, s_lo - 0.5, s_hi + 0.5, 100, 0.4, 1.6);
  frame->SetStats(0);
  frame->SetDirectory(nullptr);
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  frame->Draw();
  TGraphErrors *band = new TGraphErrors(s_hi - s_lo + 1);
  TGraph *line = new TGraph(s_hi - s_lo + 1);
  TGraph *beam = new TGraph(s_hi - s_lo + 1);
  for (Int_t s = s_lo; s <= s_hi; s++) {
    band->SetPoint(s - s_lo, s, tpl[s]);
    band->SetPointError(s - s_lo, 0.0, StripSumScatter::StripSigma(s));
    line->SetPoint(s - s_lo, s, tpl[s]);
    beam->SetPoint(s - s_lo, s, 1.0);
  }
  band->SetFillColorAlpha(kRed - 7, 0.35);
  band->Draw("3 SAME");
  beam->SetLineColor(kGray + 3);
  beam->SetLineWidth(2);
  beam->SetLineStyle(2);
  beam->Draw("L SAME");
  line->SetLineColor(kRed + 2);
  line->SetLineWidth(3);
  line->Draw("L SAME");
  TLegend *leg = PlottingUtils::AddLegend(0.40, 0.875, 0.72, 0.86);
  leg->AddEntry(line, Form("(#alpha,n) template, %.0f events", n_events), "l");
  leg->AddEntry(band, "#pm1#sigma beam width per strip", "f");
  leg->AddEntry(beam, "Beam", "l");
  leg->Draw();
  PlottingUtils::AddText(Form("reaction strip %d", reac), 0.875, 0.68);
  PlottingUtils::SaveFigure(c, Form("template_reac%d", reac), "template_match",
                            PlotSaveOptions::kLINEAR);
  delete c;
  delete frame;
}

void TemplateMatch::WriteReport(Long64_t seen, Long64_t considered,
                                const char *tag_part) {
  const StripSumScatterConfig &C = Constants::cfg.STRIP_SUM_SCATTER_CONFIG;
  TString out;
  out += Form("template-match: %s, %lld events seen, %lld past the "
              "event-level cuts, delta chi2 %.1f, %lld bootstrap trials\n",
              Paths::DatasetName().Data(), seen, considered,
              C.TEMPLATE_DELTA_CHI2, C.TEMPLATE_BOOTSTRAP_TRIALS);
  out += Form("tag: %s\n\n", tag_part);
  out += "Both methods count the same reaction. 'tagged' is the tag's count "
         "inside the (a,n) region and 'corr tag' divides it by the bootstrap "
         "tag efficiency (1 when no tag-efficiency store). 'matched' is the "
         "template classifier's count, 'both' how many of those the tag also "
         "took, 'fp' the flat-beam false-positive fraction, 'bkg' that times "
         "the beam count at the strip, and 'corr tpl' (matched - bkg) / eff. "
         "Both efficiencies come from the same Gaussian resampling of the "
         "template, so neither knows about correlated noise or the spread of "
         "reaction vertices within a strip.\n\n";
  out += Form("%4s %9s %9s %7s %10s %9s %8s %7s %7s %8s %9s %10s %7s\n", "reac",
              "tpl evts", "tagged", "eff tag", "corr tag", "matched", "both",
              "eff tpl", "mig tpl", "fp", "bkg", "corr tpl", "ratio");
  for (Int_t i = 0; i < Int_t(results_.size()); i++) {
    const TemplateStripResult &r = results_[i];
    const Double_t ratio =
        r.corrected_tag > 0.0 ? r.corrected_template / r.corrected_tag : 0.0;
    out += Form("%4d %9.0f %9lld %7.3f %10.1f %9lld %8lld %7.3f %7.3f %8.1e "
                "%9.1f %10.1f %7.3f\n",
                r.reac, r.n_template_events, r.n_tagged, r.eff_tag,
                r.corrected_tag, r.n_matched, r.n_both, r.eff_template,
                r.migrate_template, r.false_positive, r.background,
                r.corrected_template, ratio);
  }
  std::cout << out;
  const TString dir = Paths::ResultsDir() + "/plots/template_match";
  gSystem->mkdir(dir, kTRUE);
  const TString path = dir + "/template_match.txt";
  std::ofstream f(path.Data());
  if (f) {
    f << out;
    f.close();
    std::cout << "template-match: report written to " << path << std::endl;
  } else {
    std::cerr << "template-match: cannot write " << path << std::endl;
  }
  // The numbers, for anything that wants to consume them later.
  const TString store = Paths::ResultsDir() + "/root_files/template_match.root";
  TFile rf(store, "RECREATE");
  if (!rf.IsZombie()) {
    TNamed("tag_part", tag_part).Write();
    for (Int_t i = 0; i < Int_t(results_.size()); i++) {
      const TemplateStripResult &r = results_[i];
      TParameter<Long64_t>(Form("n_matched_r%d", r.reac), r.n_matched).Write();
      TParameter<Long64_t>(Form("n_both_r%d", r.reac), r.n_both).Write();
      TParameter<Double_t>(Form("eff_template_r%d", r.reac), r.eff_template)
          .Write();
      TParameter<Double_t>(Form("migrate_template_r%d", r.reac),
                           r.migrate_template)
          .Write();
      TParameter<Double_t>(Form("false_positive_r%d", r.reac), r.false_positive)
          .Write();
      TParameter<Double_t>(Form("corrected_template_r%d", r.reac),
                           r.corrected_template)
          .Write();
      TParameter<Double_t>(Form("corrected_tag_r%d", r.reac), r.corrected_tag)
          .Write();
    }
    rf.Close();
  }
}

Bool_t TemplateMatch::Run() {
  const StripSumScatterConfig &C = Constants::cfg.STRIP_SUM_SCATTER_CONFIG;
  const CrossSectionConfig &X = Constants::cfg.CROSS_SECTION_CONFIG;
  const Int_t kReacMin = C.REACTION_STRIP_MIN;
  const Int_t kReacMax = C.REACTION_STRIP_MAX;
  const Int_t nReac = kReacMax - kReacMin + 1;
  ROOT::EnableThreadSafety();

  const TString cache_path =
      IO::GetRootFilesBaseDir() + "/" + StripSumScatter::CacheName();
  TFile cache(cache_path, "READ");
  if (cache.IsZombie()) {
    std::cerr << "template-match: no scatter cache at " << cache_path
              << "; run strip-sum-scatter first" << std::endl;
    return kFALSE;
  }
  if (!LoadSigmas(cache)) {
    std::cerr << "template-match: cache carries no per-strip noise sigmas; "
                 "rebuild it"
              << std::endl;
    return kFALSE;
  }
  TTree *reservoir = static_cast<TTree *>(cache.Get("traces"));
  if (!reservoir) {
    std::cerr << "template-match: no reservoir in the cache" << std::endl;
    return kFALSE;
  }
  TString tag_part = "(unknown)";
  if (TNamed *fp = static_cast<TNamed *>(cache.Get("fingerprint"))) {
    tag_part = fp->GetTitle();
    const Ssiz_t bar = tag_part.Index(" | ");
    if (bar >= 0)
      tag_part = tag_part(0, bar);
  }

  // Templates: the mean trace of the events tagged at each strip and inside
  // its (a,n) region. Without a region for a strip, every event tagged there.
  std::vector<TCutG *> cuts(nReac, nullptr);
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++)
    cuts[reac - kReacMin] = RegionCutStore::Load("region_an", reac);
  std::vector<std::vector<Double_t>> sum(nReac, std::vector<Double_t>(18, 0.0));
  std::vector<Double_t> cnt(nReac, 0.0);
  {
    Float_t total[18];
    UInt_t mask = 0;
    reservoir->SetBranchStatus("*", 0);
    reservoir->SetBranchStatus("total", 1);
    reservoir->SetBranchStatus("reac_mask", 1);
    reservoir->SetBranchAddress("total", total);
    reservoir->SetBranchAddress("reac_mask", &mask);
    const Long64_t n = reservoir->GetEntries();
    for (Long64_t i = 0; i < n; i++) {
      reservoir->GetEntry(i);
      if (mask == 0)
        continue;
      Double_t td[18];
      for (Int_t s = 0; s < 18; s++)
        td[s] = total[s];
      for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
        const Int_t k = reac - kReacMin;
        if (!(mask & (1u << k)))
          continue;
        if (cuts[k]) {
          Double_t x, y;
          StripSumScatter::PlaneXY(td, reac, x, y);
          if (!cuts[k]->IsInside(x, y))
            continue;
        }
        for (Int_t s = 0; s < 18; s++)
          sum[k][s] += td[s];
        cnt[k] += 1.0;
      }
    }
    reservoir->SetBranchStatus("*", 1);
  }
  templates_.assign(nReac, std::vector<Double_t>());
  Int_t nTemplates = 0;
  for (Int_t k = 0; k < nReac; k++) {
    if (cnt[k] < Double_t(C.TEMPLATE_MIN_EVENTS)) {
      std::cout << Form("template-match: reac %2d: %.0f events, below the "
                        "%d needed; no template",
                        kReacMin + k, cnt[k], C.TEMPLATE_MIN_EVENTS)
                << std::endl;
      continue;
    }
    templates_[k].assign(18, 0.0);
    for (Int_t s = 0; s < 18; s++)
      templates_[k][s] = sum[k][s] / cnt[k];
    nTemplates++;
    std::cout << Form("template-match: reac %2d: template from %.0f events "
                      "(%s)",
                      kReacMin + k, cnt[k],
                      cuts[k] ? "inside region_an" : "all tagged, no region")
              << std::endl;
    DrawTemplate(kReacMin + k, templates_[k], cnt[k]);
  }
  if (nTemplates == 0) {
    std::cerr << "template-match: no strip has enough tagged events for a "
                 "template"
              << std::endl;
    return kFALSE;
  }

  // Bootstrap: each template resampled with the beam widths and classified;
  // a flat trace resampled the same way for the false-positive rate.
  const Long64_t kTrials = C.TEMPLATE_BOOTSTRAP_TRIALS;
  std::vector<Double_t> eff(nReac, 0.0), mig(nReac, 0.0), fp(nReac, 0.0);
  {
    TRandom3 rng(7);
    Double_t trace[18];
    for (Int_t k = 0; k < nReac; k++) {
      if (templates_[k].empty())
        continue;
      Long64_t here = 0, next = 0;
      for (Long64_t t = 0; t < kTrials; t++) {
        for (Int_t s = 0; s < 18; s++)
          trace[s] =
              templates_[k][s] + rng.Gaus(0.0, StripSumScatter::StripSigma(s));
        const Int_t r = StripSumScatter::TemplateClassify(
            trace, templates_, C.TEMPLATE_DELTA_CHI2);
        if (r == kReacMin + k)
          here++;
        else if (r == kReacMin + k + 1)
          next++;
      }
      eff[k] = Double_t(here) / Double_t(kTrials);
      mig[k] = Double_t(next) / Double_t(kTrials);
    }
    const Long64_t kBeamTrials = C.TEMPLATE_BEAM_TRIALS;
    std::vector<Long64_t> fp_n(nReac, 0);
    for (Long64_t t = 0; t < kBeamTrials; t++) {
      for (Int_t s = 0; s < 18; s++)
        trace[s] = 1.0 + rng.Gaus(0.0, StripSumScatter::StripSigma(s));
      const Int_t r = StripSumScatter::TemplateClassify(trace, templates_,
                                                        C.TEMPLATE_DELTA_CHI2);
      if (r >= 0)
        fp_n[r - kReacMin]++;
    }
    for (Int_t k = 0; k < nReac; k++)
      fp[k] = Double_t(fp_n[k]) / Double_t(kBeamTrials);
  }

  // The data pass.
  std::vector<Int_t> run_order;
  std::map<Int_t, TChain *> chain_by_run = FileSet::GroupEventsByRun(run_order);
  if (run_order.empty()) {
    std::cerr << "template-match: no runs found" << std::endl;
    return kFALSE;
  }
  StripSumScatter sss;
  const StripSumScatter::TemplatePassResult pass = sss.TemplateMatchPass(
      run_order, chain_by_run, templates_, C.TEMPLATE_DELTA_CHI2);

  // The comparison.
  const TString channel =
      X.CHANNELS.empty() ? TString("an") : X.CHANNELS[0].name;
  results_.clear();
  for (Int_t k = 0; k < nReac; k++) {
    if (templates_[k].empty())
      continue;
    const Int_t reac = kReacMin + k;
    TemplateStripResult r;
    r.reac = reac;
    r.n_template_events = cnt[k];
    r.n_tagged = ReadCount(cache, Form("n_tagged_r%d", reac));
    TagEfficiencyRecord er;
    if (TagEfficiencyStore::Load(channel, reac, er) && er.eff > 0.0)
      r.eff_tag = er.eff;
    r.corrected_tag = r.eff_tag > 0.0 ? cnt[k] / r.eff_tag : cnt[k];
    r.n_matched = pass.matched[k];
    r.n_both = pass.matched_and_tagged[k];
    r.eff_template = eff[k];
    r.migrate_template = mig[k];
    r.false_positive = fp[k];
    const Long64_t n_beam_at = ReadCount(cache, Form("n_normed_r%d", reac));
    r.background = fp[k] * Double_t(n_beam_at > 0 ? n_beam_at : 0);
    r.corrected_template =
        eff[k] > 0.0 ? (Double_t(r.n_matched) - r.background) / eff[k] : 0.0;
    results_.push_back(r);
  }
  for (Int_t k = 0; k < nReac; k++)
    delete cuts[k];
  WriteReport(pass.seen, pass.considered, tag_part.Data());
  return kTRUE;
}
