#include "TagEfficiency.hpp"
#include "Constants.hpp"
#include "InitUtils.hpp"
#include "Normalization.hpp"
#include "Paths.hpp"
#include "PlottingUtils.hpp"
#include "StripSumScatter.hpp"
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TTree.h>
#include <cmath>
#include <iostream>

// ---------------------------------------------------------------------------
// Store
// ---------------------------------------------------------------------------
namespace TagEfficiencyStore {

TString Path() {
  return Paths::ResultsDir() + "/root_files/tag_efficiency.root";
}

static const char *kFields[5] = {"n_counted", "eff", "eff_err", "migrate",
                                 "tag_eff"};

static TString Key(const TString &channel, const char *field, Int_t reac) {
  return Form("%s_%s_r%d", channel.Data(), field, reac);
}

void Write(const TString &channel,
           const std::vector<TagEfficiencyRecord> &records,
           const TString &method) {
  TFile f(Path(), "UPDATE");
  if (f.IsZombie()) {
    std::cerr << "tag-efficiency: cannot write " << Path() << std::endl;
    return;
  }
  f.cd();
  TNamed("method", method.Data()).Write("method", TObject::kOverwrite);
  // Drop the channel's previous records so a strip that lost its cut does
  // not keep a stale one.
  for (TIter it(f.GetListOfKeys()); TObject *k = it();)
    if (TString(k->GetName()).BeginsWith(channel + "_"))
      f.Delete(Form("%s;*", k->GetName()));
  for (Int_t i = 0; i < Int_t(records.size()); i++) {
    const TagEfficiencyRecord &r = records[i];
    const Double_t v[5] = {r.n_counted, r.eff, r.eff_err, r.migrate, r.tag_eff};
    for (Int_t k = 0; k < 5; k++) {
      const TString key = Key(channel, kFields[k], r.reac);
      TParameter<Double_t>(key, v[k]).Write(key, TObject::kOverwrite);
    }
  }
  f.Close();
}

Bool_t Load(const TString &channel, Int_t reac, TagEfficiencyRecord &record) {
  if (gSystem->AccessPathName(Path()))
    return kFALSE;
  TFile f(Path(), "READ");
  if (f.IsZombie())
    return kFALSE;
  Double_t v[5];
  for (Int_t k = 0; k < 5; k++) {
    TParameter<Double_t> *p = dynamic_cast<TParameter<Double_t> *>(
        f.Get(Key(channel, kFields[k], reac)));
    if (!p)
      return kFALSE;
    v[k] = p->GetVal();
  }
  record.reac = reac;
  record.n_counted = v[0];
  record.eff = v[1];
  record.eff_err = v[2];
  record.migrate = v[3];
  record.tag_eff = v[4];
  return kTRUE;
}

TString Method() {
  if (gSystem->AccessPathName(Path()))
    return "";
  TFile f(Path(), "READ");
  TNamed *m = f.IsZombie() ? nullptr : dynamic_cast<TNamed *>(f.Get("method"));
  return m ? TString(m->GetTitle()) : TString("");
}

} // namespace TagEfficiencyStore

// ---------------------------------------------------------------------------
// Bootstrap method
// ---------------------------------------------------------------------------
namespace {

const Long64_t kBootstrapTrials = 20000;
const Double_t kMinEventsInCut = 20.0;

} // namespace

Bool_t TagEfficiency::LoadSigmas(TFile &cache) {
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

// Strips with a hand-drawn cut for this channel, one past XS_STRIP_MAX for
// the migration out of the last strip. The fitted ellipses are not used: until
// they reproduce the drawn cuts they also sweep in the beam tail under the
// island, which outnumbers the reaction a hundred to one at the low strips and
// drags the mean onto the ridge.
void TagEfficiency::LoadDrawnCuts(const TString &region) {
  const CrossSectionConfig &X = Constants::cfg.CROSS_SECTION_CONFIG;
  strips_.clear();
  for (Int_t reac = X.XS_STRIP_MIN; reac <= X.XS_STRIP_MAX + 1; reac++) {
    Strip st;
    st.cut = RegionCutStore::LoadDrawn(region, reac);
    if (st.cut)
      strips_[reac] = st;
  }
}

// Mean trace of the events inside each strip's drawn cut, one pass over the
// reservoir.
void TagEfficiency::MeanTraces(TTree *tt) {
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  Float_t total[18];
  UInt_t mask = 0;
  tt->SetBranchStatus("*", 0);
  tt->SetBranchStatus("total", 1);
  tt->SetBranchStatus("reac_mask", 1);
  tt->SetBranchAddress("total", total);
  tt->SetBranchAddress("reac_mask", &mask);
  for (std::map<Int_t, Strip>::iterator it = strips_.begin();
       it != strips_.end(); ++it) {
    it->second.mean.assign(18, 0.0);
    it->second.count = 0.0;
  }
  const Long64_t n = tt->GetEntries();
  for (Long64_t i = 0; i < n; i++) {
    tt->GetEntry(i);
    if (mask == 0)
      continue;
    Double_t td[18];
    for (Int_t s = 0; s < 18; s++)
      td[s] = total[s];
    for (std::map<Int_t, Strip>::iterator it = strips_.begin();
         it != strips_.end(); ++it) {
      const Int_t reac = it->first;
      if (!(mask & (1u << (reac - kReacMin))))
        continue;
      Double_t x, y;
      StripSumScatter::PlaneXY(td, reac, x, y);
      if (!it->second.cut->IsInside(x, y))
        continue;
      for (Int_t s = 0; s < 18; s++)
        it->second.mean[s] += td[s];
      it->second.count += 1.0;
    }
  }
  for (std::map<Int_t, Strip>::iterator it = strips_.begin();
       it != strips_.end(); ++it)
    if (it->second.count > 0.0)
      for (Int_t s = 0; s < 18; s++)
        it->second.mean[s] /= it->second.count;
  tt->SetBranchStatus("*", 1);
}

// The mean trace resampled with the measured beam widths, each resample
// through the tag and the cut at reac, and through the next strip's for the
// migration.
TagEfficiency::Outcome TagEfficiency::Bootstrap(Int_t reac, TH2F *plane) {
  const Strip &st = strips_[reac];
  TCutG *cut_next = strips_.count(reac + 1) ? strips_[reac + 1].cut : nullptr;
  TRandom3 rng(12345 + reac);
  Outcome o;
  for (Long64_t t = 0; t < kBootstrapTrials; t++) {
    Double_t total[18];
    for (Int_t s = 0; s < 18; s++)
      total[s] = st.mean[s] + rng.Gaus(0.0, StripSumScatter::StripSigma(s));
    EnergyView ev;
    for (Int_t s = 0; s < 18; s++)
      ev.total[s] = total[s];
    o.n++;
    Double_t x, y;
    StripSumScatter::PlaneXY(total, reac, x, y);
    if (plane)
      plane->Fill(x, y);
    if (StripSumScatter::PassesReaction(ev, reac)) {
      o.tagged++;
      if (st.cut->IsInside(x, y))
        o.inside++;
    }
    if (cut_next && StripSumScatter::PassesReaction(ev, reac + 1)) {
      Double_t xn, yn;
      StripSumScatter::PlaneXY(total, reac + 1, xn, yn);
      if (cut_next->IsInside(xn, yn))
        o.next_inside++;
    }
  }
  return o;
}

void TagEfficiency::DrawMeanTrace(const TString &subdir, Int_t reac,
                                  const Strip &st) {
  TH2F *frame =
      new TH2F(Form("h_eff_frame_%s_%d", subdir.Data(), reac),
               ";Strip;#DeltaE [a.u.]", 18, -0.5, 17.5, 100, 0.6, 1.6);
  frame->SetStats(0);
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  frame->Draw();
  TGraphErrors *band = new TGraphErrors(18);
  TGraph *line = new TGraph(18);
  for (Int_t s = 0; s < 18; s++) {
    band->SetPoint(s, s, st.mean[s]);
    band->SetPointError(s, 0.0, StripSumScatter::StripSigma(s));
    line->SetPoint(s, s, st.mean[s]);
  }
  band->SetFillColorAlpha(kRed + 1, 0.25);
  band->SetLineColor(kRed + 1);
  band->Draw("3 SAME");
  line->SetLineColor(kRed + 1);
  line->SetLineWidth(2);
  line->Draw("L SAME");
  TLegend *leg = PlottingUtils::AddLegend(0.40, 0.875, 0.72, 0.86);
  leg->AddEntry(
      line, Form("Mean trace, %.0f events in the drawn cut", st.count), "l");
  leg->AddEntry(band, "#pm1#sigma beam width per strip", "f");
  leg->Draw();
  PlottingUtils::AddText(Form("reaction strip %d", reac), 0.875, 0.68);
  PlottingUtils::SaveFigure(c, Form("mean_trace_reac%d", reac), subdir,
                            PlotSaveOptions::kLINEAR);
  delete c;
  delete frame;
}

void TagEfficiency::DrawPlane(const TString &subdir, Int_t reac, TH2F *plane,
                              TCutG *cut) {
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  c->SetLogz(kTRUE);
  PlottingUtils::ConfigureAndDraw2DHistogram(
      plane, c, Form("reac %d: bootstrap on the drawn cut", reac));
  cut->SetLineColor(kBlack);
  cut->SetLineWidth(2);
  cut->SetFillStyle(0);
  cut->Draw("L SAME");
  PlottingUtils::SaveFigure(c, Form("bootstrap_plane_reac%d", reac), subdir,
                            PlotSaveOptions::kLINEAR);
  delete c;
}

Bool_t TagEfficiency::RunChannel(const CrossSectionChannel &ch,
                                 TTree *reservoir) {
  const StripSumScatterConfig &C = Constants::cfg.STRIP_SUM_SCATTER_CONFIG;
  const CrossSectionConfig &X = Constants::cfg.CROSS_SECTION_CONFIG;
  const TString region = "region_" + ch.name;
  const TString subdir = "tag_efficiency/" + ch.name;
  LoadDrawnCuts(region);
  if (strips_.empty()) {
    std::cerr << "tag-efficiency: channel " << ch.name
              << ": no hand-drawn cuts (" << region
              << "); draw them in strip-sum-scatter first" << std::endl;
    return kFALSE;
  }
  std::cout << "tag-efficiency: channel " << ch.name << ": mean traces inside "
            << region << " from " << reservoir->GetEntries()
            << " reservoir events, jump gate " << C.REAC_JUMP_NSIGMA
            << " sigma; " << kBootstrapTrials << " bootstrap trials per strip"
            << std::endl;
  MeanTraces(reservoir);

  std::vector<TagEfficiencyRecord> records;
  std::cout << Form("%5s %10s %9s %9s %9s", "strip", "in cut", "tag eff", "eff",
                    "migrate")
            << std::endl;
  for (Int_t reac = X.XS_STRIP_MIN; reac <= X.XS_STRIP_MAX; reac++) {
    if (strips_.find(reac) == strips_.end() ||
        strips_[reac].count < kMinEventsInCut) {
      std::cout << Form("%5d %10.0f   no drawn cut or too few events", reac,
                        strips_.count(reac) ? strips_[reac].count : 0.0)
                << std::endl;
      continue;
    }
    const Strip &st = strips_[reac];
    TH2F *plane = new TH2F(Form("h_boot_plane_%s_%d", ch.name.Data(), reac),
                           ";norm. #DeltaE strips 1#rightarrow16 [a.u.];norm. "
                           "#DeltaE post window [a.u.]",
                           300, C.X_DISPLAY_MIN, C.X_DISPLAY_MAX, 300,
                           C.Y_DISPLAY_MIN, C.Y_DISPLAY_MAX);
    const Outcome o = Bootstrap(reac, plane);

    TagEfficiencyRecord r;
    r.reac = reac;
    r.n_counted = st.count;
    r.eff = Double_t(o.inside) / Double_t(o.n);
    // Binomial on the trials; the sample's own error rides on n_counted.
    r.eff_err = std::sqrt(r.eff * (1.0 - r.eff) / Double_t(o.n));
    r.migrate = Double_t(o.next_inside) / Double_t(o.n);
    r.tag_eff = Double_t(o.tagged) / Double_t(o.n);
    records.push_back(r);
    std::cout << Form("%5d %10.0f %9.3f %9.3f %9.3f", reac, r.n_counted,
                      r.tag_eff, r.eff, r.migrate)
              << std::endl;

    DrawMeanTrace(subdir, reac, st);
    DrawPlane(subdir, reac, plane, st.cut);
    delete plane;
  }
  TagEfficiencyStore::Write(ch.name, records, MethodLabel());
  return !records.empty();
}

Bool_t TagEfficiency::Run() {
  const CrossSectionConfig &X = Constants::cfg.CROSS_SECTION_CONFIG;
  if (X.CHANNELS.empty()) {
    std::cerr << "tag-efficiency: this dataset declares no CROSS_SECTION_CONFIG"
                 ".CHANNELS"
              << std::endl;
    return kFALSE;
  }
  TString cache_path =
      IO::GetRootFilesBaseDir() + "/" + StripSumScatter::CacheName();
  TFile cache(cache_path, "READ");
  if (cache.IsZombie()) {
    std::cerr << "tag-efficiency: no scatter cache at " << cache_path
              << "; run strip-sum-scatter first" << std::endl;
    return kFALSE;
  }
  if (!LoadSigmas(cache)) {
    std::cerr << "tag-efficiency: cache carries no per-strip noise sigmas; "
                 "rebuild it"
              << std::endl;
    return kFALSE;
  }
  TTree *reservoir = static_cast<TTree *>(cache.Get("traces"));
  if (!reservoir) {
    std::cerr << "tag-efficiency: no reservoir in the cache" << std::endl;
    return kFALSE;
  }
  Int_t done = 0;
  for (Int_t c = 0; c < Int_t(X.CHANNELS.size()); c++)
    if (RunChannel(X.CHANNELS[c], reservoir))
      done++;
  std::cout << "tag-efficiency: " << done << " of " << X.CHANNELS.size()
            << " channels written to " << TagEfficiencyStore::Path()
            << std::endl;
  return done > 0;
}
