#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Normalization.hpp"
#include "Pipeline.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TColor.h>
#include <TCutG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TKey.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <vector>

namespace RegionTracesNS {

void RegionTracesOneFile(const TString &cut_name, const FileSpec &spec) {
  TString events_path = EventsName(spec) + ".root";
  TString file_label = FileLabel(spec);

  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");

  TFile *file = IO::OpenForWriting(events_path, "UPDATE");
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open " << events_path << std::endl;
    return;
  }

  TDirectory *cuts_dir = file->GetDirectory("cuts");
  if (!cuts_dir) {
    std::cerr << "[" << file_label << "] No cuts/ directory in " << events_path
              << std::endl;
    file->Close();
    return;
  }

  struct CutState {
    TCutG *cut;
    TString name;
    TString hist_name;
    TMultiGraph *mg;
    Double_t sum_per_strip[18];
    Long64_t n_passed;
    Long64_t n_traced;
    Int_t strip_index;
    TH1F *h1_strip_e;
  };

  std::vector<CutState> states;

  if (!cut_name.IsNull()) {
    TCutG *c = static_cast<TCutG *>(cuts_dir->Get(cut_name));
    if (!c) {
      std::cerr << "[" << file_label << "] No cut '" << cut_name << "' in "
                << events_path << ":/cuts/" << std::endl;
      file->Close();
      return;
    }
    CutState s;
    s.cut = c;
    s.name = cut_name;
    s.hist_name = c->GetVarX();
    s.mg = new TMultiGraph();
    for (Int_t k = 0; k < 18; k++)
      s.sum_per_strip[k] = 0;
    s.n_passed = 0;
    s.n_traced = 0;
    s.strip_index = -1;
    s.h1_strip_e = nullptr;
    states.push_back(s);
  } else {
    TIter next(cuts_dir->GetListOfKeys());
    while (TKey *key = static_cast<TKey *>(next())) {
      TObject *obj = key->ReadObj();
      if (!obj || !obj->InheritsFrom(TCutG::Class()))
        continue;
      CutState s;
      s.cut = static_cast<TCutG *>(obj);
      s.name = key->GetName();
      s.hist_name = s.cut->GetVarX();
      s.mg = new TMultiGraph();
      for (Int_t k = 0; k < 18; k++)
        s.sum_per_strip[k] = 0;
      s.n_passed = 0;
      s.strip_index = -1;
      s.h1_strip_e = nullptr;
      states.push_back(s);
    }
    if (states.empty()) {
      std::cerr << "[" << file_label << "] No TCutG objects in " << events_path
                << ":/cuts/" << std::endl;
      file->Close();
      return;
    }
  }

  for (Int_t i = 0; i < Int_t(states.size()); i++) {
    if (states[i].hist_name.IsNull()) {
      std::cerr << "[" << file_label << "] Cut '" << states[i].name
                << "' has no VarX tag; can't infer 2D space" << std::endl;
      for (Int_t k = 0; k < Int_t(states.size()); k++)
        delete states[k].mg;
      file->Close();
      return;
    }
  }

  TTree *tree = static_cast<TTree *>(file->Get("events"));
  if (!tree) {
    std::cerr << "[" << file_label << "] No events tree in " << events_path
              << std::endl;
    for (Int_t k = 0; k < Int_t(states.size()); k++)
      delete states[k].mg;
    file->Close();
    return;
  }

  EnergyView ev;
  ev.Attach(tree);
  if (!ev.is_mev)
    std::cerr << "[" << file_label
              << "] WARNING: no MeV branches; using uncalibrated ADC values"
              << std::endl;
  const char *unit = ev.Unit();
  Double_t strip_e_min =
      ev.is_mev ? Constants::STRIP_E_MIN_MEV : Constants::STRIP_E_MIN_ADC;
  Double_t strip_e_max =
      ev.is_mev ? Constants::STRIP_E_MAX_MEV : Constants::STRIP_E_MAX_ADC;

  TString prefix = "h2_totalE_vs_stripE_s";
  for (Int_t i = 0; i < Int_t(states.size()); i++) {
    CutState &st = states[i];
    if (!st.hist_name.BeginsWith(prefix))
      continue;
    TString num = st.hist_name(prefix.Length(), st.hist_name.Length());
    Int_t s = num.Atoi();
    if (s < 0 || s >= 18)
      continue;
    st.strip_index = s;
    TString hname =
        Form("h1_strip%d_%s_%s", s, st.name.Data(), file_label.Data());
    TString title = Form(";Strip %d #DeltaE [%s];Counts", s, unit);
    st.h1_strip_e = new TH1F(hname, title, 200, strip_e_min, strip_e_max);
    st.h1_strip_e->SetDirectory(nullptr);
  }

  Long64_t n_entries = tree->GetEntries();
  std::cout << "[" << file_label << "] Region traces for " << states.size()
            << " cut(s) over " << n_entries << " events..." << std::endl;

  Double_t strip_x[18];
  for (Int_t s = 0; s < 18; s++)
    strip_x[s] = Double_t(s);

  Int_t indiv_color = TColor::GetColorTransparent(kGray + 2, 0.05);

  for (Long64_t j = 0; j < n_entries; j++) {
    tree->GetEntry(j);
    ev.Decode();

    Double_t event_total = 0;
    for (Int_t s = 0; s < 18; s++)
      event_total += ev.total[s];

    for (Int_t i = 0; i < Int_t(states.size()); i++) {
      CutState &st = states[i];
      CutXY xy = ComputeCutXY(st.hist_name, ev.total, event_total);
      if (!st.cut->IsInside(xy.x, xy.y))
        continue;
      st.n_passed++;

      if (st.h1_strip_e && st.strip_index >= 0)
        st.h1_strip_e->Fill(ev.total[st.strip_index]);

      if (st.n_traced >= Constants::REGION_TRACES_MAX_INDIV)
        continue;

      TGraph *g = new TGraph(18, strip_x, ev.total);
      g->SetLineColor(indiv_color);
      g->SetLineWidth(1);
      g->SetMarkerStyle(0);
      st.mg->Add(g, "L");

      for (Int_t s = 0; s < 18; s++)
        st.sum_per_strip[s] += ev.total[s];
      st.n_traced++;
    }

    if ((j + 1) % 100000 == 0)
      std::cout << "  " << (j + 1) << "/" << n_entries << std::endl;
  }

  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");
  TString out_subdir = "region_traces/" + file_label;

  for (Int_t i = 0; i < Int_t(states.size()); i++) {
    CutState &st = states[i];
    if (st.n_passed == 0) {
      std::cerr << "[" << file_label << "] Cut '" << st.name
                << "' selected zero events; skipping" << std::endl;
      delete st.mg;
      continue;
    }

    Long64_t denom = (st.n_traced > 0) ? st.n_traced : 1;
    Double_t avg_per_strip[18];
    for (Int_t s = 0; s < 18; s++)
      avg_per_strip[s] = st.sum_per_strip[s] / Double_t(denom);

    std::cout << "  '" << st.name << "': " << st.n_passed << " / " << n_entries
              << " (" << (100.0 * st.n_passed / n_entries) << "%)" << std::endl;

    TGraph *avg = new TGraph(18, strip_x, avg_per_strip);
    PlottingUtils::ConfigureGraph(
        avg, kRed + 1,
        Form("Region '%s' (%lld evts);Strip Index;#DeltaE [%s]", st.name.Data(),
             st.n_passed, unit));
    avg->SetMarkerColor(kRed + 1);
    avg->SetLineWidth(3);

    st.mg->SetTitle(Form("'%s' (%lld evts);Strip Index;#DeltaE [%s]",
                         st.name.Data(), st.n_passed, unit));

    Double_t y_max =
        ev.is_mev ? Constants::STRIP_E_MAX_MEV : Constants::STRIP_E_MAX_ADC;
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    st.mg->GetYaxis()->SetRangeUser(0, y_max);
    st.mg->Draw("A");
    avg->Draw("L SAME");
    PlottingUtils::SaveFigure(c, "region_" + st.name, out_subdir,
                              PlotSaveOptions::kLINEAR);
    delete c;
    delete avg;
    delete st.mg;
  }

  Int_t color_palette[8] = {kRed + 1,    kBlue + 1, kGreen + 2,  kMagenta + 1,
                            kOrange + 7, kCyan + 2, kViolet + 1, kBlack};
  for (Int_t s_target = 0; s_target < 18; s_target++) {
    std::vector<Int_t> idx_for_strip;
    for (Int_t i = 0; i < Int_t(states.size()); i++) {
      if (states[i].strip_index == s_target && states[i].h1_strip_e &&
          states[i].n_passed > 0)
        idx_for_strip.push_back(i);
    }
    if (idx_for_strip.empty())
      continue;

    Double_t y_max = 0;
    for (Int_t k = 0; k < Int_t(idx_for_strip.size()); k++) {
      Double_t m = states[idx_for_strip[k]].h1_strip_e->GetMaximum();
      if (m > y_max)
        y_max = m;
    }

    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    TLegend *leg = PlottingUtils::AddLegend(0.65, 0.99, 0.70, 0.92);
    for (Int_t k = 0; k < Int_t(idx_for_strip.size()); k++) {
      CutState &st = states[idx_for_strip[k]];
      Int_t color = color_palette[k % 8];
      st.h1_strip_e->SetLineColor(color);
      st.h1_strip_e->SetLineWidth(2);
      st.h1_strip_e->SetMinimum(10);
      st.h1_strip_e->SetMaximum(1.15 * y_max);
      st.h1_strip_e->SetTitle(
          Form(";Strip %d #DeltaE [%s];Counts", s_target, unit));
      st.h1_strip_e->Draw(k == 0 ? "HIST" : "HIST SAME");
      leg->AddEntry(st.h1_strip_e,
                    Form("%s (%lld)", st.name.Data(), st.n_passed), "l");
    }
    leg->Draw();
    PlottingUtils::SaveFigure(c, Form("strip%d_cut_overlay", s_target),
                              out_subdir, PlotSaveOptions::kLOG);

    file->cd();
    c->Write(Form("strip%d_cut_overlay", s_target), TObject::kOverwrite);
    delete leg;
    delete c;
  }

  for (Int_t i = 0; i < Int_t(states.size()); i++)
    delete states[i].h1_strip_e;

  file->Close();
}

} // namespace RegionTracesNS

void RegionTraces(TString cut_name = "", TString file_label = "") {
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

  for (Int_t k = 0; k < Int_t(specs.size()); k++)
    RegionTracesNS::RegionTracesOneFile(cut_name, specs[k]);
}
