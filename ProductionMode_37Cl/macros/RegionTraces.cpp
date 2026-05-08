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
#include <TKey.h>
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

  TFile *file = IO::OpenForReading(events_path);
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open " << events_path << std::endl;
    return;
  }

  TDirectory *cuts_dir = file->GetDirectory("cuts");
  if (!cuts_dir) {
    std::cerr << "[" << file_label << "] No cuts/ directory in "
              << events_path << std::endl;
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
      states.push_back(s);
    }
    if (states.empty()) {
      std::cerr << "[" << file_label << "] No TCutG objects in "
                << events_path << ":/cuts/" << std::endl;
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

  Baseline baseline;
  if (!LoadBaseline(file, baseline))
    std::cerr << "[" << file_label
              << "] WARNING: no baseline tree; using raw values" << std::endl;

  TTree *tree = static_cast<TTree *>(file->Get("event"));
  if (!tree) {
    std::cerr << "[" << file_label << "] No event tree in " << events_path
              << std::endl;
    for (Int_t k = 0; k < Int_t(states.size()); k++)
      delete states[k].mg;
    file->Close();
    return;
  }

  Int_t leftdE[18], rightdE[18], totaldE[18];
  tree->SetBranchAddress("LeftdE", leftdE);
  tree->SetBranchAddress("RightdE", rightdE);
  tree->SetBranchAddress("TotaldE", totaldE);

  Long64_t n_entries = tree->GetEntries();
  std::cout << "[" << file_label << "] Region traces for " << states.size()
            << " cut(s) over " << n_entries << " events..." << std::endl;

  Double_t strip_x[18];
  for (Int_t s = 0; s < 18; s++)
    strip_x[s] = Double_t(s);

  Int_t indiv_color = TColor::GetColorTransparent(kGray + 2, 0.05);

  for (Long64_t j = 0; j < n_entries; j++) {
    tree->GetEntry(j);

    Double_t normalized[18];
    ComputeNormalized(baseline, leftdE, rightdE, totaldE, normalized);

    Double_t event_total = 0;
    for (Int_t s = 0; s < 18; s++)
      event_total += normalized[s];

    for (Int_t i = 0; i < Int_t(states.size()); i++) {
      CutState &st = states[i];
      if (st.n_passed >= Constants::REGION_TRACES_MAX_INDIV)
        continue;
      CutXY xy = ComputeCutXY(st.hist_name, normalized, event_total);
      if (!st.cut->IsInside(xy.x, xy.y))
        continue;

      TGraph *g = new TGraph(18, strip_x, normalized);
      g->SetLineColor(indiv_color);
      g->SetLineWidth(1);
      g->SetMarkerStyle(0);
      st.mg->Add(g, "L");

      for (Int_t s = 0; s < 18; s++)
        st.sum_per_strip[s] += normalized[s];
      st.n_passed++;
    }

    if ((j + 1) % 100000 == 0)
      std::cout << "  " << (j + 1) << "/" << n_entries << std::endl;
  }

  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                TString(gSystem->pwd()) + "/plots",
                                TString(gSystem->pwd()) + "/root_files");
  TString out_subdir = "region_traces/" + file_label;

  for (Int_t i = 0; i < Int_t(states.size()); i++) {
    CutState &st = states[i];
    if (st.n_passed == 0) {
      std::cerr << "[" << file_label << "] Cut '" << st.name
                << "' selected zero events; skipping" << std::endl;
      delete st.mg;
      continue;
    }

    Double_t avg_per_strip[18];
    for (Int_t s = 0; s < 18; s++)
      avg_per_strip[s] = st.sum_per_strip[s] / st.n_passed;

    std::cout << "  '" << st.name << "': " << st.n_passed << " / " << n_entries
              << " (" << (100.0 * st.n_passed / n_entries) << "%)" << std::endl;

    TGraph *avg = new TGraph(18, strip_x, avg_per_strip);
    PlottingUtils::ConfigureGraph(
        avg, kRed + 1,
        Form("Region '%s' (%lld evts);Strip Index;Normalized #DeltaE [ADC]",
             st.name.Data(), st.n_passed));
    avg->SetMarkerColor(kRed + 1);
    avg->SetLineWidth(3);

    st.mg->SetTitle(Form("'%s' (%lld evts);Strip Index;#DeltaE [ADC]",
                         st.name.Data(), st.n_passed));

    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    st.mg->GetYaxis()->SetRangeUser(0, 2000);
    st.mg->Draw("A");
    avg->Draw("L SAME");
    PlottingUtils::SaveFigure(c, "region_" + st.name, out_subdir,
                              PlotSaveOptions::kLINEAR);
    delete c;
    delete avg;
    delete st.mg;
  }

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
