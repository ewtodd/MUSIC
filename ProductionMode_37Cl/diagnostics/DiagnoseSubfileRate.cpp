R__ADD_INCLUDE_PATH(../ macros)
#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Pipeline.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TMath.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TString.h>
#include <TTree.h>
#include <algorithm>
#include <iostream>
#include <vector>

namespace SubfileRate {

struct SubfilePoint {
  Int_t index;
  Double_t rate_hz;
  Long64_t n_events;
  TString label;
};

inline Bool_t ReadGridRate(const TString &events_path, Double_t &rate_hz,
                           Long64_t &n_events) {
  rate_hz = 0.0;
  n_events = 0;
  TFile *f = IO::OpenForReading(events_path);
  if (!f || f->IsZombie()) {
    if (f)
      f->Close();
    return kFALSE;
  }
  TParameter<Double_t> *p =
      dynamic_cast<TParameter<Double_t> *>(f->Get("grid_rate_hz"));
  if (!p) {
    f->Close();
    return kFALSE;
  }
  rate_hz = p->GetVal();
  TTree *t = static_cast<TTree *>(f->Get("events"));
  if (t)
    n_events = t->GetEntries();
  f->Close();
  return kTRUE;
}

void RateOneRun(Int_t run, const std::vector<FileSpec> &specs) {
  Int_t n_sub = Int_t(specs.size());
  std::vector<SubfilePoint> pts;
  pts.reserve(n_sub);

  for (Int_t s = 0; s < n_sub; s++) {
    TString events_path = EventsName(specs[s]) + ".root";
    Double_t rate = 0.0;
    Long64_t nev = 0;
    if (!ReadGridRate(events_path, rate, nev)) {
      std::cerr << "  missing grid_rate_hz in " << events_path
                << " (rebuild events to populate)" << std::endl;
      continue;
    }
    SubfilePoint pt;
    pt.index = s;
    pt.rate_hz = rate;
    pt.n_events = nev;
    pt.label = FileLabel(specs[s]);
    pts.push_back(pt);
  }

  if (pts.empty()) {
    std::cerr << "Run " << run << ": no grid_rate_hz values; skipping plot"
              << std::endl;
    return;
  }

  TGraph *g = new TGraph();
  for (Int_t k = 0; k < Int_t(pts.size()); k++)
    g->SetPoint(k, Double_t(pts[k].index), pts[k].rate_hz);
  g->SetTitle(";Subfile index;Grid trigger rate [Hz]");
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.9);
  g->SetMarkerColor(kBlue + 1);
  g->SetLineColor(kBlue + 1);
  g->SetLineWidth(2);

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  g->Draw("ALP");
  canvas->Modified();
  canvas->Update();

  TString plot_subdir = Form("subfile_rate/run%d", run);
  if (Constants::SAVE_PLOTS)
    PlottingUtils::SaveFigure(canvas, "grid_rate_vs_subfile", plot_subdir,
                              PlotSaveOptions::kLINEAR);

  std::vector<Double_t> sorted_rates;
  sorted_rates.reserve(pts.size());
  for (Int_t k = 0; k < Int_t(pts.size()); k++)
    sorted_rates.push_back(pts[k].rate_hz);
  std::sort(sorted_rates.begin(), sorted_rates.end());
  Double_t lo = sorted_rates.front();
  Double_t hi = sorted_rates.back();
  Double_t median = sorted_rates[sorted_rates.size() / 2];

  std::cout << "Run " << run << " grid trigger rate summary (" << pts.size()
            << " subfiles)" << std::endl;
  std::cout << "  min = " << lo << " Hz, median = " << median
            << " Hz, max = " << hi << " Hz, peak-to-peak = " << (hi - lo)
            << " Hz" << std::endl;
  for (Int_t k = 0; k < Int_t(pts.size()); k++) {
    std::cout << "    " << pts[k].label << " (idx " << pts[k].index
              << "):  " << pts[k].rate_hz << " Hz   N=" << pts[k].n_events
              << std::endl;
  }

  delete canvas;
  delete g;
}

} // namespace SubfileRate

void DiagnoseSubfileRate() {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");

  std::vector<FileSpec> all_specs = BuildFileSpecs();
  for (Int_t r = 0; r < Int_t(Constants::RUN_NUMBERS.size()); r++) {
    Int_t run = Constants::RUN_NUMBERS[r];
    std::vector<FileSpec> run_specs;
    for (Int_t k = 0; k < Int_t(all_specs.size()); k++) {
      if (all_specs[k].run == run)
        run_specs.push_back(all_specs[k]);
    }
    if (run_specs.empty()) {
      std::cout << "Skipping run " << run << " (no subfiles)" << std::endl;
      continue;
    }
    std::cout << "Subfile grid rate for run " << run << " (" << run_specs.size()
              << " subfiles)" << std::endl;
    SubfileRate::RateOneRun(run, run_specs);
  }
}
