#include "DiagnoseSubfileDrift.hpp"

const Int_t HIST_BINS = 400;
const Double_t HIST_MIN = 0.0;
const Double_t HIST_MAX = 10000.0;

// Per-channel peak windows that bracket the lowest-peak feature.
const Double_t FIT_LO_S0 = 1200.0;
const Double_t FIT_HI_S0 = 6000.0;
const Double_t FIT_LO_LR = 500.0;
const Double_t FIT_HI_LR = 1200.0;
const Double_t FIT_LO_S17 = 700.0;
const Double_t FIT_HI_S17 = 1800.0;

const Long64_t MIN_ENTRIES = 1000;
const Double_t MIN_PEAK_COUNTS = 20.0;

struct Channel {
  TString name;
  Int_t strip;
  Char_t side;
};

inline std::vector<Channel> LongChannels() {
  std::vector<Channel> chans;
  chans.push_back({"Strip0", 0, 'S'});
  for (Int_t s = 1; s <= 15; s += 2)
    chans.push_back({Form("L%d", s), s, 'L'});
  for (Int_t s = 2; s <= 16; s += 2)
    chans.push_back({Form("R%d", s), s, 'R'});
  chans.push_back({"Strip17", 17, 'S'});
  return chans;
}

inline void GetFitWindow(const Channel &c, Double_t &lo, Double_t &hi) {
  if (c.side == 'S' && c.strip == 0) {
    lo = FIT_LO_S0;
    hi = FIT_HI_S0;
  } else if (c.side == 'S' && c.strip == 17) {
    lo = FIT_LO_S17;
    hi = FIT_HI_S17;
  } else {
    lo = FIT_LO_LR;
    hi = FIT_HI_LR;
  }
}

struct PeakLoc {
  Double_t mu;
  Bool_t ok;
};

inline PeakLoc WindowedPeak(TH1D *h, Double_t lo, Double_t hi) {
  PeakLoc p{0, kFALSE};
  if (h->GetEntries() < MIN_ENTRIES)
    return p;
  Int_t lo_bin = h->FindBin(lo);
  Int_t hi_bin = h->FindBin(hi);
  if (lo_bin < 2)
    lo_bin = 2;
  if (hi_bin >= h->GetNbinsX())
    hi_bin = h->GetNbinsX() - 1;

  Int_t max_bin = lo_bin;
  Double_t max_val = h->GetBinContent(lo_bin);
  for (Int_t b = lo_bin + 1; b <= hi_bin; b++) {
    if (h->GetBinContent(b) > max_val) {
      max_val = h->GetBinContent(b);
      max_bin = b;
    }
  }
  if (max_val < MIN_PEAK_COUNTS)
    return p;

  // Parabolic interpolation across 3 bins around max for sub-bin precision.
  Double_t y0 = h->GetBinContent(max_bin - 1);
  Double_t y1 = h->GetBinContent(max_bin);
  Double_t y2 = h->GetBinContent(max_bin + 1);
  Double_t denom = (y0 - 2.0 * y1 + y2);
  Double_t delta = (denom != 0.0) ? 0.5 * (y0 - y2) / denom : 0.0;
  if (delta < -1.0)
    delta = -1.0;
  if (delta > 1.0)
    delta = 1.0;
  p.mu = h->GetBinCenter(max_bin) + delta * h->GetBinWidth(max_bin);
  p.ok = kTRUE;
  return p;
}

inline void FillSubfilePeaks(const TString &events_path,
                             const std::vector<Channel> &chans,
                             std::vector<PeakLoc> &out_peaks,
                             std::vector<TH1D *> &out_chan_hists,
                             Long64_t &out_entries) {
  out_peaks.assign(chans.size(), {0, kFALSE});
  out_chan_hists.assign(chans.size(), nullptr);
  out_entries = 0;
  TFile *f = IO::OpenForReading(events_path);
  if (!f || f->IsZombie()) {
    std::cerr << "  cannot open " << events_path << std::endl;
    if (f)
      f->Close();
    return;
  }
  TTree *t = static_cast<TTree *>(f->Get("events"));
  if (!t) {
    std::cerr << "  no events tree in " << events_path << std::endl;
    f->Close();
    return;
  }

  Int_t leftdE[18], rightdE[18], totaldE[18];
  t->SetBranchAddress("LeftdE", leftdE);
  t->SetBranchAddress("RightdE", rightdE);
  t->SetBranchAddress("TotaldE", totaldE);

  std::vector<TH1D *> hists(chans.size(), nullptr);
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    hists[i] = new TH1D(PlottingUtils::GetRandomName(), "", HIST_BINS, HIST_MIN,
                        HIST_MAX);
    hists[i]->SetDirectory(nullptr);
  }

  Long64_t n = t->GetEntries();
  out_entries = n;
  for (Long64_t j = 0; j < n; j++) {
    t->GetEntry(j);
    for (Int_t i = 0; i < Int_t(chans.size()); i++) {
      const Channel &c = chans[i];
      Int_t v = (c.side == 'S')
                    ? totaldE[c.strip]
                    : (c.side == 'L' ? leftdE[c.strip] : rightdE[c.strip]);
      if (v > 0)
        hists[i]->Fill(Double_t(v));
    }
  }

  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    Double_t lo = 0, hi = 0;
    GetFitWindow(chans[i], lo, hi);
    out_peaks[i] = WindowedPeak(hists[i], lo, hi);
    out_chan_hists[i] = hists[i];
  }
  f->Close();
}

inline Double_t DisplayMax(const Channel &c) {
  if (c.side == 'S' && c.strip == 17)
    return 4000.0;
  Double_t lo = 0, hi = 0;
  GetFitWindow(c, lo, hi);
  return TMath::Min(2.0 * hi, HIST_MAX);
}

void DriftOneRun(Int_t run, const std::vector<FileSpec> &specs) {
  std::vector<Channel> chans = LongChannels();
  Int_t n_chans = Int_t(chans.size());
  Int_t n_sub = Int_t(specs.size());

  std::vector<std::vector<PeakLoc>> peaks(
      n_chans, std::vector<PeakLoc>(n_sub, {0, kFALSE}));
  // chan_hists[chan][sub] holds the per-subfile ADC histogram for every
  // channel.
  std::vector<std::vector<TH1D *>> chan_hists(
      n_chans, std::vector<TH1D *>(n_sub, nullptr));
  std::vector<TString> labels(n_sub);
  for (Int_t s = 0; s < n_sub; s++)
    labels[s] = FileSet::FileLabel(specs[s]);

  std::queue<Int_t> work;
  for (Int_t s = 0; s < n_sub; s++)
    work.push(s);
  std::mutex work_mutex;
  std::mutex result_mutex;

  Int_t n_workers =
      TMath::Min(Int_t(std::thread::hardware_concurrency()), n_sub);
  n_workers = TMath::Min(n_workers, Constants::MAX_DIAGNOSE_WORKERS);
  if (n_workers < 1)
    n_workers = 1;
  std::cout << "  reading " << n_sub << " subfiles on " << n_workers
            << " worker(s)..." << std::endl;

  std::vector<std::thread> workers;
  for (Int_t w = 0; w < n_workers; w++) {
    workers.emplace_back([&]() {
      while (true) {
        Int_t s;
        {
          std::lock_guard<std::mutex> lk(work_mutex);
          if (work.empty())
            return;
          s = work.front();
          work.pop();
        }
        TString events_path = FileSet::EventsName(specs[s]) + ".root";
        std::vector<PeakLoc> p;
        std::vector<TH1D *> sh;
        Long64_t n_entries = 0;
        FillSubfilePeaks(events_path, chans, p, sh, n_entries);
        {
          std::lock_guard<std::mutex> lk(result_mutex);
          for (Int_t i = 0; i < n_chans; i++) {
            peaks[i][s] = p[i];
            chan_hists[i][s] = sh[i];
          }
          std::cout << "    done " << events_path << " (" << n_entries
                    << " events)" << std::endl;
        }
      }
    });
  }
  for (Int_t w = 0; w < Int_t(workers.size()); w++)
    workers[w].join();

  TString plot_subdir = Form("subfile_drift/run%d", run);
  std::vector<Int_t> colors = PlottingUtils::GetDefaultColors();

  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle(";Subfile index;Peak / Peak(first valid)");

  Double_t y_min = 1.0, y_max = 1.0;
  std::vector<Int_t> legend_channel_idx;

  for (Int_t i = 0; i < n_chans; i++) {
    if (chans[i].side == 'S')
      continue; // guard strips get their own per-subfile overlay below
    Double_t mu_ref = 0;
    for (Int_t s = 0; s < n_sub; s++) {
      if (peaks[i][s].ok) {
        mu_ref = peaks[i][s].mu;
        break;
      }
    }
    if (mu_ref <= 0) {
      std::cerr << "  channel " << chans[i].name
                << " has no valid peak in any subfile" << std::endl;
      continue;
    }

    TGraph *g = new TGraph();
    Int_t pt = 0;
    for (Int_t s = 0; s < n_sub; s++) {
      if (!peaks[i][s].ok)
        continue;
      Double_t y = peaks[i][s].mu / mu_ref;
      g->SetPoint(pt, Double_t(s), y);
      if (y < y_min)
        y_min = y;
      if (y > y_max)
        y_max = y;
      pt++;
    }
    if (pt == 0) {
      delete g;
      continue;
    }
    Int_t color = colors[i % Int_t(colors.size())];
    g->SetLineColor(color);
    g->SetMarkerColor(color);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.8);
    g->SetLineWidth(2);
    g->SetName(chans[i].name);
    g->SetTitle(chans[i].name);
    mg->Add(g, "LP");
    legend_channel_idx.push_back(i);
  }

  if (legend_channel_idx.empty()) {
    std::cerr << "Run " << run << ": no valid peak data; skipping plot"
              << std::endl;
    delete mg;
    return;
  }

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  canvas->SetRightMargin(0.20);
  mg->Draw("A");
  Double_t pad = TMath::Max(0.02, 0.1 * (y_max - y_min));
  mg->GetYaxis()->SetRangeUser(y_min - pad, y_max + pad);

  TLegend *leg = PlottingUtils::AddLegend(0.81, 0.99, 0.10, 0.95);
  TList *graph_list = mg->GetListOfGraphs();
  for (Int_t k = 0; k < Int_t(legend_channel_idx.size()); k++) {
    TGraph *g = static_cast<TGraph *>(graph_list->At(k));
    leg->AddEntry(g, chans[legend_channel_idx[k]].name.Data(), "lp");
  }
  leg->Draw();
  canvas->Modified();
  canvas->Update();

  if (Constants::SAVE_PLOTS)
    PlottingUtils::SaveFigure(canvas, "normalized_peak_vs_subfile", plot_subdir,
                              PlotSaveOptions::kLINEAR);

  // Summary: peak-to-peak spread per channel.
  std::cout << "Run " << run
            << " summary (peak-to-peak normalized spread):" << std::endl;
  for (Int_t i = 0; i < n_chans; i++) {
    Double_t mu_ref = 0;
    for (Int_t s = 0; s < n_sub; s++)
      if (peaks[i][s].ok) {
        mu_ref = peaks[i][s].mu;
        break;
      }
    if (mu_ref <= 0)
      continue;
    Double_t lo = 1e9, hi = -1e9;
    Int_t pts = 0;
    for (Int_t s = 0; s < n_sub; s++) {
      if (!peaks[i][s].ok)
        continue;
      Double_t y = peaks[i][s].mu / mu_ref;
      if (y < lo)
        lo = y;
      if (y > hi)
        hi = y;
      pts++;
    }
    if (pts < 2)
      continue;
    std::cout << "  " << chans[i].name << ":  spread = " << 100.0 * (hi - lo)
              << "%  (" << pts << " pts)" << std::endl;
  }

  delete canvas;
  delete leg;
  delete mg;

  gStyle->SetPalette(kViridis);
  Int_t n_palette = TColor::GetNumberOfColors();
  for (Int_t i = 0; i < n_chans; i++) {
    Double_t x_max = DisplayMax(chans[i]);
    Double_t y_top = 0;
    Int_t n_valid = 0;
    for (Int_t s = 0; s < n_sub; s++) {
      if (!chan_hists[i][s])
        continue;
      n_valid++;
      Int_t hi_bin = chan_hists[i][s]->FindBin(x_max);
      for (Int_t b = 1; b <= hi_bin; b++) {
        Double_t v = chan_hists[i][s]->GetBinContent(b);
        if (v > y_top)
          y_top = v;
      }
    }
    if (n_valid == 0)
      continue;

    TCanvas *c_chan = PlottingUtils::GetConfiguredCanvas(kFALSE);
    for (Int_t s = 0; s < n_sub; s++) {
      if (!chan_hists[i][s])
        continue;
      Int_t idx = (n_sub > 1) ? (s * (n_palette - 1)) / (n_sub - 1) : 0;
      Int_t color = TColor::GetColorPalette(idx);
      chan_hists[i][s]->SetLineColor(color);
      chan_hists[i][s]->SetLineWidth(2);
      chan_hists[i][s]->SetTitle(";#DeltaE [ADC];Counts");
      chan_hists[i][s]->SetMaximum(1.15 * y_top);
      chan_hists[i][s]->GetXaxis()->SetRangeUser(0.0, x_max);
      chan_hists[i][s]->GetXaxis()->SetNdivisions(505);
      chan_hists[i][s]->Draw(s == 0 ? "HIST" : "HIST SAME");
    }
    PlottingUtils::AddText("purple = beginning, yellow = end", 0.88, 0.82);
    if (Constants::SAVE_PLOTS)
      PlottingUtils::SaveFigure(c_chan, "overlay_" + chans[i].name, plot_subdir,
                                PlotSaveOptions::kLINEAR);
    delete c_chan;
  }

  for (Int_t i = 0; i < n_chans; i++)
    for (Int_t s = 0; s < n_sub; s++)
      delete chan_hists[i][s];
}

void DiagnoseSubfileDrift::Run() {
  ROOT::EnableThreadSafety();
  const TString project_root = Paths::DatasetDir();
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");

  std::vector<FileSpec> all_specs = FileSet::BuildFileSpecs();
  for (Int_t r = 0; r < Int_t(Constants::RUN_NUMBERS.size()); r++) {
    Int_t run = Constants::RUN_NUMBERS[r];
    std::vector<FileSpec> run_specs;
    for (Int_t k = 0; k < Int_t(all_specs.size()); k++) {
      if (all_specs[k].run == run)
        run_specs.push_back(all_specs[k]);
    }
    if (run_specs.size() < 2) {
      std::cout << "Skipping run " << run << " (" << run_specs.size()
                << " subfile(s); need >= 2)" << std::endl;
      continue;
    }
    std::cout << "Subfile drift for run " << run << " (" << run_specs.size()
              << " subfiles)" << std::endl;
    DriftOneRun(run, run_specs);
  }
}
