#include "Constants.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Pipeline.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLine.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <mutex>
#include <ostream>
#include <vector>

namespace BaselineNorm {

// Strip0 / Strip17 are single-readout (in totaldE), tagged 'S'; long
// strips are split L/R.
struct Channel {
  TString name;
  Int_t strip;
  Char_t side;
  TH1D *hist;
};

std::vector<Channel> BuildChannels(const TString &file_label) {
  std::vector<Channel> chans;
  Channel c;

  c.name = "Strip0";
  c.strip = 0;
  c.side = 'S';
  c.hist = new TH1D(Form("h_baseline_%s_Strip0", file_label.Data()),
                    Form("Strip0 (%s);#DeltaE [ADC];Counts", file_label.Data()),
                    Constants::BASELINE_HIST_BINS, 0,
                    Constants::BASELINE_HIST_MAX_ADC);
  chans.push_back(c);

  c.name = "Strip17";
  c.strip = 17;
  c.side = 'S';
  c.hist = new TH1D(Form("h_baseline_%s_Strip17", file_label.Data()),
                    Form("Strip17 (%s);#DeltaE [ADC];Counts", file_label.Data()),
                    Constants::BASELINE_HIST_BINS, 0,
                    Constants::BASELINE_HIST_MAX_ADC);
  chans.push_back(c);

  for (Int_t s = 1; s <= 16; s++) {
    c.name = Form("L%d", s);
    c.strip = s;
    c.side = 'L';
    c.hist = new TH1D(Form("h_baseline_%s_L%d", file_label.Data(), s),
                      Form("L%d (%s);#DeltaE [ADC];Counts", s, file_label.Data()),
                      Constants::BASELINE_HIST_BINS, 0,
                      Constants::BASELINE_HIST_MAX_ADC);
    chans.push_back(c);

    c.name = Form("R%d", s);
    c.strip = s;
    c.side = 'R';
    c.hist = new TH1D(Form("h_baseline_%s_R%d", file_label.Data(), s),
                      Form("R%d (%s);#DeltaE [ADC];Counts", s, file_label.Data()),
                      Constants::BASELINE_HIST_BINS, 0,
                      Constants::BASELINE_HIST_MAX_ADC);
    chans.push_back(c);
  }
  return chans;
}

void FillChannels(TTree *tree, std::vector<Channel> &chans) {
  Int_t leftdE[18], rightdE[18], totaldE[18];
  tree->SetBranchAddress("LeftdE", leftdE);
  tree->SetBranchAddress("RightdE", rightdE);
  tree->SetBranchAddress("TotaldE", totaldE);

  Long64_t n = tree->GetEntries();
  std::cout << "Filling baseline histograms from " << n << " entries..."
            << std::endl;
  for (Long64_t j = 0; j < n; j++) {
    tree->GetEntry(j);
    for (Int_t k = 0; k < Int_t(chans.size()); k++) {
      const Channel &c = chans[k];
      Int_t v = 0;
      if (c.side == 'S')
        v = totaldE[c.strip];
      else if (c.side == 'L')
        v = leftdE[c.strip];
      else
        v = rightdE[c.strip];
      if (v > 0)
        c.hist->Fill(Double_t(v));
    }
    if ((j + 1) % 100000 == 0) {
      std::cout << "  " << (j + 1) << "/" << n << std::endl;
    }
  }
}

Float_t MaxBinCenter(TH1D *hist) {
  if (hist->GetEntries() < Constants::BASELINE_MIN_ENTRIES)
    return 0;
  return Float_t(hist->GetBinCenter(hist->GetMaximumBin()));
}

void SavePeakPlot(TH1D *hist, Double_t mu, const TString &plot_subdir,
                  const TString &filename, Double_t xmax = -1) {
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureAndDrawHistogram(hist, kBlack);
  if (xmax > 0)
    hist->GetXaxis()->SetRangeUser(0, xmax);
  TLine *line =
      new TLine(mu, 0, mu, hist->GetBinContent(hist->GetMaximumBin()));
  line->SetLineColor(kRed + 1);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw("same");
  PlottingUtils::SaveFigure(c, filename, plot_subdir, PlotSaveOptions::kLOG);
  delete line;
  delete c;
}

} // namespace BaselineNorm

void BuildBaselineNormalization(std::vector<TString> input_filenames,
                                std::vector<TString> file_labels,
                                Bool_t reprocess = kFALSE) {
  if (!reprocess)
    return;

  Int_t n_files = Int_t(input_filenames.size());
  for (Int_t i = 0; i < n_files; i++) {
    TString input_filepath = input_filenames[i] + ".root";
    TString file_label = file_labels[i];

    TFile *file = IO::OpenForWriting(input_filepath, "UPDATE");
    if (!file || file->IsZombie()) {
      std::cerr << "Cannot open " << input_filepath << std::endl;
      continue;
    }
    TTree *tree = static_cast<TTree *>(file->Get("event"));
    if (!tree) {
      std::cerr << "No event tree in " << input_filepath << std::endl;
      file->Close();
      continue;
    }

    std::vector<BaselineNorm::Channel> chans = BaselineNorm::BuildChannels(file_label);
    BaselineNorm::FillChannels(tree, chans);

    TDirectory *hist_dir =
        file->mkdir("strip_histograms", "strip_histograms", kTRUE);
    TString plot_subdir = "baseline/" + file_label;

    Float_t peak[18] = {0};
    Float_t scale[18];
    for (Int_t s = 0; s < 18; s++)
      scale[s] = 1.0f;

    {
      std::lock_guard<std::mutex> lock(g_plot_mutex);
      for (Int_t k = 0; k < Int_t(chans.size()); k++) {
        const BaselineNorm::Channel &c = chans[k];
        Float_t p = BaselineNorm::MaxBinCenter(c.hist);
        std::cout << "  " << c.name << ": peak=" << p << std::endl;

        // Per-strip "peak" used downstream: Strip0/17 from totaldE, L_odd
        // from L_s, R_even from R_s. The unused side stays at 0 and is
        // skipped in the long-strip averaging below.
        Bool_t is_long = (c.side == 'L' && (c.strip % 2 == 1)) ||
                         (c.side == 'R' && (c.strip % 2 == 0));
        if (c.side == 'S' || is_long)
          peak[c.strip] = p;

        BaselineNorm::SavePeakPlot(c.hist, p, plot_subdir,
                               "baseline_" + c.name);
        if (c.name == "L1" || c.name == "R2")
          BaselineNorm::SavePeakPlot(c.hist, p, plot_subdir,
                                 "baseline_" + c.name + "_zoom", 1500.0);

        hist_dir->cd();
        c.hist->Write(Form("h_%s", c.name.Data()), TObject::kOverwrite);
      }
    }

    Double_t long_sum = 0;
    Int_t long_count = 0;
    for (Int_t s = 1; s <= 16; s++) {
      if (peak[s] > 0) {
        long_sum += peak[s];
        long_count++;
      }
    }
    if (long_count > 0 && peak[0] > 0) {
      Double_t avg_long = long_sum / long_count;
      scale[0] = Float_t(avg_long / peak[0]);
      std::cout << "Avg long-strip peak = " << avg_long << " ADC ("
                << long_count << " strips)" << std::endl;
      std::cout << "  Strip0 peak=" << peak[0] << " => scale=" << scale[0]
                << std::endl;
    } else {
      std::cout << "Skipping Strip0 gain match (missing peaks)" << std::endl;
    }

    file->cd();
    if (TObject *old = file->Get("baseline"))
      old->Delete();
    TTree *cal = new TTree("baseline", "Strip0/17 gain match scales");
    cal->Branch("Scale", scale, "Scale[18]/F");
    cal->Branch("Peak", peak, "Peak[18]/F");
    cal->Fill();
    cal->Write("baseline", TObject::kOverwrite);

    for (Int_t k = 0; k < Int_t(chans.size()); k++)
      delete chans[k].hist;

    file->Close();
    std::cout << "Wrote baseline scales for " << file_label << std::endl;
  }
}

void BaselineNormalization() {
  std::vector<TString> input_filenames, file_labels;
  std::vector<FileSpec> specs = BuildFileSpecs();
  for (Int_t k = 0; k < Int_t(specs.size()); k++) {
    input_filenames.push_back(EventsName(specs[k]));
    file_labels.push_back(FileLabel(specs[k]));
  }
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                TString(gSystem->pwd()) + "/plots",
                                TString(gSystem->pwd()) + "/root_files");
  BuildBaselineNormalization(input_filenames, file_labels, kTRUE);
}
