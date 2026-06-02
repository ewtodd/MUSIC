#include "StripScatterOverlay.hpp"

namespace {

// Axis ranges mirror the experimental ProductionMode StripScatter so the sim
// overlays are directly comparable. The window plot spans three strips, so its
// x-axis is 3x the single-strip range (as in the experimental code).
const Double_t kStripEMin = -0.2;
const Double_t kStripEMax = 10.0;
const Double_t kTotalEMin = 10.0;
const Double_t kTotalEMax = 150.0;

} // namespace

// Common setup for the stride-sampled scatter builders below.
void StripScatterOverlay::PrepareTree(TTree *tree, Float_t *totaldE,
                                      Long64_t max_points, Long64_t &n_total,
                                      Long64_t &n_use, Long64_t &stride) {
  const char *total_name =
      tree->GetBranch("TotaldE") ? "TotaldE" : "TotaldEMeV";
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus(total_name, 1);
  tree->SetBranchAddress(total_name, totaldE);

  n_total = tree->GetEntries();
  n_use = (max_points > 0 && n_total > max_points) ? max_points : n_total;
  // Stride-sample so we cover the full file rather than just its head.
  stride = (n_use > 0) ? (n_total / n_use) : 1;
  if (stride < 1)
    stride = 1;
}

void StripScatterOverlay::StyleScatter(TGraph *g, Int_t color) {
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.3);
  g->SetMarkerColorAlpha(color, 0.35);
  g->SetLineColor(color);
}

// x = strip-s #DeltaE, y = full event total #DeltaE.
TGraph *StripScatterOverlay::BuildStripScatter(TTree *tree, Int_t color,
                                               Long64_t max_points,
                                               Int_t strip) {
  Float_t totaldE[18] = {0};
  Long64_t n_total = 0, n_use = 0, stride = 1;
  PrepareTree(tree, totaldE, max_points, n_total, n_use, stride);

  const Int_t s_lo = Remix::FirstStrip();
  const Int_t s_hi = Remix::LastStrip();

  TGraph *g = new TGraph(n_use);
  Long64_t k = 0;
  for (Long64_t j = 0; j < n_total && k < n_use; j += stride) {
    tree->GetEntry(j);
    Double_t x = Double_t(totaldE[strip]);
    Double_t y = 0.0;
    for (Int_t s = s_lo; s <= s_hi; s++)
      y += Double_t(totaldE[s]);
    g->SetPoint(k++, x, y);
  }
  g->Set(k);
  StyleScatter(g, color);
  return g;
}

// x = #DeltaE summed over {center-1, center, center+1}, y = #DeltaE summed from
// the reaction (center) strip to the downstream end of the detector.
TGraph *StripScatterOverlay::BuildReactionToEndScatter(TTree *tree, Int_t color,
                                                       Long64_t max_points,
                                                       Int_t center) {
  Float_t totaldE[18] = {0};
  Long64_t n_total = 0, n_use = 0, stride = 1;
  PrepareTree(tree, totaldE, max_points, n_total, n_use, stride);

  const Int_t s_lo = Remix::FirstStrip();
  const Int_t s_hi = Remix::LastStrip();

  TGraph *g = new TGraph(n_use);
  Long64_t k = 0;
  for (Long64_t j = 0; j < n_total && k < n_use; j += stride) {
    tree->GetEntry(j);
    Double_t x = 0.0;
    for (Int_t s = center - 1; s <= center + 1; s++)
      if (s >= s_lo && s <= s_hi)
        x += Double_t(totaldE[s]);
    Double_t y = 0.0;
    for (Int_t s = center; s <= s_hi; s++)
      y += Double_t(totaldE[s]);
    g->SetPoint(k++, x, y);
  }
  g->Set(k);
  StyleScatter(g, color);
  return g;
}

void StripScatterOverlay::MakeOverlay(
    const std::vector<TString> &filepaths, const std::vector<TString> &labels,
    ScatterBuilderFn builder, Int_t builder_arg, const TString &output_name,
    const TString &output_subdir, const TString &frame_title, Double_t xmin,
    Double_t xmax, Double_t ymin, Double_t ymax, Long64_t max_points) {
  if (filepaths.size() != labels.size() || filepaths.empty()) {
    std::cerr << "MakeOverlay: bad inputs" << std::endl;
    return;
  }

  std::vector<Int_t> palette = PlottingUtils::GetDefaultColors();
  std::vector<TFile *> files;
  std::vector<TGraph *> graphs;
  std::vector<TString> kept_labels;

  for (std::size_t i = 0; i < filepaths.size(); i++) {
    TFile *f = IO::OpenForReading(filepaths[i]);
    if (!f || f->IsZombie()) {
      std::cerr << "Error opening: " << filepaths[i] << std::endl;
      continue;
    }
    TTree *tree = static_cast<TTree *>(f->Get("events_MeV"));
    if (!tree)
      tree = static_cast<TTree *>(f->Get("events_cal"));
    if (!tree) {
      std::cerr << "No events tree in: " << filepaths[i] << std::endl;
      f->Close();
      continue;
    }

    Int_t color = palette[i % palette.size()];
    std::cout << "Building " << output_name << " scatter for " << labels[i]
              << " (" << tree->GetEntries() << " entries)..." << std::endl;
    graphs.push_back(builder(tree, color, max_points, builder_arg));
    kept_labels.push_back(labels[i]);
    files.push_back(f);
  }
  if (graphs.empty())
    return;

  // Solid-color proxy markers for the legend (the real graphs are alpha-blended
  // so their swatches would be too faint).
  std::vector<TGraph *> legend_proxies;
  for (std::size_t i = 0; i < graphs.size(); i++) {
    Int_t color = palette[i % palette.size()];
    TGraph *proxy = new TGraph(1);
    proxy->SetPoint(0, -1e9, -1e9);
    proxy->SetMarkerStyle(20);
    proxy->SetMarkerSize(1.4);
    proxy->SetMarkerColor(color);
    legend_proxies.push_back(proxy);
  }

  std::lock_guard<std::mutex> lock(g_plot_mutex);

  TString frame_name = Form("frame_%s", output_name.Data());
  TH2D *frame =
      new TH2D(frame_name, frame_title, 10, xmin, xmax, 10, ymin, ymax);
  frame->SetStats(0);
  frame->GetXaxis()->SetTitleSize(0.06);
  frame->GetYaxis()->SetTitleSize(0.06);
  frame->GetXaxis()->SetLabelSize(0.06);
  frame->GetYaxis()->SetLabelSize(0.06);
  frame->GetXaxis()->SetTitleOffset(1.2);
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->GetXaxis()->SetNdivisions(506);
  frame->GetYaxis()->SetNdivisions(506);

  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  frame->Draw();
  for (std::size_t i = 0; i < graphs.size(); i++)
    graphs[i]->Draw("P SAME");

  TLegend *legend = PlottingUtils::AddLegend(0.16, 0.40, 0.74, 0.90);
  for (std::size_t i = 0; i < graphs.size(); i++)
    legend->AddEntry(legend_proxies[i], kept_labels[i], "p");
  legend->Draw();

  PlottingUtils::SaveFigure(c, output_name, output_subdir,
                            PlotSaveOptions::kLINEAR);

  delete legend;
  delete c;
  delete frame;
  for (std::size_t i = 0; i < graphs.size(); i++)
    delete graphs[i];
  for (std::size_t i = 0; i < legend_proxies.size(); i++)
    delete legend_proxies[i];
  for (std::size_t i = 0; i < files.size(); i++)
    files[i]->Close();
}

TString StripScatterOverlay::PrettyLabel(const TString &tag) {
  TString base = RemixSim::TagWithoutStrip(tag);
  base.ReplaceAll("_eres", "");
  if (base == "aa")
    return "#alpha + #alpha";
  if (base == "an")
    return "#alpha + n";
  if (base == "beam")
    return "Beam (^{37}Cl)";
  return base;
}

void StripScatterOverlay::Run() {
  const TString project_root = Paths::DatasetDir();
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");

  std::vector<RemixSim::SimFileSpec> all_specs = RemixSim::BuildFileSpecs();

  // Split the eres datasets into reacted ones (a real reaction strip encoded in
  // the tag) and unreacted references (beam), which overlay onto every group.
  std::vector<Int_t> reacted_strip;
  std::vector<TString> reacted_files, reacted_labels;
  std::vector<TString> ref_files, ref_labels;
  for (std::size_t i = 0; i < all_specs.size(); i++) {
    const RemixSim::SimFileSpec &s = all_specs[i];
    if (!RemixSim::IsEresTag(s.tag))
      continue;
    TString file = RemixSim::SimRootPath(s);
    TString label = PrettyLabel(s.tag);
    Int_t strip = RemixSim::ReactionStripOf(s.tag);
    if (strip < 0) {
      ref_files.push_back(file);
      ref_labels.push_back(label);
    } else {
      reacted_strip.push_back(strip);
      reacted_files.push_back(file);
      reacted_labels.push_back(label);
    }
  }

  if (reacted_files.empty()) {
    std::cerr << "StripScatterOverlay: no reacted eres datasets found"
              << std::endl;
    return;
  }

  // Collect the distinct reaction strips so each gets its own overlay group;
  // mixing datasets from different reaction strips on one axis is meaningless.
  std::vector<Int_t> strips;
  for (std::size_t i = 0; i < reacted_strip.size(); i++) {
    Bool_t seen = kFALSE;
    for (std::size_t m = 0; m < strips.size(); m++)
      if (strips[m] == reacted_strip[i])
        seen = kTRUE;
    if (!seen)
      strips.push_back(reacted_strip[i]);
  }

  const Int_t s_hi = Remix::LastStrip();
  for (std::size_t g = 0; g < strips.size(); g++) {
    Int_t reaction_strip = strips[g];

    // Datasets sharing this reaction strip, plus the common beam reference.
    std::vector<TString> files, labels;
    for (std::size_t i = 0; i < reacted_strip.size(); i++)
      if (reacted_strip[i] == reaction_strip) {
        files.push_back(reacted_files[i]);
        labels.push_back(reacted_labels[i]);
      }
    for (std::size_t i = 0; i < ref_files.size(); i++) {
      files.push_back(ref_files[i]);
      labels.push_back(ref_labels[i]);
    }

    // One namespace per reaction strip; each group's datasets share that strip.
    TString subdir = Form("strip_scatter/react_s%d", reaction_strip);

    // Per-strip scatter (strip #DeltaE vs full event total) for the reaction
    // strip itself.
    MakeOverlay(
        files, labels, BuildStripScatter, reaction_strip,
        Form("stripscatter_s%d_eres", reaction_strip), subdir,
        Form(";Strip %d #DeltaE [MeV];Total #DeltaE [MeV]", reaction_strip),
        kStripEMin, kStripEMax, kTotalEMin, kTotalEMax);

    // Window sum (strip-1 + strip + strip+1) vs downstream sum (reaction strip
    // to the end of the detector).
    MakeOverlay(files, labels, BuildReactionToEndScatter, reaction_strip,
                "stripscatter_window_vs_downstream_eres", subdir,
                Form(";#DeltaE strips %d+%d+%d [MeV];#DeltaE strips "
                     "%d#rightarrow%d [MeV]",
                     reaction_strip - 1, reaction_strip, reaction_strip + 1,
                     reaction_strip, s_hi),
                3.0 * kStripEMin, 3.0 * kStripEMax, kTotalEMin, kTotalEMax);
  }
}
