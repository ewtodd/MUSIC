#include "GateCache.hpp"

namespace {

// Gate: 2D bigaus fit on Strip1 long-side (L, LeftdEMeV[1]) vs Strip2
// (R, RightdEMeV[2]). Mirrors the k=1 ellipse gate used in CalibrateBeam.
const Double_t kGateMin = 0.0;
const Double_t kGateMax = 12.0;
const Int_t kGateBins = 240;
const Double_t kGateNSigma = 2.0;
const Int_t kSeedHalfBins = 40;
const Double_t kSeedFrac = 0.30;

} // namespace

Double_t GateCache::NSigma() { return kGateNSigma; }

TString GateCache::FileSubpath(Int_t run) {
  return Form("Gated_Run%d.root", run);
}

TString GateCache::FileName(Int_t run) {
  return IO::GetRootFilesBaseDir() + "/" + FileSubpath(run);
}

Bool_t GateCache::Exists(Int_t run) {
  return !gSystem->AccessPathName(FileName(run));
}

TH2F *GateCache::BuildGateHist(TChain *chain, const TString &name) {
  Float_t leftdE[18] = {0};
  Float_t rightdE[18] = {0};
  chain->SetBranchStatus("*", 0);
  chain->SetBranchStatus("LeftdEMeV", 1);
  chain->SetBranchStatus("RightdEMeV", 1);
  chain->SetBranchAddress("LeftdEMeV", leftdE);
  chain->SetBranchAddress("RightdEMeV", rightdE);

  TH2F *h = new TH2F(name, ";#DeltaE S0 [MeV];#DeltaE L1 [MeV]", kGateBins,
                     kGateMin, kGateMax, kGateBins, kGateMin, kGateMax);
  h->SetDirectory(nullptr);

  Long64_t n = chain->GetEntries();
  for (Long64_t j = 0; j < n; j++) {
    chain->GetEntry(j);
    Double_t x = Double_t(leftdE[1]);
    Double_t y = Double_t(rightdE[2]);
    if (x > 0 && y > 0)
      h->Fill(x, y);
  }
  return h;
}

BeamFit2D GateCache::FitBigPeak(TH2F *h, const TString &tag) {
  BeamFit2D out;
  if (!h || h->GetEntries() < 100)
    return out;

  Double_t bw_x = h->GetXaxis()->GetBinWidth(1);
  Double_t bw_y = h->GetYaxis()->GetBinWidth(1);

  Int_t bx = 0, by = 0, bz = 0;
  h->GetMaximumBin(bx, by, bz);
  Double_t peak_val = h->GetBinContent(bx, by);
  if (peak_val <= 0)
    return out;

  Int_t lo_bx = std::max(1, bx - kSeedHalfBins);
  Int_t hi_bx = std::min(h->GetNbinsX(), bx + kSeedHalfBins);
  Int_t lo_by = std::max(1, by - kSeedHalfBins);
  Int_t hi_by = std::min(h->GetNbinsY(), by + kSeedHalfBins);

  Moments2D m = BeamFitUtils::ComputeMoments(h, lo_bx, hi_bx, lo_by, hi_by,
                                             kSeedFrac * peak_val, bw_x, bw_y);
  if (m.weight <= 0) {
    std::cerr << "  FitBigPeak " << tag << ": no bins above seed threshold"
              << std::endl;
    return out;
  }

  Double_t x_lo = std::max(h->GetXaxis()->GetXmin(), m.mu_x - 3.0 * m.sigma_x);
  Double_t x_hi = std::min(h->GetXaxis()->GetXmax(), m.mu_x + 3.0 * m.sigma_x);
  Double_t y_lo = std::max(h->GetYaxis()->GetXmin(), m.mu_y - 3.0 * m.sigma_y);
  Double_t y_hi = std::min(h->GetYaxis()->GetXmax(), m.mu_y + 3.0 * m.sigma_y);

  out = BeamFitUtils::FitBigausInWindow(h, x_lo, x_hi, y_lo, y_hi, m, peak_val,
                                        Form("f_gate_%s", tag.Data()),
                                        TString("FitBigPeak ") + tag);

  std::cout << "  gate fit " << tag << ": mu=(" << out.mu_x << "," << out.mu_y
            << ") sigma=(" << out.sigma_x << "," << out.sigma_y
            << ") rho=" << out.rho << std::endl;
  return out;
}

void GateCache::SaveGatePlot(TH2F *gate_hist, const BeamFit2D &gate,
                             const TString &subdir) {
  std::lock_guard<std::mutex> lock(g_plot_mutex);
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureAndDraw2DHistogram(gate_hist, c);
  if (gate.ok) {
    Double_t major, minor, theta_deg;
    BeamFitUtils::DiagonalizeCov(gate, major, minor, theta_deg);
    TEllipse *e = new TEllipse(gate.mu_x, gate.mu_y, kGateNSigma * major,
                               kGateNSigma * minor, 0.0, 360.0, theta_deg);
    e->SetFillStyle(0);
    e->SetLineColor(kRed + 1);
    e->SetLineWidth(2);
    e->Draw();
  }
  PlottingUtils::SaveFigure(c, "gate_S0_vs_L1", subdir,
                            PlotSaveOptions::kLINEAR);
  delete c;
}

void GateCache::BuildForRun(Int_t run, TChain *src_chain) {
  if (!src_chain || src_chain->GetEntries() == 0) {
    std::cerr << "Run " << run << ": empty source chain; no gate cache"
              << std::endl;
    return;
  }

  TString tag = Form("run%d", run);
  TH2F *gate_hist =
      BuildGateHist(src_chain, Form("h2_gate_S0_L1_%s", tag.Data()));
  BeamFit2D gate = FitBigPeak(gate_hist, tag);
  SaveGatePlot(gate_hist, gate, Form("gate/run%d", run));
  if (!gate.ok) {
    std::cerr << "Run " << run << ": gate fit failed; no gate cache"
              << std::endl;
    delete gate_hist;
    return;
  }

  Float_t leftdE[18] = {0};
  Float_t rightdE[18] = {0};
  Float_t totaldE[18] = {0};
  Float_t cathode = 0;
  src_chain->SetBranchStatus("*", 0);
  src_chain->SetBranchStatus("LeftdEMeV", 1);
  src_chain->SetBranchStatus("RightdEMeV", 1);
  src_chain->SetBranchStatus("TotaldEMeV", 1);
  src_chain->SetBranchStatus("CathodeMeV", 1);
  src_chain->SetBranchAddress("LeftdEMeV", leftdE);
  src_chain->SetBranchAddress("RightdEMeV", rightdE);
  src_chain->SetBranchAddress("TotaldEMeV", totaldE);
  src_chain->SetBranchAddress("CathodeMeV", &cathode);

  TFile *out = IO::OpenForWriting(FileSubpath(run), "RECREATE");
  if (!out || out->IsZombie()) {
    std::cerr << "Run " << run << ": cannot open gate cache for writing"
              << std::endl;
    if (out)
      delete out;
    delete gate_hist;
    return;
  }
  TTree *tree =
      new TTree("events_gated", Form("gate-passing events, run %d", run));
  tree->Branch("LeftdEMeV", leftdE, "LeftdEMeV[18]/F");
  tree->Branch("RightdEMeV", rightdE, "RightdEMeV[18]/F");
  tree->Branch("TotaldEMeV", totaldE, "TotaldEMeV[18]/F");
  tree->Branch("CathodeMeV", &cathode, "CathodeMeV/F");

  Long64_t n = src_chain->GetEntries();
  Long64_t kept = 0;
  std::cout << "Run " << run << ": gating " << n << " events..." << std::endl;
  for (Long64_t j = 0; j < n; j++) {
    src_chain->GetEntry(j);
    Double_t xg = Double_t(leftdE[1]);
    Double_t yg = Double_t(rightdE[2]);
    if (!(xg > 0 && yg > 0))
      continue;
    if (!BeamFitUtils::InEllipse(gate, xg, yg, kGateNSigma))
      continue;
    tree->Fill();
    kept++;
  }

  TVectorD params(8);
  params[0] = gate.ok ? 1.0 : 0.0;
  params[1] = gate.amp;
  params[2] = gate.mu_x;
  params[3] = gate.mu_y;
  params[4] = gate.sigma_x;
  params[5] = gate.sigma_y;
  params[6] = gate.rho;
  params[7] = kGateNSigma;

  out->cd();
  tree->Write();
  params.Write("gate_params");
  gate_hist->Write("h2_gate_S0_L1");
  out->Close();
  delete out;
  delete gate_hist;

  std::cout << "Run " << run << ": cached " << kept << "/" << n
            << " gate-passing events to " << FileSubpath(run) << std::endl;
}

void GateCache::EnsureForRun(Int_t run, TChain *src_chain) {
  if (Constants::SKIP_EXISTING && Exists(run)) {
    std::cout << "[gate-cache] run " << run << ": " << FileSubpath(run)
              << " exists; reusing" << std::endl;
    return;
  }
  BuildForRun(run, src_chain);
}

Bool_t GateCache::LoadGate(TFile *f, BeamFit2D &gate) {
  if (!f)
    return kFALSE;
  TVectorD *p = static_cast<TVectorD *>(f->Get("gate_params"));
  if (!p || p->GetNoElements() < 7)
    return kFALSE;
  gate.ok = ((*p)[0] != 0.0);
  gate.amp = (*p)[1];
  gate.mu_x = (*p)[2];
  gate.mu_y = (*p)[3];
  gate.sigma_x = (*p)[4];
  gate.sigma_y = (*p)[5];
  gate.rho = (*p)[6];
  return gate.ok;
}

void GateCache::Run() {
  const TString project_root = Paths::DatasetDir();
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");

  std::vector<FileSpec> all_specs = FileSet::BuildFileSpecs();
  std::map<Int_t, TChain *> chain_by_run;
  std::vector<Int_t> run_order;
  for (std::size_t i = 0; i < all_specs.size(); i++) {
    const FileSpec &s = all_specs[i];
    TString full = IO::GetRootFilesBaseDir() + "/" + FileSet::CalSidecarName(s);
    if (gSystem->AccessPathName(full)) {
      std::cerr << "Missing cal sidecar: " << full << std::endl;
      continue;
    }
    if (chain_by_run.find(s.run) == chain_by_run.end()) {
      chain_by_run[s.run] = new TChain("events_cal");
      run_order.push_back(s.run);
    }
    chain_by_run[s.run]->Add(full);
  }

  for (std::size_t i = 0; i < run_order.size(); i++)
    EnsureForRun(run_order[i], chain_by_run[run_order[i]]);

  for (std::size_t i = 0; i < run_order.size(); i++)
    delete chain_by_run[run_order[i]];
}
