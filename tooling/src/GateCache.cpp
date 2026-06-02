#include "GateCache.hpp"

Double_t GateCache::NSigma() { return BeamFitUtils::GateNSigma(); }

TString GateCache::FileSubpath(Int_t run) {
  return Form("Gated_Run%d.root", run);
}

TString GateCache::FileName(Int_t run) {
  return IO::GetRootFilesBaseDir() + "/" + FileSubpath(run);
}

Bool_t GateCache::Exists(Int_t run) {
  return !gSystem->AccessPathName(FileName(run));
}

void GateCache::SaveGatePlot(TH2F *gate_hist, const BeamFit2D &gate,
                             const TString &subdir) {
  std::lock_guard<std::mutex> lock(g_plot_mutex);
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureAndDraw2DHistogram(gate_hist, c);
  BeamFitUtils::DrawGateEllipse(gate, BeamFitUtils::GateNSigma());
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
  TH2F *gate_hist = BeamFitUtils::BuildGateHist(
      src_chain, Form("h2_gate_S0_L1_%s", tag.Data()));
  BeamFit2D gate = BeamFitUtils::FitBigPeak(gate_hist, tag);
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
    if (!BeamFitUtils::InEllipse(gate, xg, yg, BeamFitUtils::GateNSigma()))
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
  params[7] = BeamFitUtils::GateNSigma();

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

  std::vector<Int_t> run_order;
  std::map<Int_t, TChain *> chain_by_run =
      FileSet::GroupCalSidecarsByRun(run_order);

  for (std::size_t i = 0; i < run_order.size(); i++)
    EnsureForRun(run_order[i], chain_by_run[run_order[i]]);

  for (std::size_t i = 0; i < run_order.size(); i++)
    delete chain_by_run[run_order[i]];
}
