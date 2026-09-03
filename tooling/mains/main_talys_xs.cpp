// talys-xs: Hauser-Feshbach cross sections from TALYS for this dataset's
// reaction, written as ROOT graphs cross-section overlays.
//
// Runs TALYS in normal kinematics (alpha on the beam nucleus as target; the
// cross section is the same either way) over a grid of centre-of-mass
// energies spanning the strips the cross section reports and the reference
// table, once per configured model, and writes every residual-production
// channel TALYS produced as a TGraph of E_cm [MeV] vs sigma [mb], named by
// TALYS's own file stem (rp039090 for Z=39, A=90), into one directory per
// model (m0, m1, ...) of root_files/talys/talys_xs.root, with the model's
// label and exact input as TNameds. Which channels make up (a,xn) is the
// analysis's business, so cross-section sums them itself. Like srim-cache in
// the simulator, this is the one place the external code is driven;
// everything downstream only reads the file.
//
// TALYS comes from the talys-nix flake input, whose path the build compiles
// in; TALYS_BIN in the environment overrides it.
#include "BeamEnergies.hpp"
#include "Constants.hpp"
#include "InitUtils.hpp"
#include "Paths.hpp"
#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TMath.h>
#include <TNamed.h>
#include <TString.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifndef MUSIC_TALYS_BIN
#define MUSIC_TALYS_BIN ""
#endif

namespace {

// Grid step and the margin past the reported strips' energies. TALYS is
// cheap at these energies, so the grid errs on the fine side.
const Double_t kEcmStep = 0.25;
const Double_t kEcmMargin = 1.0;

// One TALYS residual-production file: E [MeV] and xs [mb] rows after the
// YANDF header, with the lab energy converted back to centre-of-mass.
TGraph *ReadResidual(const TString &path, Double_t cm_per_lab) {
  std::ifstream in(path.Data());
  if (!in)
    return nullptr;
  std::vector<Double_t> e, xs;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream ss(line);
    Double_t el = 0.0, s = 0.0;
    if (!(ss >> el >> s))
      continue;
    e.push_back(el * cm_per_lab);
    xs.push_back(s);
  }
  if (e.empty())
    return nullptr;
  return new TGraph(Int_t(e.size()), &e[0], &xs[0]);
}

// Run one model in its own work directory and write its graphs into `dir`.
// Returns the number of residual channels written, -1 on a TALYS failure.
Int_t RunModel(const TalysModel &model, const TString &work,
               const TString &talys, const TString &energies,
               Double_t cm_per_lab, TDirectory *dir) {
  const CrossSectionConfig &X = Constants::cfg.CROSS_SECTION_CONFIG;
  gSystem->mkdir(work, kTRUE);
  TString input;
  input += "# Written by talys-xs for " + Paths::DatasetName() + ": " +
           model.label + "\n";
  input += "projectile a\n";
  input += "element " + X.BEAM_ELEMENT + "\n";
  input += Form("mass %d\n", X.BEAM_A);
  input += "energy energies\n";
  for (Int_t k = 0; k < Int_t(model.keywords.size()); k++)
    input += model.keywords[k] + "\n";
  {
    std::ofstream en((work + "/energies").Data());
    en << energies;
    std::ofstream inp((work + "/talys.inp").Data());
    inp << input;
  }
  for (Int_t k = 0; k < Int_t(model.keywords.size()); k++)
    std::cout << "    " << model.keywords[k] << std::endl;
  const TString cmd = Form("cd '%s' && '%s' < talys.inp > talys.out 2>&1",
                           work.Data(), talys.Data());
  if (std::system(cmd.Data()) != 0) {
    std::cerr << "talys-xs: TALYS failed; see " << work << "/talys.out"
              << std::endl;
    return -1;
  }

  dir->cd();
  TNamed("label", model.label.Data()).Write();
  TNamed("input", input.Data()).Write();
  Int_t n_graphs = 0;
  TSystemDirectory sd("work", work);
  TList *files = sd.GetListOfFiles();
  for (TIter it(files); TObject *o = it();) {
    const TString name = o->GetName();
    Int_t z = 0, a = 0;
    if (!name.EndsWith(".tot") || sscanf(name.Data(), "rp%3d%3d", &z, &a) != 2)
      continue;
    TGraph *g = ReadResidual(work + "/" + name, cm_per_lab);
    if (!g)
      continue;
    g->SetName(Form("rp%03d%03d", z, a));
    g->SetTitle(Form("Z=%d A=%d residual production;E_{c.m.} [MeV];#sigma "
                     "[mb]",
                     z, a));
    g->Write();
    delete g;
    n_graphs++;
  }
  delete files;
  return n_graphs;
}

} // namespace

int main() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                Paths::ResultsDir() + "/plots",
                                Paths::ResultsDir() + "/root_files");
  const CrossSectionConfig &X = Constants::cfg.CROSS_SECTION_CONFIG;
  if (X.BEAM_A <= 0 || X.BEAM_Z <= 0 || X.BEAM_ELEMENT.IsNull() ||
      X.BEAM_SIM_FILE.IsNull() || X.TALYS_MODELS.empty()) {
    std::cerr << "talys-xs: this dataset's CROSS_SECTION_CONFIG needs BEAM_A, "
                 "BEAM_Z, BEAM_ELEMENT, BEAM_SIM_FILE and at least one "
                 "TALYS model"
              << std::endl;
    return 1;
  }
  const Char_t *env = gSystem->Getenv("TALYS_BIN");
  TString talys = (env && env[0] != '\0') ? TString(env) : MUSIC_TALYS_BIN;
  if (talys.IsNull())
    talys = "talys";

  // The grid spans the strips the cross section reports, with a margin, at
  // the energies it reports them at, and the reference table so the curve
  // is drawn under every point.
  Double_t dE[18];
  Double_t e_strip0 = 0.0;
  if (!BeamEnergies::Profile(BeamEnergies::SimPath(), dE, e_strip0)) {
    std::cerr << "talys-xs: cannot read the simulated beam at "
              << BeamEnergies::SimPath() << std::endl;
    return 1;
  }
  const Double_t cm_frac = BeamEnergies::CmFraction();
  // TALYS runs the reaction the other way round, alpha on the beam nucleus,
  // so its lab energy is E_cm * (A_b + A_t) / A_b.
  const Double_t lab_per_cm =
      Double_t(X.BEAM_A + TargetGasA(X.TARGET_GAS)) / Double_t(X.BEAM_A);
  Double_t e_hi =
      BeamEnergies::LabAtStrip(dE, e_strip0, X.XS_STRIP_MIN) * cm_frac +
      kEcmMargin;
  Double_t e_lo = TMath::Max(
      kEcmStep, (BeamEnergies::LabAtStrip(dE, e_strip0, X.XS_STRIP_MAX) -
                 dE[X.XS_STRIP_MAX]) *
                        cm_frac -
                    kEcmMargin);
  for (Int_t k = 0; k < Int_t(X.REFERENCE_XS.size()); k++) {
    e_lo = TMath::Min(e_lo, X.REFERENCE_XS[k][0] - kEcmMargin);
    e_hi = TMath::Max(e_hi, X.REFERENCE_XS[k][0] + kEcmMargin);
  }
  const Double_t ecm_lo =
      kEcmStep * std::floor(TMath::Max(kEcmStep, e_lo) / kEcmStep);
  const Double_t ecm_hi = kEcmStep * std::ceil(e_hi / kEcmStep);
  TString energies;
  for (Double_t ecm = ecm_lo; ecm <= ecm_hi + 1.0e-9; ecm += kEcmStep)
    energies += Form("%.4f\n", ecm * lab_per_cm);

  // Everything TALYS touches lives under root_files/talys, which git
  // ignores: one run per model under work/, the graphs beside them.
  const TString talys_dir = Paths::ResultsDir() + "/root_files/talys";
  const TString out_path = talys_dir + "/talys_xs.root";
  gSystem->mkdir(talys_dir, kTRUE);
  TFile out(out_path, "RECREATE");
  if (out.IsZombie()) {
    std::cerr << "talys-xs: cannot write " << out_path << std::endl;
    return 1;
  }
  TNamed("talys", talys.Data()).Write();
  std::cout << "talys-xs: alpha + " << X.BEAM_A << X.BEAM_ELEMENT << ", E_cm "
            << ecm_lo << ".." << ecm_hi << " MeV (strips " << X.XS_STRIP_MIN
            << ".." << X.XS_STRIP_MAX << "), " << X.TALYS_MODELS.size()
            << " model(s), " << talys << std::endl;
  for (Int_t m = 0; m < Int_t(X.TALYS_MODELS.size()); m++) {
    const TalysModel &model = X.TALYS_MODELS[m];
    std::cout << "  m" << m << ": " << model.label << std::endl;
    TDirectory *dir = out.mkdir(Form("m%d", m));
    const Int_t n = RunModel(model, talys_dir + Form("/work/m%d", m), talys,
                             energies, 1.0 / lab_per_cm, dir);
    if (n <= 0) {
      if (n == 0)
        std::cerr << "talys-xs: no residual-production files for "
                  << model.label << std::endl;
      return 1;
    }
    std::cout << "    " << n << " residual channels" << std::endl;
  }
  out.Close();
  std::cout << "talys-xs: -> " << out_path << std::endl;
  return 0;
}
