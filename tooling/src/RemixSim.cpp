#include "RemixSim.hpp"

RemixSim::Event::Event() : strip0(0.0f), strip17(0.0f) {
  for (Int_t k = 0; k < 16; k++) {
    left[k] = 0.0f;
    right[k] = 0.0f;
  }
}

Bool_t RemixSim::Event::Attach(TTree *t) {
  if (!t->GetBranch("LeftdE") || !t->GetBranch("Strip0dE")) {
    std::cerr << "RemixSim: events_MeV tree has no LeftdE/Strip0dE branches; "
                 "it predates the per-end layout and must be re-simulated"
              << std::endl;
    return kFALSE;
  }
  t->SetBranchAddress("LeftdE", left);
  t->SetBranchAddress("RightdE", right);
  t->SetBranchAddress("Strip0dE", &strip0);
  t->SetBranchAddress("Strip17dE", &strip17);
  return kTRUE;
}

Double_t RemixSim::Event::Left(Int_t s) const {
  if (s == 0)
    return strip0;
  if (s == 17)
    return strip17;
  return left[s - 1];
}

Double_t RemixSim::Event::Right(Int_t s) const {
  if (s == 0 || s == 17)
    return 0.0;
  return right[s - 1];
}

Bool_t RemixSim::SimFileSpecTagLess(const SimFileSpec &a,
                                    const SimFileSpec &b) {
  return a.tag.CompareTo(b.tag) < 0;
}

TString RemixSim::ControlDir() { return Paths::DatasetDir() + "/sim_control"; }

TString RemixSim::TagFromControlFile(const TString &filepath) {
  std::ifstream in(filepath.Data());
  if (!in.is_open())
    return "";
  TString prefix = TString("traces_") + Paths::DatasetName() + "_";
  std::string raw;
  while (std::getline(in, raw)) {
    TString line(raw.c_str());
    line = line.Strip(TString::kLeading);
    if (!line.BeginsWith("output"))
      continue;
    Ssiz_t start = line.Index(prefix);
    if (start == kNPOS)
      return "";
    TString rest = line(start, line.Length() - start);
    Ssiz_t dot = rest.Index(".root");
    if (dot == kNPOS)
      return "";
    TString base = rest(0, dot);
    base.Remove(0, prefix.Length());
    return base;
  }
  return "";
}

std::vector<RemixSim::SimFileSpec> RemixSim::BuildFileSpecs() {
  std::vector<SimFileSpec> specs;
  TString dir = ControlDir();
  void *d = gSystem->OpenDirectory(dir.Data());
  if (!d) {
    std::cerr << "BuildFileSpecs: cannot open control dir " << dir << std::endl;
    return specs;
  }
  const char *entry;
  while ((entry = gSystem->GetDirEntry(d))) {
    TString name(entry);
    if (!name.EndsWith(".toml"))
      continue;
    TString tag = TagFromControlFile(dir + "/" + name);
    if (tag.Length() > 0) {
      SimFileSpec s;
      s.tag = tag;
      specs.push_back(s);
    }
  }
  gSystem->FreeDirectory(d);
  std::sort(specs.begin(), specs.end(), SimFileSpecTagLess);
  // One simulation per class and strip: a duplicate has no defined
  // precedence, so it stops the run rather than silently picking one.
  for (Int_t i = 1; i < Int_t(specs.size()); i++) {
    const TString a = TagWithoutStrip(specs[i - 1].tag);
    const TString b = TagWithoutStrip(specs[i].tag);
    if (a == b &&
        ReactionStripOf(specs[i - 1].tag) == ReactionStripOf(specs[i].tag)) {
      std::cerr << "RemixSim: two control files in " << dir
                << " simulate the same reaction at the same strip: "
                << specs[i - 1].tag << " and " << specs[i].tag << "; keep one"
                << std::endl;
      std::exit(1);
    }
  }
  return specs;
}

TString RemixSim::TracesName(const SimFileSpec &s) {
  return Form("traces_%s_%s", Paths::DatasetName().Data(), s.tag.Data());
}

TString RemixSim::SimRootPath(const SimFileSpec &s) {
  return Paths::DatasetDir() + "/sim_root_files/" + TracesName(s) + ".root";
}

Int_t RemixSim::ReactionStripOf(const TString &tag) {
  Ssiz_t pos = tag.Last('_');
  if (pos == kNPOS)
    return -1;
  TString last = tag(pos + 1, tag.Length() - pos - 1);
  if (last.Length() < 2 || last[0] != 's')
    return -1;
  TString digits = last(1, last.Length() - 1);
  if (!digits.IsDigit())
    return -1;
  return digits.Atoi();
}

TString RemixSim::TagWithoutStrip(const TString &tag) {
  if (ReactionStripOf(tag) < 0)
    return tag;
  Ssiz_t pos = tag.Last('_');
  return tag(0, pos);
}
