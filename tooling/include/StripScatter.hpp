#ifndef STRIP_SCATTER_HPP
#define STRIP_SCATTER_HPP

#include "BeamFit2D.hpp"
#include "Constants.hpp"
#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <Rtypes.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TEllipse.h>
#include <TH2F.h>
#include <TString.h>
#include <TSystem.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <mutex>
#include <vector>

class StripScatter {
public:
  static void Run();

private:
  static void SaveH2(TH2F *h, const TString &save_name, const TString &subdir);
  static void SaveGatePlot(TH2F *gate_hist, const BeamFit2D &gate,
                           const TString &subdir);
  static void ProcessRun(Int_t run, TChain *chain, const TString &subdir);
};

#endif
