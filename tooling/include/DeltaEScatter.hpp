#ifndef DELTA_E_SCATTER_HPP
#define DELTA_E_SCATTER_HPP

#include "BeamFit2D.hpp"
#include "Constants.hpp"
#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TChain.h>
#include <TEllipse.h>
#include <TF2.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <mutex>
#include <vector>

class DeltaEScatter {
public:
  static void Run();
};

#endif
