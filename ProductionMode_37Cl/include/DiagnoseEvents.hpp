#ifndef DIAGNOSE_EVENTS_HPP
#define DIAGNOSE_EVENTS_HPP

#include "Constants.hpp"
#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

class DiagnoseEvents {
public:
  static void Run();
};

#endif
