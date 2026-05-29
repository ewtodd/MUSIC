#ifndef DIAGNOSE_TIMING_HPP
#define DIAGNOSE_TIMING_HPP

#include "Constants.hpp"
#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

class DiagnoseTiming {
public:
  static void Run();
};

#endif
