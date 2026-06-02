#ifndef SI_CALIBRATION_HPP
#define SI_CALIBRATION_HPP

#include "Calibration.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Paths.hpp"
#include "PlottingUtils.hpp"
#include "SiCalibConstants.hpp"
#include <Rtypes.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TMath.h>
#include <TObject.h>
#include <TString.h>
#include <TTree.h>
#include <cstdio>
#include <iostream>
#include <map>
#include <vector>

class SiCalibration {
public:
  static void Run();
};

#endif
