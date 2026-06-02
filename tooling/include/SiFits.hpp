#ifndef SI_FITS_HPP
#define SI_FITS_HPP

#include "FittingUtils.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Paths.hpp"
#include "PlottingUtils.hpp"
#include "SiCalibConstants.hpp"
#include <Rtypes.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>
#include <TObject.h>
#include <TString.h>
#include <TTree.h>
#include <iostream>

class SiFits {
public:
  static void Run();
};

#endif
