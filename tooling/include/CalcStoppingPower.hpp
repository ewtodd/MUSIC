#ifndef CALC_STOPPING_POWER_HPP
#define CALC_STOPPING_POWER_HPP

#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Paths.hpp"
#include "SiCalibConstants.hpp"
#include "StoppingPowerLISE.hpp"
#include <Rtypes.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <algorithm>
#include <iostream>

class CalcStoppingPower {
public:
  static void Run();
};

#endif
