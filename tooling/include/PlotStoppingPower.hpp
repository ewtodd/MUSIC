#ifndef PLOT_STOPPING_POWER_HPP
#define PLOT_STOPPING_POWER_HPP

#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Paths.hpp"
#include "PlottingUtils.hpp"
#include "SiCalibConstants.hpp"
#include <Rtypes.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TTree.h>
#include <iostream>
#include <vector>

class PlotStoppingPower {
public:
  static void Run();
};

#endif
