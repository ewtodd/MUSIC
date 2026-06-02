#ifndef PIPELINE_HPP
#define PIPELINE_HPP

#include "BinaryToRoot.hpp"
#include "BinaryUtils.hpp"
#include "CalibrateBeam.hpp"
#include "Constants.hpp"
#include "EventBuilder.hpp"
#include "FileSet.hpp"
#include "GpuAccel.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Timing.hpp"
#include "TraceCreator.hpp"
#include <TError.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <queue>
#include <set>
#include <streambuf>
#include <thread>
#include <utility>
#include <vector>

class Pipeline {
public:
  static void Run();
};

#endif
