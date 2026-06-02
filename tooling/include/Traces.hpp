#ifndef TRACES_HPP
#define TRACES_HPP

#include "Constants.hpp"
#include "FileSet.hpp"
#include "InitUtils.hpp"
#include "TraceCreator.hpp"
#include <TString.h>
#include <iostream>
#include <vector>

class Traces {
public:
  static void Run(const TString &file_label = "");
};

#endif
