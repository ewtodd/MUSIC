#include "Constants.hpp"
#include "StripSumScatter.hpp"

int main() {
  // The one epoch the chain works on, when the config names one: its file
  // tag, channel map and decode flags are what the events files carry.
  Constants::ActivateAnalysisEpoch();
  StripSumScatter scatter;
  scatter.Run();
  return 0;
}
