#include "BinaryUtils.hpp"
#include <iostream>
#include <sstream>

int main(int argc, char *argv[]) {
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0]
              << " <input.sol> <output_dir> <chunk_seconds>" << std::endl;
    return 1;
  }

  const char *inputFile = argv[1];
  const char *outputDir = argv[2];
  Double_t chunkSeconds = std::stod(argv[3]);

  std::cout << "Splitting " << inputFile << " into " << chunkSeconds
            << "s chunks..." << std::endl;

  Int_t totalBlocks = 0;
  Int_t totalChunks = 0;
  std::vector<TString> outputFiles = SOLReader::SplitSolFileByTime(
      inputFile, outputDir, chunkSeconds, totalBlocks, totalChunks);

  if (outputFiles.empty()) {
    std::cerr << "ERROR: Failed to split file" << std::endl;
    return 1;
  }

  std::cout << "Output files:" << std::endl;
  for (std::vector<TString>::const_iterator it = outputFiles.begin();
       it != outputFiles.end(); ++it) {
    std::cout << "  " << *it << std::endl;
  }

  return 0;
}
