#ifndef PIPELINE_HPP
#define PIPELINE_HPP

#include <TString.h>
#include <mutex>
#include <vector>

inline std::mutex g_plot_mutex;

// One sim file. `tag` becomes both the on-disk basename
// (root_files/traces_37Cl_<tag>.root) and the plot subdirectory leaf.
struct FileSpec {
  TString tag;
};

inline std::vector<FileSpec> BuildFileSpecs() {
  return {{"aa"}, {"an"}, {"beam"}, {"aa_eres"}, {"an_eres"}, {"beam_eres"}};
}

inline TString TracesName(const FileSpec &s) {
  return Form("traces_37Cl_%s", s.tag.Data());
}

inline TString FileLabel(const FileSpec &s) {
  return Form("37Cl_%s", s.tag.Data());
}

#endif
