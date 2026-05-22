#include "Constants.hpp"
#include "InitUtils.hpp"
#include "Pipeline.hpp"
#include "TraceCreator.hpp"
#include <TString.h>
#include <iostream>
#include <vector>

void Traces(TString file_label = "") {
  const TString project_root = Paths::ProjectRootOf(__FILE__);
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");

  std::vector<FileSpec> specs;
  if (file_label.IsNull()) {
    specs = BuildFileSpecs();
    if (specs.empty()) {
      std::cerr << "No file specs from BuildFileSpecs()" << std::endl;
      return;
    }
  } else {
    FileSpec s = ResolveFileSpec(file_label);
    if (s.run < 0) {
      std::cerr << "Could not resolve file label '" << file_label << "'"
                << std::endl;
      return;
    }
    specs.push_back(s);
  }

  std::vector<TString> events_names;
  std::vector<TString> file_labels;
  for (Int_t k = 0; k < Int_t(specs.size()); k++) {
    events_names.push_back(EventsName(specs[k]));
    file_labels.push_back(FileLabel(specs[k]));
  }

  BuildTraces(events_names, file_labels, /*save_plots=*/kTRUE,
              /*reprocess=*/kTRUE);
}
