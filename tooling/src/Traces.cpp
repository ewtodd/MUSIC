#include "Traces.hpp"

void Traces::Run(const TString &file_label) {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                Paths::ResultsDir() + "/plots",
                                Paths::ResultsDir() + "/root_files");

  std::vector<FileSpec> specs;
  if (file_label.IsNull()) {
    specs = FileSet::BuildProcessedFileSpecs();
    if (specs.empty()) {
      std::cerr << "No file specs from FileSet::BuildProcessedFileSpecs()"
                << std::endl;
      return;
    }
  } else {
    FileSpec s = FileSet::ResolveFileSpec(file_label);
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
    events_names.push_back(FileSet::EventsName(specs[k]));
    file_labels.push_back(FileSet::FileLabel(specs[k]));
  }

  TraceCreator::BuildTraces(events_names, file_labels, /*save_plots=*/kTRUE,
                            /*reprocess=*/kTRUE);
}
