#ifndef SELECTION_DIAGRAM_HPP
#define SELECTION_DIAGRAM_HPP

#include "StripSumScatter.hpp"
#include <TString.h>
#include <vector>

/**
 * @file SelectionDiagram.hpp
 * @brief Block diagram of the event selection, generated from the
 *        configuration.
 *
 * Serializes the steps StripSumScatter::DescribeSelection() lists, in their
 * order, as a Graphviz DOT graph and as a Mermaid flowchart: one box per
 * step with its short name and the condition to pass, the numbers filled in
 * from the configuration; every cut has a reject exit. Only what the
 * configuration applies is drawn: a step that is off is left out. The
 * stages (input, event-level cuts, per reaction strip, outcome) are
 * clusters, and the reaction-strip range is in the loop marker that opens
 * the per-strip stage.
 *
 * Both are text and always written. The DOT lays the stages out in three
 * columns (input and event-level cuts, per strip, outcome) so the figure is
 * a landscape page rather than a strip, and is rendered to PNG and PDF when
 * `dot` (Graphviz) is on the PATH, which it is in the dataset shells. The
 * Mermaid flowchart is top to bottom and renders wherever Markdown does
 * (GitHub, GitLab, VS Code, Obsidian).
 */
namespace SelectionDiagram {

/// @brief The steps as a Graphviz DOT graph.
/// @param steps From StripSumScatter::DescribeSelection().
/// @param title Graph title, e.g. "37Cl event selection".
TString Dot(const std::vector<SelectionStep> &steps, const TString &title);

/// @brief The steps as a Mermaid flowchart.
/// @param steps From StripSumScatter::DescribeSelection().
/// @param title Chart title (a `---` front-matter block).
TString Mermaid(const std::vector<SelectionStep> &steps, const TString &title);

/**
 * @brief Write `selection.dot` and `selection.mmd` under `dir`, and render
 *        the DOT to `selection.png` and `selection.pdf` when `dot` is
 *        available.
 * @param steps From StripSumScatter::DescribeSelection().
 * @param dir   Output directory, created if missing.
 * @param title Diagram title.
 */
void Write(const std::vector<SelectionStep> &steps, const TString &dir,
           const TString &title);

} // namespace SelectionDiagram

#endif
