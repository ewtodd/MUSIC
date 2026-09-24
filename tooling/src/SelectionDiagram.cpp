#include "SelectionDiagram.hpp"
#include <TSystem.h>
#include <fstream>
#include <iostream>
#include <sstream>

namespace {

// Greedy word wrap to about `width` characters per line, lines joined by
// `br`. Node labels are one short line otherwise and the boxes grow wide.
TString Wrap(const TString &text, Int_t width, const char *br) {
  TString out, line, tok;
  Ssiz_t from = 0;
  while (text.Tokenize(tok, from, " ")) {
    if (line.Length() > 0 && line.Length() + 1 + tok.Length() > width) {
      out += line + br;
      line = "";
    }
    if (line.Length() > 0)
      line += " ";
    line += tok;
  }
  return out + line;
}

TString EscapeHtml(const TString &s) {
  TString e(s);
  e.ReplaceAll("&", "&amp;");
  e.ReplaceAll("<", "&lt;");
  e.ReplaceAll(">", "&gt;");
  e.ReplaceAll("\"", "&quot;");
  return e;
}

// Mermaid labels sit in double quotes; entity codes stand in for the
// characters its parser would otherwise read.
TString EscapeMermaid(const TString &s) {
  TString e(s);
  e.ReplaceAll("&", "&amp;");
  e.ReplaceAll("<", "&lt;");
  e.ReplaceAll(">", "&gt;");
  e.ReplaceAll("\"", "#quot;");
  return e;
}

// Only what the configuration applies is drawn: a step that is off is not
// in the diagram at all.
std::vector<SelectionStep> Live(const std::vector<SelectionStep> &steps) {
  std::vector<SelectionStep> live;
  for (size_t i = 0; i < steps.size(); i++)
    if (steps[i].on)
      live.push_back(steps[i]);
  return live;
}

// A step that filters: it has a reject exit. Input, the loop marker and the
// outcome nodes do not.
Bool_t IsCut(const SelectionStep &st) {
  return (st.pre_cut >= 0 || (st.tag_cut >= 0 && st.tag_cut != kTagPass));
}

const char *RejectLabel(const SelectionStep &st) {
  return st.stage == SelectionStep::kEventLevel ? "drop event"
                                                : "no tag at reac";
}

const char *StageTitle(SelectionStep::Stage s) {
  switch (s) {
  case SelectionStep::kInput:
    return "input";
  case SelectionStep::kEventLevel:
    return "event-level cuts";
  case SelectionStep::kPerStrip:
    return "per reaction strip";
  case SelectionStep::kOutcome:
    return "outcome";
  }
  return "";
}

// The DOT layout runs in three columns, so the figure is a landscape page
// rather than a strip: input and the event-level cuts, the per-strip
// conditions, the outcome.
Int_t Column(SelectionStep::Stage s) {
  switch (s) {
  case SelectionStep::kInput:
  case SelectionStep::kEventLevel:
    return 0;
  case SelectionStep::kPerStrip:
    return 1;
  case SelectionStep::kOutcome:
    return 2;
  }
  return 0;
}

// The chain: every step feeds the next within the input, event-level and
// per-strip stages and across their boundaries; the outcome nodes are wired
// by id: "tagged" from the last per-strip step, "xs" from "tagged" and
// "n_beam". The beam count's own origin is in its text: an edge from the
// event-level column to it would cross the per-strip column.
struct Edge {
  TString from, to;
  Bool_t across; // between columns
};

std::vector<Edge> Chain(const std::vector<SelectionStep> &steps) {
  std::vector<Edge> edges;
  TString prev, last_strip;
  Int_t prev_col = -1;
  for (size_t i = 0; i < steps.size(); i++) {
    const SelectionStep &st = steps[i];
    const Int_t col = Column(st.stage);
    if (st.stage == SelectionStep::kOutcome) {
      if (st.id == "tagged" && last_strip.Length())
        edges.push_back({last_strip, st.id, kTRUE});
      else if (st.id == "xs") {
        edges.push_back({"tagged", st.id, kFALSE});
        edges.push_back({"n_beam", st.id, kFALSE});
      }
      continue;
    }
    if (prev.Length())
      edges.push_back({prev, st.id, col != prev_col});
    prev = st.id;
    prev_col = col;
    if (st.stage == SelectionStep::kPerStrip)
      last_strip = st.id;
  }
  return edges;
}

} // namespace

TString SelectionDiagram::Dot(const std::vector<SelectionStep> &all,
                              const TString &title) {
  const std::vector<SelectionStep> steps = Live(all);
  const Int_t kWrap = 44;
  std::ostringstream s;
  s << "digraph selection {" << std::endl;
  // Orthogonal edges: the column connectors run as clean L shapes through
  // the gaps between the columns.
  s << Form("  graph [rankdir=TB, newrank=true, splines=ortho, nodesep=0.25, "
            "ranksep=0.3, "
            "fontname=\"Helvetica\", fontsize=13, labelloc=t, "
            "label=<<B>%s</B>>];",
            EscapeHtml(title).Data())
    << std::endl;
  s << "  node [shape=box, style=\"rounded,filled\", fillcolor=\"#f4f4f4\", "
       "color=\"#555555\", fontname=\"Helvetica\", fontsize=10, "
       "margin=\"0.16,0.08\"];"
    << std::endl;
  s << "  edge [color=\"#555555\", fontname=\"Helvetica\", fontsize=9, "
       "arrowsize=0.7];"
    << std::endl;
  for (Int_t stage = SelectionStep::kInput; stage <= SelectionStep::kOutcome;
       stage++) {
    s << Form("  subgraph cluster_%d {", stage) << std::endl;
    s << Form("    label=\"%s\"; labeljust=\"l\"; fontname=\"Helvetica-Bold\"; "
              "fontsize=11; color=\"#aaaaaa\"; style=\"rounded\";",
              StageTitle(SelectionStep::Stage(stage)))
      << std::endl;
    for (size_t i = 0; i < steps.size(); i++) {
      const SelectionStep &st = steps[i];
      if (st.stage != stage)
        continue;
      const TString name = EscapeHtml(st.name);
      const TString detail = Wrap(EscapeHtml(st.detail), kWrap, "<BR/>");
      const char *fill = IsCut(st)           ? "#f4f4f4"
                         : st.id == "tagged" ? "#e3f2e1"
                                             : "#e8eef7";
      s << Form("    %s [fillcolor=\"%s\", label=<<B>%s</B><BR/>%s>];",
                st.id.Data(), fill, name.Data(), detail.Data())
        << std::endl;
      if (IsCut(st)) {
        // The reject exit beside the box, on the same rank and to its left,
        // so the gap to the right of each column stays free for the
        // connector to the next one (a flat edge keeps its tail on the
        // left; dir=back puts the arrowhead on the exit).
        s << Form("    rej_%s [shape=plaintext, style=\"\", fontsize=9, "
                  "fontcolor=\"#b00020\", label=\"%s\"];",
                  st.id.Data(), RejectLabel(st))
          << std::endl;
        s << Form("    rej_%s -> %s [dir=back, color=\"#b00020\", "
                  "arrowsize=0.6];",
                  st.id.Data(), st.id.Data())
          << std::endl;
        s << Form("    {rank=same; rej_%s; %s;}", st.id.Data(), st.id.Data())
          << std::endl;
      }
    }
    s << "  }" << std::endl;
  }
  // The columns: the first node of each on one rank, ordered left to right
  // by an invisible flat edge; the edges between columns carry no rank.
  TString heads;
  Int_t last_col = -1;
  for (size_t i = 0; i < steps.size(); i++) {
    const Int_t col = Column(steps[i].stage);
    if (col == last_col)
      continue;
    heads += (heads.Length() ? " -> " : "") + steps[i].id;
    last_col = col;
  }
  s << Form("  {rank=same; %s [style=invis];}", heads.Data()) << std::endl;
  const std::vector<Edge> edges = Chain(steps);
  for (size_t i = 0; i < edges.size(); i++)
    s << Form("  %s -> %s%s;", edges[i].from.Data(), edges[i].to.Data(),
              edges[i].across ? " [constraint=false]" : "")
      << std::endl;
  s << "}" << std::endl;
  return s.str();
}

TString SelectionDiagram::Mermaid(const std::vector<SelectionStep> &all,
                                  const TString &title) {
  const std::vector<SelectionStep> steps = Live(all);
  const Int_t kWrap = 44;
  std::ostringstream s;
  s << "---" << std::endl;
  s << "title: " << title << std::endl;
  s << "---" << std::endl;
  s << "flowchart TB" << std::endl;
  s << "  classDef cut fill:#f4f4f4,stroke:#555,color:#000;" << std::endl;
  s << "  classDef note fill:#e8eef7,stroke:#555,color:#000;" << std::endl;
  s << "  classDef tagged fill:#e3f2e1,stroke:#555,color:#000;" << std::endl;
  s << "  classDef rej fill:none,stroke:none,color:#b00020;" << std::endl;
  for (Int_t stage = SelectionStep::kInput; stage <= SelectionStep::kOutcome;
       stage++) {
    s << Form("  subgraph stage_%d[\"%s\"]", stage,
              StageTitle(SelectionStep::Stage(stage)))
      << std::endl;
    s << "    direction TB" << std::endl;
    for (size_t i = 0; i < steps.size(); i++) {
      const SelectionStep &st = steps[i];
      if (st.stage != stage)
        continue;
      const TString name = EscapeMermaid(st.name);
      const TString detail = Wrap(EscapeMermaid(st.detail), kWrap, "<br/>");
      const char *cls = IsCut(st)           ? "cut"
                        : st.id == "tagged" ? "tagged"
                                            : "note";
      s << Form("    %s[\"<b>%s</b><br/>%s\"]:::%s", st.id.Data(), name.Data(),
                detail.Data(), cls)
        << std::endl;
      if (IsCut(st))
        s << Form("    rej_%s[\"%s\"]:::rej", st.id.Data(), RejectLabel(st))
          << std::endl;
    }
    s << "  end" << std::endl;
  }
  const std::vector<Edge> edges = Chain(steps);
  for (size_t i = 0; i < edges.size(); i++)
    s << Form("  %s --> %s", edges[i].from.Data(), edges[i].to.Data())
      << std::endl;
  for (size_t i = 0; i < steps.size(); i++)
    if (IsCut(steps[i]))
      s << Form("  %s -.-> rej_%s", steps[i].id.Data(), steps[i].id.Data())
        << std::endl;
  return s.str();
}

void SelectionDiagram::Write(const std::vector<SelectionStep> &steps,
                             const TString &dir, const TString &title) {
  gSystem->mkdir(dir, kTRUE);
  const TString base = dir + "/selection";
  {
    std::ofstream f((base + ".dot").Data());
    if (!f) {
      std::cerr << "selection diagram: cannot write " << base << ".dot"
                << std::endl;
      return;
    }
    f << Dot(steps, title);
  }
  {
    std::ofstream f((base + ".mmd").Data());
    if (f)
      f << Mermaid(steps, title);
  }
  // Graphviz renders the DOT when it is about; the text files stand alone
  // otherwise.
  char *dot = gSystem->Which(gSystem->Getenv("PATH"), "dot");
  if (!dot) {
    std::cout << "selection diagram: " << base << ".dot and .mmd written "
              << "(no `dot` on the PATH, not rendered)" << std::endl;
    return;
  }
  const char *fmt[2] = {"png", "pdf"};
  Bool_t ok = kTRUE;
  for (Int_t k = 0; k < 2; k++) {
    const TString cmd =
        Form("\"%s\" -T%s%s -o \"%s.%s\" \"%s.dot\"", dot, fmt[k],
             k == 0 ? " -Gdpi=150" : "", base.Data(), fmt[k], base.Data());
    if (gSystem->Exec(cmd) != 0)
      ok = kFALSE;
  }
  delete[] dot;
  std::cout << "selection diagram: " << base << ".{dot,mmd"
            << (ok ? ",png,pdf}" : "} (dot failed to render)") << std::endl;
}
