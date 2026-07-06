// Newick parsing, instance construction, and solution output.
#pragma once
#include "core.hpp"
#include "instance.hpp"
#include <string>
#include <vector>
#include <istream>
#include <ostream>
#include <functional>
#include <stdexcept>
#include <cctype>
#include <utility>

namespace maf {

// Recursive-descent Newick parser. `p` is the read cursor into `s`; `nleaves` is incremented for
// every leaf encountered.
inline RecTree parse_newick(const std::string& s, size_t& p, int& nleaves) {
  if (p < s.size() && std::isdigit((unsigned char)s[p])) {
    size_t start = p;
    while (p < s.size() && std::isdigit((unsigned char)s[p])) ++p;
    int label = std::stoi(s.substr(start, p - start));
    nleaves += 1;
    return RecTree::leaf(label);
  } else if (p < s.size() && s[p] == '(') {
    ++p;  // consume '('
    std::vector<RecTree> children;
    while (p < s.size() && s[p] != ')') {
      if (s[p] == ',') ++p;
      children.push_back(parse_newick(s, p, nleaves));
    }
    if (p >= s.size()) throw std::runtime_error("unterminated '(' in Newick");
    ++p;  // consume ')'
    return RecTree::node(std::move(children));
  }
  throw std::runtime_error("expected number or '(' in Newick");
}

// Convert a RecTree into a DiGraph with 2n-1 vertices: leaves keep their labels 1..n, internal nodes
// are numbered in post-order starting at n+1, so the root is 2n-1. Returns (vertex id of this
// subtree's root, next free number).
inline std::pair<int, int> build_edges(DiGraph& g, const RecTree& t, int i) {
  if (t.is_leaf) return {t.label, i};
  std::vector<int> child_ids;
  for (const RecTree& c : t.children) {
    auto pr = build_edges(g, c, i);
    child_ids.push_back(pr.first);
    i = pr.second;
  }
  for (int cid : child_ids) g.add_edge(i, cid);
  return {i, i + 1};
}

inline DiGraph rec_to_graph(const RecTree& tree, int n) {
  DiGraph g(2 * n - 1);
  auto pr = build_edges(g, tree, n + 1);
  if (pr.first != 2 * n - 1 || pr.second != 2 * n)
    throw std::runtime_error("root must be the last node");
  g.finalize();
  return g;
}

inline std::string strip(const std::string& s) {
  size_t a = 0, b = s.size();
  while (a < b && std::isspace((unsigned char)s[a])) ++a;
  while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
  return s.substr(a, b - a);
}

// Read a MAFP instance: lines starting with '#' (and blank lines) are skipped; each remaining line
// is one Newick tree terminated by ';'.  All trees must share the same leaf count.
inline MAFP read_instance(std::istream& in) {
  std::vector<RecTree> trees;
  int n = -1;
  std::string line;
  while (std::getline(in, line)) {
    std::string l = strip(line);
    if (l.empty() || l[0] == '#') continue;
    size_t p = 0; int ntree = 0;
    RecTree t = parse_newick(l, p, ntree);
    if (p >= l.size() || l[p] != ';')
      throw std::runtime_error("missing semicolon in input");
    if (n != -1 && n != ntree) throw std::runtime_error("trees have differently sized leaf-set");
    trees.push_back(std::move(t));
    n = ntree;
  }
  MAFP maf;
  maf.n = (n < 0 ? 0 : n);
  for (const RecTree& t : trees) maf.T.push_back(rec_to_graph(t, maf.n));
  return maf;
}

// ---- solution output (PACE format) -----------------------------------------------------------
// Write a RecTree as a Newick string (no edge labels, no trailing ';').
inline void show_newick(std::ostream& os, const RecTree& t) {
  if (t.is_leaf) { os << t.label; return; }
  os << '(';
  for (size_t i = 0; i < t.children.size(); ++i) { if (i) os << ','; show_newick(os, t.children[i]); }
  os << ')';
}

// Induced rooted subtree of `tree` on leaf set S, with degree-2 (unary) nodes suppressed, as a
// RecTree carrying ORIGINAL leaf labels. Build down from the root, keep only vertices on S's
// root-paths, and collapse any node left with a single kept child. Tree topology is unordered for the
// checker.
inline RecTree induced_rectree(const DiGraph& tree, const BitSet& S, int n) {
  std::vector<char> keep(tree.nv() + 1, 0);
  for (int x : S) { int v = x; keep[v] = 1; while (tree.in_degree(v) == 1) { v = tree.in_neighbors(v)[0]; keep[v] = 1; } }
  std::function<RecTree(int)> build = [&](int v) -> RecTree {
    if (tree.out_degree(v) == 0) return RecTree::leaf(v);   // original leaf label
    std::vector<RecTree> ch;
    for (int w : tree.out_neighbors(v)) if (keep[w]) ch.push_back(build(w));
    if (ch.size() == 1) return std::move(ch[0]);            // suppress unary node
    return RecTree::node(std::move(ch));
  };
  return build(2 * n - 1);                                  // root = 2n-1
}

// Write one Newick component per line, induced in tree T[1].
inline void write_forest(std::ostream& os, const MAFP& maf, const std::vector<BitSet>& sol) {
  for (const BitSet& comp : sol) {
    show_newick(os, induced_rectree(maf.T[0], comp, maf.n));
    os << ";\n";
  }
}

} // namespace maf
