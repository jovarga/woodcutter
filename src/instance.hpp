// MAFP (a MAF instance) and RecTree (recursive tree used during Newick parsing).
#pragma once
#include "core.hpp"
#include <vector>

namespace maf {

// Recursive tree: leaves are integer labels; internal nodes hold children.
struct RecTree {
  bool is_leaf = false;
  int label = 0;                       // valid when is_leaf
  std::vector<RecTree> children;       // valid when !is_leaf

  static RecTree leaf(int l) { RecTree t; t.is_leaf = true; t.label = l; return t; }
  static RecTree node(std::vector<RecTree> ch) { RecTree t; t.is_leaf = false; t.children = std::move(ch); return t; }
};

// A Maximum Agreement Forest instance with `n` leaves and trees `T`.
// Each tree is a directed graph; vertices 1..n are the leaves (indexed by label) and vertices are
// in reverse topological order, so vertex 2n-1 is the root.
struct MAFP {
  int n = 0;
  std::vector<DiGraph> T;
};

} // namespace maf
