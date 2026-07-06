// Graph/topology helpers on the input trees, plus solution validation.
#pragma once
#include "core.hpp"
#include "instance.hpp"
#include <vector>
#include <algorithm>

namespace maf {

// Path from `v` up to the root, returned root-first.
inline std::vector<int> root_path(const DiGraph& g, int v) {
  std::vector<int> path{v};
  while (g.in_degree(v) >= 1) { v = g.in_neighbors(v)[0]; path.push_back(v); }
  std::reverse(path.begin(), path.end());
  return path;
}

// Inner nodes of the subtree spanned by leaves `xs` (returns ascending vertex list).
// Works for any range of ints (BitSet, vector, initializer_list).
template <class Range>
inline std::vector<int> inner_vs(const DiGraph& g, const Range& xs) {
  BitSet vs;
  for (int x : xs) for (int u : root_path(g, x)) vs.insert(u);
  int v = g.nv();
  // Peel away the path above the spanning subtree's root: while exactly one child of v is in vs.
  while (true) {
    int cnt = 0, child = -1;
    for (int c : g.out_neighbors(v)) if (vs.contains(c)) { ++cnt; child = c; }
    if (cnt != 1) break;
    vs.erase(v); v = child;
  }
  for (int x : xs) vs.erase(x);
  return vs.to_vec();
}
inline std::vector<int> inner_vs(const DiGraph& g, std::initializer_list<int> xs) {
  return inner_vs<std::initializer_list<int>>(g, xs);
}

// Depth (1-based) at which two root-paths first diverge. Identical prefixes return min(len)+1.
inline int divergence_depth(const std::vector<int>& a, const std::vector<int>& b) {
  size_t m = std::min(a.size(), b.size());
  for (size_t i = 0; i < m; ++i) if (a[i] != b[i]) return (int)(i + 1);
  return (int)m + 1;
}

// The "odd one out" of a rooted triple (the leaf that splits off highest), or 0 if unresolved.
inline int higher_node(const DiGraph& g, int x1, int x2, int x3) {
  std::vector<int> rp1 = root_path(g, x1), rp2 = root_path(g, x2), rp3 = root_path(g, x3);
  int d23 = divergence_depth(rp2, rp3), d13 = divergence_depth(rp1, rp3), d12 = divergence_depth(rp1, rp2);
  if (d23 > d13 && d13 == d12) return x1;
  if (d13 > d12 && d12 == d23) return x2;
  if (d12 > d23 && d23 == d13) return x3;
  return 0;  // all diverge at the same depth (binary tree => genuine 3-way / unresolved)
}

// A rooted triple is compatible iff all trees agree on its topology.
inline bool is_compatible(const MAFP& maf, int x1, int x2, int x3) {
  int h0 = higher_node(maf.T[0], x1, x2, x3);
  for (size_t i = 1; i < maf.T.size(); ++i)
    if (higher_node(maf.T[i], x1, x2, x3) != h0) return false;
  return true;
}

// desc[v] = leaves below v.  Post-order vertex numbering => children have smaller ids than their
// parent, so ascending iteration is reverse topological order.
inline std::vector<BitSet> compute_descendants(const DiGraph& tree) {
  int N = tree.nv();
  std::vector<BitSet> desc(N + 1);
  for (int v = 1; v <= N; ++v) {
    if (tree.out_degree(v) == 0) { desc[v].insert(v); }
    else for (int w : tree.out_neighbors(v)) desc[v].union_with(desc[w]);
  }
  return desc;
}

// Lowest common ancestor of a,b given precomputed descendant sets; 0 == "nothing" (different
// subtrees / not comparable).
inline int lca(const DiGraph& tree, const std::vector<BitSet>& desc, int a, int b) {
  int x = a;
  if (desc[x].contains(b)) return x;
  while (tree.in_degree(x) == 1) {
    x = tree.in_neighbors(x)[0];
    if (!desc[x].contains(a)) return 0;
    if (desc[x].contains(b)) return x;
  }
  return 0;
}

// ---- solution validation ---------------------------------------------------------------------
inline bool is_partition_of_leaves(const std::vector<BitSet>& sol, int n) {
  BitSet covered;
  for (const BitSet& comp : sol) {
    if (comp.empty()) return false;
    for (int leaf : comp) {
      if (leaf < 1 || leaf > n) return false;
      if (covered.contains(leaf)) return false;
      covered.insert(leaf);
    }
  }
  BitSet all; for (int i = 1; i <= n; ++i) all.insert(i);
  return covered == all;
}

inline bool is_valid_component(const MAFP& maf, const BitSet& comp) {
  std::vector<int> v = comp.to_vec();  // ascending
  int m = (int)v.size();
  for (int i = 0; i < m - 2; ++i)
    for (int j = i + 1; j < m - 1; ++j)
      for (int k = j + 1; k < m; ++k)
        if (!is_compatible(maf, v[i], v[j], v[k])) return false;
  return true;
}

inline bool has_disjoint_steiner_vertices(const DiGraph& tree, const std::vector<BitSet>& sol) {
  BitSet used;
  for (const BitSet& comp : sol) {
    std::vector<int> iv = inner_vs(tree, comp);
    BitSet inner(iv.begin(), iv.end());
    if (!isdisjoint(used, inner)) return false;
    used.union_with(inner);
  }
  return true;
}

inline bool is_valid_solution(const MAFP& maf, const std::vector<BitSet>& sol) {
  if (!is_partition_of_leaves(sol, maf.n)) return false;
  for (const BitSet& c : sol) if (!is_valid_component(maf, c)) return false;
  for (const DiGraph& t : maf.T) if (!has_disjoint_steiner_vertices(t, sol)) return false;
  return true;
}

} // namespace maf
