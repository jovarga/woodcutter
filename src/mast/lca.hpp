// Binary-lifting LCA over a canonical rooted binary tree.  Templated on the tree type so it serves
// both PhyloTree (multi path) and TreeView (2-tree path); both expose root, left/right/parent (1-based)
// and isleaf(v). Flat ancestor table `up[(level-1)*nv + v]`.
#pragma once
#include "trees.hpp"
#include <vector>

namespace maf::mast {

struct LCAIndex {
  std::vector<int> depth;   // 1-based, root depth 0
  std::vector<NodeId> up;   // flat, 1-based: up[(level-1)*nv + v]; level 1 = immediate parent
  int levels;
};

template <class Tree>
inline int _nverts(const Tree& t) { return (int)t.parent.size() - 1; }   // 1-based arrays

template <class Tree>
inline LCAIndex build_lca(const Tree& t) {
  const int nv = _nverts(t);
  int levels = 1;
  while ((1 << levels) <= nv) ++levels;

  LCAIndex idx;
  idx.levels = levels;
  idx.depth.assign(nv + 1, 0);
  idx.up.assign((size_t)nv * levels + 1, 0);   // 1-based; index (level-1)*nv+v in [1, levels*nv]

  std::vector<NodeId> stack;
  stack.reserve(nv);
  stack.push_back(t.root);
  while (!stack.empty()) {
    NodeId v = stack.back(); stack.pop_back();
    if (!t.isleaf(v)) {
      NodeId lc = t.left[v], rc = t.right[v];
      idx.depth[lc] = idx.depth[v] + 1;
      idx.depth[rc] = idx.depth[v] + 1;
      idx.up[lc] = v;     // level 1 == (1-1)*nv + lc == lc
      idx.up[rc] = v;
      stack.push_back(rc);   // push right then left -> left processed first
      stack.push_back(lc);
    }
  }
  for (int level = 2; level <= levels; ++level) {
    size_t off = (size_t)(level - 1) * nv, prev = (size_t)(level - 2) * nv;
    for (NodeId v = 1; v <= nv; ++v) {
      NodeId p = idx.up[prev + v];
      idx.up[off + v] = (p == 0) ? 0 : idx.up[prev + p];
    }
  }
  return idx;
}

template <class Tree>
inline NodeId lca_query(const Tree& t, const LCAIndex& idx, NodeId a, NodeId b) {
  const int nv = (int)idx.depth.size() - 1;   // == nverts
  if (idx.depth[a] < idx.depth[b]) std::swap(a, b);
  int diff = idx.depth[a] - idx.depth[b];
  int level = 1;
  while (diff != 0) {
    if (diff & 1) a = idx.up[(size_t)(level - 1) * nv + a];
    diff >>= 1;
    ++level;
  }
  if (a == b) return a;
  for (int lv = idx.levels; lv >= 1; --lv) {
    size_t off = (size_t)(lv - 1) * nv;
    NodeId ua = idx.up[off + a], ub = idx.up[off + b];
    if (ua != ub && ua != 0 && ub != 0) { a = ua; b = ub; }
  }
  return idx.up[a];   // level 1 == parent of a
}

} // namespace maf::mast
