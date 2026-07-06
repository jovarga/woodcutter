// Canonical rooted binary tree (PhyloTree) + the small taxa/weight helpers shared by both the
// 2-tree centroid solver and the multi-tree DP.  Everything here is in namespace maf::mast to keep
// it separate from the input-graph helpers in maf_util.hpp (which has its own `lca`, etc.).
//
// Canonical layout: leaves = ids 1..ntaxa (leaf id == taxon id), internals = ntaxa+1..2*ntaxa-1,
// root = 2*ntaxa-1, every child id < parent id. This is exactly the MAFP reverse-topological
// numbering, so building a PhyloTree from a maf::DiGraph is validate + cache-fill. 1-based indexing
// throughout (index 0 unused).
#pragma once
#include "../core.hpp"
#include <vector>
#include <cstdint>
#include <algorithm>
#include <stdexcept>
#include <cmath>

namespace maf::mast {

using NodeId = int;
using TaxonId = int;
using W = double;   // production fixes the weight type to Float64

// ---- taxa helpers ----------------------------------------------------------------------------
// Element-wise <, shorter-is-smaller on a prefix tie. Compares SORTED taxa.
inline bool lexless(const std::vector<TaxonId>& a, const std::vector<TaxonId>& b) {
  size_t lim = std::min(a.size(), b.size());
  for (size_t i = 0; i < lim; ++i) if (a[i] != b[i]) return a[i] < b[i];
  return a.size() < b.size();
}
// _sorted_taxa_union (TR:36): strictly-increasing dedup merge of two sorted taxa vectors.
inline std::vector<TaxonId> sorted_taxa_union(const std::vector<TaxonId>& l, const std::vector<TaxonId>& r) {
  std::vector<TaxonId> m; m.reserve(l.size() + r.size());
  size_t i = 0, j = 0;
  while (i < l.size() || j < r.size()) {
    TaxonId v;
    if (j >= r.size() || (i < l.size() && l[i] < r[j])) v = l[i++];
    else if (i >= l.size() || r[j] < l[i]) v = r[j++];
    else { v = l[i]; ++i; ++j; }
    if (m.empty() || m.back() != v) m.push_back(v);
  }
  return m;
}

// ---- solution reference / score unit (TR:102-112, SOL:54-67) ---------------------------------
namespace ref_kind {
  constexpr uint8_t EMPTY = 0x00, RED_GREEN_PAIR = 0x53, WHITE_WITH_LOWER = 0x56,
                    MATCHING_LEAF = 0x61, STORED_TAXA = 0x80, INVALID = 0xff;
}
struct SolutionRef {
  uint8_t kind = ref_kind::EMPTY;
  int a = 0, b = 0, c = 0;
};
inline const SolutionRef EMPTY_REF{ref_kind::EMPTY, 0, 0, 0};
inline const SolutionRef INVALID_REF{ref_kind::INVALID, 0, 0, 0};

struct ScoreRef { W score; SolutionRef ref; };

inline W invalid_score() { return -std::numeric_limits<W>::infinity(); }   // _invalid_score (SOL:65)
inline bool is_valid_score(W s) { return std::isfinite(s); }               // _is_valid_score (SOL:67)

// ---- reconstructed-solution structs (TR:134, 173) --------------------------------------------
struct MASTSolution {                 // ranked candidate (multi Branch B)
  W score;
  std::vector<TaxonId> taxa;          // sorted ascending
  SolutionRef root_ref;
};
struct RootedSolution {               // best agreement subtree rooted at one vertex
  NodeId root;
  W score;
  std::vector<TaxonId> taxa;          // sorted ascending (empty if !reconstruct)
  SolutionRef root_ref;
};

// ---- PhyloTree (TR:76) -----------------------------------------------------------------------
struct PhyloTree {
  NodeId root;
  int ntaxa;
  // 1-based arrays (index 0 unused), length nv+1 where nv = 2*ntaxa-1.
  std::vector<NodeId> parent, left, right;
  std::vector<W> internal_weight;

  int nv() const { return 2 * ntaxa - 1; }
  bool isleaf(NodeId v) const { return left[v] == 0 && right[v] == 0; }   // TR:869
};

// Build a canonical PhyloTree from a maf::DiGraph whose numbering already satisfies the canonical
// contract (the MAFP convention).  `iw` is the internal-weight vector in the "n-1 form" (length
// ntaxa-1, one per internal vertex, iw[v-ntaxa-1 (0-based)] == weight of internal vertex v) or the
// full "nv form" (length nv).  Mirrors _canonical_phylo_tree_from_graph_with_mapping (TR:674) +
// validate_binary_phylo_tree (TR:753).
inline PhyloTree build_phylo_tree(const DiGraph& g, const std::vector<W>& iw, int ntaxa) {
  const int nv = g.nv();
  if (ntaxa < 1) throw std::runtime_error("ntaxa must be positive");
  if (nv != 2 * ntaxa - 1) throw std::runtime_error("canonical trees must have 2*ntaxa-1 vertices");

  PhyloTree t;
  t.root = nv;            // canonical root id == nv
  t.ntaxa = ntaxa;
  t.parent.assign(nv + 1, 0);
  t.left.assign(nv + 1, 0);
  t.right.assign(nv + 1, 0);
  t.internal_weight.assign(nv + 1, 0.0);

  const int ninternal = nv - ntaxa;
  const bool full_form = ((int)iw.size() == nv);
  if (!full_form && (int)iw.size() != ninternal)
    throw std::runtime_error("internal_weight must have length nv or ninternal");

  for (NodeId v = 1; v <= nv; ++v) {
    const std::vector<int>& ch = g.out_neighbors(v);   // sorted ascending (DiGraph::finalize)
    if (ch.empty()) continue;                          // leaf
    if (ch.size() != 2) throw std::runtime_error("canonical internal vertices must have two children");
    t.left[v] = ch[0];
    t.right[v] = ch[1];
    t.parent[ch[0]] = v;
    t.parent[ch[1]] = v;
  }
  for (NodeId v = ntaxa + 1; v <= nv; ++v)
    t.internal_weight[v] = full_form ? iw[v - 1] : iw[v - ntaxa - 1];

  // Lightweight invariant checks; the heavy ones are implied by the MAFP numbering.
  for (NodeId v = 1; v <= nv; ++v) {
    if (v == t.root) { if (t.parent[v] != 0) throw std::runtime_error("parent[root] must be 0"); }
    else if (t.parent[v] <= v) throw std::runtime_error("children must have smaller ids than parent");
    if (v <= ntaxa) { if (!t.isleaf(v)) throw std::runtime_error("canonical taxa must be leaves"); }
    else if (t.isleaf(v)) throw std::runtime_error("canonical internals must have two children");
  }
  return t;
}

// ---- traversal-derived caches (order-independent given canonical id ordering) ----------------
// Number of leaves under each vertex. Children have smaller ids => ascending id order visits children
// before parents.
inline std::vector<int> leaf_count(const PhyloTree& t) {
  const int nv = t.nv();
  std::vector<int> c(nv + 1, 0);
  for (NodeId v = 1; v <= nv; ++v)
    c[v] = t.isleaf(v) ? 1 : c[t.left[v]] + c[t.right[v]];
  return c;
}
// Inclusive root->v prefix sum of internal_weight, accumulated in preorder. Parents have larger ids
// => descending id order visits parents before children, and the prefix value is independent of
// sibling order.
inline std::vector<W> internal_prefix(const PhyloTree& t) {
  const int nv = t.nv();
  std::vector<W> p(nv + 1, 0.0);
  for (NodeId v = nv; v >= 1; --v) {
    W base = (t.parent[v] == 0) ? 0.0 : p[t.parent[v]];
    p[v] = base + t.internal_weight[v];
  }
  return p;
}
// _prefix_before (TR:993): prefix at v's parent (0 if root).
inline W prefix_before(const std::vector<W>& prefix, const PhyloTree& t, NodeId v) {
  return t.parent[v] == 0 ? 0.0 : prefix[t.parent[v]];
}
// attached_descent_weight (TR:1009): sum on the open path strictly between ancestor and descendant.
inline W attached_descent_weight(const PhyloTree& t, const std::vector<W>& prefix, NodeId anc, NodeId desc) {
  return prefix_before(prefix, t, desc) - prefix_before(prefix, t, anc);
}
// strict_connector_weight (TR:1038): includes the ancestor vertex itself.
inline W strict_connector_weight(const PhyloTree& t, const std::vector<W>& prefix, NodeId anc, NodeId desc) {
  if (anc == desc) return 0.0;
  return prefix_before(prefix, t, desc) - prefix[anc];
}

} // namespace maf::mast
