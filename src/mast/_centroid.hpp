// HEURISTIC TRACK ONLY.  This is the O(n log n) 2-tree centroid MAST engine (centroid_mast_solve),
// recovered from git (removed from the exact track in b7df81d in favour of the simpler O(n^2)
// quad_dp).  It is included ONLY in heuristic-track builds (from heuristics.hpp) and drives the
// heuristic's MAST-peel.  The exact track deliberately keeps quad_dp (checkable/deterministic).  A
// CMake flag can override the heuristic back to quad_dp (see MAF_HEUR_MAST / MAF_HEUR_QUAD_MAST);
// by default the heuristic uses THIS centroid engine.  centroid_mast_solve is interface-identical
// to quad_mast_solve, so the two are drop-in interchangeable.
//
// Returns the per-rooted-node best agreement subtrees as (score, taxa) pairs, filtered to |taxa|>1,
// in ascending vertex id, tree1 before tree2.
//
// PRODUCTION SCOPE only: allow_empty=false, solution_limit=1, partition=nothing, blacklist=nothing,
// reconstruct=true.  Every materialized solution stores taxa (STORED_TAXA index into recon store).
// The attached engine is also used when all internal weights are 0; lifting by 0 is a no-op.
//
// See docs/CENTROID_SPEC.md for design details.
#pragma once
#include "trees.hpp"
#include "lca.hpp"
#include <vector>
#include <utility>
#include <cstdint>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <functional>
#include <memory>
#include <cmath>

namespace maf::mast {

// ============================================================================================
//  Local value carriers (production: reconstruct=true => taxa stored by index in a recon store)
// ============================================================================================
// A ScoreRef bundles a score with a 1-based index into the reconstruction store (0 == "no taxa":
// used for the invalid sentinel and the empty ref). Copied refs keep their taxa stable even if the
// source scratch vector is recycled. The score==-inf marks invalid.
struct CScoreRef {
  W score;
  int taxa = 0;   // 1-based recon index; 0 means empty/invalid taxa
};

// ScoreStates: length-4 array indexed by the 2-bit attachment mask (ATT:1-9, TR:447).
// bit0 (value 1) = tree1 attached upward, bit1 (value 2) = tree2 attached (ATT:11-12).
struct CScoreStates {
  CScoreRef s[4];
  CScoreRef& operator[](int m) { return s[m]; }
  const CScoreRef& operator[](int m) const { return s[m]; }
};

// ---- TreeView (TR:92): a (sub-)tree over LOCAL node ids, 1-based (index 0 unused). ----
struct TreeView {
  uint8_t tree_id = 0;          // 1 or 2 = which original tree
  NodeId root = 0;              // local
  std::vector<NodeId> orig;     // local -> original id (1-based; orig[0] unused)
  std::vector<NodeId> parent, left, right;  // local
  std::vector<TaxonId> taxon;   // local -> taxon (0 for internal)

  int nlocal() const { return (int)parent.size() - 1; }   // 1-based arrays
  bool isleaf(NodeId v) const { return left[v] == 0 && right[v] == 0; }   // TR:869
};

// ---- Centroid types (CEN:1-13) ----
struct CentroidPath {
  NodeId head = 0;
  std::vector<NodeId> nodes;        // heavy path top->bottom
  std::vector<NodeId> side_roots;   // light child per path node (terminal: own id)
};
struct CentroidInfo {
  std::vector<int> path_id;             // per local vertex -> its path (1-based)
  std::vector<int> path_pos;            // per local vertex -> 1-based position on its path
  std::vector<CentroidPath> paths;      // paths[0] unused; path 1 = root path => paths index 1..
  std::vector<NodeId> path_head_of;     // per local vertex -> head of its path
  std::vector<int> postorder_path_ids;  // paths children-before-parent
};

// ---- Weighted search tree over right positions (AM:1-10) ----
struct WeightedSearchTree {
  int root = 0;
  std::vector<int> left, right, parent;       // post-order ids 1..nnodes; 0 == null
  std::vector<int> leaf_for_right_pos;        // right_pos (1-based) -> leaf node id
  std::vector<int> leaf_lo, leaf_hi;          // covered position range per node, inclusive
  bool is_leaf(int node) const { return left[node] == 0 && right[node] == 0; }  // AM:186
};

// ============================================================================================
//  The solver context: holds instance topology, prefix sums, LCAs, recon store, accumulators.
// ============================================================================================
struct CentroidSolver {
  const PhyloTree& t1;
  const PhyloTree& t2;
  // taxon_weight is all 1.0 (production); leaf_score(taxon) == 1.0.  We keep it implicit.

  std::vector<W> iprefix1, iprefix2;   // internal_prefix per original tree (preorder), 1-based
  LCAIndex lca1, lca2;

  // Reconstruction store (SOL:87-132): recon[idx-1] is the sorted-deduped taxa for STORED_TAXA(idx).
  std::vector<std::vector<TaxonId>> recon;

  // Rooted accumulators (SOL:50-51/620), indexed by ORIGINAL vertex id (1-based, [0] unused).
  std::vector<CScoreRef> rooted_tree1, rooted_tree2;

  // Epoch-stamped scratch (SOL:771-801) — taxon membership + taxon->b-leaf.
  std::vector<int> taxon_marks;        // size ntaxa+1
  int taxon_mark_epoch = 0;
  std::vector<NodeId> leaf_for_taxon;  // taxon -> b-local leaf
  std::vector<int> leaf_for_taxon_epoch;
  int leaf_for_taxon_mark_epoch = 0;

  // Segment-tree cache keyed by nright (ATT:558-589); lookup/insert only.
  std::unordered_map<int, WeightedSearchTree> tree_cache;

  CentroidSolver(const PhyloTree& a, const PhyloTree& b) : t1(a), t2(b) {
    iprefix1 = internal_prefix(t1);
    iprefix2 = internal_prefix(t2);
    lca1 = build_lca(t1);
    lca2 = build_lca(t2);
    int n1 = t1.nv(), n2 = t2.nv();
    CScoreRef inv{invalid_score(), 0};
    rooted_tree1.assign(n1 + 1, inv);
    rooted_tree2.assign(n2 + 1, inv);
    int nt = t1.ntaxa;
    taxon_marks.assign(nt + 1, 0);
    leaf_for_taxon.assign(nt + 1, 0);
    leaf_for_taxon_epoch.assign(nt + 1, 0);
  }

  const PhyloTree& tree_for(const TreeView& v) const { return v.tree_id == 1 ? t1 : t2; }
  const std::vector<W>& prefix_for(const TreeView& v) const { return v.tree_id == 1 ? iprefix1 : iprefix2; }
  const LCAIndex& lca_for(const TreeView& v) const { return v.tree_id == 1 ? lca1 : lca2; }

  // ---- recon store helpers (SOL:87-132) ----
  CScoreRef store_taxa(W score, std::vector<TaxonId>&& taxa) {
    recon.push_back(std::move(taxa));
    return CScoreRef{score, (int)recon.size()};   // 1-based index
  }
  CScoreRef score_ref_from_taxon(W score, TaxonId taxon) {
    return store_taxa(score, std::vector<TaxonId>{taxon});
  }
  // taxa for a ref (REC:1): STORED_TAXA -> stored vector; index 0 (empty/invalid) -> empty.
  const std::vector<TaxonId>& taxa_view(const CScoreRef& r) const {
    static const std::vector<TaxonId> empty;
    return r.taxa == 0 ? empty : recon[r.taxa - 1];
  }
  // _score_ref_union (SOL:117): merge taxa of refs (sorted union), store, return STORED_TAXA.
  CScoreRef score_ref_union2(W score, const CScoreRef& a, const CScoreRef& b) {
    std::vector<TaxonId> merged = sorted_taxa_union(taxa_view(a), taxa_view(b));
    return store_taxa(score, std::move(merged));
  }

  // ---- _better_ref (SOL:947): max score, then lexless of sorted taxa on tie. ----
  bool better_ref(const CScoreRef& cand, const CScoreRef& cur) const {
    if (!is_valid_score(cur.score)) return is_valid_score(cand.score);
    if (!is_valid_score(cand.score)) return false;
    if (cand.score != cur.score) return cand.score > cur.score;
    return lexless(taxa_view(cand), taxa_view(cur));   // reconstruct=true
  }
  // _choose_ref (ATT:27): keep LEFT on tie (better_ref false on tie).
  CScoreRef choose_ref(const CScoreRef& l, const CScoreRef& r) const {
    return better_ref(l, r) ? l : r;
  }
  CScoreStates choose_states(const CScoreStates& l, const CScoreStates& r) const {
    CScoreStates out;
    for (int m = 0; m < 4; ++m) out.s[m] = choose_ref(l.s[m], r.s[m]);
    return out;
  }

  // ---- rooted-slot writers (SOL:635/648) ----
  void update_rooted_tree1(NodeId root, const CScoreRef& cand) {
    if (!is_valid_score(cand.score)) return;
    if (better_ref(cand, rooted_tree1[root])) rooted_tree1[root] = cand;
  }
  void update_rooted_tree2(NodeId root, const CScoreRef& cand) {
    if (!is_valid_score(cand.score)) return;
    if (better_ref(cand, rooted_tree2[root])) rooted_tree2[root] = cand;
  }

  // ---- epoch helpers (SOL:771-801) ----
  int next_taxon_mark_epoch() {
    if (taxon_mark_epoch == std::numeric_limits<int>::max()) { std::fill(taxon_marks.begin(), taxon_marks.end(), 0); taxon_mark_epoch = 1; }
    else ++taxon_mark_epoch;
    return taxon_mark_epoch;
  }
  int next_leaf_for_taxon_epoch() {
    if (leaf_for_taxon_mark_epoch == std::numeric_limits<int>::max()) { std::fill(leaf_for_taxon_epoch.begin(), leaf_for_taxon_epoch.end(), 0); leaf_for_taxon_mark_epoch = 1; }
    else ++leaf_for_taxon_mark_epoch;
    return leaf_for_taxon_mark_epoch;
  }

  // ---- connector weights via view (SOL:847) ----
  W attached_descent(const TreeView& view, NodeId parent_local, NodeId child_local) const {
    if (parent_local == child_local) return 0.0;   // _attached_descent_weight (ATT:379) guard
    const PhyloTree& t = tree_for(view); const std::vector<W>& p = prefix_for(view);
    return attached_descent_weight(t, p, view.orig[parent_local], view.orig[child_local]);
  }
  W strict_connector(const TreeView& view, NodeId parent_local, NodeId child_local) const {
    const PhyloTree& t = tree_for(view); const std::vector<W>& p = prefix_for(view);
    return strict_connector_weight(t, p, view.orig[parent_local], view.orig[child_local]);
  }
  // prefix_before of a local node's original id (for path bases) (ATT:117).
  W prefix_before_local(const TreeView& view, NodeId local) const {
    const PhyloTree& t = tree_for(view); const std::vector<W>& p = prefix_for(view);
    return prefix_before(p, t, view.orig[local]);
  }

  // The engine entry (recursive); see below.
  struct AttachedSubproblemResult;
  AttachedSubproblemResult solve_attached_agree(const TreeView& a, const TreeView& b);
};

// ============================================================================================
//  Mask / score helpers (ATT:14-115)  — free functions taking the solver for tie-breaks.
// ============================================================================================
inline CScoreRef invalid_ref() { return CScoreRef{invalid_score(), 0}; }
inline CScoreStates invalid_states() { CScoreStates s; CScoreRef iv = invalid_ref(); for (int m=0;m<4;++m) s.s[m]=iv; return s; }
inline CScoreStates same_states(const CScoreRef& r) { CScoreStates s; for (int m=0;m<4;++m) s.s[m]=r; return s; }

// _add_score (ATT:31): invalid in -> invalid out unchanged; else score+delta, ref preserved.
inline CScoreRef add_score(const CScoreRef& r, W delta) {
  if (!is_valid_score(r.score)) return r;
  return CScoreRef{r.score + delta, r.taxa};
}
inline bool has_valid_state(const CScoreStates& s) {
  return is_valid_score(s.s[0].score) || is_valid_score(s.s[1].score) ||
         is_valid_score(s.s[2].score) || is_valid_score(s.s[3].score);
}
// _lift_tree1 (ATT:60): add delta to masks {1,3}.  _lift_tree2 (ATT:69): masks {2,3}.
inline CScoreStates lift_tree1(const CScoreStates& st, W d) {
  CScoreStates o = st; o.s[1] = add_score(st.s[1], d); o.s[3] = add_score(st.s[3], d); return o;
}
inline CScoreStates lift_tree2(const CScoreStates& st, W d) {
  CScoreStates o = st; o.s[2] = add_score(st.s[2], d); o.s[3] = add_score(st.s[3], d); return o;
}
// _adjust_states (ATT:109): lift_tree2(lift_tree1(states,l), r) — tree1 first.
inline CScoreStates adjust_states(const CScoreStates& st, W l, W r) {
  return lift_tree2(lift_tree1(st, l), r);
}
inline CScoreRef adjust_ref(const CScoreRef& r, W d) { return add_score(r, d); }

// _lift_between_tree1/2 (ATT:78/89).
inline CScoreStates lift_between_tree1(const CentroidSolver& cs, const TreeView& v, const CScoreStates& st, NodeId anc_local, NodeId desc_local) {
  const PhyloTree& t = cs.tree_for(v); const std::vector<W>& p = cs.prefix_for(v);
  W d = attached_descent_weight(t, p, v.orig[anc_local], v.orig[desc_local]);
  return lift_tree1(st, d);
}
inline CScoreStates lift_between_tree2(const CentroidSolver& cs, const TreeView& v, const CScoreStates& st, NodeId anc_local, NodeId desc_local) {
  const PhyloTree& t = cs.tree_for(v); const std::vector<W>& p = cs.prefix_for(v);
  W d = attached_descent_weight(t, p, v.orig[anc_local], v.orig[desc_local]);
  return lift_tree2(st, d);
}
// _attached_join_weight (ATT:152): leaf->0; else (node_weight + connL) + connR, left-assoc.
inline W attached_join_weight(const CentroidSolver& cs, const TreeView& view, NodeId node) {
  if (view.isleaf(node)) return 0.0;
  NodeId lc = view.left[node], rc = view.right[node];
  const PhyloTree& t = cs.tree_for(view);
  W nw = t.internal_weight[view.orig[node]];   // node_weight (TR:963)
  return nw + cs.strict_connector(view, node, lc) + cs.strict_connector(view, node, rc);
}

// _attached_leaf_states (ATT:184): reconstruct -> score_ref_from_taxon; wrapped via same_states.
inline CScoreStates attached_leaf_states(CentroidSolver& cs, W score, TaxonId taxon) {
  return same_states(cs.score_ref_from_taxon(score, taxon));
}

// ============================================================================================
//  internal_prefix-derived "node base" path bases (ATT:117/863)
// ============================================================================================
inline std::vector<W> path_prefix_bases(const CentroidSolver& cs, const TreeView& view, const std::vector<NodeId>& nodes) {
  std::vector<W> bases(nodes.size());
  for (size_t i = 0; i < nodes.size(); ++i) bases[i] = cs.prefix_before_local(view, nodes[i]);
  return bases;
}

// ============================================================================================
//  leaf counts / postorder / preorder ranges / depths
// ============================================================================================
inline std::vector<int> tv_leaf_count(const TreeView& t) {
  // TR:920 postorder fold; for a TreeView, local ids are post-order ranks (children before parents
  // in side/induced views) but the value is sibling-order-independent => any order with children
  // first works.  Use an explicit postorder to be safe for arbitrary local id orderings.
  int n = t.nlocal();
  std::vector<int> c(n + 1, 0);
  std::vector<NodeId> stack; stack.reserve(n);
  stack.push_back(t.root);
  // two-pass via sign marker
  std::vector<NodeId> order; order.reserve(n);
  while (!stack.empty()) {
    NodeId v = stack.back(); stack.pop_back();
    if (v < 0) { NodeId node = -v; c[node] = t.isleaf(node) ? 1 : c[t.left[node]] + c[t.right[node]]; continue; }
    stack.push_back(-v);
    if (!t.isleaf(v)) { stack.push_back(t.right[v]); stack.push_back(t.left[v]); }
  }
  return c;
}

// _postorder! (TR:873): children-before-parent, left before right.
inline void tv_postorder(std::vector<NodeId>& order, const TreeView& t, NodeId root) {
  order.clear();
  std::vector<NodeId> stack; stack.push_back(root);
  while (!stack.empty()) {
    NodeId v = stack.back(); stack.pop_back();
    if (v < 0) { order.push_back(-v); continue; }
    stack.push_back(-v);
    if (!t.isleaf(v)) { stack.push_back(t.right[v]); stack.push_back(t.left[v]); }
  }
}

// _tree_preorder_ranges! (IND:65): preorder_index (entry time), subtree_last_index (max in subtree).
inline void tv_preorder_ranges(std::vector<int>& pre, std::vector<int>& last, const TreeView& t) {
  int n = t.nlocal();
  pre.assign(n + 1, 0); last.assign(n + 1, 0);
  int timer = 0;
  std::vector<NodeId> stack; stack.push_back(t.root);
  while (!stack.empty()) {
    NodeId v = stack.back(); stack.pop_back();
    if (v < 0) { last[-v] = timer; continue; }
    ++timer; pre[v] = timer;
    stack.push_back(-v);
    if (!t.isleaf(v)) { stack.push_back(t.right[v]); stack.push_back(t.left[v]); }
  }
}

// _depths (MG:48): per-local depth, root 0, left-first DFS.
inline std::vector<int> tv_depths(const TreeView& t) {
  int n = t.nlocal();
  std::vector<int> depth(n + 1, 0);
  std::vector<NodeId> stack; stack.push_back(t.root);
  while (!stack.empty()) {
    NodeId v = stack.back(); stack.pop_back();
    if (!t.isleaf(v)) {
      NodeId lc = t.left[v], rc = t.right[v];
      depth[lc] = depth[v] + 1; depth[rc] = depth[v] + 1;
      stack.push_back(rc); stack.push_back(lc);
    }
  }
  return depth;
}

// ============================================================================================
//  Centroid decomposition (CEN:43)
// ============================================================================================
inline CentroidInfo centroid_decomposition(const TreeView& t) {
  std::vector<int> counts = tv_leaf_count(t);
  int n = t.nlocal();
  CentroidInfo info;
  info.path_id.assign(n + 1, 0);
  info.path_pos.assign(n + 1, 0);
  info.path_head_of.assign(n + 1, 0);
  info.paths.clear();
  info.paths.push_back(CentroidPath{});   // paths[0] unused; 1-based path ids
  info.postorder_path_ids.clear();

  // iterative DFS over side-roots to mirror the recursion (decompose), preserving discovery order.
  std::vector<NodeId> work; work.push_back(t.root);
  // We must replicate the recursive structure exactly: assign path id in discovery order, then
  // recurse into side_roots[0..end-2] (in order), THEN push the path id to postorder.  A simple
  // explicit stack with a "needs postorder push" marker reproduces this.
  struct Frame { NodeId subroot; int pid; std::vector<NodeId> side_roots; size_t next_child; };
  std::vector<Frame> stack;

  auto build_path = [&](NodeId subroot) -> Frame {
    int pid = (int)info.paths.size();   // 1-based: paths currently has size = next id; size()==pid means new index pid
    // Actually info.paths size is pid_so_far+1 (index 0 dummy). New pid = current size (== count+1? no).
    // We compute pid = info.paths.size() (since index 0 is dummy, the first real path gets size()==1 => pid 1).
    pid = (int)info.paths.size();
    std::vector<NodeId> nodes, side_roots;
    NodeId current = subroot;
    while (true) {
      nodes.push_back(current);
      if (t.isleaf(current)) { side_roots.push_back(current); break; }
      NodeId lc = t.left[current], rc = t.right[current];
      int left_count = counts[lc], right_count = counts[rc];
      NodeId heavy, side;
      if (right_count > left_count) { heavy = rc; side = lc; }   // ties -> left heavy (strict >)
      else { heavy = lc; side = rc; }
      side_roots.push_back(side);
      current = heavy;
    }
    NodeId head = nodes[0];
    for (size_t pos = 0; pos < nodes.size(); ++pos) {
      NodeId node = nodes[pos];
      info.path_id[node] = pid;
      info.path_pos[node] = (int)pos + 1;   // 1-based
      info.path_head_of[node] = head;
    }
    CentroidPath cp; cp.head = head; cp.nodes = nodes; cp.side_roots = side_roots;
    info.paths.push_back(std::move(cp));
    Frame f; f.subroot = subroot; f.pid = pid; f.side_roots = side_roots; f.next_child = 0;
    return f;
  };

  stack.push_back(build_path(t.root));
  while (!stack.empty()) {
    Frame& f = stack.back();
    // recurse into side_roots[0 .. size-2] in order
    if (f.next_child + 1 < f.side_roots.size()) {
      NodeId sr = f.side_roots[f.next_child];
      ++f.next_child;
      stack.push_back(build_path(sr));
    } else {
      info.postorder_path_ids.push_back(f.pid);
      stack.pop_back();
    }
  }
  (void)work;
  return info;
}

inline const CentroidPath& root_centroid_path(const CentroidInfo& info) { return info.paths[1]; }  // path 1 (CEN:111)

// ============================================================================================
//  Side-tree view (CEN:127): postorder renumber a subtree of `t` rooted at side_root.
// ============================================================================================
inline TreeView side_tree_view(const TreeView& t, NodeId side_root) {
  std::vector<NodeId> post;
  tv_postorder(post, t, side_root);
  int nlocal = (int)post.size();
  std::vector<NodeId> mapping(t.parent.size(), 0);   // old local -> new local (1-based)
  for (int i = 0; i < nlocal; ++i) mapping[post[i]] = i + 1;

  TreeView v;
  v.tree_id = t.tree_id;
  v.orig.assign(nlocal + 1, 0);
  v.parent.assign(nlocal + 1, 0);
  v.left.assign(nlocal + 1, 0);
  v.right.assign(nlocal + 1, 0);
  v.taxon.assign(nlocal + 1, 0);
  for (NodeId old_id : post) {
    NodeId nid = mapping[old_id];
    v.orig[nid] = t.orig[old_id];
    v.taxon[nid] = t.taxon[old_id];
    if (!t.isleaf(old_id)) {
      NodeId nl = mapping[t.left[old_id]], nr = mapping[t.right[old_id]];
      v.left[nid] = nl; v.right[nid] = nr;
      v.parent[nl] = nid; v.parent[nr] = nid;
    }
  }
  v.root = mapping[side_root];
  return v;
}

// ============================================================================================
//  Induced subtree (IND:115): contract b onto `leaves` (preorder) + adjacent-pair LCAs.
// ============================================================================================
// Returns whether an induced view was built (false == empty leaves => nothing, ATT:1003).
inline bool induced_subtree_from_ordered_leaves(
    TreeView& out, const TreeView& base,
    const std::vector<NodeId>& leaves, const std::vector<NodeId>& adjacent_lcas,
    const std::vector<int>& pre, const std::vector<int>& last) {
  int nleaves = (int)leaves.size();
  // important = interleave [leaf1, lca(1,2), leaf2, lca(2,3), ...]  (IND:131-150)
  std::vector<NodeId> important;
  important.reserve(std::max(1, 2 * nleaves - 1));
  int previous_first = 0;
  for (int idx = 0; idx < nleaves; ++idx) {
    NodeId leaf = leaves[idx];
    int first = pre[leaf];
    // first > previous_first guaranteed by construction (strictly increasing preorder).
    (void)first; (void)previous_first; previous_first = first;
    important.push_back(leaf);
    if (idx > 0) important.push_back(adjacent_lcas[idx - 1]);
  }
  // sort by preorder index, then adjacent dedup (IND:152-160).
  std::sort(important.begin(), important.end(), [&](NodeId x, NodeId y){ return pre[x] < pre[y]; });
  {
    size_t write = 0;
    for (size_t read = 1; read < important.size(); ++read)
      if (important[read] != important[write]) important[++write] = important[read];
    important.resize(important.empty() ? 0 : write + 1);
  }
  int m = (int)important.size();

  // build virtual parent/child via a preorder ancestor stack (IND:162-198).
  std::vector<int> first_child(m, 0), second_child(m, 0);
  std::vector<uint8_t> child_count(m, 0);
  std::vector<int> stack;
  for (int posi = 0; posi < m; ++posi) {
    NodeId node = important[posi];
    while (!stack.empty()) {
      int ancestor_pos = stack.back();
      NodeId ancestor = important[ancestor_pos];
      if (pre[ancestor] <= pre[node] && pre[node] <= last[ancestor]) break;   // inclusive both ends
      stack.pop_back();
    }
    if (!stack.empty()) {
      int parent_pos = stack.back();
      uint8_t cnt = (uint8_t)(child_count[parent_pos] + 1);
      child_count[parent_pos] = cnt;
      if (cnt == 1) first_child[parent_pos] = posi;
      else if (cnt == 2) second_child[parent_pos] = posi;
      // (third child is a structural error; binary trees only — not expected here)
    }
    stack.push_back(posi);
  }

  // post-order emit (IND:200-222): root = position 0; push second then first => first popped first.
  std::vector<int> emit_order; emit_order.reserve(m);
  std::vector<int> emit_stack; emit_stack.push_back(0);   // position index
  while (!emit_stack.empty()) {
    int posi = emit_stack.back(); emit_stack.pop_back();
    if (posi < 0) { emit_order.push_back(-posi - 1); continue; }   // encode close as -(pos+1)
    emit_stack.push_back(-(posi + 1));
    uint8_t cnt = child_count[posi];
    if (cnt == 0) { /* leaf: base.isleaf(important[posi]) assumed */ }
    else { /* cnt==2 */ emit_stack.push_back(second_child[posi]); emit_stack.push_back(first_child[posi]); }
  }

  // assign local ids in emit (post) order (IND:224-252).
  std::vector<NodeId> local_for_pos(m, 0);
  out.tree_id = base.tree_id;
  out.orig.assign(m + 1, 0);
  out.parent.assign(m + 1, 0);
  out.left.assign(m + 1, 0);
  out.right.assign(m + 1, 0);
  out.taxon.assign(m + 1, 0);
  NodeId next_local = 0;
  for (int posi : emit_order) {
    ++next_local;
    NodeId node = important[posi];
    local_for_pos[posi] = next_local;
    out.orig[next_local] = base.orig[node];
    if (child_count[posi] == 0) {
      out.taxon[next_local] = base.taxon[node];
    } else {
      NodeId lc = local_for_pos[first_child[posi]], rc = local_for_pos[second_child[posi]];
      out.left[next_local] = lc; out.right[next_local] = rc;
      out.parent[lc] = next_local; out.parent[rc] = next_local;
    }
  }
  out.root = local_for_pos[0];
  return true;
}

// ============================================================================================
//  Candidate buckets (MG)
// ============================================================================================
struct Candidate { int right_pos; NodeId mapped; int depth; };
struct PathCandidateBucket { int path_id; std::vector<Candidate> candidates; };

// linear find-or-create by path_id (MG:157).
inline std::vector<Candidate>* candidate_bucket(std::vector<PathCandidateBucket>& buckets, int path_id, bool create) {
  for (auto& b : buckets) if (b.path_id == path_id) return &b.candidates;
  if (!create) return nullptr;
  buckets.push_back(PathCandidateBucket{path_id, {}});
  return &buckets.back().candidates;
}

// _emit_candidate! (MG:195): dedup by right_pos, keep shallowest depth (strict <, first-wins).
inline void emit_candidate(std::vector<PathCandidateBucket>& path_buckets, const CentroidInfo& cb,
                           NodeId right_local, NodeId mapped, int mapped_depth) {
  int path_id = cb.path_id[right_local];
  int right_pos = cb.path_pos[right_local];
  std::vector<Candidate>* bucket = candidate_bucket(path_buckets, path_id, true);
  for (auto& c : *bucket) {
    if (c.right_pos == right_pos) {
      if (mapped_depth < c.depth) c = Candidate{right_pos, mapped, mapped_depth};
      return;
    }
  }
  bucket->push_back(Candidate{right_pos, mapped, mapped_depth});
}

// _emit_ancestor_path_candidates! (MG:217): walk up b's centroid paths.
inline void emit_ancestor_path_candidates(std::vector<PathCandidateBucket>& path_buckets, const CentroidInfo& cb,
                                          const TreeView& b, NodeId current, int stop_path_id, NodeId mapped, int mapped_depth) {
  while (current != 0) {
    emit_candidate(path_buckets, cb, current, mapped, mapped_depth);
    if (stop_path_id != 0 && cb.path_id[current] == stop_path_id) break;
    NodeId head = cb.path_head_of[current];
    current = b.parent[head];
  }
}

// _orig_is_descendant (MG:67).
inline bool orig_is_descendant(const std::vector<NodeId>& parent, NodeId node, NodeId ancestor) {
  NodeId current = node;
  while (current != 0) { if (current == ancestor) return true; current = parent[current]; }
  return false;
}

// ============================================================================================
//  Recursive side info + per-mapped helpers (ATT:229-267)
// ============================================================================================

// forward declare result type used by side info
struct CentroidSolver::AttachedSubproblemResult {
  std::vector<CScoreStates> by_b_vertex;   // indexed by right-view LOCAL vertex id (1-based)
};

struct AttachedRecursiveSideInfo {
  int left_pos;
  TreeView induced_view;
  std::vector<int> induced_depth;
  CentroidSolver::AttachedSubproblemResult induced_result;
};

// _fill_mapped_local_vertices! (MG:76): for each induced-local vertex, walk its ORIGINAL-ancestor
// chain; record `mapped` for each ancestor present in b iff strictly shallower (first/lowest wins).
inline void fill_mapped_local_vertices(
    std::vector<NodeId>& mapped_by_b_vertex, std::vector<int>& mapped_depth_by_b_vertex,
    const TreeView& induced_view, const std::vector<int>& induced_depth,
    const TreeView& b, const std::vector<NodeId>& original_parent, const std::vector<NodeId>& b_orig_to_local) {
  std::fill(mapped_by_b_vertex.begin(), mapped_by_b_vertex.end(), 0);
  std::fill(mapped_depth_by_b_vertex.begin(), mapped_depth_by_b_vertex.end(), std::numeric_limits<int>::max());
  NodeId b_root_orig = b.orig[b.root];
  int nind = induced_view.nlocal();
  for (NodeId mapped = 1; mapped <= nind; ++mapped) {
    int depth = induced_depth[mapped];
    NodeId current_orig = induced_view.orig[mapped];
    while (current_orig != 0) {
      if (current_orig < (NodeId)b_orig_to_local.size()) {
        NodeId b_local = b_orig_to_local[current_orig];
        if (b_local != 0 && depth < mapped_depth_by_b_vertex[b_local]) {
          mapped_depth_by_b_vertex[b_local] = depth;
          mapped_by_b_vertex[b_local] = mapped;
        }
      }
      if (current_orig == b_root_orig) break;
      current_orig = original_parent[current_orig];
    }
  }
}

// ============================================================================================
//  Weighted search tree build (AM:109-182)
// ============================================================================================
inline int weighted_split_index(const std::vector<W>& prefix, int lo, int hi) {
  // unit weights; prefix is 1-based cumulative.  (AM:163)
  W base = (lo == 1) ? 0.0 : prefix[lo - 1];
  W total = prefix[hi] - base;
  W target = base + total / 2.0;
  // searchsortedfirst over prefix[lo..hi-1] for `target`, returns first index with prefix>=target.
  // std::lower_bound on the [lo, hi-1] subrange (1-based).
  int first_ge;
  {
    int L = lo, R = hi - 1;   // inclusive range of candidate indices
    int idx = hi;             // default: not found within [lo,hi-1] => hi (==(hi-1)+1)
    int loi = L, hii = R + 1; // binary search lower_bound
    while (loi < hii) {
      int mid = loi + (hii - loi) / 2;
      if (prefix[mid] < target) loi = mid + 1; else hii = mid;
    }
    idx = loi;   // first index in [lo, hi] with prefix>=target; if none in [lo,hi-1] gives hi
    first_ge = idx;
  }
  if (first_ge > hi - 1) return hi - 1;
  if (first_ge == lo) return lo;
  int prev = first_ge - 1;
  W prev_left = prefix[prev] - base;
  W curr_left = prefix[first_ge] - base;
  W prev_gap = std::abs(total - (prev_left + prev_left));
  W curr_gap = std::abs(total - (curr_left + curr_left));
  return curr_gap < prev_gap ? first_ge : prev;   // ties -> smaller index prev
}

inline WeightedSearchTree build_unit_weighted_search_tree(int nleaves) {
  WeightedSearchTree t;
  int nnodes = 2 * nleaves - 1;
  t.left.assign(nnodes + 1, 0);
  t.right.assign(nnodes + 1, 0);
  t.parent.assign(nnodes + 1, 0);
  t.leaf_for_right_pos.assign(nleaves + 1, 0);
  t.leaf_lo.assign(nnodes + 1, 0);
  t.leaf_hi.assign(nnodes + 1, 0);
  std::vector<W> prefix(nleaves + 1, 0.0);
  for (int i = 1; i <= nleaves; ++i) prefix[i] = prefix[i - 1] + 1.0;   // unit weights
  int next_id = 1;
  // recursive build (left first) — use explicit recursion via lambda.
  std::function<int(int,int)> build = [&](int lo, int hi) -> int {
    if (lo == hi) {
      int node = next_id++;
      t.leaf_for_right_pos[lo] = node;
      t.leaf_lo[node] = lo; t.leaf_hi[node] = hi;
      return node;
    }
    int split = weighted_split_index(prefix, lo, hi);
    int lc = build(lo, split);
    int rc = build(split + 1, hi);
    int node = next_id++;
    t.left[node] = lc; t.right[node] = rc;
    t.parent[lc] = node; t.parent[rc] = node;
    t.leaf_lo[node] = lo; t.leaf_hi[node] = hi;
    return node;
  };
  t.root = build(1, nleaves);
  return t;
}

// ============================================================================================
//  The attached matching DP (ATT:591) — the core matcher.
// ============================================================================================
// Bundles the per-(left,right path) inputs the edge emitter needs (ATT:361).
struct AttachedEdgeSource {
  CentroidSolver* cs;
  const TreeView* a;
  const TreeView* b;
  const CentroidInfo* cb;
  const CentroidPath* left_path;
  const std::vector<NodeId>* right_parent;   // original right parent
  const std::vector<TaxonId>* left_leaf_taxon;
  const std::vector<const AttachedRecursiveSideInfo*>* recursive_by_left;  // per left_pos (nullptr if none)
  std::vector<std::vector<PathCandidateBucket>>* candidates_by_left;       // per left_pos
  std::vector<std::vector<CScoreStates>*>* full_left_suffix;               // per right path id (nullptr if none)
  int path_id;
  const CentroidPath* right_path;
};

struct AttachedMultiEdge {
  int left_pos, right_pos;
  CScoreRef white_branch;
  CScoreStates red_single, green_single;
  CScoreRef red_cross, green_cross;
  W join_weight;
};

// _attached_edge_from_parts (ATT:389).
inline AttachedMultiEdge attached_edge_from_parts(int left_pos, int right_pos, const CScoreRef& white,
    const CScoreStates& red_base, const CScoreStates& green_base, W left_descent, W right_descent, W left_join, W right_join) {
  AttachedMultiEdge e;
  e.left_pos = left_pos; e.right_pos = right_pos;
  e.white_branch = white;
  e.red_single = lift_tree2(red_base, right_descent);
  e.green_single = lift_tree1(green_base, left_descent);
  e.red_cross = add_score(red_base[3], right_join);
  e.green_cross = add_score(green_base[3], left_join);
  e.join_weight = left_join + right_join;
  return e;
}

// _emit_attached_edges_for_left! (ATT:416).
template <class HandleEdge>
inline void emit_attached_edges_for_left(AttachedEdgeSource& src, int left_pos, HandleEdge&& handle_edge) {
  CentroidSolver& cs = *src.cs;
  const TreeView& a = *src.a;
  const TreeView& b = *src.b;
  const CentroidInfo& cb = *src.cb;
  const CentroidPath& left_path = *src.left_path;
  const CentroidPath& right_path = *src.right_path;

  std::vector<Candidate>* path_candidates = candidate_bucket((*src.candidates_by_left)[left_pos], src.path_id, false);
  if (path_candidates == nullptr) return;
  std::sort(path_candidates->begin(), path_candidates->end(), [](const Candidate& x, const Candidate& y){ return x.right_pos < y.right_pos; });

  NodeId left_node = left_path.nodes[left_pos - 1];
  NodeId left_side = left_path.side_roots[left_pos - 1];
  W left_descent = cs.attached_descent(a, left_node, left_side);
  W left_join = attached_join_weight(cs, a, left_node);
  bool is_boundary_left = (left_pos == (int)left_path.nodes.size());
  TaxonId left_taxon = (*src.left_leaf_taxon)[left_pos];
  bool is_singleton_left = (left_taxon != 0);

  CScoreRef inv_ref = invalid_ref();
  CScoreStates inv_states = invalid_states();

  for (const Candidate& cand : *path_candidates) {
    int right_pos = cand.right_pos;
    NodeId mapped = cand.mapped;
    NodeId right_node = right_path.nodes[right_pos - 1];
    NodeId right_side = right_path.side_roots[right_pos - 1];
    W right_descent = cs.attached_descent(b, right_node, right_side);
    W right_join = attached_join_weight(cs, b, right_node);
    bool is_boundary_right = (right_pos == (int)right_path.nodes.size());

    CScoreRef white = inv_ref;
    CScoreStates red_base = inv_states;
    CScoreStates green_base = inv_states;

    if (is_boundary_left || is_singleton_left) {
      // leaf-left (ATT:456). leaf_score == 1.0 (alpha).
      CScoreStates leaf_states = attached_leaf_states(cs, 1.0, left_taxon);
      // side_states = _lift_between_tree2(b, leaf_states, b.orig[right_side], left_taxon)
      // here desc_orig == NodeId(left_taxon) which is the original tree-2 leaf id for that taxon.
      {
        const PhyloTree& t = cs.tree_for(b); const std::vector<W>& p = cs.prefix_for(b);
        W d = attached_descent_weight(t, p, b.orig[right_side], (NodeId)left_taxon);
        CScoreStates side_states = lift_tree2(leaf_states, d);
        white = side_states[3];
        if (is_boundary_left) {
          red_base = side_states;
        } else if (is_boundary_right) {
          red_base = lift_tree1(side_states, left_descent);
        } else {
          int child_path_id = cb.path_id[right_side];
          std::vector<CScoreStates>* child_left = (*src.full_left_suffix)[child_path_id];
          red_base = (child_left == nullptr) ? inv_states : (*child_left)[left_pos];
        }
        green_base = lift_tree2(side_states, right_descent);
      }
    } else {
      const AttachedRecursiveSideInfo* side_info = (*src.recursive_by_left)[left_pos];
      if (side_info == nullptr) continue;
      NodeId right_orig = b.orig[right_node];
      NodeId side_root_orig = b.orig[right_side];
      // green_base = _attached_states_from_mapped(side_info, b, right_orig, mapped)
      //   = lift_between_tree2(b, side_info.by_b_vertex[mapped].nonempty, right_orig, induced_view.orig[mapped])
      {
        const TreeView& iv = side_info->induced_view;
        const CScoreStates& states = side_info->induced_result.by_b_vertex[mapped];
        const PhyloTree& t = cs.tree_for(b); const std::vector<W>& p = cs.prefix_for(b);
        W d = attached_descent_weight(t, p, right_orig, iv.orig[mapped]);
        green_base = lift_tree2(states, d);
      }
      if (is_boundary_right) {
        white = green_base[3];
        red_base = lift_tree1(green_base, left_descent);
      } else {
        // recursive_white = _attached_white_from_mapped(...) (ATT:241)
        const TreeView& iv = side_info->induced_view;
        NodeId target = 0;
        if (iv.orig[mapped] != right_orig) {
          target = mapped;
        } else {
          for (NodeId child : {iv.left[mapped], iv.right[mapped]}) {
            if (child == 0) continue;
            if (orig_is_descendant(*src.right_parent, iv.orig[child], side_root_orig)) { target = child; break; }
          }
        }
        if (target == 0) {
          white = inv_ref;
        } else {
          const CScoreStates& states = side_info->induced_result.by_b_vertex[target];
          const PhyloTree& t = cs.tree_for(b); const std::vector<W>& p = cs.prefix_for(b);
          W d = attached_descent_weight(t, p, side_root_orig, iv.orig[target]);
          CScoreStates lifted = lift_tree2(states, d);
          white = lifted[3];
        }
        int child_path_id = cb.path_id[right_side];
        std::vector<CScoreStates>* child_left = (*src.full_left_suffix)[child_path_id];
        red_base = (child_left == nullptr) ? inv_states : (*child_left)[left_pos];
      }
    }

    handle_edge(left_pos, attached_edge_from_parts(left_pos, right_pos, white, red_base, green_base,
                                                   left_descent, right_descent, left_join, right_join));
  }
}

// Result of the matching DP (ATT:471): per-left_pos best ScoreStates over the right path's subtree.
struct AttachedGraphMatchingResult {
  std::vector<CScoreStates> best_from_left_suffix;   // per left_pos (1-based index, [0] unused)
};

// _dp_attached_agreement_matching_from_edges (ATT:591).
// right_suffix_cb(right_node, states); exact_candidate_cb(left_node, right_node, ref).
template <class RightSuffixCb, class ExactCb>
inline AttachedGraphMatchingResult dp_attached_agreement_matching_from_edges(
    CentroidSolver& cs, const std::vector<NodeId>& left_nodes, const std::vector<NodeId>& right_nodes,
    const std::vector<W>& left_bases, const std::vector<W>& right_bases, AttachedEdgeSource& edge_source,
    bool have_right_suffix_cb, RightSuffixCb&& right_suffix_cb,
    std::vector<CScoreStates>& best_from_left_buffer, ExactCb&& exact_cb) {
  CScoreRef inv_ref = invalid_ref();
  CScoreStates inv_states = invalid_states();
  int nright = (int)right_nodes.size();
  int nleft = (int)left_nodes.size();

  // segment tree (cached by nright).
  auto it = cs.tree_cache.find(nright);
  if (it == cs.tree_cache.end()) it = cs.tree_cache.emplace(nright, build_unit_weighted_search_tree(nright)).first;
  const WeightedSearchTree& tree = it->second;
  int nnodes = (int)tree.parent.size() - 1;

  // best_from_left buffer (1-based; index 0 unused), filled invalid.
  best_from_left_buffer.assign(nleft + 1, inv_states);

  // DP arrays (1-based by tree node).
  std::vector<CScoreRef> red_self(nnodes + 1, inv_ref), red_max(nnodes + 1, inv_ref);
  std::vector<CScoreStates> proper_self(nnodes + 1, inv_states), proper_max(nnodes + 1, inv_states);
  std::vector<CScoreStates> white_self(nnodes + 1, inv_states), white_max(nnodes + 1, inv_states);
  std::vector<CScoreStates> best_max(nnodes + 1, inv_states);

  auto pull = [&](int node) {
    int lc = tree.left[node], rc = tree.right[node];
    CScoreRef red = red_self[node];
    CScoreStates proper = proper_self[node];
    CScoreStates white = white_self[node];
    if (lc != 0) {
      red = cs.choose_ref(red_max[lc], red);
      proper = cs.choose_states(proper_max[lc], proper);
      white = cs.choose_states(white_max[lc], white);
    }
    if (rc != 0) {
      red = cs.choose_ref(red_max[rc], red);
      proper = cs.choose_states(proper_max[rc], proper);
      white = cs.choose_states(white_max[rc], white);
    }
    red_max[node] = red; proper_max[node] = proper; white_max[node] = white;
    best_max[node] = cs.choose_states(proper, white);   // proper wins ties
  };
  auto recompute_to_root = [&](int node) {
    int current = node;
    while (current != 0) { pull(current); current = tree.parent[current]; }
  };
  auto update_leaf_ref = [&](std::vector<CScoreRef>& dest, int right_pos, const CScoreRef& cand) {
    if (!is_valid_score(cand.score)) return;
    int leaf = tree.leaf_for_right_pos[right_pos];
    if (cs.better_ref(cand, dest[leaf])) { dest[leaf] = cand; recompute_to_root(leaf); }
  };
  auto states_eq = [](const CScoreStates& x, const CScoreStates& y) {
    for (int m = 0; m < 4; ++m) if (x.s[m].score != y.s[m].score || x.s[m].taxa != y.s[m].taxa) return false;
    return true;
  };
  auto update_leaf_states = [&](std::vector<CScoreStates>& dest, int right_pos, const CScoreStates& cand) {
    if (!has_valid_state(cand)) return;
    int leaf = tree.leaf_for_right_pos[right_pos];
    CScoreStates updated = cs.choose_states(cand, dest[leaf]);
    if (!states_eq(updated, dest[leaf])) { dest[leaf] = updated; recompute_to_root(leaf); }
  };

  std::vector<int> range_stack;
  auto query_range = [&](int lo, int hi, W left_base, W right_base) -> CScoreStates {
    if (lo > hi) return inv_states;
    CScoreStates best = inv_states;
    range_stack.clear();
    range_stack.push_back(tree.root);
    while (!range_stack.empty()) {
      int node = range_stack.back(); range_stack.pop_back();
      if (tree.leaf_hi[node] < lo) continue;
      if (hi < tree.leaf_lo[node]) continue;
      if (lo <= tree.leaf_lo[node] && tree.leaf_hi[node] <= hi) {
        best = cs.choose_states(best_max[node], best);
        continue;
      }
      int lc = tree.left[node], rc = tree.right[node];
      if (lc != 0) range_stack.push_back(lc);
      if (rc != 0) range_stack.push_back(rc);
    }
    return adjust_states(best, -left_base, -right_base);
  };

  auto update_pair_range = [&](int lo, int hi, const CScoreRef& green, int left_pos) {
    if (lo > hi) return;
    if (!is_valid_score(green.score)) return;
    if (!(left_pos < (int)left_bases.size())) return;
    W red_left_base = left_bases[left_pos];
    range_stack.clear();
    range_stack.push_back(tree.root);
    while (!range_stack.empty()) {
      int node = range_stack.back(); range_stack.pop_back();
      if (tree.leaf_hi[node] < lo) continue;
      if (hi < tree.leaf_lo[node]) continue;
      if (!is_valid_score(red_max[node].score)) continue;
      if (tree.is_leaf(node)) {
        CScoreRef red = red_self[node];
        if (!is_valid_score(red.score)) continue;
        int red_right_pos = tree.leaf_lo[node];
        W green_score = green.score - right_bases[red_right_pos];   // right_bases[rp+1] 1-based == index rp here
        W red_score = red.score - red_left_base;
        W score = green_score + red_score;
        CScoreRef candidate = cs.score_ref_union2(score, red, green);
        exact_cb(left_nodes[left_pos - 1], right_nodes[red_right_pos - 1], candidate);
        CScoreStates adjusted = adjust_states(same_states(candidate), left_bases[left_pos - 1], right_bases[red_right_pos - 1]);
        update_leaf_states(proper_self, red_right_pos, adjusted);
        continue;
      }
      int lc = tree.left[node], rc = tree.right[node];
      if (lc != 0) range_stack.push_back(lc);
      if (rc != 0) range_stack.push_back(rc);
    }
  };

  // pending vectors (per left iteration).
  std::vector<std::pair<int,CScoreStates>> pending_proper, pending_white;
  std::vector<std::pair<int,CScoreRef>> pending_red, pending_pairs;

  auto handle_edge = [&](int left_pos, const AttachedMultiEdge& edge) {
    int right_pos = edge.right_pos;
    W left_base = left_bases[left_pos - 1];
    W right_base = right_bases[right_pos - 1];
    if (has_valid_state(edge.red_single))
      pending_proper.emplace_back(right_pos, adjust_states(edge.red_single, left_base, right_base));
    if (has_valid_state(edge.green_single))
      pending_proper.emplace_back(right_pos, adjust_states(edge.green_single, left_base, right_base));
    if (is_valid_score(edge.red_cross.score))
      pending_red.emplace_back(right_pos, adjust_ref(edge.red_cross, left_base));
    if (is_valid_score(edge.green_cross.score) && right_pos > 1)
      pending_pairs.emplace_back(right_pos - 1, adjust_ref(edge.green_cross, right_base));
    if (is_valid_score(edge.white_branch.score) && right_pos < nright && left_pos < nleft) {
      // lower = query_range(right_pos+1, nright, left_bases[left_pos+1], right_bases[right_pos+1])[3]
      CScoreStates lower_states = query_range(right_pos + 1, nright, left_bases[left_pos], right_bases[right_pos]);
      CScoreRef lower = lower_states[3];
      if (is_valid_score(lower.score)) {
        W score = edge.white_branch.score + lower.score + edge.join_weight;
        CScoreRef candidate = cs.score_ref_union2(score, edge.white_branch, lower);
        exact_cb(left_nodes[left_pos - 1], right_nodes[right_pos - 1], candidate);
        pending_white.emplace_back(right_pos, adjust_states(same_states(candidate), left_base, right_base));
      }
    }
  };

  // Main left loop: LAST -> FIRST (ATT:792).
  for (int left_pos = nleft; left_pos >= 1; --left_pos) {
    pending_proper.clear(); pending_white.clear(); pending_red.clear(); pending_pairs.clear();
    emit_attached_edges_for_left(edge_source, left_pos, handle_edge);
    // commit pairs -> proper -> white -> red (ATT:798).
    for (auto& pr : pending_pairs) update_pair_range(1, pr.first, pr.second, left_pos);
    for (auto& pr : pending_proper) update_leaf_states(proper_self, pr.first, pr.second);
    for (auto& pr : pending_white) update_leaf_states(white_self, pr.first, pr.second);
    for (auto& pr : pending_red) update_leaf_ref(red_self, pr.first, pr.second);
    best_from_left_buffer[left_pos] = adjust_states(best_max[tree.root], -left_bases[left_pos - 1], -right_bases[0]);
  }

  // Right-suffix pass (ATT:813).
  for (int right_pos = 1; right_pos <= nright; ++right_pos) {
    CScoreStates suffix = query_range(right_pos, nright, left_bases[0], right_bases[right_pos - 1]);
    if (have_right_suffix_cb) right_suffix_cb(right_nodes[right_pos - 1], suffix);
  }

  AttachedGraphMatchingResult res;
  res.best_from_left_suffix = std::move(best_from_left_buffer);
  return res;
}

// ============================================================================================
//  Matching phase _build_attached_matching_graphs_and_results! (ATT:825)
// ============================================================================================
inline void build_attached_matching_graphs_and_results(
    CentroidSolver& cs, const TreeView& a, const TreeView& b, const CentroidInfo& ca, const CentroidInfo& cb,
    const std::vector<AttachedRecursiveSideInfo>& recursive_results, const std::vector<NodeId>& b_orig_to_local,
    std::vector<CScoreStates>& by_b_vertex) {
  const CentroidPath& left_path = root_centroid_path(ca);
  const PhyloTree& right_tree = cs.tree_for(b);
  const std::vector<NodeId>& right_parent = right_tree.parent;
  int npaths = (int)cb.paths.size() - 1;   // 1-based; real paths 1..npaths

  // _mark_b_leaves_by_taxon! (MG:267).
  int leaf_epoch = cs.next_leaf_for_taxon_epoch();
  for (NodeId v = 1; v <= b.nlocal(); ++v) {
    if (b.isleaf(v)) {
      TaxonId taxon = b.taxon[v];
      cs.leaf_for_taxon[taxon] = v;
      cs.leaf_for_taxon_epoch[taxon] = leaf_epoch;
    }
  }

  // _left_leaf_taxa (MG:281): per left_pos, taxon of side_root if it's a leaf (1-based index).
  int nleft = (int)left_path.nodes.size();
  std::vector<TaxonId> left_leaf_taxon(nleft + 1, 0);
  for (int i = 1; i <= nleft; ++i) {
    NodeId side_root = left_path.side_roots[i - 1];
    if (a.isleaf(side_root)) left_leaf_taxon[i] = a.taxon[side_root];
  }

  // recursive_by_left (ATT:211): scatter recursive_results by left_pos.
  std::vector<const AttachedRecursiveSideInfo*> recursive_by_left(nleft + 1, nullptr);
  for (const auto& item : recursive_results) recursive_by_left[item.left_pos] = &item;

  // candidate buckets per left_pos (MG:174 + _populate_candidate_buckets! MG:305).
  std::vector<std::vector<PathCandidateBucket>> candidates_by_left(nleft + 1);
  for (int left_pos = 1; left_pos <= nleft; ++left_pos) {
    std::vector<PathCandidateBucket>& path_buckets = candidates_by_left[left_pos];
    if (left_pos == nleft || left_leaf_taxon[left_pos] != 0) {
      TaxonId taxon = left_leaf_taxon[left_pos];
      NodeId twin_local = (taxon == 0 || cs.leaf_for_taxon_epoch[taxon] != leaf_epoch) ? 0 : cs.leaf_for_taxon[taxon];
      if (twin_local != 0)
        emit_ancestor_path_candidates(path_buckets, cb, b, twin_local, 0, 1, 0);
    } else {
      const AttachedRecursiveSideInfo* side_info = recursive_by_left[left_pos];
      if (side_info != nullptr) {
        // _emit_induced_candidates! (MG:235).
        const TreeView& view = side_info->induced_view;
        for (NodeId mapped = 1; mapped <= view.nlocal(); ++mapped) {
          NodeId parent_local = view.parent[mapped];
          int stop_path = (parent_local == 0) ? 0 : cb.path_id[b_orig_to_local[view.orig[parent_local]]];
          NodeId current = b_orig_to_local[view.orig[mapped]];
          int mapped_depth = side_info->induced_depth[mapped];
          emit_ancestor_path_candidates(path_buckets, cb, b, current, stop_path, mapped, mapped_depth);
        }
      }
    }
  }

  // full_left_suffix per right path id (nullptr until filled). We own the storage via a pool.
  std::vector<std::vector<CScoreStates>*> full_left_suffix(npaths + 1, nullptr);
  std::vector<std::unique_ptr<std::vector<CScoreStates>>> suffix_storage;  // keeps vectors alive
  std::vector<std::unique_ptr<std::vector<CScoreStates>>> suffix_pool;     // reusable buffers

  std::vector<W> left_bases = path_prefix_bases(cs, a, left_path.nodes);

  // dummy callbacks for the rooted writers (exact_cb) and update_vertex (right_suffix_cb).
  auto update_vertex = [&](NodeId v, const CScoreStates& candidate) {
    by_b_vertex[v] = cs.choose_states(candidate, by_b_vertex[v]);
  };
  auto right_suffix_cb = [&](NodeId right_node, const CScoreStates& suffix) {
    update_vertex(right_node, suffix);
  };
  auto exact_cb = [&](NodeId left_root, NodeId right_root, const CScoreRef& candidate) {
    cs.update_rooted_tree1(a.orig[left_root], candidate);
    cs.update_rooted_tree2(b.orig[right_root], candidate);
  };

  for (int path_id : cb.postorder_path_ids) {
    const CentroidPath& right_path = cb.paths[path_id];
    std::vector<W> right_bases = path_prefix_bases(cs, b, right_path.nodes);

    // acquire a left_suffix buffer (pooled).
    std::unique_ptr<std::vector<CScoreStates>> left_suffix_buffer;
    if (!suffix_pool.empty()) { left_suffix_buffer = std::move(suffix_pool.back()); suffix_pool.pop_back(); }
    else left_suffix_buffer = std::make_unique<std::vector<CScoreStates>>();

    AttachedEdgeSource edge_source;
    edge_source.cs = &cs; edge_source.a = &a; edge_source.b = &b; edge_source.cb = &cb;
    edge_source.left_path = &left_path; edge_source.right_parent = &right_parent;
    edge_source.left_leaf_taxon = &left_leaf_taxon; edge_source.recursive_by_left = &recursive_by_left;
    edge_source.candidates_by_left = &candidates_by_left; edge_source.full_left_suffix = &full_left_suffix;
    edge_source.path_id = path_id; edge_source.right_path = &right_path;

    AttachedGraphMatchingResult graph_result = dp_attached_agreement_matching_from_edges(
        cs, left_path.nodes, right_path.nodes, left_bases, right_bases, edge_source,
        true, right_suffix_cb, *left_suffix_buffer, exact_cb);

    // move result into storage so full_left_suffix can point at it.
    suffix_storage.push_back(std::make_unique<std::vector<CScoreStates>>(std::move(graph_result.best_from_left_suffix)));
    std::vector<CScoreStates>* combined_ptr = suffix_storage.back().get();
    // the buffer we passed was consumed (moved-from); recycle the now-stale unique_ptr shell.
    suffix_pool.push_back(std::move(left_suffix_buffer));

    if ((int)right_path.side_roots.size() == 1) {
      full_left_suffix[path_id] = combined_ptr;
    } else {
      std::vector<CScoreStates>& combined_left = *combined_ptr;   // 1-based, index 0 unused
      for (size_t idx = 0; idx + 1 < right_path.side_roots.size(); ++idx) {
        NodeId side_root = right_path.side_roots[idx];
        int child_path_id = cb.path_id[side_root];
        std::vector<CScoreStates>* child_left = full_left_suffix[child_path_id];
        if (child_left == nullptr) continue;
        W child_delta = cs.attached_descent(b, right_path.head, side_root);
        for (size_t left_pos = 1; left_pos < combined_left.size(); ++left_pos) {
          CScoreStates lifted_child = lift_tree2((*child_left)[left_pos], child_delta);
          combined_left[left_pos] = cs.choose_states(lifted_child, combined_left[left_pos]);
        }
        full_left_suffix[child_path_id] = nullptr;
        // child_left buffer freed when its suffix_storage entry is later destroyed; we don't pool it
        // (it lives in suffix_storage). Harmless to keep.
      }
      full_left_suffix[path_id] = combined_ptr;
    }
  }
}

// ============================================================================================
//  The engine impl _solve_attached_agree_impl! (ATT:936)
// ============================================================================================
inline CentroidSolver::AttachedSubproblemResult
CentroidSolver::solve_attached_agree(const TreeView& a, const TreeView& b) {
  CentroidSolver& cs = *this;
  CScoreStates inv_states = invalid_states();
  AttachedSubproblemResult result;
  result.by_b_vertex.assign(b.nlocal() + 1, inv_states);   // 1-based; [0] unused

  CentroidInfo ca = centroid_decomposition(a);
  CentroidInfo cb = centroid_decomposition(b);
  const PhyloTree& right_tree = cs.tree_for(b);
  const LCAIndex& right_lca = cs.lca_for(b);

  // b_orig_to_local (SOL:868): indexed by ORIGINAL id, size = max(b.orig)+1.
  NodeId max_orig = 0;
  for (NodeId v = 1; v <= b.nlocal(); ++v) max_orig = std::max(max_orig, b.orig[v]);
  std::vector<NodeId> b_orig_to_local(max_orig + 1, 0);
  for (NodeId v = 1; v <= b.nlocal(); ++v) b_orig_to_local[b.orig[v]] = v;

  std::vector<int> b_pre, b_last;
  tv_preorder_ranges(b_pre, b_last, b);

  const CentroidPath& left_path = root_centroid_path(ca);

  // ordered_b_leaves: b's leaves left-to-right (stack push right then left) (ATT:979).
  std::vector<NodeId> ordered_b_leaves;
  {
    std::vector<NodeId> stack; stack.push_back(b.root);
    while (!stack.empty()) {
      NodeId v = stack.back(); stack.pop_back();
      if (b.isleaf(v)) ordered_b_leaves.push_back(v);
      else { stack.push_back(b.right[v]); stack.push_back(b.left[v]); }
    }
  }

  // induced_view_for_side_view (ATT:992).
  auto induced_view_for_side_view = [&](const TreeView& side_view, TreeView& out) -> bool {
    int epoch = cs.next_taxon_mark_epoch();
    for (NodeId v = 1; v <= side_view.nlocal(); ++v)
      if (side_view.isleaf(v)) cs.taxon_marks[side_view.taxon[v]] = epoch;
    std::vector<NodeId> leaves;
    for (NodeId leaf : ordered_b_leaves)
      if (cs.taxon_marks[b.taxon[leaf]] == epoch) leaves.push_back(leaf);
    if (leaves.empty()) return false;
    std::vector<NodeId> adjacent_lcas(std::max(0, (int)leaves.size() - 1));
    for (size_t i = 1; i < leaves.size(); ++i) {
      NodeId original_lca = lca_query(right_tree, right_lca, b.orig[leaves[i - 1]], b.orig[leaves[i]]);
      NodeId local_lca = (original_lca < (NodeId)b_orig_to_local.size()) ? b_orig_to_local[original_lca] : 0;
      // require local_lca != 0 (invariant); if it ever fails the view is malformed.
      adjacent_lcas[i - 1] = local_lca;
    }
    return induced_subtree_from_ordered_leaves(out, b, leaves, adjacent_lcas, b_pre, b_last);
  };

  // recursion over left side trees (ATT:1032).
  std::vector<AttachedRecursiveSideInfo> recursive_results;
  std::vector<NodeId> mapped_by_b_vertex(b.nlocal() + 1, 0);
  std::vector<int> mapped_depth_by_b_vertex(b.nlocal() + 1, std::numeric_limits<int>::max());

  for (int left_pos = 1; left_pos + 1 <= (int)left_path.side_roots.size(); ++left_pos) {
    NodeId side_root = left_path.side_roots[left_pos - 1];
    if (a.isleaf(side_root)) continue;

    TreeView side_view = side_tree_view(a, side_root);
    TreeView induced_view;
    if (!induced_view_for_side_view(side_view, induced_view)) continue;

    AttachedSubproblemResult induced_result = cs.solve_attached_agree(side_view, induced_view);
    std::vector<int> induced_depth = tv_depths(induced_view);
    fill_mapped_local_vertices(mapped_by_b_vertex, mapped_depth_by_b_vertex, induced_view, induced_depth,
                               b, right_tree.parent, b_orig_to_local);

    // We need stable storage for induced_view/result so AttachedRecursiveSideInfo pointers used later
    // (recursive_by_left) remain valid; push first, then fold using the stored copy.
    recursive_results.push_back(AttachedRecursiveSideInfo{left_pos, std::move(induced_view), std::move(induced_depth), std::move(induced_result)});
    AttachedRecursiveSideInfo& stored = recursive_results.back();

    for (NodeId v = 1; v <= b.nlocal(); ++v) {
      NodeId mapped = mapped_by_b_vertex[v];
      if (mapped == 0) continue;
      CScoreStates states = stored.induced_result.by_b_vertex[mapped];
      states = lift_between_tree1(cs, a, states, a.root, side_root);   // a-root -> side-root (uses orig internally)
      // _lift_between_tree2(b, states, b.orig[v], induced_view.orig[mapped])
      {
        const PhyloTree& t = cs.tree_for(b); const std::vector<W>& p = cs.prefix_for(b);
        W d = attached_descent_weight(t, p, b.orig[v], stored.induced_view.orig[mapped]);
        states = lift_tree2(states, d);
      }
      result.by_b_vertex[v] = cs.choose_states(states, result.by_b_vertex[v]);
    }
  }

  build_attached_matching_graphs_and_results(cs, a, b, ca, cb, recursive_results, b_orig_to_local, result.by_b_vertex);

  return result;
}

// ============================================================================================
//  root_view (TR:865) for an original PhyloTree.
// ============================================================================================
inline TreeView root_view(const PhyloTree& t, uint8_t tree_id) {
  int nv = t.nv();
  TreeView v;
  v.tree_id = tree_id;
  v.root = t.root;
  v.orig.assign(nv + 1, 0);
  for (NodeId i = 1; i <= nv; ++i) v.orig[i] = i;   // identity (local == original)
  v.parent = t.parent; v.left = t.left; v.right = t.right;
  v.taxon.assign(nv + 1, 0);
  for (NodeId i = 1; i <= t.ntaxa; ++i) v.taxon[i] = i;
  return v;
}

// ============================================================================================
//  Top driver + emission (SOL:736, M:267).
// ============================================================================================
inline std::vector<std::pair<W, std::vector<int>>>
centroid_mast_solve(const PhyloTree& t1, const PhyloTree& t2,
                    const std::vector<W>* leaf_w = nullptr) {
  CentroidSolver cs(t1, t2);

  // Seed leaf accumulators (SOL:698): each leaf slot v=taxon in both trees gets the singleton.  The
  // seed IS the leaf weight (alpha): 1.0 uniform => size-MAST; a per-leaf leaf_w[taxon-1] lets callers
  // bias the agreement subtree toward/away from individual taxa (composes with the internal weights).
  for (TaxonId taxon = 1; taxon <= t1.ntaxa; ++taxon) {
    W score = leaf_w ? (*leaf_w)[taxon - 1] : (W)1.0;   // alpha (per-leaf if supplied)
    CScoreRef ref = cs.score_ref_from_taxon(score, taxon);
    cs.update_rooted_tree1((NodeId)taxon, ref);
    cs.update_rooted_tree2((NodeId)taxon, ref);
  }

  // Run the engine on the full root views.
  TreeView va = root_view(t1, 1);
  TreeView vb = root_view(t2, 2);
  cs.solve_attached_agree(va, vb);

  // Emission: vcat(tree1, tree2) ascending vertex id, filter |taxa|>1 (M:275, §4).
  std::vector<std::pair<W, std::vector<int>>> out;
  auto emit_tree = [&](const std::vector<CScoreRef>& slots) {
    for (size_t root = 1; root < slots.size(); ++root) {
      const CScoreRef& cand = slots[root];
      if (!is_valid_score(cand.score)) continue;
      const std::vector<TaxonId>& taxa = cs.taxa_view(cand);
      if ((int)taxa.size() > 1) out.emplace_back(cand.score, taxa);   // taxa already sorted
    }
  };
  emit_tree(cs.rooted_tree1);
  emit_tree(cs.rooted_tree2);
  return out;
}

} // namespace maf::mast
