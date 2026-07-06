// Multi-tree rooted MAST DP production implementation.
// The classic O(n^3) two-leaf-anchored agreement DP generalized to t trees, producing for every
// vertex of every tree the best non-blacklisted agreement subtree rooted there.
//
// Indexing: trees are 0-based (ti in 0..t-1); vertices/taxa are 1-based (PhyloTree arrays); the
// ordered-pair index is 1-based, _pair_index(n,i,j)=i+(j-1)*n in 1..n*n (arrays sized n*n+1, index
// 0 unused). Ranks are 1-based (converted to 0-based at array access).
#pragma once
#include "trees.hpp"
#include "lca.hpp"
#include "ranked.hpp"
#include <vector>
#include <optional>
#include <algorithm>

namespace maf::mast {

struct MultiMASTInstance {
  std::vector<PhyloTree> trees;     // 0-based
  std::vector<W> taxon_weight;      // 1-based, size ntaxa+1
};

inline int _pair_index(int n, int i, int j) { return i + (j - 1) * n; }   // 1-based, MD:57

// _multi_preorder_ranges (MD:87): (entry,exit) preorder timestamps; push right then left.
inline void _preorder_ranges(const PhyloTree& t, std::vector<int>& first, std::vector<int>& last) {
  const int nv = t.nv();
  first.assign(nv + 1, 0);
  last.assign(nv + 1, 0);
  std::vector<NodeId> stack{t.root};
  int timer = 0;
  while (!stack.empty()) {
    NodeId v = stack.back(); stack.pop_back();
    if (v < 0) { last[-v] = timer; continue; }
    ++timer; first[v] = timer;
    stack.push_back(-v);
    if (!t.isleaf(v)) { stack.push_back(t.right[v]); stack.push_back(t.left[v]); }
  }
}
// _multi_child_below_ancestor (MD:111): climb from leaf up to the child of `ancestor`.
inline NodeId _child_below_ancestor(const PhyloTree& t, NodeId ancestor, TaxonId leaf) {
  NodeId child = leaf, parent = t.parent[child];
  while (parent != ancestor) { child = parent; parent = t.parent[child]; }
  return child;
}
inline bool _is_descendant(const std::vector<int>& pre, const std::vector<int>& last, NodeId node, NodeId anc) {
  return pre[anc] <= pre[node] && pre[node] <= last[anc];   // MD:121
}

struct MultiDPData {
  const MultiMASTInstance* inst;
  int ntaxa;
  bool has_internal_weights;
  std::vector<LCAIndex> lca_indices;                  // [ti]
  std::vector<std::vector<W>> internal_prefixes;      // [ti] 1-based (empty if no iw)
  std::vector<std::vector<int>> leaf_counts;          // [ti] 1-based
  std::vector<std::vector<int>> preorder_index;       // [ti] 1-based
  std::vector<std::vector<int>> subtree_last_index;   // [ti] 1-based
  std::vector<std::vector<int>> pair_roots;           // [ti][pair] pair 1-based
  std::vector<std::vector<int>> pair_anchor_side_roots; // [ti][pair]
  std::vector<int> pair_measure;                      // [pair]
  std::vector<W> pair_root_weight;                    // [pair] or empty
  std::vector<W> singleton_attach;                    // [pair] or empty
};

inline bool _multi_has_nonzero_iw(const std::vector<PhyloTree>& trees) {   // MD:66
  for (const PhyloTree& t : trees)
    for (NodeId v = t.ntaxa + 1; v <= t.nv(); ++v)
      if (t.internal_weight[v] != 0.0) return true;
  return false;
}

// _multi_dp_data (MD:131): the immutable precompute.
inline MultiDPData build_multi_dp_data(const MultiMASTInstance& inst) {
  MultiDPData d;
  d.inst = &inst;
  const int n = (int)inst.taxon_weight.size() - 1;   // taxon_weight 1-based
  d.ntaxa = n;
  const int total = n * n;
  const int T = (int)inst.trees.size();
  d.has_internal_weights = _multi_has_nonzero_iw(inst.trees);

  d.lca_indices.resize(T);
  d.internal_prefixes.assign(T, {});
  d.leaf_counts.resize(T);
  d.preorder_index.resize(T);
  d.subtree_last_index.resize(T);
  d.pair_roots.assign(T, std::vector<int>(total + 1, 0));
  d.pair_anchor_side_roots.assign(T, std::vector<int>(total + 1, 0));
  d.pair_measure.assign(total + 1, 0);
  d.pair_root_weight = d.has_internal_weights ? std::vector<W>(total + 1, 0.0) : std::vector<W>{};
  d.singleton_attach = d.has_internal_weights ? std::vector<W>(total + 1, 0.0) : std::vector<W>{};

  for (int ti = 0; ti < T; ++ti) {
    const PhyloTree& tree = inst.trees[ti];
    d.lca_indices[ti] = build_lca(tree);
    if (d.has_internal_weights) d.internal_prefixes[ti] = internal_prefix(tree);
    d.leaf_counts[ti] = leaf_count(tree);
    _preorder_ranges(tree, d.preorder_index[ti], d.subtree_last_index[ti]);
  }

  for (int i = 1; i <= n - 1; ++i)
    for (int j = i + 1; j <= n; ++j) {
      int idx = _pair_index(n, i, j), ridx = _pair_index(n, j, i);
      int measure = 0; W root_weight = 0.0;
      for (int ti = 0; ti < T; ++ti) {
        const PhyloTree& tree = inst.trees[ti];
        NodeId root = lca_query(tree, d.lca_indices[ti], i, j);
        d.pair_roots[ti][idx] = root;
        d.pair_roots[ti][ridx] = root;
        d.pair_anchor_side_roots[ti][idx] = _child_below_ancestor(tree, root, i);
        d.pair_anchor_side_roots[ti][ridx] = _child_below_ancestor(tree, root, j);
        measure += d.leaf_counts[ti][root];
        if (d.has_internal_weights) root_weight += tree.internal_weight[root];
      }
      d.pair_measure[idx] = measure;
      d.pair_measure[ridx] = measure;
      if (d.has_internal_weights) { d.pair_root_weight[idx] = root_weight; d.pair_root_weight[ridx] = root_weight; }
    }

  if (d.has_internal_weights)
    for (int i = 1; i <= n; ++i)
      for (int j = 1; j <= n; ++j) {
        if (i == j) continue;
        int idx = _pair_index(n, i, j);
        W delta = 0.0;
        for (int ti = 0; ti < T; ++ti) {
          const PhyloTree& tree = inst.trees[ti];
          NodeId root = d.pair_roots[ti][idx];
          delta += strict_connector_weight(tree, d.internal_prefixes[ti], root, i);
        }
        d.singleton_attach[idx] = delta;
      }
  return d;
}

inline bool _pair_below_all(const MultiDPData& d, TaxonId anchor, TaxonId other, TaxonId candidate) {  // MD:243
  int pidx = _pair_index(d.ntaxa, anchor, other);
  for (int ti = 0; ti < (int)d.inst->trees.size(); ++ti)
    if (!_is_descendant(d.preorder_index[ti], d.subtree_last_index[ti], candidate, d.pair_anchor_side_roots[ti][pidx]))
      return false;
  return true;
}
inline W _pair_attach_delta(const MultiDPData& d, TaxonId anchor, TaxonId other, TaxonId candidate) {  // MD:263
  if (!d.has_internal_weights) return 0.0;
  int pidx = _pair_index(d.ntaxa, anchor, other), cidx = _pair_index(d.ntaxa, anchor, candidate);
  W delta = 0.0;
  for (int ti = 0; ti < (int)d.inst->trees.size(); ++ti) {
    const PhyloTree& tree = d.inst->trees[ti];
    delta += strict_connector_weight(tree, d.internal_prefixes[ti], d.pair_roots[ti][pidx], d.pair_roots[ti][cidx]);
  }
  return delta;
}

// ---- lazy ranked workspace (per block) -------------------------------------------------------
struct MultiLazyWorkspace {
  std::vector<std::vector<MASTSolution>> pair_solutions;   // [pair] 1-based (size n*n+1)
  std::vector<std::vector<MASTSolution>> child_solutions;  // [pair]
  std::vector<TaxonId> block;                              // sorted taxa
  std::vector<int> pairs;                                  // pair indices, by measure
  std::vector<int> pair_limit, child_limit;
  std::vector<char> pair_in_progress, child_in_progress;
};

inline int _lazy_target_limit(int current, int requested) {   // MD:565
  if (current >= requested) return current;
  if (current <= 0) return requested;
  return std::max(requested, current << 1);
}

// mutual recursion
inline void _ensure_pair_ranked(MultiLazyWorkspace& lz, const MultiDPData& d, int pair_idx, int requested);

inline void _ensure_child_ranked(MultiLazyWorkspace& lz, const MultiDPData& d, TaxonId anchor, TaxonId other, int requested) {  // MD:571
  if (requested < 1) return;
  int n = d.ntaxa, idx = _pair_index(n, anchor, other);
  if (lz.child_limit[idx] >= requested) return;
  int limit = _lazy_target_limit(lz.child_limit[idx], requested);
  lz.child_in_progress[idx] = 1;
  const SolutionRef ref{ref_kind::STORED_TAXA, 0, 0, 0};
  std::vector<MASTSolution>& sols = lz.child_solutions[idx];
  sols.clear();
  W attach = d.has_internal_weights ? d.singleton_attach[idx] : 0.0;
  W singleton = d.inst->taxon_weight[anchor] + attach;
  insert_unique_solution_score_bounded(sols, MASTSolution{singleton, {anchor}, ref}, limit);
  for (TaxonId candidate : lz.block) {
    if (candidate == anchor || candidate == other) continue;
    if (!(candidate > anchor)) continue;
    if (!_pair_below_all(d, anchor, other, candidate)) continue;
    int cpidx = _pair_index(n, anchor, candidate);
    _ensure_pair_ranked(lz, d, cpidx, limit);
    const std::vector<MASTSolution>& cps = lz.pair_solutions[cpidx];
    if (cps.empty()) continue;
    W delta = _pair_attach_delta(d, anchor, other, candidate);
    for (const MASTSolution& cs : cps) {
      W score = cs.score + delta;
      if ((int)sols.size() >= limit && score <= sols.back().score) break;
      insert_unique_solution_score_bounded(sols, MASTSolution{score, cs.taxa, ref}, limit);
    }
  }
  lz.child_limit[idx] = limit;
  lz.child_in_progress[idx] = 0;
}

inline void _ensure_pair_ranked(MultiLazyWorkspace& lz, const MultiDPData& d, int pair_idx, int requested) {  // MD:623
  if (requested < 1) return;
  int n = d.ntaxa;
  int i = ((pair_idx - 1) % n) + 1, j = ((pair_idx - 1) / n) + 1;
  if (i == j) return;
  if (i > j) { std::swap(i, j); pair_idx = _pair_index(n, i, j); }
  int rev = _pair_index(n, j, i);
  int current = std::max(lz.pair_limit[pair_idx], lz.pair_limit[rev]);
  if (current >= requested) return;
  int limit = _lazy_target_limit(current, requested);
  lz.pair_in_progress[pair_idx] = 1; lz.pair_in_progress[rev] = 1;
  _ensure_child_ranked(lz, d, i, j, limit);
  _ensure_child_ranked(lz, d, j, i, limit);
  std::vector<MASTSolution>& ps = lz.pair_solutions[pair_idx];
  ps.clear();
  W root_score = d.has_internal_weights ? d.pair_root_weight[pair_idx] : 0.0;
  const SolutionRef ref{ref_kind::STORED_TAXA, 0, 0, 0};
  combine_unique_solution_lists_score_bounded(ps, lz.child_solutions[_pair_index(n, i, j)],
                                              lz.child_solutions[_pair_index(n, j, i)], root_score, limit, ref);
  lz.pair_solutions[rev] = lz.pair_solutions[pair_idx];   // share-by-copy (kept in sync each rebuild)
  lz.pair_limit[pair_idx] = limit; lz.pair_limit[rev] = limit;
  lz.pair_in_progress[pair_idx] = 0; lz.pair_in_progress[rev] = 0;
}

// ---- _MultiLazySource heap (MD:48, 672-733) --------------------------------------------------
enum : uint8_t { LAZY_EMPTY = 0, LAZY_SINGLETON = 1, LAZY_PAIR = 2 };
struct MultiLazySource { W score; int block_idx; uint8_t kind; TaxonId taxon; int pair_idx; int rank; };
// 6-key strict total order (MD:672): higher score, then EMPTY<SINGLETON<PAIR, then smaller
// block_idx, pair_idx, taxon, rank.
inline bool _lazy_precedes(const MultiLazySource& a, const MultiLazySource& b) {
  if (a.score != b.score) return a.score > b.score;
  if (a.kind != b.kind) return a.kind < b.kind;
  if (a.block_idx != b.block_idx) return a.block_idx < b.block_idx;
  if (a.pair_idx != b.pair_idx) return a.pair_idx < b.pair_idx;
  if (a.taxon != b.taxon) return a.taxon < b.taxon;
  return a.rank < b.rank;
}
inline void _lazy_sift_up(std::vector<MultiLazySource>& h, int idx) {
  while (idx > 0) { int p = (idx - 1) / 2; if (!_lazy_precedes(h[idx], h[p])) break; std::swap(h[idx], h[p]); idx = p; }
}
inline void _lazy_sift_down(std::vector<MultiLazySource>& h, int idx) {
  int n = (int)h.size();
  while (true) {
    int left = 2 * idx + 1; if (left >= n) break;
    int right = left + 1;
    int best = (right < n && _lazy_precedes(h[right], h[left])) ? right : left;
    if (!_lazy_precedes(h[best], h[idx])) break;
    std::swap(h[idx], h[best]); idx = best;
  }
}
inline void _lazy_push(std::vector<MultiLazySource>& h, MultiLazySource it) { h.push_back(it); _lazy_sift_up(h, (int)h.size() - 1); }
inline MultiLazySource _lazy_pop(std::vector<MultiLazySource>& h) {
  MultiLazySource top = h[0], last = h.back(); h.pop_back();
  if (!h.empty()) { h[0] = last; _lazy_sift_down(h, 0); }
  return top;
}

inline MASTSolution _lazy_source_solution(const MultiLazySource& it, std::vector<MultiLazyWorkspace>& wss) {  // MD:735
  if (it.kind == LAZY_EMPTY) return MASTSolution{it.score, {}, EMPTY_REF};
  if (it.kind == LAZY_SINGLETON) return MASTSolution{it.score, {it.taxon}, SolutionRef{ref_kind::STORED_TAXA, 0, 0, 0}};
  return wss[it.block_idx].pair_solutions[it.pair_idx][it.rank - 1];   // rank 1-based
}

inline bool _lazy_push_pair_rank(std::vector<MultiLazySource>& heap, std::vector<MultiLazyWorkspace>& wss,
                                 const MultiDPData& d, int block_idx, int pair_idx, int rank) {  // MD:752
  MultiLazyWorkspace& lz = wss[block_idx];
  _ensure_pair_ranked(lz, d, pair_idx, rank);
  const std::vector<MASTSolution>& ps = lz.pair_solutions[pair_idx];
  if (rank > (int)ps.size()) return false;
  _lazy_push(heap, MultiLazySource{ps[rank - 1].score, block_idx, LAZY_PAIR, 0, pair_idx, rank});
  return true;
}

inline std::vector<int> _block_pairs(const MultiDPData& d, const std::vector<TaxonId>& block) {  // MD:771
  int n = d.ntaxa;
  std::vector<int> pairs;
  for (int a = 0; a < (int)block.size() - 1; ++a)
    for (int b = a + 1; b < (int)block.size(); ++b)
      pairs.push_back(_pair_index(n, block[a], block[b]));
  std::stable_sort(pairs.begin(), pairs.end(), [&](int x, int y) { return d.pair_measure[x] < d.pair_measure[y]; });
  return pairs;
}

inline MultiLazyWorkspace _make_lazy_workspace(const MultiDPData& d, const std::vector<TaxonId>& block) {  // MD:225
  MultiLazyWorkspace lz;
  int total = d.ntaxa * d.ntaxa;
  lz.pair_solutions.assign(total + 1, {});
  lz.child_solutions.assign(total + 1, {});
  lz.block = block; std::sort(lz.block.begin(), lz.block.end());
  lz.pairs = _block_pairs(d, lz.block);
  lz.pair_limit.assign(total + 1, 0);
  lz.child_limit.assign(total + 1, 0);
  lz.pair_in_progress.assign(total + 1, 0);
  lz.child_in_progress.assign(total + 1, 0);
  return lz;
}

// Branch B (MD:1060): fill rooted_slots[ti][v] with the best non-blacklisted agreement subtree
// rooted at v.  rooted_slots[ti][v] empty optional == "nothing".
inline void _solve_blocks_rooted_blacklist_lazy(
    std::vector<std::vector<std::optional<RootedSolution>>>& rooted_slots, const MultiDPData& d,
    const std::vector<std::vector<TaxonId>>& blocks, const BlacklistSet& blacklist) {
  const int T = (int)d.inst->trees.size();
  std::vector<MultiLazyWorkspace> wss;
  wss.reserve(blocks.size());
  for (const auto& b : blocks) wss.push_back(_make_lazy_workspace(d, b));

  // root_heaps[ti][v]
  std::vector<std::vector<std::vector<MultiLazySource>>> root_heaps(T);
  for (int ti = 0; ti < T; ++ti) root_heaps[ti].assign(d.inst->trees[ti].nv() + 1, {});

  for (int block_idx = 0; block_idx < (int)wss.size(); ++block_idx) {
    MultiLazyWorkspace& lz = wss[block_idx];
    for (TaxonId taxon : lz.block) {
      W score = d.inst->taxon_weight[taxon];
      for (int ti = 0; ti < T; ++ti)
        _lazy_push(root_heaps[ti][taxon], MultiLazySource{score, block_idx, LAZY_SINGLETON, taxon, 0, 1});
    }
    for (int pair_idx : lz.pairs) {
      _ensure_pair_ranked(lz, d, pair_idx, 1);
      const std::vector<MASTSolution>& ps = lz.pair_solutions[pair_idx];
      if (ps.empty()) continue;
      W score = ps[0].score;
      for (int ti = 0; ti < T; ++ti) {
        int root = d.pair_roots[ti][pair_idx];
        _lazy_push(root_heaps[ti][root], MultiLazySource{score, block_idx, LAZY_PAIR, 0, pair_idx, 1});
      }
    }
  }

  for (int ti = 0; ti < T; ++ti) {
    auto& tree_heaps = root_heaps[ti];
    auto& tree_slots = rooted_slots[ti];
    for (int root = 1; root < (int)tree_heaps.size(); ++root) {
      auto& heap = tree_heaps[root];
      if (heap.empty()) continue;
      std::unordered_set<std::vector<TaxonId>, VecHash> seen;
      while (!heap.empty()) {
        MultiLazySource item = _lazy_pop(heap);
        MASTSolution solution = _lazy_source_solution(item, wss);
        if (seen.find(solution.taxa) == seen.end()) {
          seen.insert(solution.taxa);
          if (!is_blacklisted(blacklist, solution.taxa)) {
            tree_slots[root] = RootedSolution{root, solution.score, solution.taxa,
                                              SolutionRef{ref_kind::STORED_TAXA, 0, 0, 0}};
            break;
          }
        }
        if (item.kind == LAZY_PAIR)
          _lazy_push_pair_rank(heap, wss, d, item.block_idx, item.pair_idx, item.rank + 1);
      }
    }
  }
}

// Top-level multi entry, matching mast_solve's else-branch output: the flattened list of
// (score, sorted taxa) over tree1 (asc vertex id), tree2, ..., filtered to |taxa| > 1.
// `blacklist` is a list of taxa sets (each sorted); production passes the master's columns.
inline std::vector<std::pair<W, std::vector<TaxonId>>>
multi_mast_solve(const std::vector<PhyloTree>& trees, const std::vector<W>& taxon_weight_1based,
                 const std::vector<std::vector<int>>& blacklist) {
  MultiMASTInstance inst{trees, taxon_weight_1based};
  MultiDPData d = build_multi_dp_data(inst);
  const int T = (int)trees.size();
  std::vector<std::vector<std::optional<RootedSolution>>> rooted_slots(T);
  for (int ti = 0; ti < T; ++ti) rooted_slots[ti].assign(trees[ti].nv() + 1, std::nullopt);

  int n = (int)taxon_weight_1based.size() - 1;
  std::vector<std::vector<TaxonId>> blocks(1);
  for (int x = 1; x <= n; ++x) blocks[0].push_back(x);
  BlacklistSet bl = normalize_blacklist(blacklist);
  _solve_blocks_rooted_blacklist_lazy(rooted_slots, d, blocks, bl);

  std::vector<std::pair<W, std::vector<TaxonId>>> out;
  for (int ti = 0; ti < T; ++ti)
    for (int v = 1; v < (int)rooted_slots[ti].size(); ++v) {
      const auto& slot = rooted_slots[ti][v];
      if (slot && (int)slot->taxa.size() > 1) out.emplace_back(slot->score, slot->taxa);
    }
  return out;
}

} // namespace maf::mast
