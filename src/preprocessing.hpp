// The LS11 (Linz & Semple 2011) cluster decomposition and its graph-editing primitives.  This is
// the PRODUCTION wrapper around solve_bcp: solve_clusters(maf, solver) recursively splits the
// instance at common clusters, solves each piece with `solver`, and recombines.
//
// PRODUCTION PATH ONLY.
#pragma once
#include "core.hpp"
#include "instance.hpp"
#include "maf_util.hpp"
#include "budget.hpp"        // check_budget
#include "io.hpp"            // rec_to_graph / RecTree
#include <vector>
#include <set>
#include <map>
#include <tuple>
#include <cmath>
#include <functional>
#include <algorithm>
#include <stdexcept>

namespace maf {

// ---- layout / partition assertions -----------------------------------------------------------
// Kept as active checks because cluster decomposition is not a hot path.
inline void assert_maf_layout(const MAFP& maf) {
  if (maf.n < 1) throw std::runtime_error("assert_maf_layout: n<1");
  int root = 2 * maf.n - 1;
  for (const DiGraph& tree : maf.T) {
    if (tree.nv() != root) throw std::runtime_error("assert_maf_layout: nv != 2n-1");
    if (tree.in_degree(root) != 0) throw std::runtime_error("assert_maf_layout: root has parent");
    for (int v = 1; v < root; ++v)
      if (tree.in_degree(v) != 1) throw std::runtime_error("assert_maf_layout: non-root indegree != 1");
    for (int v = 1; v <= maf.n; ++v)
      if (tree.out_degree(v) != 0) throw std::runtime_error("assert_maf_layout: leaf has children");
  }
}

inline void assert_partition(const std::vector<BitSet>& sol, int n) {
  BitSet covered;
  for (const BitSet& comp : sol) {
    if (comp.empty()) throw std::runtime_error("assert_partition: empty component");
    for (int x : comp) {
      if (x < 1 || x > n) throw std::runtime_error("assert_partition: label out of range");
      if (covered.contains(x)) throw std::runtime_error("assert_partition: overlap");
      covered.insert(x);
    }
  }
  BitSet all; for (int i = 1; i <= n; ++i) all.insert(i);
  if (!(covered == all)) throw std::runtime_error("assert_partition: not a cover");
}

// Subgraph on the vertices in `perm`, RELABELED 1..|perm| in perm order, keeping every edge with
// both endpoints in perm. `perm` is a list of OLD vertex ids.
inline DiGraph induced_subgraph(const DiGraph& g, const std::vector<int>& perm) {
  std::vector<int> pos(g.nv() + 1, 0);             // old vertex -> new index (1-based), 0 = absent
  for (int j = 0; j < (int)perm.size(); ++j) pos[perm[j]] = j + 1;
  DiGraph h((int)perm.size());
  for (int j = 0; j < (int)perm.size(); ++j) {
    int u = perm[j];
    for (int v : g.out_neighbors(u)) if (pos[v]) h.add_edge(pos[u], pos[v]);
  }
  h.finalize();
  return h;
}

// Sorted leaf-index vectors of all minimal non-trivial clusters common to all trees (a leaf set that
// is a node's descendant set in EVERY tree; minimal = no proper subset is also such a cluster).
inline std::vector<std::vector<int>> minimal_clusters(const MAFP& maf) {
  assert_maf_layout(maf);
  std::vector<std::vector<BitSet>> desc;
  for (const DiGraph& t : maf.T) desc.push_back(compute_descendants(t));

  // all_desc = intersect over trees of the SET of descendant BitSets.
  std::set<BitSet> common;
  for (const BitSet& b : desc[0]) common.insert(b);
  for (size_t i = 1; i < desc.size(); ++i) {
    std::set<BitSet> cur(desc[i].begin(), desc[i].end());
    std::set<BitSet> next;
    for (const BitSet& b : common) if (cur.count(b)) next.insert(b);
    common.swap(next);
  }
  // Drop singletons and the full leaf set.
  std::vector<BitSet> nontrivial;
  for (const BitSet& c : common) { int s = c.size(); if (1 < s && s < maf.n) nontrivial.push_back(c); }
  // Keep only minimal: no proper subset d ⊊ c with d also non-trivial.
  std::vector<std::vector<int>> minimal;
  for (const BitSet& c : nontrivial) {
    bool has_proper_subset = false;
    for (const BitSet& d : nontrivial) {
      if (&d == &c) continue;
      // d ⊊ c : d subset of c and d != c
      if (d != c && set_diff(d, c).empty()) { has_proper_subset = true; break; }
    }
    if (!has_proper_subset) minimal.push_back(c.to_vec());   // to_vec is ascending (== sort)
  }
  std::sort(minimal.begin(), minimal.end());                 // sort([sort(collect(c)) ...])
  return minimal;
}

// Rooted restriction maf|leaves with all degree-2 vertices contracted; returns (sub_maf, leaf_map)
// where leaf_map[k] = leaves[k].
// `leaves` must be sorted, distinct, non-empty, ⊆ 1..n.
inline std::pair<MAFP, std::vector<int>> restrict_to_leaves(const MAFP& maf, const std::vector<int>& leaves) {
  assert_maf_layout(maf);
  if (leaves.empty()) throw std::runtime_error("restrict_to_leaves: empty");
  int nc = (int)leaves.size();
  if (nc == 1) {
    MAFP sub; sub.n = 1;
    for (size_t i = 0; i < maf.T.size(); ++i) sub.T.emplace_back(1);   // single isolated vertex
    assert_maf_layout(sub);
    return {sub, leaves};
  }

  std::vector<int> leaf_pos(maf.n + 1, 0);          // old leaf label -> new position (1..nc)
  for (int k = 0; k < nc; ++k) leaf_pos[leaves[k]] = k + 1;

  MAFP sub; sub.n = nc;
  for (const DiGraph& tree : maf.T) {
    std::vector<char> keep(tree.nv() + 1, 0);
    for (int x : leaves) { int v = x; keep[v] = 1; while (tree.in_degree(v) == 1) { v = tree.in_neighbors(v)[0]; keep[v] = 1; } }

    // build(v): contract degree-1 (single kept child) nodes; leaves carry their new position.
    std::function<RecTree(int)> build = [&](int v) -> RecTree {
      if (tree.out_degree(v) == 0) return RecTree::leaf(leaf_pos[v]);   // kept leaf
      std::vector<RecTree> childs;
      for (int w : tree.out_neighbors(v)) { if (!keep[w]) continue; childs.push_back(build(w)); }
      if (childs.empty()) throw std::runtime_error("restrict_to_leaves: empty subtree");
      if (childs.size() == 1) return std::move(childs[0]);
      return RecTree::node(std::move(childs));
    };
    RecTree rec = build(2 * maf.n - 1);
    sub.T.push_back(rec_to_graph(rec, nc));
  }
  assert_maf_layout(sub);
  return {sub, leaves};
}

// Replace the subtree spanned by the common cluster `group` with a single representative leaf at
// position n_noncl+1. Returns
// (reduced_maf, leaf_map): leaf_map[k] = original label for surviving leaves, -1 for the rep.
inline std::pair<MAFP, std::vector<int>> collapse_to_leaf(const MAFP& maf, const std::vector<int>& group) {
  assert_maf_layout(maf);
  if (group.empty() || (int)group.size() >= maf.n) throw std::runtime_error("collapse_to_leaf: bad group");
  BitSet group_set(group.begin(), group.end());
  std::vector<int> noncl_leaves;
  for (int x = 1; x <= maf.n; ++x) if (!group_set.contains(x)) noncl_leaves.push_back(x);  // sorted
  int n_noncl = (int)noncl_leaves.size();
  int orig_root = 2 * maf.n - 1;

  MAFP reduced; reduced.n = n_noncl + 1;
  for (const DiGraph& tree : maf.T) {
    std::vector<int> iv = inner_vs(tree, group);            // ascending
    BitSet iv_set(iv.begin(), iv.end());
    int r = iv.back();                                      // Steiner root = highest-indexed inner
    std::vector<int> noncl_inner;
    for (int x = maf.n + 1; x <= orig_root - 1; ++x) if (!iv_set.contains(x)) noncl_inner.push_back(x);
    std::vector<int> perm;
    perm.insert(perm.end(), noncl_leaves.begin(), noncl_leaves.end());
    perm.push_back(r);
    perm.insert(perm.end(), noncl_inner.begin(), noncl_inner.end());
    perm.push_back(orig_root);
    reduced.T.push_back(induced_subgraph(tree, perm));
  }
  std::vector<int> leaf_map = noncl_leaves;
  leaf_map.push_back(-1);
  assert_maf_layout(reduced);
  return {reduced, leaf_map};
}

// Add leaf n+1 as a sibling of the root in all trees.
// Old inner nodes n+1..2n-1 shift up by one to free index n+1 for the new leaf.
inline MAFP add_pendant_leaf(const MAFP& maf) {
  assert_maf_layout(maf);
  int n = maf.n, n_new = n + 1;
  auto shift = [n](int k) { return k <= n ? k : k + 1; };
  MAFP out; out.n = n_new;
  for (const DiGraph& tree : maf.T) {
    DiGraph nt(2 * n_new - 1);
    for (int u = 1; u <= tree.nv(); ++u)
      for (int v : tree.out_neighbors(u)) nt.add_edge(shift(u), shift(v));
    nt.add_edge(2 * n_new - 1, 2 * n);     // new root -> old root (shifted to 2n)
    nt.add_edge(2 * n_new - 1, n_new);     // new root -> new pendant leaf
    nt.finalize();
    out.T.push_back(std::move(nt));
  }
  assert_maf_layout(out);
  return out;
}

// ---- cluster sequence + solve ----------------------------------------------------------------
// Encoded labels: positive l = original leaf l; negative -j = representative a_j of cluster seq[j].
struct ClusterEntry { MAFP sub_maf; std::vector<int> leaf_map; int rep_label; };

inline void cluster_sequence(const MAFP& maf,
                             std::vector<ClusterEntry>& seq, MAFP& main_maf, std::vector<int>& main_label_map) {
  assert_maf_layout(maf);
  std::vector<int> label_map(maf.n);
  for (int k = 0; k < maf.n; ++k) label_map[k] = k + 1;
  MAFP current = maf;
  seq.clear();

  while (true) {
    std::vector<std::vector<int>> clusters = minimal_clusters(current);
    if (clusters.empty()) break;

    // Snapshot encoded labels for each cluster leaf before any collapse this round.
    std::vector<std::vector<int>> clusters_labels;
    for (const auto& c : clusters) {
      std::vector<int> lc;
      for (int k : c) lc.push_back(label_map[k - 1]);   // c holds 1-based positions in `current`
      clusters_labels.push_back(std::move(lc));
    }
    int seq_base = (int)seq.size();

    // Extract cluster subinstances (read-only from current).
    for (int ci = 0; ci < (int)clusters.size(); ++ci) {
      int i = seq_base + ci + 1;                          // 1-based cluster id
      auto [sub, _] = restrict_to_leaves(current, clusters[ci]);
      MAFP sub_rho = add_pendant_leaf(sub);
      std::vector<int> lmap = clusters_labels[ci];
      lmap.push_back(0);                                  // rho_i encoded as 0 (last position)
      seq.push_back(ClusterEntry{std::move(sub_rho), std::move(lmap), -i});
    }

    // Collapse clusters one by one, updating label_map after each.
    for (int ci = 0; ci < (int)clusters_labels.size(); ++ci) {
      int i = seq_base + ci + 1;
      // Encoded labels can be negative (cluster representatives), so use std::set<int> not BitSet.
      std::set<int> label_set(clusters_labels[ci].begin(), clusters_labels[ci].end());
      std::vector<int> current_indices;                       // k ascending => already sorted
      for (int k = 1; k <= current.n; ++k)
        if (label_set.count(label_map[k - 1])) current_indices.push_back(k);
      auto [collapsed, red_lmap] = collapse_to_leaf(current, current_indices);
      std::vector<int> new_lmap(red_lmap.size());
      for (int k = 0; k < (int)red_lmap.size(); ++k) {
        int lk = red_lmap[k];
        new_lmap[k] = (lk == -1) ? -i : label_map[lk - 1];
      }
      current = std::move(collapsed);
      label_map = std::move(new_lmap);
    }
  }
  main_maf = std::move(current);
  main_label_map = std::move(label_map);
}

// Reconstruct a solution from cluster choices and the reduced main solution.
struct ClusterChoice { std::vector<BitSet> sol; std::vector<int> leaf_map; };

inline std::vector<BitSet> reconstruct_solution(const std::vector<ClusterChoice>& cluster_choices,
                                                const std::vector<BitSet>& sol_main,
                                                const std::vector<int>& main_label_map, int n) {
  std::vector<BitSet> result;
  std::vector<BitSet> expanded_rho(cluster_choices.size());   // ρ_i component (ρ_i excluded)

  auto expand = [&](BitSet& bs, int label) {
    if (label == 0) throw std::runtime_error("reconstruct_solution: expand label 0");
    if (label > 0) bs.insert(label);
    else bs.union_with(expanded_rho[-label - 1]);             // -label is 1-based cluster id
  };

  for (int i = 0; i < (int)cluster_choices.size(); ++i) {
    const std::vector<BitSet>& cl_sol = cluster_choices[i].sol;
    const std::vector<int>& lmap = cluster_choices[i].leaf_map;
    assert_partition(cl_sol, (int)lmap.size());

    int rho_pos = 0;                                          // 1-based position of label 0, or 0
    for (int k = 0; k < (int)lmap.size(); ++k) if (lmap[k] == 0) { rho_pos = k + 1; break; }

    if (rho_pos == 0) {
      expanded_rho[i] = BitSet();
    } else {
      int rho_comp = -1;
      for (int c = 0; c < (int)cl_sol.size(); ++c) if (cl_sol[c].contains(rho_pos)) { rho_comp = c; break; }
      if (rho_comp < 0) throw std::runtime_error("reconstruct_solution: rho not in any component");
      BitSet rho_bs;
      for (int k : cl_sol[rho_comp]) { if (k == rho_pos) continue; expand(rho_bs, lmap[k - 1]); }
      expanded_rho[i] = rho_bs;
    }

    for (const BitSet& comp : cl_sol) {
      if (rho_pos != 0 && comp.contains(rho_pos)) continue;   // ρ_i component never survives
      BitSet bs;
      for (int k : comp) expand(bs, lmap[k - 1]);
      if (!bs.empty()) result.push_back(bs);
    }
  }

  assert_partition(sol_main, (int)main_label_map.size());
  for (const BitSet& comp : sol_main) {
    BitSet bs;
    for (int k : comp) expand(bs, main_label_map[k - 1]);
    if (!bs.empty()) result.push_back(bs);
  }

  assert_partition(result, n);
  return result;
}

// Active leaf positions for LS11 Step 2'.
inline std::vector<int> active_leaf_positions(const std::vector<int>& leaf_map,
                                              const std::set<int>& rho_singleton, bool keep_rho) {
  std::vector<int> positions;
  for (int k = 0; k < (int)leaf_map.size(); ++k) {
    int label = leaf_map[k];
    if (label == 0) { if (keep_rho) positions.push_back(k + 1); }
    else if (!(label < 0 && rho_singleton.count(-label))) positions.push_back(k + 1);
  }
  return positions;
}

// The solver interface used by solve_clusters: solver(sub_maf) -> (sol partition, dual LB).
using ClusterSolver = std::function<std::pair<std::vector<BitSet>, double>(const MAFP&)>;

// Restrict to active positions and solve.
// Returns (sol, dual, enc_map) where enc_map[k] = encoded label of sub leaf k.
inline std::tuple<std::vector<BitSet>, double, std::vector<int>>
solve_restricted_instance(const MAFP& maf, const std::vector<int>& full_leaf_map,
                          const std::vector<int>& active_positions, const ClusterSolver& solver) {
  if (active_positions.empty()) return {std::vector<BitSet>{}, 0.0, std::vector<int>{}};
  auto [sub_maf, pos_map] = restrict_to_leaves(maf, active_positions);
  std::vector<int> enc_map;
  for (int k : pos_map) enc_map.push_back(full_leaf_map[k - 1]);
  if (sub_maf.n == 1) return {std::vector<BitSet>{BitSet{1}}, 1.0, enc_map};
  auto [sol, dual] = solver(sub_maf);
  assert_partition(sol, sub_maf.n);
  return {sol, dual, enc_map};
}

// LS11 cluster decomposition.
inline std::pair<std::vector<BitSet>, double> solve_clusters(const MAFP& maf, const ClusterSolver& solver) {
  assert_maf_layout(maf);
  std::vector<ClusterEntry> seq; MAFP main_maf; std::vector<int> main_label_map;
  cluster_sequence(maf, seq, main_maf, main_label_map);
  if (seq.empty()) return solver(maf);

  std::set<int> rho_singleton;
  std::vector<ClusterChoice> cluster_choices(seq.size());
  std::vector<double> cluster_duals(seq.size());
  std::vector<int> cluster_costs(seq.size());

  for (int i = 0; i < (int)seq.size(); ++i) {
    check_budget();                 // bail out on timeout/SIGINT (per cluster sub-solve)
    const MAFP& sub_maf = seq[i].sub_maf;
    const std::vector<int>& leaf_map = seq[i].leaf_map;

    // Option A: keep ρ_i.
    std::vector<int> active_A = active_leaf_positions(leaf_map, rho_singleton, /*keep_rho=*/true);
    auto [sol_A, dual_A, map_A] = solve_restricted_instance(sub_maf, leaf_map, active_A, solver);
    int cost_A = (int)sol_A.size() - 1;
    double dual_cost_A = std::max(std::ceil(dual_A - 1e-8) - 1.0, 0.0);

    // Option B: remove ρ_i (force it to be a singleton).
    std::vector<int> active_B = active_leaf_positions(leaf_map, rho_singleton, /*keep_rho=*/false);
    auto [sol_B, dual_B, map_B] = solve_restricted_instance(sub_maf, leaf_map, active_B, solver);
    int cost_B = (int)sol_B.size();
    double dual_cost_B = std::max(std::ceil(dual_B - 1e-8), 0.0);

    if (cost_B <= cost_A) {                                   // ties -> B (ρ_i singleton)
      cluster_choices[i] = ClusterChoice{sol_B, map_B};
      cluster_duals[i] = dual_cost_B; cluster_costs[i] = cost_B;
      rho_singleton.insert(i + 1);                            // 1-based cluster id
    } else {
      cluster_choices[i] = ClusterChoice{sol_A, map_A};
      cluster_duals[i] = dual_cost_A; cluster_costs[i] = cost_A;
    }
  }

  std::vector<int> active_main = active_leaf_positions(main_label_map, rho_singleton, /*keep_rho=*/false);
  auto [sol_main, dual_main, main_map] = solve_restricted_instance(main_maf, main_label_map, active_main, solver);

  std::vector<BitSet> sol = reconstruct_solution(cluster_choices, sol_main, main_map, maf.n);
  int primal = 0; for (int c : cluster_costs) primal += c; primal += (int)sol_main.size();
  double dual = 0.0; for (double d : cluster_duals) dual += d;
  dual += std::max(std::ceil(dual_main - 1e-8), 0.0);

  if ((int)sol.size() != primal) throw std::runtime_error("solve_clusters: size != primal");
  if (!(dual <= (double)sol.size() + 1e-8)) throw std::runtime_error("solve_clusters: dual > primal");
  return {sol, std::max(dual, 0.0)};
}

} // namespace maf
