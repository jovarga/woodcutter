// Used by ctx_price! once cuts exist (the no-cuts path uses multi_dp/centroid). A self-contained
// DFS branch-and-bound over taxon insertion points; reuses only PhyloTree/LCA/connector primitives
// (trees.hpp, lca.hpp). See docs/PRICING_SPEC.md. 1-based vertices/taxa; pair index 1-based
// column-major a+(b-1)*n.
//
// STATUS: foundation (structs + input normalization + precompute helpers). The B&B core
// (_pricing_pair_conflicts_direct, _pricing_precompute assembly, cut activation, child-state,
// search, :rooted_best output) is added next — see docs/PRICING_SPEC.md §3 + §7.
#pragma once
#include "trees.hpp"
#include "lca.hpp"
#include "../core.hpp"
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <cstdint>
#include <optional>
#include <set>
#include <tuple>
#include <array>
#include <limits>

namespace maf::mast {

// ---- output carriers (TR:295/349) ------------------------------------------------------------
// Downstream (ctx_price!) reads only .objective and .taxa of each solution.
struct MASTPricingSolution {
  W objective = 0.0;
  std::vector<TaxonId> taxa;        // sorted ascending
  int base_size = 0;               // == taxa.size()
  W penalty = 0.0;                 // current_internal_weight + cut_penalty; objective == base_size + penalty
};
struct MASTPricingResult {
  std::vector<MASTPricingSolution> solutions;   // the only field the caller consumes
  W objective = 0.0;
  bool improving = false;
};

// ---- B&B state (PR:23/28/41/71) --------------------------------------------------------------
struct PricingInsertionPoint {
  TaxonId representative_descendant;
  std::vector<TaxonId> taxa;
};
struct PricingSearchState {
  TaxonId initial_a, initial_b;
  std::vector<TaxonId> selected_taxa;           // sorted
  std::vector<PricingInsertionPoint> insertion_points;
  std::vector<int> insertion_point_of;          // length n+1 (1-based taxon -> slot, 0 if not addable)
  W current_internal_weight = 0.0;
  W mast_upper_bound = 0.0;
  BitSet unconflicted_variables;
  BitSet inactive_cuts;
  W cut_penalty = 0.0;
};
struct PricingIncumbent { MASTPricingSolution solution; };

struct PricingData {
  std::vector<PhyloTree> trees;
  std::vector<std::vector<std::vector<NodeId>>> lca_tables;   // [tree][a][b], a,b in 1..nv (dense)
  std::vector<std::vector<int>> preorder_first, preorder_last; // [tree][v] 1-based
  std::vector<std::vector<W>> internal_prefixes;              // [tree][v]
  std::vector<int> rank;                                      // rank[x] = x (taxon_order = 1..n)
  std::vector<std::pair<TaxonId, TaxonId>> canonical_pairs;   // lexicographic (a<b)
  std::vector<std::vector<TaxonId>> initial_side_taxa;        // [pair index] (1-based, size n*n+1)
  std::vector<W> weighted_pair_bounds;                        // [pair index]
  int variable_count = 0;
  std::vector<BitSet> pair_conflicting_variables;            // [pair index]
  std::vector<BitSet> required_variables;                    // [cut] member var indices
  std::vector<W> cut_weights;                                 // [cut], <= 0
  int ntaxa = 0;
};

inline int pricing_pair_index(int n, int a, int b) { return a + (b - 1) * n; }   // PR:75, 1-based

// _normalized_taxa (TR:17): sort+unique, bounds-checked.
inline std::vector<TaxonId> normalized_taxa(std::vector<TaxonId> taxa, int max_taxon) {
  std::sort(taxa.begin(), taxa.end());
  taxa.erase(std::unique(taxa.begin(), taxa.end()), taxa.end());
  for (TaxonId t : taxa) if (t < 1 || t > max_taxon) throw std::runtime_error("taxon id out of bounds");
  return taxa;
}

// ---- input normalization (PR:121-169, 1043) --------------------------------------------------
// _pricing_flatten_cut_sets (PR:132): every member occurrence becomes a variable (NO cross-cut
// dedup); each cut owns the 1-based indices of its members in push order.
inline void pricing_flatten_cut_sets(const std::vector<std::vector<std::vector<TaxonId>>>& cut_sets, int ntaxa,
                                     std::vector<std::vector<TaxonId>>& variables, std::vector<std::vector<int>>& cuts) {
  variables.clear(); cuts.clear();
  for (const auto& cut : cut_sets) {
    if (cut.empty()) throw std::runtime_error("pricing cuts must be nonempty");
    std::vector<int> ids;
    for (const auto& var : cut) {
      std::vector<TaxonId> taxa = normalized_taxa(var, ntaxa);
      if (taxa.size() < 2) throw std::runtime_error("pricing variables need >=2 taxa");
      variables.push_back(std::move(taxa));
      ids.push_back((int)variables.size());        // 1-based
    }
    cuts.push_back(std::move(ids));
  }
}

// _pricing_normalize_cuts (PR:149): bounds + no-dup-per-cut checks; clamp weights in (0,atol] -> 0;
// assert weight <= atol.
inline void pricing_normalize_cuts(const std::vector<std::vector<int>>& cuts, const std::vector<W>& weights,
                                   int nvariables, W atol, std::vector<std::vector<int>>& out_cuts,
                                   std::vector<W>& out_weights) {
  if (cuts.size() != weights.size()) throw std::runtime_error("cuts/weights length mismatch");
  out_cuts.clear(); out_weights.clear();
  for (size_t c = 0; c < cuts.size(); ++c) {
    const std::vector<int>& ids = cuts[c];
    if (ids.empty()) throw std::runtime_error("pricing cuts must be nonempty");
    BitSet seen;
    for (int id : ids) {
      if (id < 1 || id > nvariables) throw std::runtime_error("cut var index out of bounds");
      if (seen.contains(id)) throw std::runtime_error("duplicate variable index in cut");
      seen.insert(id);
    }
    W w = weights[c];
    if (!std::isfinite(w)) throw std::runtime_error("cut weight must be finite");
    if (w > atol) throw std::runtime_error("cut weight must be non-positive up to atol");
    out_cuts.push_back(ids);
    out_weights.push_back((0.0 < w && w <= atol) ? 0.0 : w);
  }
}

// _pricing_graph_internal_weights (PR:1043): beta Dict{(i,v)->w} -> per-tree full-form (length nv)
// weight vector (0-based here: vec[ti][v-1] = weight of vertex v; leaves stay 0). Order-independent.
inline std::vector<std::vector<W>>
pricing_graph_internal_weights(const std::vector<DiGraph>& graphs,
                               const std::unordered_map<long long, W>& beta /*key = i*BIG+v*/, long long BIG) {
  std::vector<std::vector<W>> vecs(graphs.size());
  for (size_t ti = 0; ti < graphs.size(); ++ti) vecs[ti].assign(graphs[ti].nv(), 0.0);
  for (const auto& kv : beta) {
    int i = (int)(kv.first / BIG), v = (int)(kv.first % BIG);
    vecs[i][v - 1] = kv.second;                    // full-form, 0-based by vertex
  }
  return vecs;
}

// Assert internal weights finite ONLY. The old non-positivity check + (0,atol]->0 clamp were removed
// because branching adds branch-constraint duals to beta, which can make internal weights POSITIVE;
// the B&B objective/bound DP handles positive weights directly. (atol kept for signature
// compatibility; no longer used here.)
inline void pricing_clamp_internal_weights(std::vector<PhyloTree>& trees, W atol) {
  (void)atol;
  for (PhyloTree& t : trees)
    for (NodeId v = t.ntaxa + 1; v <= t.nv(); ++v) {
      W w = t.internal_weight[v];
      if (!std::isfinite(w)) throw std::runtime_error("internal weights must be finite");
    }
}

// ---- precompute helpers (PR:171-290) ---------------------------------------------------------
// dense nv x nv LCA matrix (PR:171). Indexed [a][b], 1-based; entry 0/row 0 unused.
inline std::vector<std::vector<NodeId>> pricing_lca_table(const PhyloTree& tree) {
  int n = tree.nv();
  LCAIndex idx = build_lca(tree);
  std::vector<std::vector<NodeId>> table(n + 1, std::vector<NodeId>(n + 1, 0));
  for (NodeId a = 1; a <= n; ++a) {
    table[a][a] = a;
    for (NodeId b = a + 1; b <= n; ++b) { NodeId q = lca_query(tree, idx, a, b); table[a][b] = q; table[b][a] = q; }
  }
  return table;
}
// preorder in/out timestamps (PR:186), push right then left.
inline void pricing_preorder_ranges(const PhyloTree& tree, std::vector<int>& first, std::vector<int>& last) {
  int n = tree.nv();
  first.assign(n + 1, 0); last.assign(n + 1, 0);
  int timer = 0;
  std::vector<NodeId> stack{tree.root};
  while (!stack.empty()) {
    NodeId v = stack.back(); stack.pop_back();
    if (v < 0) { last[-v] = timer; continue; }
    ++timer; first[v] = timer;
    stack.push_back(-v);
    if (!tree.isleaf(v)) { stack.push_back(tree.right[v]); stack.push_back(tree.left[v]); }
  }
}
inline bool pricing_is_descendant(const std::vector<int>& first, const std::vector<int>& last, NodeId node, NodeId anc) {
  return first[anc] <= first[node] && first[node] <= last[anc];      // PR:209
}
inline NodeId pricing_child_containing(const PhyloTree& tree, const std::vector<int>& first, const std::vector<int>& last,
                                       NodeId ancestor, NodeId node) {                     // PR:213
  NodeId left = tree.left[ancestor];
  return pricing_is_descendant(first, last, node, left) ? left : tree.right[ancestor];
}

// _pricing_initial_sides (PR:218): for ordered (a,b), taxa x (rank[x]>rank[a]) on a's side of
// lca(a,b) in EVERY tree.
inline std::vector<std::vector<TaxonId>>
pricing_initial_sides(const std::vector<PhyloTree>& trees, const std::vector<std::vector<std::vector<NodeId>>>& lca_tables,
                      const std::vector<std::vector<int>>& firsts, const std::vector<std::vector<int>>& lasts,
                      const std::vector<int>& rank) {
  int n = trees[0].ntaxa, T = (int)trees.size();
  std::vector<std::vector<TaxonId>> sides(n * n + 1);
  for (TaxonId a = 1; a <= n; ++a)
    for (TaxonId b = 1; b <= n; ++b) {
      if (a == b) continue;
      std::vector<TaxonId>& out = sides[pricing_pair_index(n, a, b)];
      for (TaxonId x = 1; x <= n; ++x) {
        if (x == a || x == b || rank[x] <= rank[a]) continue;
        bool valid = true;
        for (int i = 0; i < T; ++i) {
          NodeId q = lca_tables[i][a][b];
          NodeId side = pricing_child_containing(trees[i], firsts[i], lasts[i], q, a);
          if (!pricing_is_descendant(firsts[i], lasts[i], x, side)) { valid = false; break; }
        }
        if (valid) out.push_back(x);
      }
    }
  return sides;
}

// Admissible per-pair UB DP. Pairs processed in increasing total-leaf-count-at-LCA order; ties keep
// lex (a,b) push order. Reads earlier pairs' bounds, so this order is reproducibility-critical.
inline std::vector<W>
pricing_weighted_pair_bounds(const std::vector<PhyloTree>& trees, const std::vector<std::vector<std::vector<NodeId>>>& lca_tables,
                             const std::vector<std::vector<W>>& prefixes, const std::vector<std::vector<TaxonId>>& initial_sides) {
  int n = trees[0].ntaxa, T = (int)trees.size();
  int total = n * n;
  std::vector<W> pair_bounds(total + 1, -std::numeric_limits<W>::infinity());
  std::vector<W> child_bounds(total + 1, -std::numeric_limits<W>::infinity());
  std::vector<std::vector<int>> counts(T);
  for (int i = 0; i < T; ++i) counts[i] = leaf_count(trees[i]);

  std::vector<std::pair<TaxonId, TaxonId>> pairs;
  for (TaxonId a = 1; a <= n - 1; ++a) for (TaxonId b = a + 1; b <= n; ++b) pairs.emplace_back(a, b);
  auto key = [&](const std::pair<TaxonId, TaxonId>& p) {
    long long s = 0; for (int i = 0; i < T; ++i) s += counts[i][lca_tables[i][p.first][p.second]]; return s;
  };
  std::stable_sort(pairs.begin(), pairs.end(), [&](const auto& x, const auto& y) { return key(x) < key(y); });

  for (const auto& pr : pairs) {
    TaxonId a = pr.first, b = pr.second;
    for (int side_sel = 0; side_sel < 2; ++side_sel) {
      TaxonId anchor = side_sel == 0 ? a : b, other = side_sel == 0 ? b : a;
      int idx = pricing_pair_index(n, anchor, other);
      W score = 1.0;
      for (int i = 0; i < T; ++i) {
        NodeId q = lca_tables[i][anchor][other];
        score += strict_connector_weight(trees[i], prefixes[i], q, anchor);
      }
      for (TaxonId x : initial_sides[idx]) {
        int child_idx = pricing_pair_index(n, anchor, x);
        W candidate = pair_bounds[child_idx];
        if (!std::isfinite(candidate)) continue;
        for (int i = 0; i < T; ++i) {
          NodeId parent_root = lca_tables[i][anchor][other];
          NodeId child_root = lca_tables[i][anchor][x];
          candidate += strict_connector_weight(trees[i], prefixes[i], parent_root, child_root);
        }
        score = std::max(score, candidate);
      }
      child_bounds[idx] = score;
    }
    int idx = pricing_pair_index(n, a, b);
    W score = child_bounds[idx] + child_bounds[pricing_pair_index(n, b, a)];
    for (int i = 0; i < T; ++i) score += trees[i].internal_weight[lca_tables[i][a][b]];
    pair_bounds[idx] = score;
    pair_bounds[pricing_pair_index(n, b, a)] = score;
  }
  return pair_bounds;
}

// ---- B&B core (PR:361-1005) ------------------------------------------------------------------
// _pricing_taxa_lca (PR:361): fold pairwise LCA over a variable's taxa via the dense table.
inline NodeId pricing_taxa_lca(const std::vector<std::vector<NodeId>>& table, const std::vector<TaxonId>& taxa) {
  NodeId root = taxa[0];
  for (size_t i = 1; i < taxa.size(); ++i) root = table[root][taxa[i]];
  return root;
}

// _pricing_pair_conflicts_direct (PR:400): per ordered pair, the variables whose induced partition
// separates a from b (union over trees of path-conflicts at the pair LCA).
inline std::vector<BitSet>
pricing_pair_conflicts_direct(const std::vector<PhyloTree>& trees, const std::vector<std::vector<TaxonId>>& variables,
                              const std::vector<std::vector<std::vector<NodeId>>>& lca_tables,
                              const std::vector<std::pair<TaxonId, TaxonId>>& canonical_pairs) {
  int n = trees[0].ntaxa, T = (int)trees.size();
  std::vector<BitSet> pair_conflicts(n * n + 1);
  if (variables.empty()) return pair_conflicts;

  for (int ti = 0; ti < T; ++ti) {
    const PhyloTree& tree = trees[ti];
    const auto& lca_table = lca_tables[ti];
    int nv = tree.nv();
    std::vector<BitSet> covering(nv + 1);
    std::vector<int> last_variable(nv + 1, 0);
    for (int vidx = 1; vidx <= (int)variables.size(); ++vidx) {
      const std::vector<TaxonId>& variable = variables[vidx - 1];
      NodeId root = pricing_taxa_lca(lca_table, variable);
      for (TaxonId taxon : variable) {
        NodeId vertex = tree.parent[taxon];
        while (true) {
          if (last_variable[vertex] == vidx) break;
          last_variable[vertex] = vidx;
          covering[vertex].insert(vidx);
          if (vertex == root) break;
          vertex = tree.parent[vertex];
        }
      }
    }
    // path_conflicts[taxon][vertex] = union of covering[w] for w on taxon's path up to vertex.
    std::vector<std::vector<BitSet>> path_conflicts(n + 1, std::vector<BitSet>(nv + 1));
    for (TaxonId taxon = 1; taxon <= n; ++taxon) {
      BitSet conflicts;
      NodeId vertex = tree.parent[taxon];
      while (vertex != 0) {
        conflicts.union_with(covering[vertex]);
        path_conflicts[taxon][vertex] = conflicts;
        vertex = tree.parent[vertex];
      }
    }
    for (const auto& pr : canonical_pairs) {
      TaxonId a = pr.first, b = pr.second;
      NodeId root = lca_table[a][b];
      BitSet& c = pair_conflicts[pricing_pair_index(n, a, b)];
      c.union_with(path_conflicts[a][root]);
      c.union_with(path_conflicts[b][root]);
    }
  }
  for (const auto& pr : canonical_pairs)
    pair_conflicts[pricing_pair_index(n, pr.second, pr.first)] = pair_conflicts[pricing_pair_index(n, pr.first, pr.second)];
  return pair_conflicts;
}

// _pricing_precompute (PR:458): assemble PricingData.
inline PricingData pricing_precompute(std::vector<PhyloTree> trees, const std::vector<std::vector<TaxonId>>& variables,
                                      const std::vector<std::vector<int>>& cuts, const std::vector<W>& cut_weights) {
  PricingData d;
  int n = trees[0].ntaxa;
  d.ntaxa = n;
  d.rank.assign(n + 1, 0); for (int x = 1; x <= n; ++x) d.rank[x] = x;
  for (const PhyloTree& t : trees) d.lca_tables.push_back(pricing_lca_table(t));
  d.preorder_first.resize(trees.size()); d.preorder_last.resize(trees.size());
  for (size_t i = 0; i < trees.size(); ++i) pricing_preorder_ranges(trees[i], d.preorder_first[i], d.preorder_last[i]);
  for (const PhyloTree& t : trees) d.internal_prefixes.push_back(internal_prefix(t));
  d.initial_side_taxa = pricing_initial_sides(trees, d.lca_tables, d.preorder_first, d.preorder_last, d.rank);
  d.weighted_pair_bounds = pricing_weighted_pair_bounds(trees, d.lca_tables, d.internal_prefixes, d.initial_side_taxa);
  for (TaxonId a = 1; a <= n - 1; ++a) for (TaxonId b = a + 1; b <= n; ++b) d.canonical_pairs.emplace_back(a, b);
  d.pair_conflicting_variables = pricing_pair_conflicts_direct(trees, variables, d.lca_tables, d.canonical_pairs);
  d.variable_count = (int)variables.size();
  d.required_variables.resize(cuts.size());
  for (size_t c = 0; c < cuts.size(); ++c) for (int id : cuts[c]) d.required_variables[c].insert(id);
  d.cut_weights = cut_weights;
  d.trees = std::move(trees);
  return d;
}

// _pricing_activate_completed_cuts! (PR:543): activate (add weight to penalty) every still-inactive
// cut whose required variables are disjoint from `unconflicted`. Cut indices iterated ascending.
inline W pricing_activate_completed_cuts(BitSet& inactive, const BitSet& unconflicted, const PricingData& d, W penalty) {
  for (int cut = 0; cut < (int)d.required_variables.size(); ++cut)
    if (inactive.contains(cut) && isdisjoint(unconflicted, d.required_variables[cut])) {
      inactive.erase(cut);
      penalty += d.cut_weights[cut];
    }
  return penalty;
}
inline void pricing_initial_cut_state(TaxonId a, TaxonId b, const PricingData& d, BitSet& unconflicted, BitSet& inactive, W& penalty) {
  int n = d.ntaxa;
  unconflicted = BitSet();
  for (int v = 1; v <= d.variable_count; ++v) unconflicted.insert(v);
  unconflicted.setdiff_with(d.pair_conflicting_variables[pricing_pair_index(n, a, b)]);
  inactive = BitSet();
  for (int c = 0; c < (int)d.required_variables.size(); ++c) inactive.insert(c);
  penalty = pricing_activate_completed_cuts(inactive, unconflicted, d, 0.0);
}
inline void pricing_updated_cut_state(const PricingSearchState& state, TaxonId new_taxon, const PricingData& d,
                                      BitSet& unconflicted, BitSet& inactive, W& penalty) {
  int n = d.ntaxa;
  unconflicted = state.unconflicted_variables;
  unconflicted.setdiff_with(d.pair_conflicting_variables[pricing_pair_index(n, state.initial_a, new_taxon)]);
  inactive = state.inactive_cuts;
  penalty = pricing_activate_completed_cuts(inactive, unconflicted, d, state.cut_penalty);
}

// _pricing_rebuild_lookup (PR:517): taxon -> 1-based slot index (0 if not addable).
inline std::vector<int> pricing_rebuild_lookup(const std::vector<PricingInsertionPoint>& points, int ntaxa) {
  std::vector<int> lookup(ntaxa + 1, 0);
  for (int pi = 1; pi <= (int)points.size(); ++pi)
    for (TaxonId taxon : points[pi - 1].taxa) {
      if (lookup[taxon] != 0) throw std::runtime_error("insertion-point taxa must be disjoint");
      lookup[taxon] = pi;
    }
  return lookup;
}

// _pricing_initial_state (PR:586): seed the DFS for pair (a,b).
inline PricingSearchState pricing_initial_state(TaxonId a, TaxonId b, W bound, const PricingData& d) {
  int n = d.ntaxa, T = (int)d.trees.size();
  PricingSearchState s;
  s.initial_a = a; s.initial_b = b;
  const std::vector<TaxonId>& side_a = d.initial_side_taxa[pricing_pair_index(n, a, b)];
  const std::vector<TaxonId>& side_b = d.initial_side_taxa[pricing_pair_index(n, b, a)];
  if (!side_a.empty()) s.insertion_points.push_back({a, side_a});
  if (!side_b.empty()) s.insertion_points.push_back({b, side_b});
  W internal = 0.0;
  for (int i = 0; i < T; ++i) {
    const PhyloTree& tree = d.trees[i];
    NodeId root = d.lca_tables[i][a][b];
    internal += tree.internal_weight[root];
    internal += strict_connector_weight(tree, d.internal_prefixes[i], root, a);
    internal += strict_connector_weight(tree, d.internal_prefixes[i], root, b);
  }
  s.selected_taxa = {a, b}; std::sort(s.selected_taxa.begin(), s.selected_taxa.end());
  s.insertion_point_of = pricing_rebuild_lookup(s.insertion_points, n);
  s.current_internal_weight = internal;
  s.mast_upper_bound = bound;
  pricing_initial_cut_state(a, b, d, s.unconflicted_variables, s.inactive_cuts, s.cut_penalty);
  return s;
}

// _pricing_component_label (PR:616): 1/2/3 = candidate is below the new taxon's side / the
// representative's side / neither (i.e. on the connector itself), consistent across trees.
inline uint8_t pricing_component_label(const PricingData& d, int ti, NodeId connection, TaxonId candidate,
                                       TaxonId new_taxon, TaxonId descendant) {
  const PhyloTree& tree = d.trees[ti];
  const std::vector<int>& first = d.preorder_first[ti];
  const std::vector<int>& last = d.preorder_last[ti];
  NodeId new_child = pricing_child_containing(tree, first, last, connection, new_taxon);
  NodeId lower_child = pricing_child_containing(tree, first, last, connection, descendant);
  if (new_child == lower_child) throw std::runtime_error("insertion-point representative does not define two sides");
  if (pricing_is_descendant(first, last, candidate, new_child)) return 1;
  if (pricing_is_descendant(first, last, candidate, lower_child)) return 2;
  return 3;
}

// _pricing_child_state (PR:690): add `new_taxon`; re-derive addable frontier (rank-filtered for
// symmetry breaking); early-prune on the constant pair bound + cut penalty.
inline std::optional<PricingSearchState>
pricing_child_state(const PricingSearchState& state, TaxonId new_taxon, W incumbent_objective, const PricingData& d, W atol) {
  int n = d.ntaxa, T = (int)d.trees.size();
  BitSet unconflicted, inactive; W cut_penalty;
  pricing_updated_cut_state(state, new_taxon, d, unconflicted, inactive, cut_penalty);
  if (state.mast_upper_bound + cut_penalty <= incumbent_objective + atol) return std::nullopt;

  int old_idx = state.insertion_point_of[new_taxon];
  if (old_idx == 0) throw std::runtime_error("selected taxon is not currently addable");
  TaxonId descendant = state.insertion_points[old_idx - 1].representative_descendant;
  std::vector<PricingInsertionPoint> points;
  std::vector<TaxonId> old_taxa;
  for (int idx = 1; idx <= (int)state.insertion_points.size(); ++idx) {
    const PricingInsertionPoint& point = state.insertion_points[idx - 1];
    std::vector<TaxonId> later;
    for (TaxonId taxon : point.taxa) if (d.rank[taxon] > d.rank[new_taxon]) later.push_back(taxon);
    if (idx == old_idx) old_taxa = std::move(later);
    else if (!later.empty()) points.push_back({point.representative_descendant, std::move(later)});
  }
  std::vector<NodeId> connections(T);
  for (int i = 0; i < T; ++i) connections[i] = d.lca_tables[i][new_taxon][descendant];
  std::array<std::vector<TaxonId>, 3> groups;
  for (TaxonId candidate : old_taxa) {
    uint8_t label = pricing_component_label(d, 0, connections[0], candidate, new_taxon, descendant);
    bool consistent = true;
    for (int i = 1; i < T; ++i)
      if (pricing_component_label(d, i, connections[i], candidate, new_taxon, descendant) != label) { consistent = false; break; }
    if (consistent) groups[label - 1].push_back(candidate);
  }
  if (!groups[0].empty()) points.push_back({new_taxon, std::move(groups[0])});
  if (!groups[1].empty()) points.push_back({descendant, std::move(groups[1])});
  if (!groups[2].empty()) points.push_back({descendant, std::move(groups[2])});

  std::vector<TaxonId> selected = state.selected_taxa;
  selected.push_back(new_taxon); std::sort(selected.begin(), selected.end());
  W internal = state.current_internal_weight;
  for (int i = 0; i < T; ++i)
    internal += strict_connector_weight(d.trees[i], d.internal_prefixes[i], connections[i], new_taxon);

  PricingSearchState child;
  child.initial_a = state.initial_a; child.initial_b = state.initial_b;
  child.insertion_point_of = pricing_rebuild_lookup(points, n);   // before moving points
  child.insertion_points = std::move(points);
  child.selected_taxa = std::move(selected);
  child.current_internal_weight = internal;
  child.mast_upper_bound = state.mast_upper_bound;                // :none mode => parent bound
  child.unconflicted_variables = std::move(unconflicted);
  child.inactive_cuts = std::move(inactive);
  child.cut_penalty = cut_penalty;
  return child;
}

inline std::vector<TaxonId> pricing_ordered_addable(const PricingSearchState& state) {   // PR:774
  std::vector<TaxonId> taxa;
  for (const auto& point : state.insertion_points) taxa.insert(taxa.end(), point.taxa.begin(), point.taxa.end());
  std::sort(taxa.begin(), taxa.end());
  return taxa;
}
// objective: (|S| + internal) + cut_penalty.
inline W pricing_objective(const PricingSearchState& s) {
  return ((W)s.selected_taxa.size() + s.current_internal_weight) + s.cut_penalty;
}
// solution (PR:844): penalty = internal + cut_penalty; objective = |S| + penalty  (NOTE: different
// FP association than pricing_objective — replicated exactly).
inline MASTPricingSolution pricing_solution(const PricingSearchState& s) {
  W penalty = s.current_internal_weight + s.cut_penalty;
  MASTPricingSolution out;
  out.objective = (W)s.selected_taxa.size() + penalty;
  out.taxa = s.selected_taxa;
  out.base_size = (int)s.selected_taxa.size();
  out.penalty = penalty;
  return out;
}

// _pricing_unique_solutions (PR:854): first-seen dedup by key (objective, penalty, sorted taxa).
inline std::vector<MASTPricingSolution> pricing_unique_solutions(const std::vector<MASTPricingSolution>& solutions) {
  std::vector<MASTPricingSolution> result;
  std::set<std::tuple<W, W, std::vector<TaxonId>>> seen;
  for (const auto& s : solutions) {
    auto key = std::make_tuple(s.objective, s.penalty, s.taxa);
    if (seen.count(key)) continue;
    seen.insert(key);
    result.push_back(s);
  }
  return result;
}

// search (PR:901), specialized to mode=:rooted_best (no early-stop): DFS that updates `incumbent`.
inline void pricing_search(const PricingSearchState& state, PricingIncumbent& incumbent, const PricingData& d, W atol) {
  if (state.mast_upper_bound + state.cut_penalty <= incumbent.solution.objective + atol) return;
  W objective = pricing_objective(state);
  if (objective > incumbent.solution.objective + atol) incumbent.solution = pricing_solution(state);
  std::vector<TaxonId> candidates = pricing_ordered_addable(state);
  for (TaxonId taxon : candidates) {
    std::optional<PricingSearchState> child = pricing_child_state(state, taxon, incumbent.solution.objective, d, atol);
    if (child) pricing_search(*child, incumbent, d, atol);
  }
}

// _pricing_solve (PR:866), mode=:rooted_best, state_bound=:none, debug/stats off.
inline MASTPricingResult pricing_solve(std::vector<PhyloTree> trees, const std::vector<std::vector<TaxonId>>& variables_in,
                                       const std::vector<std::vector<int>>& cuts_in, const std::vector<W>& cut_weights_in, W atol = 1e-8) {
  pricing_clamp_internal_weights(trees, atol);              // _pricing_validate_trees
  int n = trees[0].ntaxa;
  std::vector<std::vector<TaxonId>> variables;             // _pricing_normalize_variables
  for (const auto& var : variables_in) {
    std::vector<TaxonId> t = normalized_taxa(var, n);
    if (t.size() < 2) throw std::runtime_error("pricing variables need >=2 taxa");
    variables.push_back(std::move(t));
  }
  std::vector<std::vector<int>> cuts; std::vector<W> cut_weights;
  pricing_normalize_cuts(cuts_in, cut_weights_in, (int)variables.size(), atol, cuts, cut_weights);

  PricingData d = pricing_precompute(std::move(trees), variables, cuts, cut_weights);
  MASTPricingSolution singleton{1.0, {1}, 1, 0.0};
  std::vector<MASTPricingSolution> pair_incumbents(n * n + 1, singleton);   // 1-based by pair index

  for (const auto& pr : d.canonical_pairs) {
    TaxonId a = pr.first, b = pr.second;
    int pair_idx = pricing_pair_index(n, a, b);
    PricingIncumbent incumbent{pair_incumbents[pair_idx]};
    W bound = d.weighted_pair_bounds[pair_idx];
    if (bound <= incumbent.solution.objective + atol) continue;
    PricingSearchState state = pricing_initial_state(a, b, bound, d);
    pricing_search(state, incumbent, d, atol);
    pair_incumbents[pair_idx] = incumbent.solution;
  }

  std::vector<MASTPricingSolution> pair_solutions;
  for (const auto& pr : d.canonical_pairs) {
    const MASTPricingSolution& s = pair_incumbents[pricing_pair_index(n, pr.first, pr.second)];
    if (s.objective > 1.0 + atol) pair_solutions.push_back(s);
  }
  MASTPricingResult res;
  res.solutions = pricing_unique_solutions(pair_solutions);
  if (res.solutions.empty()) { res.solutions.push_back(singleton); res.objective = 1.0; }
  else {
    res.objective = res.solutions[0].objective;
    for (size_t i = 1; i < res.solutions.size(); ++i) if (res.solutions[i].objective > res.objective + atol) res.objective = res.solutions[i].objective;
  }
  res.improving = res.objective > 1.0 + atol;
  return res;
}

// Top-level production entry: trees already built (PhyloTrees with internal_weight = beta duals),
// cut_sets = per-cut list of taxa-sets (the clique members), cut_weights <= 0.  Returns the per-pair
// best solutions in lex pair order (== ctx_price!'s consume order).
inline std::vector<MASTPricingSolution>
pricing_mast_solve(std::vector<PhyloTree> trees, const std::vector<std::vector<std::vector<TaxonId>>>& cut_sets,
                   const std::vector<W>& cut_weights, W atol = 1e-8) {
  int ntaxa = trees[0].ntaxa;
  std::vector<std::vector<TaxonId>> variables; std::vector<std::vector<int>> cuts;
  pricing_flatten_cut_sets(cut_sets, ntaxa, variables, cuts);
  return pricing_solve(std::move(trees), variables, cuts, cut_weights, atol).solutions;
}

} // namespace maf::mast
