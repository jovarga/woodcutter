// Hand-written greedy clique separator over the set-packing conflict graph ("glpk clique cuts" —
// there is NO GLPK library).
//
// All separator state is indexed by COLUMN INDEX j into the caller's `vars` vector (0-based here).
// Relative order is what the tie-breaks use, and the caller fixes the order.
// Determinism notes (CUTS_SPEC §5): every sort here uses a strict total order (descending weight,
// ties by lexicographically-smaller support); _bitset_weight folds in ASCENDING column order; dedup
// is by exact set value (std::set<BitSet>).
#pragma once
#include "../core.hpp"
#include <vector>
#include <set>
#include <utility>
#include <algorithm>
#include <cmath>
#include <limits>

namespace maf {

// _bitset_weight (SP:2230): sum lp/weights over a set, folded in ASCENDING element order.
inline double _bitset_weight(const BitSet& s, const std::vector<double>& w) {
  double total = 0.0;
  for (int j : s) total += w[j];      // BitSet iterates ascending
  return total;
}

// _lexicographically_smaller (SP:1195): compare ascending element lists positionally; on a common
// prefix the SHORTER set is smaller.  (NB: this is NOT BitSet::operator<, which is set-membership
// order — they disagree, e.g. {1,3} vs {1,2}.)
inline bool _lex_smaller(const BitSet& a, const BitSet& b) {
  std::vector<int> va = a.to_vec(), vb = b.to_vec();
  size_t m = std::min(va.size(), vb.size());
  for (size_t i = 0; i < m; ++i) if (va[i] != vb[i]) return va[i] < vb[i];
  return va.size() < vb.size();
}

// _keep_best_clique_supports! (SP:1212): sort descending by weight, ties by lex-smaller support,
// then truncate to max_cuts.  Total order ⇒ std::sort (stability irrelevant).
inline void _keep_best_clique_supports(std::vector<std::pair<double, BitSet>>& found, int max_cuts) {
  std::sort(found.begin(), found.end(), [](const std::pair<double, BitSet>& a, const std::pair<double, BitSet>& b) {
    if (a.first == b.first) return _lex_smaller(a.second, b.second);
    return a.first > b.first;
  });
  if ((int)found.size() > max_cuts) found.resize(max_cuts);
}

// _build_full_conflict_graph (SP:500): neighbors[v] = all columns sharing a (tree,vertex) row with
// v.  Rows of size < 2 are skipped; no self-loops.
inline std::vector<BitSet> _build_full_conflict_graph(int nvars, const std::vector<BitSet>& cliques) {
  std::vector<BitSet> neighbors(nvars);
  for (const BitSet& row : cliques) {
    if (row.size() < 2) continue;
    for (int v : row) neighbors[v].union_with(row);
  }
  for (int v = 0; v < nvars; ++v) neighbors[v].erase(v);
  return neighbors;
}

// _best_glpk_clique_candidate (SP:1144): pick the candidate maximizing weight + intra-candidate
// neighbor weight; ties go to the smaller index.  Returns -1 if `candidates` is empty.
inline int _best_glpk_clique_candidate(const BitSet& candidates, const std::vector<BitSet>& neighbors,
                                       const std::vector<double>& weights) {
  int best = -1; double best_score = -std::numeric_limits<double>::infinity();
  for (int c : candidates) {                                   // ascending
    double score = weights[c] + _bitset_weight(set_intersect(candidates, neighbors[c]), weights);
    if (score > best_score || (score == best_score && (best == -1 || c < best))) { best = c; best_score = score; }
  }
  return best;
}

// _glpk_greedy_clique_candidates (SP:1090): per-seed greedy clique growth with relative bound prune
// and `marked` absorption.  Returns up to max_candidates deduped (weight, clique) pairs.
inline std::vector<std::pair<double, BitSet>>
_glpk_greedy_clique_candidates(const std::vector<double>& weights, const std::vector<BitSet>& neighbors,
                               const std::vector<double>& cumulative, const BitSet& promising,
                               double bound_prune_rel_tol, int max_candidates, bool use_bound_prune) {
  std::vector<int> seeds = promising.to_vec();                 // ascending
  std::sort(seeds.begin(), seeds.end(), [&](int a, int b) {
    if (cumulative[a] == cumulative[b]) return a < b;
    return cumulative[a] > cumulative[b];                      // descending cumulative
  });
  BitSet marked;
  std::vector<std::pair<double, BitSet>> out;
  std::set<BitSet> seen;
  double best_weight = 0.0;
  for (int seed : seeds) {
    if (marked.contains(seed)) continue;
    BitSet clique; clique.insert(seed);
    double clique_weight = weights[seed];
    BitSet candidates = set_intersect(promising, neighbors[seed]);
    double candidate_weight = _bitset_weight(candidates, weights);
    while (!candidates.empty()) {
      if (use_bound_prune) {
        double prune_gap = bound_prune_rel_tol * (1.0 + std::abs(best_weight));
        if (clique_weight + candidate_weight < best_weight + prune_gap) break;
      }
      int next_vertex = _best_glpk_clique_candidate(candidates, neighbors, weights);
      if (next_vertex == -1) break;
      clique.insert(next_vertex); clique_weight += weights[next_vertex];
      candidates = set_intersect(candidates, neighbors[next_vertex]);
      candidate_weight = _bitset_weight(candidates, weights);
    }
    marked.union_with(clique);                                 // absorb the whole clique
    if (clique.size() >= 2) {
      best_weight = std::max(best_weight, clique_weight);
      if (seen.insert(clique).second) {
        out.emplace_back(clique_weight, clique);
        if ((int)out.size() > max_candidates) _keep_best_clique_supports(out, max_candidates);
      }
    }
  }
  _keep_best_clique_supports(out, std::min((int)out.size(), max_candidates));
  return out;
}

// _maximalize_clique_glpk_order (SP:1018): extend to a maximal clique in the FULL conflict graph by
// repeatedly adding the smallest-id common neighbor (allowed = all variables).
inline BitSet _maximalize_clique_glpk_order(const BitSet& clique, const std::vector<BitSet>& neighbors) {
  BitSet maximal = clique;
  BitSet candidates;
  for (int v = 0; v < (int)neighbors.size(); ++v) candidates.insert(v);   // allowed = 1:nvars
  for (int v : clique) candidates.intersect_with(neighbors[v]);
  candidates.setdiff_with(maximal);
  while (!candidates.empty()) {
    int v = candidates.first();                                // smallest id
    maximal.insert(v);
    candidates.intersect_with(neighbors[v]);
    candidates.setdiff_with(maximal);
  }
  return maximal;
}

// _separate_glpk_clique_supports (SP:921): the engine (selection=:best, maximalize=:full_glpk_order).
inline std::vector<BitSet>
_separate_glpk_clique_supports(const std::vector<double>& lp, const std::vector<BitSet>& full_neighbors,
                               int max_cuts, int candidate_multiplier, double active_tol,
                               double promising_value_tol, double min_cumulative_weight,
                               double min_cut_weight, double bound_prune_rel_tol) {
  int nvars = (int)lp.size();
  // (A) active restriction (lp[j] > active_tol), build local <-> orig index maps.
  std::vector<int> local_to_orig; std::vector<int> orig_to_local(nvars, -1);
  for (int j = 0; j < nvars; ++j)
    if (lp[j] > active_tol) { orig_to_local[j] = (int)local_to_orig.size(); local_to_orig.push_back(j); }
  int nactive = (int)local_to_orig.size();
  if (nactive < 2) return {};
  // (B) local weights + local conflict adjacency.
  std::vector<double> weights(nactive);
  for (int i = 0; i < nactive; ++i) weights[i] = lp[local_to_orig[i]];
  std::vector<BitSet> active_neighbors(nactive);
  for (int li = 0; li < nactive; ++li)
    for (int no : full_neighbors[local_to_orig[li]]) { int nl = orig_to_local[no]; if (nl != -1) active_neighbors[li].insert(nl); }
  // (C) cumulative weight + promising seed set.
  std::vector<double> cumulative(nactive, 0.0); BitSet promising;
  for (int li = 0; li < nactive; ++li) {
    double total = weights[li] + _bitset_weight(active_neighbors[li], weights);
    cumulative[li] = total;
    if (weights[li] >= promising_value_tol && total >= min_cumulative_weight) promising.insert(li);
  }
  if (promising.size() < 2) return {};
  // (D) greedy candidate cliques (use_bound_prune = true for selection=:best).
  int candidate_limit = std::max(max_cuts, candidate_multiplier * max_cuts);
  auto local_candidates = _glpk_greedy_clique_candidates(weights, active_neighbors, cumulative,
                            promising, bound_prune_rel_tol, candidate_limit, /*use_bound_prune=*/true);
  // (E) map to orig ids, maximalize in full graph, dedup, weight-filter.
  std::set<BitSet> seen; std::vector<std::pair<double, BitSet>> found;
  for (auto& [weight, local_clique] : local_candidates) {
    if (local_clique.size() < 2) continue;
    if (weight < min_cut_weight) continue;
    BitSet original_clique;
    for (int i : local_clique) original_clique.insert(local_to_orig[i]);
    BitSet support = _maximalize_clique_glpk_order(original_clique, full_neighbors);
    if (!seen.insert(support).second) continue;
    double activity = _bitset_weight(support, lp);
    if (activity < min_cut_weight) continue;
    found.emplace_back(activity, support);
    if ((int)found.size() > candidate_limit) _keep_best_clique_supports(found, candidate_limit);
  }
  // (F) final top-max_cuts (selection=:best).
  _keep_best_clique_supports(found, max_cuts);
  std::vector<BitSet> supports;
  for (auto& fs : found) supports.push_back(fs.second);
  return supports;
}

// separate_glpk_clique_cuts (SP:240, production binding: selection=:best, seed=:single_vertex):
// build the conflict graph from the set-packing rows, separate, then drop supports set-equal to an
// existing row (skip_existing_rows=true).  Returns the cliques as column-index BitSets.
inline std::vector<BitSet>
separate_glpk_clique_cuts(const std::vector<double>& lp, const std::vector<BitSet>& cliques, int max_cuts) {
  int nvars = (int)lp.size();
  std::vector<BitSet> full_neighbors = _build_full_conflict_graph(nvars, cliques);
  std::vector<BitSet> supports = _separate_glpk_clique_supports(
      lp, full_neighbors, max_cuts, /*candidate_multiplier=*/4, /*active_tol=*/1e-9,
      /*promising_value_tol=*/1e-3, /*min_cumulative_weight=*/1.010, /*min_cut_weight=*/1.07,
      /*bound_prune_rel_tol=*/1e-5);
  std::vector<BitSet> kept;
  for (const BitSet& support : supports) {
    bool exists = false;
    for (const BitSet& row : cliques) if (row == support) { exists = true; break; }
    if (!exists) kept.push_back(support);
  }
  return kept;
}

} // namespace maf
