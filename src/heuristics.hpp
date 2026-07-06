// heuristics.hpp — HEURISTIC-track algorithms.
//
// No entry point here: this header is #included by harness_heuristic.hpp (only in heuristic-track
// builds, i.e. when -DMAF_HEURISTIC is set), which provides the competition time/SIGTERM harness.
// It ships the most successful heuristics:
//   * greedy_disjoint_partition / weighted_repair — disjoint set-packing (valid by construction),
//   * feasibility_pump — LP-guided round -> repair -> pump matheuristic (the workhorse),
//   * the FP + MILP (fpmip) finish — a binary set-partition MILP over the enriched pool,
// orchestrated by run_heuristic(). The column pool + master LP/MILP are reused from the exact
// solver (BCPContext root column generation, MPModel, solve_rmp).
#pragma once

// If an instance (or, since run_heuristic is the per-block leaf, a reduced/cluster block) has at most
// this many taxa, run_heuristic mimics the EXACT solver (proven optimal when it finishes in budget).
#define MAX_TAXA_COUNT__EXACT_APPROACH 2000 //0
// also further apporaches towards a real portfolio solver could be combined, e.g., for very large
// clusters, it would be possible to just recycle the solution of one cluster part, indead of solving
// both cluster versions; the value of MAX_TAXA_COUNT__EXACT_APPROACH needs to be tuned

#include "solve.hpp"       // BCPContext, ctx_process_root, solve_rmp, MPModel, solve_maf (exact entry)
#include "io.hpp"          // write_forest
#include "budget.hpp"      // remaining_seconds, set_time_limit_seconds, g_deadline, Interrupted
#include "maf_util.hpp"    // inner_vs, is_partition_of_leaves, has_disjoint_steiner_vertices
#include "heur_params.hpp" // heur_params() — env-overridable tuning knobs (defaults = the old literals)
#include <vector>
#include <algorithm>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <cmath>
#include <cstdint>
#include <map>
#include <iostream>

// Heuristic 2-tree MAST engine. By default the heuristic uses the O(n log n) centroid engine
// (recovered, heuristic-only); the CMake flag -DMAF_HEUR_MAST=quad defines MAF_HEUR_QUAD_MAST to
// override it back to the exact track's O(n^2) quad_dp. The centroid header is compiled ONLY when
// not overridden.
#if !defined(MAF_HEUR_QUAD_MAST)
#include "mast/_centroid.hpp"
#endif

namespace maf { namespace heur {

using Sol = std::vector<BitSet>;

inline Sol singletons(const MAFP& maf) {
  Sol s; s.reserve(maf.n);
  for (int x = 1; x <= maf.n; ++x) s.push_back(BitSet{x});
  return s;
}
// AHU canonical id of a rooted tree's root.  A leaf's id is its label (leaves keep the same labels
// across all trees of a restricted sub-instance); an internal node's id comes from a dictionary keyed
// by its children's sorted ids.  Sharing dict/next_id across the trees gives leaf-labelED isomorphism:
// two trees are isomorphic (same induced topology) iff their root ids are equal.  Iterative post-order.
inline long long ahu_root(const DiGraph& t, int nleaves,
                          std::map<std::vector<long long>, long long>& dict, long long& next_id) {
  const int N = t.nv();
  int root = N; for (int v = 1; v <= N; ++v) if (t.in_degree(v) == 0) { root = v; break; }
  std::vector<long long> cid(N + 1, -1);
  std::vector<int> st{root}; std::vector<char> down(N + 1, 0);
  while (!st.empty()) {
    int v = st.back();
    const std::vector<int>& ch = t.out_neighbors(v);
    if (!down[v]) { down[v] = 1; for (int c : ch) st.push_back(c); continue; }
    st.pop_back();
    if (ch.empty()) { cid[v] = v; continue; }                 // leaf: canonical id == label (<= nleaves)
    std::vector<long long> key; key.reserve(ch.size());
    for (int c : ch) key.push_back(cid[c]);
    std::sort(key.begin(), key.end());
    auto it = dict.find(key);
    cid[v] = (it != dict.end()) ? it->second : (dict[key] = (long long)nleaves + (next_id++));
  }
  return cid[root];
}
inline bool all_trees_agree(const MAFP& sub) {                // all induced trees isomorphic (leaf-labeled)?
  std::map<std::vector<long long>, long long> dict; long long next_id = 0;
  long long r0 = ahu_root(sub.T[0], sub.n, dict, next_id);
  for (size_t i = 1; i < sub.T.size(); ++i)
    if (ahu_root(sub.T[i], sub.n, dict, next_id) != r0) return false;
  return true;
}
// Is `comp` an agreement subtree (same induced topology in every tree)?  1-/2-leaf sets trivially agree.
// Small components use the O(k^3) rooted-triple test (early-exits, never scans the whole tree); large
// (and rare) ones build the induced sub-instance once (O(t*n)) and test tree isomorphism via AHU — so a
// big VALID agreement component costs O(t*n), not O(k^3).
inline bool is_agreement_component(const MAFP& maf, const BitSet& comp) {
  const int k = (int)comp.size();
  if (k <= 2)   return true;
  if (k <= 128) return is_valid_component(maf, comp);
  std::vector<int> leaves = comp.to_vec();                    // ascending
  try { auto sr = restrict_to_leaves(maf, leaves); return all_trees_agree(sr.first); }
  catch (const std::exception&) { return is_valid_component(maf, comp); }
}

// Full validity: a leaf partition, disjoint Steiner trees in every input tree, AND every component an
// agreement subtree (induced rooted topology identical across all trees).  The topology check was once
// skipped ("pool columns are agreement subtrees by construction"), but NOT every incumbent comes from
// the priced pool — the cluster-reduction reconstruction can emit a component whose topology disagrees
// between the trees.  Such a forest passes partition + Steiner-disjointness yet is REJECTED by the PACE
// checker (Wrong Answer = DISQUALIFICATION).  guard() falls any invalid candidate back to singletons,
// so an infeasible forest can NEVER be published.
inline bool cheap_valid(const MAFP& maf, const Sol& sol) {
  if (!is_partition_of_leaves(sol, maf.n)) return false;
  for (const auto& t : maf.T) if (!has_disjoint_steiner_vertices(t, sol)) return false;
  for (const BitSet& comp : sol) if (!is_agreement_component(maf, comp)) return false;
  return true;
}
inline Sol guard(const MAFP& maf, const Sol& sol) {
  return cheap_valid(maf, sol) ? sol : singletons(maf);
}

// Precompute inner_vs(T[i], col) for the whole pool once (reused across every repair).
inline std::vector<std::vector<BitSet>> pool_inners(const MAFP& maf, const std::vector<BitSet>& pool) {
  std::vector<std::vector<BitSet>> inn(pool.size(), std::vector<BitSet>(maf.T.size()));
  for (size_t k = 0; k < pool.size(); ++k)
    for (size_t i = 0; i < maf.T.size(); ++i) {
      auto iv = inner_vs(maf.T[i], pool[k]);
      inn[k][i] = BitSet(iv.begin(), iv.end());
    }
  return inn;
}

// Weight-ordered disjoint repair (_repair / greedy_disjoint_partition): accept pool columns by
// descending (weight, size) while leaves AND per-tree Steiner vertices stay disjoint; fill uncovered
// leaves with singletons. Valid by construction.
inline Sol weighted_repair(const MAFP& maf, const std::vector<BitSet>& pool,
                           const std::vector<std::vector<BitSet>>& inn,
                           const std::vector<double>& w) {
  std::vector<int> order(pool.size());
  for (size_t k = 0; k < pool.size(); ++k) order[k] = (int)k;
  std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
    if (w[a] != w[b]) return w[a] > w[b];
    return pool[a].size() > pool[b].size();
  });
  std::vector<BitSet> used(maf.T.size());
  BitSet covered;
  Sol chosen;
  for (int k : order) {
    const BitSet& col = pool[k];
    if (col.size() < 2) continue;
    if (!isdisjoint(col, covered)) continue;
    bool ok = true;
    for (size_t i = 0; i < maf.T.size(); ++i)
      if (!isdisjoint(used[i], inn[k][i])) { ok = false; break; }
    if (!ok) continue;
    chosen.push_back(col); covered.union_with(col);
    for (size_t i = 0; i < maf.T.size(); ++i) used[i].union_with(inn[k][i]);
  }
  for (int x = 1; x <= maf.n; ++x) if (!covered.contains(x)) chosen.push_back(BitSet{x});
  return chosen;
}

inline Sol greedy_disjoint_partition(const MAFP& maf, const std::vector<BitSet>& pool,
                                     const std::vector<std::vector<BitSet>>& inn) {
  std::vector<double> w(pool.size());
  for (size_t k = 0; k < pool.size(); ++k) w[k] = (double)pool[k].size();
  return weighted_repair(maf, pool, inn, w);
}

// LP-guided FEASIBILITY PUMP (feasibility_pump): build the pool LP, then iterate
//   round (LP fractional values) -> repair to a valid forest -> pump (penalise the columns just
//   used, re-solve the LP for a different fractional point) -> repeat, keeping the best valid forest.
// Deadline-bounded; reports improvements through `consider`.
inline Sol feasibility_pump(const MAFP& maf, const std::vector<BitSet>& pool,
                            const std::vector<std::vector<BitSet>>& inn,
                            const std::function<void(const Sol&)>& consider,
                            int stall = 25, double pump_penalty = 1.0,
                            double size_weight = 0.0, double jitter = 0.4, uint64_t seed = 1) {
  if (pool.empty()) return singletons(maf);
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> U(-0.5, 0.5);
  std::unordered_map<BitSet, int> idx;
  for (size_t k = 0; k < pool.size(); ++k) idx.emplace(pool[k], (int)k);

  // size baseline (always available even if the LP never solves)
  std::vector<double> sizew(pool.size());
  for (size_t k = 0; k < pool.size(); ++k) sizew[k] = (double)pool[k].size();
  Sol best = guard(maf, weighted_repair(maf, pool, inn, sizew));
  consider(best);

  MPModel mp(maf, pool, /*binary=*/false);
  double lpobj;
  try { lpobj = mp.solve(); } catch (...) { return best; }
  int lb = (int)std::ceil(lpobj - 1e-6); if (lb < 1) lb = 1;

  std::vector<double> pcoef(pool.size(), 0.0);
  std::vector<double> w(pool.size());
  auto lpweight = [&](bool jit) {
    for (size_t k = 0; k < pool.size(); ++k)
      w[k] = mp.sol_value((int)k) + size_weight * (double)((int)pool[k].size() - 1)
           + (jit && jitter > 0 ? jitter * U(rng) : 0.0);
  };

  int it = 0, since = 0;
  while (remaining_seconds() > 0.3 && (int)best.size() > lb && since < stall) {
    ++it;
    lpweight((it % 2) != 0);                       // jitter on odd iterations
    Sol P = guard(maf, weighted_repair(maf, pool, inn, w));
    if (P.size() < best.size()) { best = P; consider(best); since = 0; }
    else                        { ++since; }

    // PUMP: penalise every (size>=2) column used in P, then re-solve the perturbed LP.
    for (const BitSet& c : P) {
      if (c.size() < 2) continue;
      auto it2 = idx.find(c);
      if (it2 == idx.end()) continue;
      int k = it2->second;
      pcoef[k] += pump_penalty;
      mp.hs->changeColCost((HighsInt)(mp.real_base + k), -(double)((int)pool[k].size() - 1) + pcoef[k]);
    }
    try { mp.solve(); } catch (...) { break; }
  }
  return best;
}

// Heuristic-track MAST dispatch: for 2 trees use the O(n log n) centroid engine by default; the
// exact track's O(n^2) quad_dp is used instead when -DMAF_HEUR_MAST=quad (MAF_HEUR_QUAD_MAST) is set.
// For t != 2 both fall through to the multi-tree DP (via the shared mast_solve facade).
inline std::vector<std::pair<W, BitSet>>
heur_mast_solve(const std::vector<mast::PhyloTree>& trees, int n,
                const std::vector<W>* leaf_w = nullptr) {      // optional per-leaf alpha (null => uniform)
#if defined(MAF_HEUR_QUAD_MAST)
  (void)leaf_w;
  return mast_solve(trees, n, {});                             // quadratic override (exact-track engine)
#else
  if (trees.size() == 2) {
    auto raw = mast::centroid_mast_solve(trees[0], trees[1], leaf_w);  // O(n log n) 2-tree centroid
    std::vector<std::pair<W, BitSet>> out; out.reserve(raw.size());
    for (auto& p : raw) out.emplace_back(p.first, BitSet(p.second.begin(), p.second.end()));
    return out;
  }
  return mast_solve(trees, n, {});                             // t != 2 -> multi-tree DP
#endif
}

// Time-efficient MAST-PEEL: repeatedly compute the
// unit-weight maximum agreement subtree of the remaining leaves (reusing the exact solver's pricer,
// mast_solve, with all internal weights = 1 so it maximises agreement-subtree SIZE), collect every
// returned agreement subtree as a candidate column, then peel the largest to make progress. This
// enriches the pool with large, structurally-agreeing components the dual-priced root CG may never
// generate. Bounded by max_cols AND a wall-clock cap (one MAST is O(n^2) for t=2, so large n must be
// capped).
// `leaf_weight` (per ORIGINAL leaf 1..n; null => uniform) is an optional per-leaf alpha threaded into
// the centroid MAST — it biases which agreement subtree wins, composing with the internal weights below.
inline std::vector<BitSet> mast_peel(const MAFP& maf, size_t max_cols, double time_cap,
                                     const std::vector<double>* leaf_weight = nullptr) {
  std::vector<BitSet> cols;
  std::unordered_set<BitSet> seen;
  std::vector<int> remaining(maf.n);
  for (int x = 0; x < maf.n; ++x) remaining[x] = x + 1;      // sorted 1..n
  const int T = (int)maf.T.size();
  const double start_left = remaining_seconds();

  while (remaining.size() >= 2 && cols.size() < max_cols
         && remaining_seconds() > 0.3 && (start_left - remaining_seconds()) < time_cap) {
    auto sub_pair = restrict_to_leaves(maf, remaining);      // (sub_maf, lmap): lmap[k]=orig of sub-leaf k+1
    const MAFP& sub = sub_pair.first;
    const std::vector<int>& lmap = sub_pair.second;

    std::vector<mast::PhyloTree> trees; trees.reserve(T);
    std::vector<W> ones(std::max(sub.n - 1, 0), (W)1.0);     // unit internal weights => size-MAST
    for (int i = 0; i < T; ++i) trees.push_back(mast::build_phylo_tree(sub.T[i], ones, sub.n));

    // Map the per-original-leaf alpha to THIS sub-instance (sub-leaf k+1 -> original lmap[k]).
    std::vector<W> subw;
    if (leaf_weight) { subw.resize(sub.n); for (int k = 0; k < sub.n; ++k) subw[k] = (W)(*leaf_weight)[lmap[k] - 1]; }
    std::vector<std::pair<W, BitSet>> priced = heur_mast_solve(trees, sub.n, leaf_weight ? &subw : nullptr);
    if (priced.empty()) break;

    BitSet peel; int peelsz = 0;                             // remember the largest to peel
    for (auto& pr : priced) {
      std::vector<int> ov;
      for (int x : pr.second.to_vec()) ov.push_back(lmap[x - 1]);   // sub label -> original
      BitSet orig(ov.begin(), ov.end());
      if ((int)orig.size() >= 2 && seen.insert(orig).second) cols.push_back(orig);
      if ((int)pr.second.size() > peelsz) { peelsz = (int)pr.second.size(); peel = orig; }
    }
    if (peelsz < 2) break;

    std::vector<int> nxt;
    for (int x : remaining) if (!peel.contains(x)) nxt.push_back(x);
    if (nxt.size() == remaining.size()) break;               // no progress
    remaining.swap(nxt);
  }
  return cols;
}

// Orchestrator: a fast MAST-peel + greedy pack gives an immediate non-trivial floor (so an early
// SIGTERM on large n never falls back to singletons), then the exact solver's root column generation
// enriches the pool, and finally a feasibility-pump restart loop spends the rest of the budget.
// Reports every improvement through `consider` (kept SIGTERM-safe by the harness).
//
// `single_pass`: when true, skip the randomized-restart loop (do one greedy + one feasibility pump).
// Used as the per-block leaf solver inside the reduction/cluster pipeline, where the OUTER restart
// loop owns diversification and each of the many small blocks must return quickly.
inline Sol run_heuristic(const MAFP& maf, const std::function<void(const Sol&)>& consider,
                         uint64_t seed = 1, bool single_pass = false, bool allow_exact = true) {
  Sol best = singletons(maf);
  auto take = [&](const Sol& s) { Sol g = guard(maf, s); if (g.size() < best.size()) { best = g; consider(best); } };

  // Preliminary: within the exact-approach size threshold, MIMIC THE EXACT SOLVER (proven optimal
  // when it finishes in budget). solve_maf is header-only inline (solve.hpp), so it is compiled even
  // in heuristic-track builds — no copy of the exact solver is needed. It honours the global deadline
  // (check_budget/time-capped LPs); on budget exhaustion it throws Interrupted, and we fall through to
  // the heuristic pipeline so a valid incumbent is always available.
  if (allow_exact && maf.n <= heur_params().max_taxa_exact) {
    try {
      auto [sol, dual] = solve_maf(maf);
      (void)dual;
      take(sol);
      return best;                       // exact (optimal) solution
    } catch (const Interrupted&) {        // out of budget mid-exact-solve -> heuristic fallback below
    } catch (const std::exception&) {}    // any exact-solver failure -> heuristic fallback below
  }

  std::vector<BitSet> pool;
  std::unordered_set<BitSet> poolset;
  std::vector<std::vector<BitSet>> inn;                       // per-column per-tree Steiner vertices
  auto add_col = [&](const BitSet& c) {
    if ((int)c.size() < 2 || !poolset.insert(c).second) return;
    pool.push_back(c);
    std::vector<BitSet> iv(maf.T.size());
    for (size_t i = 0; i < maf.T.size(); ++i) { auto v = inner_vs(maf.T[i], c); iv[i] = BitSet(v.begin(), v.end()); }
    inn.push_back(std::move(iv));
  };

  // (1) FAST FLOOR — time-bounded MAST-peel columns + greedy pack, published immediately so an early
  //     SIGTERM on large n yields a non-trivial forest (not singletons).
  if (remaining_seconds() > 1.0) {
    const auto& hp = heur_params();
    double cap = std::min(remaining_seconds() * hp.floor_cap_frac, hp.floor_cap_max);
    for (const BitSet& c : mast_peel(maf, /*max_cols=*/(size_t)hp.floor_max_cols, cap)) add_col(c);
    if (!pool.empty()) take(greedy_disjoint_partition(maf, pool, inn));
  }

  // (2) price-and-branch at the root: enrich the pool + a valid RMP forest (ctx.incumbent).
  //     Time-boxed so the pump restarts still get a share of the budget.
  auto saved_deadline = g_deadline; bool saved_has = g_has_deadline;
  double full = remaining_seconds();
  set_time_limit_seconds(std::max(std::min(full * heur_params().root_cg_frac, full), 0.5));
  BCPContext ctx(maf);
  try { ctx_process_root(ctx); } catch (const Interrupted&) {} catch (const std::exception&) {}
  g_deadline = saved_deadline; g_has_deadline = saved_has;    // restore the full deadline
  take(ctx.incumbent);                                        // MILP-over-pool forest (fpmip's MILP)
  size_t pool_after_peel = pool.size();
  for (const auto& c : ctx.mip_vars) add_col(c);
  size_t pool_after_mipvars = pool.size();

  // Experimental enrichment #1: the heuristic used to pass only ctx.mip_vars (columns that were
  // positive in an LP solution) to the FP pool.  The root master often contains many more valid
  // dual-priced columns.  They may be LP-zero but still useful for the final integer packing.
  const auto& hp = heur_params();
  if (hp.cg_all_cols && !ctx.nodes.empty() && ctx.nodes[0]) {
    std::vector<const BitSet*> candidates;
    candidates.reserve(ctx.nodes[0]->mp.col_sets.size());
    for (const BitSet& c : ctx.nodes[0]->mp.col_sets)
      if ((int)c.size() >= hp.cg_all_min_size) candidates.push_back(&c);
    if (hp.cg_all_sort_size)
      std::stable_sort(candidates.begin(), candidates.end(), [](const BitSet* a, const BitSet* b) {
        if (a->size() != b->size()) return a->size() > b->size();
        return a->to_vec() < b->to_vec();
      });
    int added_all = 0;
    for (const BitSet* cp : candidates) {
      if (hp.cg_all_cap > 0 && added_all >= hp.cg_all_cap) break;
      size_t before = pool.size();
      add_col(*cp);
      if (pool.size() != before) ++added_all;
    }
  }
  size_t pool_after_allcg = pool.size();

  // Experimental enrichment #2: do one extra pricing pass at the final root duals, but harvest
  // near-improving columns into the heuristic pool instead of adding them to the LP master.
  if (hp.near_rc_min_score > 0.0 && !ctx.nodes.empty() && ctx.nodes[0] && remaining_seconds() > 0.5) {
    try {
      std::vector<mast::PhyloTree> trees = build_pricing_trees(ctx, *ctx.nodes[0]);
      std::vector<std::pair<W, BitSet>> priced = mast_solve(trees, maf.n, {});
      std::stable_sort(priced.begin(), priced.end(), [](const auto& a, const auto& b) {
        if (a.first != b.first) return a.first > b.first;
        return a.second.size() > b.second.size();
      });
      int added_near = 0;
      for (const auto& p : priced) {
        if (p.first < hp.near_rc_min_score) break;
        if (hp.near_rc_cap > 0 && added_near >= hp.near_rc_cap) break;
        size_t before = pool.size();
        add_col(p.second);
        if (pool.size() != before) ++added_near;
      }
    } catch (...) {}
  }
  size_t pool_after_near = pool.size();

  // Experimental enrichment #3: ask nearby dual pricing problems for alternative tied/near-tied
  // subtrees.  This keeps the perturbation LP-guided: only internal weights from the final root duals
  // are jittered, and the resulting valid MAST columns are merely added to the heuristic pool.
  if (hp.dual_jitter_passes > 0 && !ctx.nodes.empty() && ctx.nodes[0] && remaining_seconds() > 0.5) {
    try {
      std::vector<mast::PhyloTree> base = build_pricing_trees(ctx, *ctx.nodes[0]);
      std::mt19937_64 jrng(seed ^ 0x9e3779b97f4a7c15ULL);
      std::uniform_real_distribution<double> J(-hp.dual_jitter_mag, hp.dual_jitter_mag);
      int added_jit = 0;
      for (int pass = 0; pass < hp.dual_jitter_passes && remaining_seconds() > 0.5; ++pass) {
        std::vector<mast::PhyloTree> trees = base;
        for (auto& tr : trees)
          for (int v = tr.ntaxa + 1; v <= tr.nv(); ++v) tr.internal_weight[v] += J(jrng);
        std::vector<std::pair<W, BitSet>> priced = mast_solve(trees, maf.n, {});
        std::stable_sort(priced.begin(), priced.end(), [](const auto& a, const auto& b) {
          if (a.first != b.first) return a.first > b.first;
          return a.second.size() > b.second.size();
        });
        for (const auto& p : priced) {
          if (p.first < hp.dual_jitter_min_score) break;
          if (hp.dual_jitter_cap > 0 && added_jit >= hp.dual_jitter_cap) break;
          size_t before = pool.size();
          add_col(p.second);
          if (pool.size() != before) ++added_jit;
        }
        if (hp.dual_jitter_cap > 0 && added_jit >= hp.dual_jitter_cap) break;
      }
    } catch (...) {}
  }
  if (hp.pool_log) {
    std::cerr << "MAF_HEUR_POOL n=" << maf.n
              << " peel=" << pool_after_peel
              << " mipvars=" << pool_after_mipvars
              << " allcg=" << pool_after_allcg
              << " near=" << pool_after_near
              << " jitter=" << pool.size()
              << " ctx_mip_vars=" << ctx.mip_vars.size()
              << " ctx_cols=" << ((!ctx.nodes.empty() && ctx.nodes[0]) ? ctx.nodes[0]->mp.col_sets.size() : 0)
              << "\n";
  }

  if (pool.empty()) return best;

  // (3) greedy disjoint pack over the full (peel + CG) pool.
  take(greedy_disjoint_partition(maf, pool, inn));

  // Single-pass leaf mode: one feasibility pump, no restart loop (the outer pipeline diversifies).
  if (single_pass) {
    const auto& hp = heur_params();
    if (remaining_seconds() > 0.3)
      try { take(feasibility_pump(maf, pool, inn, consider, hp.fp_stall, hp.fp_pump_penalty,
                                  hp.fp_size_weight, hp.fp_jitter, seed)); }
      catch (const Interrupted&) {} catch (const std::exception&) {}
    return best;
  }

  // (4) feasibility pump with RANDOMIZED RESTARTS — spend the WHOLE remaining budget diversifying.
  //     A single pump pass stalls out well before the 5-minute limit; each restart re-runs the pump
  //     over the same pool with a fresh seed and varied steering (LP-only vs size-aware weight, and
  //     different jitter), so different LP roundings are explored. The incumbent only ever improves
  //     (reported live via `consider`, so SIGTERM always flushes the best). Loops until the budget
  //     is (nearly) gone or an Interrupted is thrown.
  // Pool LP lower bound (#components over this pool): restarts cannot beat it, so once the incumbent
  // reaches it we stop rather than spin — it is proven best for this pool.
  int lb = 1;
  if (remaining_seconds() > 0.5) {
    try { MPModel lpm(maf, pool, /*binary=*/false); double o = lpm.solve();
          lb = std::max(1, (int)std::ceil(o - 1e-6)); } catch (...) {}
  }

  uint64_t s = seed;
  for (int restart = 0; remaining_seconds() > 0.5 && (int)best.size() > lb; ++restart) {
    double size_weight = (restart % 2 == 0) ? hp.fp_size_weight : hp.fp_size_weight_alt;  // pure LP vs size-aware
    double jitter      = hp.fp_jitter + hp.fp_jitter_step * (double)(restart % 3);        // vary tie-breaking magnitude
    try {
      take(feasibility_pump(maf, pool, inn, consider, hp.fp_stall, hp.fp_pump_penalty,
                            size_weight, jitter, /*seed=*/s));
    } catch (const Interrupted&) { break; } catch (const std::exception&) { break; }
    s = s * 6364136223846793005ULL + 1442695040888963407ULL; // advance the seed (LCG)
  }

  return best;
}

}} // namespace maf::heur
