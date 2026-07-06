// _lineartime_reductions_clusters.hpp — HEURISTIC-track orchestration (uses the linear-time infra).
//
// Wraps the heuristic solver (heuristics.hpp::run_heuristic) in kernelization + cluster
// decomposition preprocessing, ALWAYS applied on the heuristic track (no CMake switch):
//   * solve_with_reductions_linear — cherry + (4,3)-chain + 3-2 rules via a mutable-forest worklist
//                                    (Day-linear O(t*n), _lineartime_day.hpp; lossless),
//   * solve_clusters               — LS11 cluster decomposition (the QUALITY-PRESERVING shared one
//                                    from preprocessing.hpp; solves both cluster parts — option A
//                                    with rho, option B without rho — INDEPENDENTLY).
// The linear cluster_sequence_linear / solve_clusters_linear live in the Day header but are NOT used
// here: the Day decomposition is coarser and loses quality (measured 766->779 on exact100), so the
// shared solve_clusters drives the decomposition. All steps preserve validity for ANY valid block
// solution (decompose + relabel + recombine + production reconstruction), so a heuristic — possibly
// suboptimal — leaf stays feasible.
//
// SIGTERM safety: a fast top-level MAST-peel floor is published FIRST (so an early SIGTERM on large
// n yields a non-trivial forest, never singletons); only then does the reduce/cluster pipeline run,
// with each block solved by the single-pass heuristic leaf. An outer restart loop re-runs the whole
// reduced pipeline with fresh seeds until the budget is spent (each improvement reported live).
#pragma once
#include "_lineartime_day.hpp"             // maf::lineartime::solve_with_reductions_linear (Day-linear, O(t*n))
#include "preprocessing.hpp"               // solve_clusters (quality-preserving decomposition), ClusterSolver
#include "heuristics.hpp"                  // maf::heur::run_heuristic, mast_peel, greedy_disjoint_partition, guard
#include "budget.hpp"
#include <vector>
#include <functional>
#include <cstdint>
#include <cstdlib>

namespace maf { namespace heur {

// Run the heuristic behind the exact solver's reduction + cluster decomposition, once, with a given
// seed. The per-block leaf is the single-pass heuristic (no inner restart loop; the OUTER loop in
// run_heuristic_reduced owns diversification). `track` receives every reconstructed full-instance
// solution (validity-guarded) so the caller can keep the global incumbent + SIGTERM buffer current.
inline Sol run_reduced_once(const MAFP& maf, uint64_t seed, int defer_mode,
                            const std::function<void(const Sol&)>& progress = {}) {
  static const int CLUSTER_EXACT_MAX = [](){
    const char* e = std::getenv("MAF_CLUSTER_EXACT_MAX");
    return (e && *e) ? std::atoi(e) : -1;   // -1 preserves the historical run_heuristic leaf policy
  }();
  // Per-block leaf: single-pass heuristic returning (partition, dual = #components) so the
  // solve_clusters primal/dual guards (dual <= primal) always hold.  The leaf forwards every block-level
  // improvement to the cluster progress bridge (set per-part by solve_clusters_linear) for PROGRESSIVE
  // full-instance registration; when no progress consumer is wired the bridge is empty (a no-op).
  ClusterSolver leaf = [seed](const MAFP& block) -> std::pair<Sol, double> {
    Sol s;
    if (CLUSTER_EXACT_MAX >= 0 && block.n <= CLUSTER_EXACT_MAX) {
      try {
        auto [exact_sol, dual] = solve_maf(block);
        (void)dual;
        s = exact_sol;
      } catch (const Interrupted&) {
        throw;
      } catch (const std::exception&) {
        s = run_heuristic(block, [](const Sol& b){
              auto& br = lineartime::cluster_progress_bridge(); if (br) br(b);
            }, seed, /*single_pass=*/true, /*allow_exact=*/false);
      }
    } else {
      s = run_heuristic(block, [](const Sol& b){
            auto& br = lineartime::cluster_progress_bridge(); if (br) br(b);
          }, seed, /*single_pass=*/true, /*allow_exact=*/CLUSTER_EXACT_MAX < 0);
    }
    return {s, (double)s.size()};
  };

  // Inline of solve_with_reductions_linear so progressive registration threads through BOTH mappings:
  // the cluster composition (solve_clusters_linear's bridge) AND the reduction un-map (reconstruct).
  // The reduced-space consider un-maps each provisional forest back to ORIGINAL labels before `progress`
  // (which validates against the original instance).  reconstruct(state, .) is non-mutating (the exact
  // hook already calls it repeatedly), so re-invoking it per improvement is safe.
  lineartime::LinearReduction R = lineartime::reduce_linear(maf);
  SolverState state(R.reduced_maf);
  state.leaf_map = R.leaf_map; state.pending_actions = R.actions; state.k_adjustment = R.k_adjustment;

  std::function<void(const std::vector<BitSet>&)> reduced_consider;
  if (progress)
    reduced_consider = [&](const std::vector<BitSet>& reduced_sol) {
      try { progress(reconstruct(state, reduced_sol)); } catch (...) {}
    };

  std::vector<BitSet> sol_reduced;
  if      (R.reduced_maf.n == 0) sol_reduced = {};
  else if (R.reduced_maf.n == 1) sol_reduced = {BitSet{1}};
  else sol_reduced = lineartime::solve_clusters_linear(R.reduced_maf, leaf, defer_mode, reduced_consider).first;

  return reconstruct(state, sol_reduced);
}

// Top-level heuristic entry for the heuristic track: fast floor -> reduce/cluster pipeline with
// randomized restarts. Always applies the reductions + linear cluster decomposition.
inline Sol run_heuristic_reduced(const MAFP& maf, const std::function<void(const Sol&)>& consider,
                                 uint64_t seed = 1) {
  Sol best = singletons(maf);
  auto take = [&](const Sol& s) { Sol g = guard(maf, s); if (g.size() < best.size()) { best = g; consider(best); } };

  // Fast top-level floor (SIGTERM safety on large n): a MAST-peel + greedy pack over the WHOLE
  // instance, before the reduce/cluster pipeline runs.
  if (remaining_seconds() > 1.0) {
    double cap = std::min(remaining_seconds() * 0.2, 45.0);
    std::vector<BitSet> cols = mast_peel(maf, /*max_cols=*/8000, cap);
    if (!cols.empty()) {
      std::vector<std::vector<BitSet>> inn(cols.size(), std::vector<BitSet>(maf.T.size()));
      for (size_t k = 0; k < cols.size(); ++k)
        for (size_t i = 0; i < maf.T.size(); ++i) {
          auto v = inner_vs(maf.T[i], cols[k]); inn[k][i] = BitSet(v.begin(), v.end());
        }
      take(greedy_disjoint_partition(maf, cols, inn));
    }
  }

  // Adaptive EXACT phase (before the heuristic restart loop).  The branch-and-price's reach is governed
  // by the IRREDUCIBLE CORE, not raw n: measured finish-frontier is exact-finishes for reduced size
  // <~450 (in <15 s) and does-not-finish by ~1300.  So on instances whose LINEARLY-reduced size is
  // within reach, run the exact BCP with LIVE incumbent registration: the hook reconstructs EVERY
  // improving branch-and-price incumbent to the full instance and registers it.  If the BCP proves
  // optimality -> that is the best possible answer (return it).  Otherwise its partial progress is kept
  // (unlike the old code, which discarded a non-finishing exact solve down to the MAST-peel floor), and
  // we fall through to the heuristic restart loop starting from that strong incumbent.  Time-boxed so
  // the restart loop keeps a share on non-finishing instances.  Reduced size (reduction only, no
  // clustering) is the gate — decomposable instances (reduced still large) are left to the restart loop.
  // TODO(retune-300s): 3000 was tuned at a 120 s proxy budget; the crossover moves up at the real 300 s
  // budget, so this can likely be raised.  See docs/TODO.md.
  static const int EXACT_INCUMBENT_MAX = [](){ const char* e=std::getenv("MAF_EXACT_MAX"); return (e&&*e)?std::atoi(e):3000; }();
  if (remaining_seconds() > 2.0) {
    try {
      lineartime::LinearReduction R = lineartime::reduce_linear(maf);
      if (R.reduced_maf.n >= 2 && (int)R.reduced_maf.n <= EXACT_INCUMBENT_MAX) {
        SolverState st(R.reduced_maf);
        st.leaf_map = R.leaf_map; st.pending_actions = R.actions; st.k_adjustment = R.k_adjustment;
        auto saved_hook = g_incumbent_hook;
        g_incumbent_hook = [&](const std::vector<BitSet>& inc, double) -> bool {
          take(reconstruct(st, inc)); return false;              // register live, never stop early
        };
        bool proven = false;
        auto sv = g_deadline; bool sh = g_has_deadline;          // time-box: leave a share for heuristics
        set_time_limit_seconds(std::min(remaining_seconds() * 0.8, remaining_seconds()));
        try { auto [sol, dual] = solve_bcp(R.reduced_maf); (void)dual; take(reconstruct(st, sol)); proven = true; }
        catch (const Interrupted&) {} catch (const std::exception&) {}
        g_deadline = sv; g_has_deadline = sh;
        g_incumbent_hook = saved_hook;
        if (proven) return best;                                 // exact optimum -> best possible
      }
    } catch (const Interrupted&) {} catch (const std::exception&) {}
  }

  // Reduce/cluster pipeline with randomized restarts until the budget is (nearly) gone.  The FIRST pass
  // is A-only (defer_mode=1): guaranteed to complete fast, registering a strong incumbent early even on
  // huge instances (where a full A+B pass would exceed the budget and be discarded to the floor).  Every
  // later pass uses the adaptive window (defer_mode=2): it dedicates both options for quality but drops
  // option B once the budget can no longer cover finishing the pass, so a near-complete refinement pass
  // still lands an improvement instead of being lost.
  // `take` doubles as the PROGRESSIVE registration sink: run_reduced_once forwards every composed
  // provisional full-instance forest here (guarded), so an early SIGTERM mid-pass keeps the best partial
  // work instead of only whole-pass results.  Registering a valid provisional incumbent can only help.
  static const int PROG = [](){ const char* e=std::getenv("MAF_PROG"); return (e&&*e)?std::atoi(e):1; }();
  std::function<void(const Sol&)> progress;
  if (PROG) progress = [&](const Sol& s2){ take(s2); };   // MAF_PROG=0 => whole-pass registration only (old behaviour)
  static const int FIRST_DEFER_OVERRIDE = [](){
    const char* e = std::getenv("MAF_FIRST_DEFER");
    return (e && *e) ? std::atoi(e) : -1;
  }();
  static const int RESTART_DEFER_OVERRIDE = [](){
    const char* e = std::getenv("MAF_RESTART_DEFER");
    return (e && *e) ? std::atoi(e) : -1;
  }();
  uint64_t s = seed;
  for (int restart = 0; remaining_seconds() > 0.5; ++restart) {
    int defer_mode = (restart == 0) ? 1 : 2;
    if (restart == 0 && FIRST_DEFER_OVERRIDE >= 0) defer_mode = FIRST_DEFER_OVERRIDE;
    if (restart > 0 && RESTART_DEFER_OVERRIDE >= 0) defer_mode = RESTART_DEFER_OVERRIDE;
    try { take(run_reduced_once(maf, s, defer_mode, progress)); }
    catch (const Interrupted&) { break; } catch (const std::exception&) { break; }
    s = s * 6364136223846793005ULL + 1442695040888963407ULL;   // advance the seed (LCG)
  }
  return best;
}

}} // namespace maf::heur
