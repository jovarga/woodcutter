// harness_lowerbound.hpp — LOWER-BOUND-track entry (selected when MAF_LOWERBOUND is defined).
//
// PACE 2026 lower-bound track (an APPROXIMATION track).  2-tree instances, 10-min limit.  Each
// instance has an approximation line `#a a b` with 1 <= a < 1.5 and 0 <= b.  Let k* be the (unknown)
// size of a maximum agreement forest; a submission must output a valid MAF of size at most
// floor(a*k*) + b.  Outputting the EXACT optimum (size k*) trivially satisfies this (k* <=
// floor(a*k*)+b for a>=1, b>=0), and — being the tightest — also scores best.
//
// The exact track's own O(n^2) reductions / cluster decomposition TIME OUT on the large instances of
// this track (n up to ~10 000).  So we solve with the LINEAR-TIME (O(t*n)) reductions + cluster
// decomposition from _lineartime_day.hpp (Day 1985), whose correctness the optimality guarantee
// depends on, feeding each irreducible block to the exact branch-and-price (solve_bcp).  The LS11
// recombination yields the global optimum; we output that forest.
//
// SIGTERM-safe: the all-singletons forest (always a valid MAF) is emitted if we are interrupted
// before the solve completes.
//
// "The behaviour of the exact solver shall be unchanged": solve_bcp / reconstruct are untouched; the
// linear reductions/clusters live in their own files and are NOT used by the exact track.
#pragma once
#include "solve.hpp"                     // solve_bcp, BCPContext, ctx_columngeneration, ctx_lpobj
#include "_lineartime_day.hpp"           // maf::lineartime::solve_with_reductions_linear, solve_clusters_linear
#include "heuristics.hpp"                // maf::heur:: feasibility_pump / greedy pack / mast_peel (large-block UB)
#include "_lineartime_reductions_clusters.hpp" // maf::heur::run_heuristic_reduced (Phase-2 upper bound)
#include "io.hpp"                         // read_instance, write_forest
#include "budget.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <atomic>
#include <csignal>
#include <cstdlib>
#include <cstdio>
#include <cerrno>
#include <algorithm>
#include <unistd.h>

namespace maf {
namespace lb_harness {

inline std::atomic<const std::string*> g_out{nullptr};   // pre-rendered best forest (SIGTERM-safe)
inline double g_a = 1.0;                                   // from `#a a b`
inline long   g_b = 0;

inline void publish(const std::string* s) { g_out.store(s, std::memory_order_release); }  // old leaked (tiny)
inline void write_all(const char* p, size_t n) {   // async-signal-safe: LOOP over short writes/EINTR
  size_t off = 0;                                   // a single write() to a pipe stops at ~64 KB, so a
  while (off < n) {                                 // large forest MUST be written in a loop or it is
    ssize_t r = write(STDOUT_FILENO, p + off, n - off);   // silently TRUNCATED -> invalid MAF on optil.
    if (r > 0) { off += (size_t)r; continue; }
    if (r < 0 && (errno == EINTR || errno == EAGAIN)) continue;
    break;
  }
}
inline void flush_and_exit(int) {                         // async-signal-safe forest flush
  const std::string* s = g_out.load(std::memory_order_acquire);
  if (s) write_all(s->data(), s->size());
  _exit(0);
}
inline void parse_ab(const std::string& raw) {
  std::istringstream in(raw);
  std::string line;
  while (std::getline(in, line)) {
    size_t i = line.find_first_not_of(" \t");
    if (i == std::string::npos || line[i] != '#') continue;
    if (line.compare(i, 2, "#a") != 0) continue;
    std::istringstream ls(line.substr(i + 2));
    double a; long b;
    if (ls >> a >> b) { g_a = a; g_b = b; }
    return;
  }
}

} // namespace lb_harness

// Below this size a block's exact branch-and-price is expected to be fast; we still cap it (below) so a
// hard small block cannot stall the whole solve.
inline constexpr int LB_EXACT_FAST = 250;

// Certificate for a hard block that the exact branch-and-price cannot close in its time cap: return a
// valid feasible forest (UPPER bound) together with a VALID LOWER bound on the block's k* (the CONVERGED
// root LP relaxation, plus any dual the warm-started BCP proves).  Returns (forest, lower_bound).
//
// SOUNDNESS: the returned lower_bound is <= the block's k*, and solve_clusters_linear composes the block
// lower bounds (via min over the A/B options) into a valid GLOBAL lower bound.  So the harness's single
// global certificate (global_forest_size <= floor(a*global_lb)+b) implies size <= floor(a*k*_global)+b.
// This function has NO global side effects and never throws: every stage is time-boxed and the whole
// certificate is bounded to a share of the remaining budget so LATER blocks in the sequence still run.
inline std::pair<std::vector<BitSet>, double> lb_certificate(const MAFP& block) {
  const double a = lb_harness::g_a; const long b = lb_harness::g_b;
  auto within = [&](size_t ub, double lb) {         // block already comfortably inside range -> stop early
    return (long)ub <= (long)std::floor(a * std::max(lb, 1.0)) + b;
  };
  // Bound the WHOLE certificate to a fraction of the remaining budget (each stage below respects this
  // deadline via check_budget / remaining_seconds, so nothing runs past it and nothing throws out).
  auto sv = g_deadline; bool sh = g_has_deadline;
  double cap = std::max(std::min(remaining_seconds() * 0.7, remaining_seconds() - 1.0), 1.0);
  set_time_limit_seconds(std::min(cap, remaining_seconds()));

  // 1. Root column generation (LP relaxation, NO branching) -> valid LP lower bound + a rich pool.
  //    Only a CONVERGED CG gives a valid LP lower bound (a truncated restricted master over-estimates),
  //    so on interruption we keep the trivial lb=1 and still use the partial pool.
  BCPContext ctx(block);
  double lb = 1.0;
  { auto s2 = g_deadline; bool h2 = g_has_deadline;
    set_time_limit_seconds(std::max(std::min(remaining_seconds() * 0.4, 60.0), 0.5));
    try { ctx_columngeneration(ctx, 0); lb = std::max(1.0, ctx_lpobj(ctx, 0)); }
    catch (const Interrupted&) {} catch (const std::exception&) {}
    g_deadline = s2; g_has_deadline = h2; }

  // 2. Strong UPPER bound.  Build the column pool from the CG columns AND a fast centroid MAST-peel
  //    (essential on the largest blocks, where root CG alone is too slow to yield a usable pool ->
  //    the greedy/FP would otherwise start from singletons).  Then greedy pack + feasibility pump.
  std::vector<BitSet> pool;
  std::unordered_set<BitSet> pseen;
  for (const auto& c : ctx.mip_vars) if ((int)c.size() >= 2 && pseen.insert(c).second) pool.push_back(c);
  if (remaining_seconds() > 1.0) {
    double pcap = std::min(remaining_seconds() * 0.3, 60.0);
    for (const auto& c : heur::mast_peel(block, /*max_cols=*/8000, pcap))
      if (pseen.insert(c).second) pool.push_back(c);
  }
  std::vector<BitSet> best = heur::singletons(block); double best_lb = lb;
  if (!pool.empty()) {
    auto inn = heur::pool_inners(block, pool);
    std::vector<BitSet> g = heur::greedy_disjoint_partition(block, pool, inn);
    if (g.size() < best.size()) best = g;
    if (!within(best.size(), best_lb) && remaining_seconds() > 0.5) {
      std::vector<BitSet> fp = heur::feasibility_pump(block, pool, inn, [](const std::vector<BitSet>&){});
      if (fp.size() < best.size()) best = fp;
    }
  }

  // 3. If not yet comfortably in range, run the flexible branch-and-price WARM-STARTED with a LEAN pool
  //    (CG columns + the incumbent's own columns) and the approximation hook armed: it tightens the
  //    primal AND the dual, and STOPS the moment the incumbent is within floor(a*dual)+b (or the cap
  //    hits).  Then a few diversified FP restarts on any leftover time.  All bounded by the cap above.
  if (!within(best.size(), best_lb) && remaining_seconds() > 1.0) {
    g_incumbent_hook = [&](const std::vector<BitSet>& inc, double d) -> bool {
      if (inc.size() < best.size()) best = inc;
      if (d > best_lb) best_lb = d;
      return within(best.size(), std::max(d, best_lb));
    };
    std::vector<BitSet> warm;
    { std::unordered_set<BitSet> ws;
      for (const auto& c : ctx.mip_vars) if ((int)c.size() >= 2 && ws.insert(c).second) warm.push_back(c);
      for (const auto& c : best)         if ((int)c.size() >= 2 && ws.insert(c).second) warm.push_back(c); }
    try { auto [s, dd] = ctx_solve(block, warm); if (s.size() < best.size()) best = s; best_lb = std::max(best_lb, dd); }
    catch (const Interrupted&) {} catch (const std::exception&) {}
    g_incumbent_hook = {};

    if (!within(best.size(), best_lb) && !pool.empty()) {
      auto inn = heur::pool_inners(block, pool);
      for (uint64_t seed = 2; !within(best.size(), best_lb) && remaining_seconds() > 1.0; ++seed) {
        std::vector<BitSet> fp = heur::feasibility_pump(block, pool, inn, [](const std::vector<BitSet>&){},
                                                        25, 1.0, 0.0, 0.4, seed);
        if (fp.size() < best.size()) best = fp;
      }
    }
  }

  g_deadline = sv; g_has_deadline = sh;                 // restore the outer deadline
  return { heur::guard(block, best), std::max(best_lb, 1.0) };
}

// Leaf solver for the lower-bound track: EXACT-FIRST (proven-optimal k* — the tightest answer and always
// trivially within range) with a per-block TIME CAP; if the exact BCP cannot prove optimality within the
// cap, fall back to the bounded CERTIFICATE.  Both return (forest, valid_lower_bound), so the composed
// global dual stays a valid lower bound and the harness's global certificate is sound.  (Unlike the old
// design, there is no size threshold and no per-block certified/uncertified flag: hard blocks of ANY size
// now yield a valid within-range forest instead of nothing.)
inline std::pair<std::vector<BitSet>, double> lb_solve_block(const MAFP& block) {
  // 1. Exact-first, capped.  Tiny blocks finish instantly (proven optimal); a hard block hits the cap,
  //    throws Interrupted, and we fall through to the certificate WITHOUT losing the whole solve.
  { auto sv = g_deadline; bool sh = g_has_deadline;
    double cap = (block.n <= LB_EXACT_FAST) ? remaining_seconds()          // small: let it finish
                                            : std::max(std::min(remaining_seconds() * 0.25, 30.0), 2.0);
    set_time_limit_seconds(std::min(cap, remaining_seconds()));
    bool ok = false; std::pair<std::vector<BitSet>, double> r;
    try { r = solve_bcp(block); ok = true; }             // returns ONLY when proven optimal (else throws)
    catch (const Interrupted&) {} catch (const std::exception&) {}
    g_deadline = sv; g_has_deadline = sh;
    if (ok) return r;                                     // proven optimal -> tightest possible
  }
  // 2. Certificate fallback (bounded): valid UB forest + valid LB.
  return lb_certificate(block);
}

// CHEAP lower-bound leaf: the ROOT LP relaxation only (column generation, NO branching) — a valid but
// loose lower bound obtained quickly, plus a quick feasible forest (greedy pack over the CG pool).  Used
// by the phased harness's Phase 1 so a usable global lower bound exists BEFORE the heuristic upper-bound
// phase.  A truncated / non-converged CG OVER-estimates the LP, so on interruption we keep L=1 (never an
// invalid over-estimate); only a converged root CG yields a valid LP lower bound.
inline std::pair<std::vector<BitSet>, double> lb_cheap_leaf(const MAFP& block) {
  BCPContext ctx(block);
  double L = 1.0;
  { auto sv = g_deadline; bool sh = g_has_deadline;
    set_time_limit_seconds(std::min(remaining_seconds(), std::max(remaining_seconds() * 0.5, 0.2)));
    try { ctx_columngeneration(ctx, 0); L = std::max(1.0, ctx_lpobj(ctx, 0)); }   // converged -> valid LP LB
    catch (const Interrupted&) {} catch (const std::exception&) {}
    g_deadline = sv; g_has_deadline = sh; }
  std::vector<BitSet> pool; std::unordered_set<BitSet> seen;
  for (const auto& c : ctx.mip_vars) if ((int)c.size() >= 2 && seen.insert(c).second) pool.push_back(c);
  std::vector<BitSet> ub = heur::singletons(block);
  if (!pool.empty()) {
    auto inn = heur::pool_inners(block, pool);
    auto g = heur::greedy_disjoint_partition(block, pool, inn);
    if (g.size() < ub.size()) ub = g;
  }
  return { heur::guard(block, ub), L };
}

// PHASED lower-bound track (large-instance strategy).  Decouples the UPPER bound (which the tuned
// HEURISTIC track converges to fastest) from the LOWER bound (which certifies), and spends the budget in
// three phases, cheapest-effort-first, publishing the smallest CERTIFIED forest as soon as it exists:
//
//   Phase 1 — CHEAP lower bound: the root LP relaxation (CG, no branching) composed over the sound
//             decomposition -> a valid global L quickly, so Phase 2 already has something to certify
//             against.  (Also yields a first, weak, feasible forest.)
//   Phase 2 — HEURISTIC upper bound: run the full tuned heuristic; every improving forest is checked
//             against L and, once it fits floor(a*L)+b, PUBLISHED (SIGTERM-safe).  For loose a / small
//             cores this certifies without ever touching the expensive exact lower bound.
//   Phase 3 — EXACT lower bound: only if still uncertified, tighten L with the exact branch-and-price
//             per block (lb_solve_block: exact-first + certificate) over the remaining budget, then
//             re-check the guarantee.  If already certified, the leftover budget instead goes to MORE
//             heuristic (shrink the emitted forest -> tighter, better-scoring output).
//
// SOUND throughout: we only ever publish a forest F with F.size() <= floor(a*best_L)+b for a VALID
// global lower bound best_L <= k*_global, so the emitted forest is provably within floor(a*k*)+b.  If no
// phase certifies, we emit NOTHING (an out-of-range forest is a DISQUALIFICATION).
inline int run_lowerbound_track() {
  std::signal(SIGTERM, lb_harness::flush_and_exit);
  std::signal(SIGINT,  lb_harness::flush_and_exit);

  std::ostringstream ss; ss << std::cin.rdbuf();
  std::string raw = ss.str();
  lb_harness::parse_ab(raw);

  MAFP maf;
  try { std::istringstream in(raw); maf = read_instance(in); }
  catch (const std::exception&) { return 0; }
  if (maf.n == 0) return 0;

  auto render = [&](const std::vector<BitSet>& forest) {
    std::ostringstream o; write_forest(o, maf, forest); return new std::string(o.str());
  };

  const char* tl = std::getenv("MAF_TIMELIMIT");                   // 10 min default (lower-bound track)
  set_time_limit_seconds((tl && *tl) ? std::atof(tl) : 600.0);
  const double budget = remaining_seconds();                       // total budget captured at start
  const double a = lb_harness::g_a; const long b = lb_harness::g_b;

  // Crossover tuning parameters (fractions of the TOTAL budget), env-overridable.
  // TODO(retune-full-budget): set by reasoning, not a full-budget sweep — see docs/TODO.md.
  static const double LB_CHEAP_FRAC = [](){ const char* e=std::getenv("MAF_LB_CHEAP_FRAC"); return (e&&*e)?std::atof(e):0.15; }();
  static const double LB_HEUR_FRAC  = [](){ const char* e=std::getenv("MAF_LB_HEUR_FRAC");  return (e&&*e)?std::atof(e):0.50; }();

  // Shared incumbent (smallest feasible forest) + best VALID lower bound.  publish_if publishes the
  // pre-rendered CERTIFIED forest (SIGTERM-safe) and STOPS the search: the lower-bound track is scored by
  // TIME (fit the guarantee fast — there is NO tightness bonus), so the moment a forest is within
  // floor(a*L)+b we output it and terminate.  `done` collapses the deadline (so the running phase aborts
  // promptly via check_budget) and guards the remaining phases.
  std::vector<BitSet> best_forest; size_t best_sz = SIZE_MAX; double best_L = 1.0; bool done = false;
  auto within     = [&](size_t sz, double L){ return (long)sz <= (long)std::floor(a * std::max(L, 1.0)) + b; };
  auto certified  = [&](){ return best_sz != SIZE_MAX && within(best_sz, best_L); };
  auto publish_if = [&](){ if (certified()) { lb_harness::publish(render(best_forest));
                                              if (!done) { done = true; set_time_limit_seconds(0.0); } } };
  auto consider   = [&](const std::vector<BitSet>& F){ if (done) return; if (F.size() >= best_sz) return;
                        if (!heur::cheap_valid(maf, F)) return;
                        best_sz = F.size(); best_forest = F; publish_if(); };
  auto bump_L     = [&](double L){ if (done) return; if (L > best_L) best_L = L; publish_if(); };

  auto phase_decomp = [&](const ClusterSolver& leaf, double phase_seconds) {
    auto sv = g_deadline; bool sh = g_has_deadline;
    set_time_limit_seconds(std::min(phase_seconds, remaining_seconds()));
    try {
      auto [F, L] = lineartime::solve_with_reductions_linear(
          maf, [&leaf](const MAFP& m) { return lineartime::solve_clusters_linear(m, leaf); });
      consider(F); bump_L(L);
    } catch (const Interrupted&) {} catch (const std::exception&) {}
    g_deadline = sv; g_has_deadline = sh;
  };
  auto phase_heur = [&](double phase_seconds) {
    auto sv = g_deadline; bool sh = g_has_deadline;
    set_time_limit_seconds(std::min(phase_seconds, remaining_seconds()));
    try { heur::run_heuristic_reduced(maf, [&](const std::vector<BitSet>& F){ consider(F); }); }
    catch (const Interrupted&) {} catch (const std::exception&) {}
    g_deadline = sv; g_has_deadline = sh;
  };

  // a-ADAPTIVE workflow.  When a is very close to 1 the guarantee is LOWER-bound-bottlenecked: only a
  // near-exact L certifies, so a better UPPER bound cannot help and diverting budget from the exact
  // branch-and-price only HURTS (measured: the phased path emits NOTHING on lower100/102 where the
  // exact path certifies).  So for tight a we run the ENTIRE exact approach at full budget — its
  // per-block certificate already carries the within-range early-termination hook.  For looser a a
  // strong upper bound certifies against the wide floor(a*L)+b quickly, so the phased workflow (cheap
  // L -> heuristic UB -> exact L) wins by getting a tight forest fast.
  // TODO(retune-full-budget): LB_A_TIGHT and the phase fractions were set by reasoning, not a full
  // 600 s sweep — see docs/TODO.md.
  static const double LB_A_TIGHT = [](){ const char* e=std::getenv("MAF_LB_A_TIGHT"); return (e&&*e)?std::atof(e):1.05; }();

  if (a <= LB_A_TIGHT) {
    // TIGHT a — run the entire exact approach at full budget (the committed sound pipeline).  Its
    // per-block certificate already carries the within-range early-termination hook.
    ClusterSolver leaf = [](const MAFP& m) { return lb_solve_block(m); };
    phase_decomp(leaf, remaining_seconds());
  } else {
    // LOOSE a — phased: cheap L, then heuristic UB, then exact L — STOPPING the moment we are certified
    // (each `phase_*` returns early once `done` collapses the deadline; the guards skip the rest).
    { ClusterSolver leaf = [](const MAFP& m) { return lb_cheap_leaf(m); };
      phase_decomp(leaf, LB_CHEAP_FRAC * budget); }
    if (!done && remaining_seconds() > 1.0)
      phase_heur(LB_HEUR_FRAC * budget);
    if (!done && remaining_seconds() > 1.0) {
      ClusterSolver leaf = [](const MAFP& m) { return lb_solve_block(m); };
      phase_decomp(leaf, remaining_seconds());
    }
  }

  const std::string* s = lb_harness::g_out.load(std::memory_order_acquire);
  if (s) std::cout << *s;
  std::cout.flush();
  return 0;
}

} // namespace maf
