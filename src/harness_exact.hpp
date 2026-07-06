// harness_exact.hpp — EXACT-track entry (selected when MAF_HEURISTIC is NOT defined).
//
// Reads a MAF instance from stdin, solves to PROVEN optimality (solve_maf), writes the agreement
// forest. Competition contract: if interrupted before a provably-optimal forest is in hand (time
// budget / SIGTERM) or on any error, emit NOTHING — a suboptimal answer claimed optimal would be a
// DISQUALIFICATION. (Timing out with no output is allowed; a wrong answer is not.)
#pragma once
#include "solve.hpp"
#include "io.hpp"
#include "budget.hpp"
#include "maf_util.hpp"      // is_valid_solution (final safety guard before output)
#include <iostream>
#include <cstdio>

namespace maf {

inline int run_exact_track() {
  install_budget_handlers();                     // SIGINT/SIGTERM + MAF_TIMELIMIT -> Interrupted

  MAFP maf;
  try {
    maf = read_instance(std::cin);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "maf: parse error: %s\n", e.what());
    return 0;                                     // no valid instance -> emit nothing
  }
  if (maf.n == 0) return 0;

  try {
    auto [sol, dual] = solve_maf(maf);
    (void)dual;
    // EXACT-TRACK SAFETY GUARD (2026-07-04): solve_maf returns a forest it believes is proven-optimal,
    // but a defect in reduction/cluster RECONSTRUCTION or pricing could return a MALFORMED forest.
    // write_forest is the ONLY thing this binary prints to stdout, and it emits O(sum|component|):
    // for a valid partition that is O(n), but for a non-partition (overlapping/duplicated components)
    // it is up to O(n^2) — which blows past OPTIL's per-instance output cap on large instances
    // (Output-Limit-Exceeded), and any topology-inconsistent component is a Wrong Answer. Both are
    // DISQUALIFICATIONS. Validate against the exact structural checks the PACE checker applies and, if
    // the forest is invalid, emit NOTHING (indistinguishable from a timeout — the sanctioned safe
    // outcome). Mirrors the heuristic track's guard() and the CSP solver's validate_result().
    if (is_valid_solution(maf, sol)) {
      write_forest(std::cout, maf, sol);
      std::cout.flush();
    } else {
      std::fprintf(stderr, "maf: internal solution failed structural validation (%zu components) — "
                           "suppressed to avoid OLE/WA disqualification\n", sol.size());
    }
  } catch (const Interrupted&) {
    // Budget exhausted / signalled before optimality: emit nothing (exact-track safe).
  } catch (const std::exception& e) {
    std::fprintf(stderr, "maf: solve error: %s\n", e.what());
  }
  return 0;
}

} // namespace maf
