// main.cpp — PACE 2026 track dispatcher.
//
// Compile-time track selection (set by CMake -DMAF_TRACK=exact|heuristic|lowerbound):
//   * default (no macro)  -> EXACT track: prove optimality, emit nothing on timeout.
//   * -DMAF_HEURISTIC     -> HEURISTIC track: 5 min, flush best incumbent on SIGTERM.
//   * -DMAF_LOWERBOUND    -> LOWER-BOUND track: 10 min, output an integer lower bound (a,b early exit).
// The tracks differ only in the harness; the solver/heuristic algorithms are shared code.
#if defined(MAF_LOWERBOUND)
  #include "harness_lowerbound.hpp"
#elif defined(MAF_HEURISTIC)
  #include "harness_heuristic.hpp"
#else
  #include "harness_exact.hpp"
#endif

int main() {
#if defined(MAF_LOWERBOUND)
  return maf::run_lowerbound_track();
#elif defined(MAF_HEURISTIC)
  return maf::run_heuristic_track();
#else
  return maf::run_exact_track();
#endif
}
