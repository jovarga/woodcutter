// harness_heuristic.hpp — HEURISTIC-track entry (selected when MAF_HEURISTIC is defined).
//
// PACE 2026 heuristic track: 2-tree, 5-minute per-instance limit (MAF_TIMELIMIT, default 300 s),
// terminates on SIGTERM (10 s grace) and must output its BEST VALID solution — an empty/infeasible
// forest is a DISQUALIFICATION; suboptimal is fine, and extra time is not penalised so the whole
// budget is spent improving. Realised here: a valid incumbent is kept PRE-RENDERED at all times
// (singletons floor, then every improvement), and SIGTERM/SIGINT flush it with a single async-
// signal-safe write() before exiting.
#pragma once
#include "heuristics.hpp"                          // maf::heur helpers
#include "_lineartime_reductions_clusters.hpp" // maf::heur::run_heuristic_reduced (always-on preprocessing)
#include "io.hpp"
#include "budget.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <atomic>
#include <functional>
#include <csignal>
#include <cstdlib>
#include <cerrno>
#include <unistd.h>

namespace maf {
namespace heur_harness {   // internal helpers for the SIGTERM-safe incumbent buffer

inline std::atomic<const std::string*> g_out{nullptr};

inline std::string render(const MAFP& maf, const heur::Sol& sol) {
  std::ostringstream oss; write_forest(oss, maf, sol); return oss.str();
}
inline void publish(const MAFP& maf, const heur::Sol& sol) {
  g_out.store(new std::string(render(maf, sol)), std::memory_order_release);  // old buffer leaked (rare, tiny)
}
inline void write_all(const char* p, size_t n) {   // async-signal-safe: LOOP over short writes/EINTR
  size_t off = 0;                                   // a single write() to a pipe stops at ~64 KB, so a
  while (off < n) {                                 // large forest MUST be written in a loop or it is
    ssize_t r = write(STDOUT_FILENO, p + off, n - off);   // silently TRUNCATED -> invalid MAF on optil.
    if (r > 0) { off += (size_t)r; continue; }
    if (r < 0 && (errno == EINTR || errno == EAGAIN)) continue;
    break;                                          // unrecoverable error -> give up (nothing else to do)
  }
}
inline void flush_and_exit(int) {                                            // async-signal-safe
  const std::string* s = g_out.load(std::memory_order_acquire);
  if (s) write_all(s->data(), s->size());
  _exit(0);
}

} // namespace heur_harness

inline int run_heuristic_track() {
  std::signal(SIGTERM, heur_harness::flush_and_exit);
  std::signal(SIGINT,  heur_harness::flush_and_exit);

  MAFP maf;
  try { maf = read_instance(std::cin); }
  catch (const std::exception&) { return 0; }
  if (maf.n == 0) return 0;

  heur::Sol incumbent = heur::singletons(maf);        // immediate valid floor
  heur_harness::publish(maf, incumbent);

  const char* tl = std::getenv("MAF_TIMELIMIT");      // 5 min default (heuristic track)
  set_time_limit_seconds((tl && *tl) ? std::atof(tl) : 300.0);

  std::function<void(const heur::Sol&)> consider = [&](const heur::Sol& cand) {
    heur::Sol g = heur::guard(maf, cand);
    if (g.size() < incumbent.size()) { incumbent = g; heur_harness::publish(maf, incumbent); }
  };

  try { heur::run_heuristic_reduced(maf, consider); }   // always applies reductions + cluster decomposition
  catch (const Interrupted&) {}
  catch (const std::exception&) {}

  const std::string* s = heur_harness::g_out.load(std::memory_order_acquire);
  if (s) std::cout << *s;
  std::cout.flush();
  return 0;
}

} // namespace maf
