// budget.hpp — time-budget + SIGINT/SIGTERM handling.
//
// Competition requirement #6 (qualification-critical): the exact-track solver must terminate within
// the time limit and on SIGINT/SIGTERM, and — because a suboptimal answer claimed optimal is a
// DISQUALIFICATION — it must emit NOTHING if it is interrupted before a provably-optimal solution
// is in hand.  The solver checks `check_budget()` at loop heads (and the master caps each HiGHS solve
// at the remaining time); on budget exhaustion it throws `Interrupted`, which main() catches and
// then prints no solution.
#pragma once
#include <atomic>
#include <chrono>
#include <stdexcept>
#include <csignal>
#include <cstdlib>

namespace maf {

struct Interrupted : std::exception {
  const char* what() const noexcept override { return "interrupted or time budget exceeded"; }
};

// Thrown by MPModel::solve() when the master is INFEASIBLE (as opposed to a time-limit interruption).
// In branch-and-price this means the node's branch constraints admit no feasible point, so the node
// has no integer solution and must be PRUNED by the caller (standard B&P) — NOT treated as a timeout.
struct NodeInfeasible : std::exception {
  const char* what() const noexcept override { return "branch-and-price node is infeasible (prune)"; }
};

inline std::atomic<bool> g_interrupted{false};
inline std::chrono::steady_clock::time_point g_deadline{};
inline bool g_has_deadline = false;

inline void set_time_limit_seconds(double seconds) {
  g_deadline = std::chrono::steady_clock::now()
             + std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(seconds));
  g_has_deadline = true;
}

// Seconds left until the deadline (huge if none); never negative-reported below 0.
inline double remaining_seconds() {
  if (!g_has_deadline) return 1e18;
  return std::chrono::duration<double>(g_deadline - std::chrono::steady_clock::now()).count();
}

inline bool budget_exceeded() {
  return g_interrupted.load(std::memory_order_relaxed)
      || (g_has_deadline && std::chrono::steady_clock::now() >= g_deadline);
}

inline void check_budget() { if (budget_exceeded()) throw Interrupted{}; }

inline void _on_signal(int) { g_interrupted.store(true, std::memory_order_relaxed); }

// Install SIGINT/SIGTERM handlers and (optionally) a deadline from the MAF_TIMELIMIT env var
// (seconds).  Call once at startup.
inline void install_budget_handlers() {
  std::signal(SIGINT, _on_signal);
  std::signal(SIGTERM, _on_signal);
  if (const char* tl = std::getenv("MAF_TIMELIMIT")) {
    double s = std::atof(tl);
    if (s > 0) set_time_limit_seconds(s);
  }
}

} // namespace maf
