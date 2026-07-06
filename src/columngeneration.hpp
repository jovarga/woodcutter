// Only the two functions that are LIVE on the production BCP path are implemented here:
//   * solve_rmp        — solve the (binary) restricted master problem, return obj + partition
//   * select_variables — pick which priced columns to add to the master each CG iteration
//
// The other experimental column-generation entry points are not on the production path.
#pragma once
#include "core.hpp"
#include "instance.hpp"
#include "mpmodel.hpp"
#include <vector>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <limits>

namespace maf {

// Solve the restricted master problem over the given columns as a binary program, then return the
// objective together with the induced partition: every selected column, plus a singleton {x} for
// each leaf left uncovered.
//
// NOTE on determinism: a set-packing MIP can have several optimal partitions; HiGHS (single-thread)
// is deterministic, but its choice among optima can depend on the column order.  The MIP *objective*
// (= number of components = k*) is order-independent, so the value returned here — and hence the BCP
// control flow that depends on it. The specific partition may differ between two equally-optimal
// solutions; both are valid optima (the exact track accepts any optimum).
//
// `mip_start` (3cb299e): optional MIP warm-start — a known feasible incumbent partition seeded into
// HiGHS as a starting solution.
// Speedup only (the MIP optimum is start-independent); nullptr reproduces the pre-3cb299e cold-start.
inline std::pair<double, std::vector<BitSet>>
solve_rmp(const MAFP& maf, const std::vector<BitSet>& var_sets,
          const std::vector<BitSet>* mip_start = nullptr) {
  MPModel mp(maf, var_sets, /*binary=*/true);
  if (mip_start) mp.set_mip_start(*mip_start);
  double obj = mp.solve();

  // Selected nonsingleton columns.
  std::vector<BitSet> sol;
  BitSet covered;
  for (int r = 0; r < (int)mp.col_sets.size(); ++r)
    if (mp.sol_value(r) >= 1 - 1e-8) {
      sol.push_back(mp.col_sets[r]);
      covered.union_with(mp.col_sets[r]);
    }
  // Add uncovered leaves as singletons.
  for (int x = 1; x <= maf.n; ++x)
    if (!covered.contains(x)) sol.push_back(BitSet{x});

  return {obj, sol};
}

// Of the priced columns (score `best_ws[k]`, taxa set `best_sets[k]`), keep those with reduced cost
// above 1 (score > 1+1e-8) that are not already in the master (`vars` = mp.vars keys), deduplicate
// by set value keeping the first occurrence, and return the first `lim` of them.
//
// Selection and `lim` truncation use the original priced order, not score order. In production,
// `lim` is effectively unbounded, so all qualifying distinct columns are returned.
inline std::vector<BitSet>
select_variables(const std::vector<W>& best_ws, const std::vector<BitSet>& best_sets,
                 const std::unordered_map<BitSet, int>& vars,
                 int lim = std::numeric_limits<int>::max()) {
  std::vector<BitSet> out;
  std::unordered_set<BitSet> seen;            // dedup by set value (unique!, first-seen)
  for (size_t k = 0; k < best_ws.size(); ++k) {
    if (best_ws[k] > 1 + 1e-8 && vars.find(best_sets[k]) == vars.end()
        && seen.find(best_sets[k]) == seen.end()) {
      seen.insert(best_sets[k]);
      out.push_back(best_sets[k]);
      if ((int)out.size() >= lim) break;
    }
  }
  return out;
}

} // namespace maf
