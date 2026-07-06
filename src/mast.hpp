// The pricing subproblem ("MAST solve"): given the master's dual prices, find columns (agreement
// subtrees) of positive reduced cost.  This file is the dispatcher mast_solve; the heavy lifting is
// in the sub-algorithm headers it includes.
//
// PRODUCTION PATH ONLY. The dispatcher returns score/taxa vectors (one entry per rooted subtree).
#pragma once
#include "core.hpp"
#include "instance.hpp"
#include "mpmodel.hpp"
#include "mast/trees.hpp"
#include "mast/multi_dp.hpp"
#include "mast/quad_dp.hpp"      // plain O(n^2) 2-tree pricer (checkable; opt-in via MAF_2TREE_PRICER)
#include "mast/_centroid.hpp"    // O(n log n) 2-tree centroid pricer
#include <vector>
#include <utility>
#include <cstdlib>
#include <string>

namespace maf {

// 2-tree MAST pricer selection. BOTH engines are compiled into every build and compute the SAME
// per-vertex max MAST score => the SAME master LP optimum and the SAME k*. They can differ only in
// WHICH equal-score subtree they reconstruct on a tie, which at most changes which columns get
// priced — never the optimum. Therefore the exact-track ANSWER is identical either way.
//   * Centroid (default): O(n log n), and the only viable choice on the huge 2-tree lower-bound /
//     heuristic instances (quad is O(n^2) SPACE).
//   * Quad: plain O(n^2) DP (mast/quad_dp.hpp) — the simple "checkable" reference engine.  Opt-in.
// Pick quad with env MAF_2TREE_PRICER=quad, or programmatically `two_tree_pricer_ref() = ...::Quad`.
// A build compiled with -DMAF_2TREE_CENTROID FORCES centroid.
enum class TwoTreePricer { Centroid, Quad };
inline TwoTreePricer& two_tree_pricer_ref() {
  static TwoTreePricer p = []{
#if defined(MAF_2TREE_CENTROID)
    return TwoTreePricer::Centroid;
#else
    if (const char* e = std::getenv("MAF_2TREE_PRICER")) {
      std::string s(e);
      if (s == "quad" || s == "quadratic") return TwoTreePricer::Quad;
    }
    return TwoTreePricer::Centroid;
#endif
  }();
  return p;
}

// Returns the priced columns as (score, taxa) pairs: the concatenation over trees of the
// per-rooted-node best subtrees, filtered to |taxa| > 1.
//
// Production specifics:
//   * alpha (taxon weights) is all 1.
//   * beta (inner-vertex prices) comes from the last master LP solve: ctx.beta[(i,v)] == mp.beta(i,v).
//     The PhyloTree internal weight of vertex v is iw[v] = -mp.beta(i, v).
//   * blacklist is ALWAYS empty on this branch: ctx_price! only reaches the no-cut branch when
//     ctx.cuts is empty (the cut branch uses pricing_mast instead).  We keep the parameter for
//     faithfulness but production never passes a non-empty one.
//
// Dispatch: two trees with empty blacklist -> the 2-tree pricer (centroid by default, quad if selected
// via two_tree_pricer_ref); otherwise the multi-tree O(n^3)-style DP.
//
// The PhyloTrees are built by the CALLER (solve.hpp ctx_price), because under branching the internal
// weights are not simply -beta: each branch constraint adds its dual to the internal weight of the
// branched vertex. `n` is maf.n (for the all-ones alpha vector).
inline std::vector<std::pair<W, BitSet>>
mast_solve(const std::vector<mast::PhyloTree>& trees, int n,
           const std::vector<std::vector<int>>& blacklist = {}) {
  const int T = (int)trees.size();

  std::vector<std::pair<W, std::vector<int>>> raw;
  if (T == 2 && blacklist.empty()) {
    // 2-tree case: centroid (default) or the plain O(n^2) quad DP (opt-in). Both give the SAME
    // per-vertex max score / LP optimum / k*.
    raw = (two_tree_pricer_ref() == TwoTreePricer::Quad)
              ? mast::quad_mast_solve(trees[0], trees[1])
              : mast::centroid_mast_solve(trees[0], trees[1]);
  } else {
    // t >= 3 (or non-empty blacklist): multi-tree DP.
    std::vector<W> taxon_weight(n + 1, 1.0);   // alpha, 1-based
    raw = mast::multi_mast_solve(trees, taxon_weight, blacklist);
  }

  // Return as (score, BitSet) preserving order. Downstream select_variables needs BitSet keys.
  std::vector<std::pair<W, BitSet>> out;
  out.reserve(raw.size());
  for (auto& s : raw) out.emplace_back(s.first, BitSet(s.second.begin(), s.second.end()));
  return out;
}

} // namespace maf
