// compute_cuts builds the set-packing conflict structure from the master LP solution and runs the
// clique separator. The production path realizes subset-row cuts through MPModel::add_cut_row.
#pragma once
#include "core.hpp"
#include "instance.hpp"
#include "maf_util.hpp"            // inner_vs
#include "cuts/set_packing_cuts.hpp"
#include <vector>
#include <utility>
#include <cstdlib>
#include <cstdio>

namespace maf {

// Env-gated diagnostic (MAF_CUT_DUMP): dump the separator's column order and emitted cuts per call,
// for trace comparison. No effect on the solver.
inline void _cut_dump(const std::vector<BitSet>& vars, const std::vector<std::vector<BitSet>>& cuts) {
  if (!std::getenv("MAF_CUT_DUMP")) return;
  static int call = 0; ++call;
  std::fprintf(stderr, "CCALL %d nvars=%d VARS", call, (int)vars.size());
  for (const BitSet& c : vars) { std::fprintf(stderr, " "); auto v = c.to_vec(); for (size_t k = 0; k < v.size(); ++k) std::fprintf(stderr, "%s%d", k ? "," : "", v[k]); }
  std::fprintf(stderr, "\n");
  for (const auto& cut : cuts) {
    std::fprintf(stderr, "  CCUT");
    for (const BitSet& col : cut) { std::fprintf(stderr, " {"); auto v = col.to_vec(); for (size_t k = 0; k < v.size(); ++k) std::fprintf(stderr, "%s%d", k ? "," : "", v[k]); std::fprintf(stderr, "}"); }
    std::fprintf(stderr, "\n");
  }
}

// Compute clique cuts from an ordered master LP solution.
// `lpsol` is the master LP solution as an ordered (column, value) list — column order defines the
// separator's variable index j (CUTS_SPEC HAZARD #1; we use the master's deterministic column order
// and trace-verify the emitted cuts).
// Returns (cuts, rhss): each cut is the list of clique-member columns (BitSets of leaves), rhs = 1.
inline std::pair<std::vector<std::vector<BitSet>>, std::vector<double>>
compute_cuts(const MAFP& maf, const std::vector<std::pair<BitSet, double>>& lpsol, int max_cuts) {
  int nvars = (int)lpsol.size();
  // CANONICAL COLUMN ORDER (lex by ascending leaves). This makes separator output independent of
  // hash/order noise and keeps cut generation deterministic.
  std::vector<std::pair<BitSet, double>> ordered(lpsol.begin(), lpsol.end());
  std::sort(ordered.begin(), ordered.end(),
            [](const std::pair<BitSet, double>& a, const std::pair<BitSet, double>& b) { return _lex_smaller(a.first, b.first); });
  std::vector<BitSet> vars(nvars);
  std::vector<double> lp(nvars);
  for (int j = 0; j < nvars; ++j) { vars[j] = ordered[j].first; lp[j] = ordered[j].second; }

  // ivs[i][j] = inner vertices of column j's Steiner tree in tree i.
  std::vector<std::vector<BitSet>> ivs(maf.T.size(), std::vector<BitSet>(nvars));
  for (size_t i = 0; i < maf.T.size(); ++i)
    for (int j = 0; j < nvars; ++j) {
      std::vector<int> iv = inner_vs(maf.T[i], vars[j]);
      ivs[i][j] = BitSet(iv.begin(), iv.end());
    }

  // cliques (set-packing rows): one per (tree i, internal vertex v), tree-major then vertex-ascending.
  // cliques[row] = { j : v in ivs[i][j] }  (columns whose Steiner tree in tree i uses vertex v).
  std::vector<BitSet> cliques;
  for (size_t i = 0; i < maf.T.size(); ++i)
    for (int v = maf.n + 1; v <= maf.T[i].nv(); ++v) {
      BitSet row;
      for (int j = 0; j < nvars; ++j) if (ivs[i][j].contains(v)) row.insert(j);
      cliques.push_back(std::move(row));
    }

  std::vector<BitSet> supports = separate_glpk_clique_cuts(lp, cliques, max_cuts);

  // Map each support (column indices, ascending) back to its member columns; rhs = 1.0.
  std::vector<std::vector<BitSet>> cuts;
  std::vector<double> rhss;
  for (const BitSet& sup : supports) {
    std::vector<BitSet> members;
    for (int j : sup) members.push_back(vars[j]);     // ascending j
    cuts.push_back(std::move(members));
    rhss.push_back(1.0);
  }
  _cut_dump(vars, cuts);
  return {cuts, rhss};
}

} // namespace maf
