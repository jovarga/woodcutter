// Undoable kernelization rules (Kelk-Linz-Meuwese order) wrapped in a SolverState, plus the
// reconstruction that lifts a reduced-instance partition back to the original leaf labels.
//
// Rules (all on the working MAFP):
//   PendantSubtree — common cherry: smaller-indexed leaf absorbs the larger (free, any t)
//   Chain43        — (4,3) common chain truncation (free, t=2 only)
//   ThreeTwo       — 3-2 near-cherry removal (k+=1, t=2 only)
//   Dummy          — never fires
#pragma once
#include "core.hpp"
#include "instance.hpp"
#include "preprocessing.hpp"     // restrict_to_leaves
#include <vector>
#include <tuple>
#include <array>
#include <optional>
#include <algorithm>
#include <stdexcept>

namespace maf {

enum class Rule { PendantSubtree, Chain43, ThreeTwo, Dummy };

struct ReductionFrame {
  MAFP maf;
  std::vector<int> leaf_map;
  size_t n_actions;
  int k_adjustment;
};

// SolverState. leaf_map[k-1] = original label of current position k.
// pending_actions: (is_singleton, orig_b, orig_a) — b removed; merged into a's comp (is_singleton
// false) or made its own singleton (true, orig_a unused = 0).  k_adjustment = #ThreeTwo firings.
struct SolverState {
  MAFP maf;
  std::vector<int> leaf_map;
  std::vector<std::tuple<bool, int, int>> pending_actions;
  int k_adjustment = 0;
  std::vector<ReductionFrame> undo_stack;

  explicit SolverState(const MAFP& m) : maf(m), leaf_map(m.n) {
    for (int k = 0; k < m.n; ++k) leaf_map[k] = k + 1;
  }
};

// ---- internal helpers ------------------------------------------------------------------------
inline void _push_frame(SolverState& s) {
  s.undo_stack.push_back(ReductionFrame{s.maf, s.leaf_map, s.pending_actions.size(), s.k_adjustment});
}

inline void _remove_leaf(SolverState& s, int pos) {
  std::vector<int> leaves;
  for (int k = 1; k <= s.maf.n; ++k) if (k != pos) leaves.push_back(k);
  auto [new_maf, lmap] = restrict_to_leaves(s.maf, leaves);   // lmap == leaves
  std::vector<int> new_leaf_map(new_maf.n);
  for (int k = 0; k < new_maf.n; ++k) new_leaf_map[k] = s.leaf_map[lmap[k] - 1];
  s.leaf_map = std::move(new_leaf_map);
  s.maf = std::move(new_maf);
}

// Sibling leaf of v if v is in a cherry, else 0.
inline int _leaf_sibling(const DiGraph& t, int v, int n) {
  if (t.in_degree(v) == 0) return 0;
  int par = t.in_neighbors(v)[0];
  const std::vector<int>& sibs = t.out_neighbors(par);
  if ((int)sibs.size() != 2) return 0;
  int sib = (sibs[0] == v) ? sibs[1] : sibs[0];
  return sib <= n ? sib : 0;
}

inline void undo_reduction(SolverState& s) {
  if (s.undo_stack.empty()) throw std::runtime_error("undo stack is empty");
  ReductionFrame f = std::move(s.undo_stack.back());
  s.undo_stack.pop_back();
  s.maf = std::move(f.maf);
  s.leaf_map = std::move(f.leaf_map);
  s.pending_actions.resize(f.n_actions);
  s.k_adjustment = f.k_adjustment;
}

// ---- PendantSubtreeRule ----------------------------------------------------------------------
inline std::optional<std::pair<int, int>> _find_pendant(const MAFP& maf) {
  if (maf.n < 2) return std::nullopt;
  for (int a = 1; a <= maf.n; ++a) {
    int sib = _leaf_sibling(maf.T[0], a, maf.n);
    if (sib == 0) continue;
    bool ok = true;
    for (size_t i = 1; i < maf.T.size(); ++i) if (_leaf_sibling(maf.T[i], a, maf.n) != sib) { ok = false; break; }
    if (!ok) continue;
    return std::make_pair(std::min(a, sib), std::max(a, sib));
  }
  return std::nullopt;
}

// ---- Chain43Rule -----------------------------------------------------------------------------
inline std::optional<std::array<int, 4>> _find_chain43(const MAFP& maf) {
  if (maf.T.size() != 2) return std::nullopt;
  if (maf.n < 4) return std::nullopt;
  const DiGraph& T1 = maf.T[0]; const DiGraph& T2 = maf.T[1];
  const int n = maf.n;
  auto lc = [n](const DiGraph& t, int v) { std::vector<int> r; for (int w : t.out_neighbors(v)) if (w <= n) r.push_back(w); return r; };
  auto only_inner = [n](const DiGraph& t, int v) -> int { for (int w : t.out_neighbors(v)) if (w > n) return w; return 0; };

  for (int x = 1; x <= n; ++x) {
    // Level 1: parent of x has exactly one leaf child (x) in both trees.
    int p1 = T1.in_neighbors(x)[0], p2 = T2.in_neighbors(x)[0];
    if ((int)lc(T1, p1).size() != 1) continue;
    if ((int)lc(T2, p2).size() != 1) continue;
    int at1 = only_inner(T1, p1), at2 = only_inner(T2, p2);
    if (at1 == 0 || at2 == 0) continue;

    // Level 2: one leaf child each, same in both trees.
    std::vector<int> lc2_1 = lc(T1, at1), lc2_2 = lc(T2, at2);
    if ((int)lc2_1.size() != 1 || (int)lc2_2.size() != 1) continue;
    if (lc2_1[0] != lc2_2[0]) continue;
    int chain2 = lc2_1[0];
    at1 = only_inner(T1, at1); at2 = only_inner(T2, at2);
    if (at1 == 0 || at2 == 0) continue;

    // Level 3: handle (1,1), (1,2), (2,1) leaf-child counts.
    std::vector<int> lc3_1 = lc(T1, at1), lc3_2 = lc(T2, at2);
    int nl1 = (int)lc3_1.size(), nl2 = (int)lc3_2.size();
    if (nl1 == 0 || nl2 == 0) continue;

    int chain3 = 0, next1 = 0, next2 = 0;
    if (nl1 == 1 && nl2 == 1) {
      if (lc3_1[0] != lc3_2[0]) continue;
      chain3 = lc3_1[0];
      next1 = only_inner(T1, at1); next2 = only_inner(T2, at2);
      if (next1 == 0 || next2 == 0) continue;
    } else if (nl1 == 1 && nl2 == 2) {                       // T2 cherry at level 3
      if (std::find(lc3_2.begin(), lc3_2.end(), lc3_1[0]) == lc3_2.end()) continue;
      chain3 = lc3_1[0];
      next1 = only_inner(T1, at1); next2 = at2;
      if (next1 == 0) continue;
    } else if (nl1 == 2 && nl2 == 1) {                       // T1 cherry at level 3
      if (std::find(lc3_1.begin(), lc3_1.end(), lc3_2[0]) == lc3_1.end()) continue;
      chain3 = lc3_2[0];
      next1 = at1; next2 = only_inner(T2, at2);
      if (next2 == 0) continue;
    } else {
      continue;
    }

    // Level 4: any leaf child of both next1 and next2.
    std::vector<int> lc4_1 = lc(T1, next1), lc4_2v = lc(T2, next2);
    BitSet lc4_2(lc4_2v.begin(), lc4_2v.end());
    for (int y : lc4_1) if (lc4_2.contains(y)) return std::array<int, 4>{x, chain2, chain3, y};
  }
  return std::nullopt;
}

// ---- ThreeTwoRule ----------------------------------------------------------------------------
inline std::optional<std::array<int, 3>> _find_three_two(const MAFP& maf) {
  if (maf.T.size() != 2) return std::nullopt;
  if (maf.n < 3) return std::nullopt;
  const DiGraph& T1 = maf.T[0]; const DiGraph& T2 = maf.T[1];
  for (int p = 1; p <= maf.n; ++p) {
    int sib1 = _leaf_sibling(T1, p, maf.n);
    int sib2 = _leaf_sibling(T2, p, maf.n);
    if (sib1 == 0 || sib2 == 0) continue;
    if (sib1 == sib2) continue;                              // common cherry -> PendantSubtree

    // Attempt A: (p,sib1) cherry in T1; T2 has ((p,sib2),sib1) -> delete sib2.
    int ppar2 = T2.in_neighbors(p)[0];
    if (T2.in_degree(ppar2) >= 1) {
      int pgpar2 = T2.in_neighbors(ppar2)[0];
      int sib1par2 = T2.in_neighbors(sib1)[0];
      if (pgpar2 == sib1par2) return std::array<int, 3>{p, sib1, sib2};
    }
    // Attempt B: (p,sib2) cherry in T2; T1 has ((p,sib1),sib2) -> delete sib1.
    int ppar1 = T1.in_neighbors(p)[0];
    if (T1.in_degree(ppar1) >= 1) {
      int pgpar1 = T1.in_neighbors(ppar1)[0];
      int sib2par1 = T1.in_neighbors(sib2)[0];
      if (pgpar1 == sib2par1) return std::array<int, 3>{p, sib2, sib1};
    }
  }
  return std::nullopt;
}

// apply_reduction!(state, rule) -> fired?
inline bool apply_reduction(SolverState& s, Rule rule) {
  switch (rule) {
    case Rule::PendantSubtree: {
      auto m = _find_pendant(s.maf);
      if (!m) return false;
      int a = m->first, b = m->second;                       // a<b; b removed, merged into a
      _push_frame(s);
      s.pending_actions.emplace_back(false, s.leaf_map[b - 1], s.leaf_map[a - 1]);
      _remove_leaf(s, b);
      return true;
    }
    case Rule::Chain43: {
      if (s.maf.T.size() != 2) return false;
      auto m = _find_chain43(s.maf);
      if (!m) return false;
      int x = (*m)[0], y = (*m)[3];                          // y (4th chain elt) merged into x
      _push_frame(s);
      s.pending_actions.emplace_back(false, s.leaf_map[y - 1], s.leaf_map[x - 1]);
      _remove_leaf(s, y);
      return true;
    }
    case Rule::ThreeTwo: {
      if (s.maf.T.size() != 2) return false;
      auto m = _find_three_two(s.maf);
      if (!m) return false;
      int r = (*m)[2];                                       // r removed -> singleton component
      _push_frame(s);
      s.pending_actions.emplace_back(true, s.leaf_map[r - 1], 0);
      s.k_adjustment += 1;
      _remove_leaf(s, r);
      return true;
    }
    case Rule::Dummy:
    default:
      return false;
  }
}

// reduce_exhaustive!(state, sequence) -> #reductions.  Restart from the front of the sequence on
// every firing, until a full pass produces nothing.
inline int reduce_exhaustive(SolverState& s,
    const std::vector<Rule>& sequence = {Rule::PendantSubtree, Rule::Chain43, Rule::ThreeTwo, Rule::PendantSubtree}) {
  int total = 0;
  while (true) {
    bool fired = false;
    for (Rule rule : sequence) {
      if (apply_reduction(s, rule)) { fired = true; ++total; break; }
    }
    if (!fired) break;
  }
  return total;
}

// Lift a reduced partition back to original labels. Re-insert pending actions in REVERSE order.
inline std::vector<BitSet> reconstruct(const SolverState& s, const std::vector<BitSet>& sol) {
  std::vector<BitSet> expanded;
  for (const BitSet& comp : sol) {
    BitSet bs;
    for (int k : comp) bs.insert(s.leaf_map[k - 1]);
    expanded.push_back(std::move(bs));
  }
  for (auto it = s.pending_actions.rbegin(); it != s.pending_actions.rend(); ++it) {
    auto [is_singleton, orig_b, orig_a] = *it;
    if (is_singleton) {
      expanded.push_back(BitSet{orig_b});
    } else {
      int idx = -1;
      for (int c = 0; c < (int)expanded.size(); ++c) if (expanded[c].contains(orig_a)) { idx = c; break; }
      if (idx < 0) throw std::runtime_error("pending merge: merge target not found in solution");
      expanded[idx].insert(orig_b);
    }
  }
  return expanded;
}

} // namespace maf
