// mast/quad_dp.hpp  —  plain O(n^2) 2-tree MAST pricer (replaces the centroid engine).
//
// The 2-tree pricing subproblem: given two canonical rooted binary PhyloTrees with internal weights
// (iw = -beta [+ branch dual], possibly POSITIVE after branching), return for every rooted vertex the
// best agreement subtree rooted (by LCA) there, as (score, sorted-taxa) pairs — T1 slots ascending
// vertex id, then T2 slots ascending vertex id, filtered to |taxa| > 1.  Same output contract and
// same per-vertex max SCORE as multi_dp/centroid (=> same LP optimum, same k*); the tie
// REPRESENTATIVE may differ but is fully deterministic (structural branch priority + fixed FP
// association) for deterministic floating-point behavior. Design: docs/QUAD_DP_SPEC.md.
//
// Three tables over pairs (u in T1, v in T2):
//   E[u][v]   = best score of an agreement subtree S with LCA_1(S)=u AND LCA_2(S)=v exactly.
//   HR2[u][v] = best over v'<=v of attach2(v,v') + E[u][v']         (peel tree 2, tree 1 fixed).
//   H[u][v]   = best over u'<=u,v'<=v of attach1(u,u')+attach2(v,v')+E[u'][v']  (attach below both).
// Score formula: |S| + sum over BOTH trees of internal_weight on each Steiner-internal vertex (once).
#pragma once
#include "trees.hpp"
#include <vector>
#include <utility>

namespace maf::mast {

struct QuadDP {
  const PhyloTree& t1;
  const PhyloTree& t2;
  int nv1, nv2, W2;
  std::vector<W> E, H, HR2;
  std::vector<signed char> Eb, Hb, HR2b;   // backpointers (branch index that won)

  QuadDP(const PhyloTree& a, const PhyloTree& b)
      : t1(a), t2(b), nv1(a.nv()), nv2(b.nv()), W2(b.nv() + 1) {
    const W NEG = invalid_score();
    const size_t sz = (size_t)(nv1 + 1) * W2;
    E.assign(sz, NEG); H.assign(sz, NEG); HR2.assign(sz, NEG);
    Eb.assign(sz, 0); Hb.assign(sz, 0); HR2b.assign(sz, 0);
    forward();
  }

  inline size_t idx(int u, int v) const { return (size_t)u * W2 + v; }

  void forward() {
    const W NEG = invalid_score();
    const std::vector<W>& iw1 = t1.internal_weight;
    const std::vector<W>& iw2 = t2.internal_weight;
    for (int u = 1; u <= nv1; ++u) {
      const bool uleaf = t1.isleaf(u);
      const int uL = t1.left[u], uR = t1.right[u];
      for (int v = 1; v <= nv2; ++v) {
        const bool vleaf = t2.isleaf(v);
        const int vL = t2.left[v], vR = t2.right[v];

        // --- E[u][v] ---  (branch 0 = parallel, 1 = cross; lowest index wins ties, strict >)
        W e = NEG; signed char eb = 0;
        if (uleaf && vleaf) {
          if (u == v) e = 1.0;                          // singleton {u}, alpha = 1
        } else if (!uleaf && !vleaf) {
          const W hLL = H[idx(uL, vL)], hRR = H[idx(uR, vR)];
          const W hLR = H[idx(uL, vR)], hRL = H[idx(uR, vL)];
          const W par = (is_valid_score(hLL) && is_valid_score(hRR)) ? hLL + hRR : NEG;
          const W cro = (is_valid_score(hLR) && is_valid_score(hRL)) ? hLR + hRL : NEG;
          W best = NEG; signed char bb = 0;
          if (is_valid_score(par) && par > best) { best = par; bb = 0; }
          if (is_valid_score(cro) && cro > best) { best = cro; bb = 1; }
          if (is_valid_score(best)) { e = (iw1[u] + iw2[v]) + best; eb = bb; }
        }
        E[idx(u, v)] = e; Eb[idx(u, v)] = eb;

        // --- HR2[u][v] ---  (branch 0 = E[u][v], 1 = descend vL, 2 = descend vR)
        W hr2 = e; signed char hr2b = 0;
        if (!vleaf) {
          const W d1 = HR2[idx(u, vL)], d2 = HR2[idx(u, vR)];
          if (is_valid_score(d1)) { W c = iw2[v] + d1; if (!is_valid_score(hr2) || c > hr2) { hr2 = c; hr2b = 1; } }
          if (is_valid_score(d2)) { W c = iw2[v] + d2; if (!is_valid_score(hr2) || c > hr2) { hr2 = c; hr2b = 2; } }
        }
        HR2[idx(u, v)] = hr2; HR2b[idx(u, v)] = hr2b;

        // --- H[u][v] ---  (branch 0 = HR2[u][v], 1 = descend uL, 2 = descend uR)
        W h = hr2; signed char hb = 0;
        if (!uleaf) {
          const W d1 = H[idx(uL, v)], d2 = H[idx(uR, v)];
          if (is_valid_score(d1)) { W c = iw1[u] + d1; if (!is_valid_score(h) || c > h) { h = c; hb = 1; } }
          if (is_valid_score(d2)) { W c = iw1[u] + d2; if (!is_valid_score(h) || c > h) { h = c; hb = 2; } }
        }
        H[idx(u, v)] = h; Hb[idx(u, v)] = hb;
      }
    }
  }

  // ---- taxa reconstruction (walk backpointers; O(|S|) per emitted slot) ----------------------
  std::vector<int> reconE(int u, int v) const {
    if (t1.isleaf(u) && t2.isleaf(v)) return {u};                       // singleton {taxon u}
    const int uL = t1.left[u], uR = t1.right[u], vL = t2.left[v], vR = t2.right[v];
    if (Eb[idx(u, v)] == 0) return sorted_taxa_union(reconH(uL, vL), reconH(uR, vR));  // parallel
    else                    return sorted_taxa_union(reconH(uL, vR), reconH(uR, vL));  // cross
  }
  std::vector<int> reconH(int u, int v) const {
    if (t1.isleaf(u)) return reconHR2(u, v);
    switch (Hb[idx(u, v)]) {
      case 1:  return reconH(t1.left[u], v);
      case 2:  return reconH(t1.right[u], v);
      default: return reconHR2(u, v);
    }
  }
  std::vector<int> reconHR2(int u, int v) const {
    if (t2.isleaf(v)) return reconE(u, v);
    switch (HR2b[idx(u, v)]) {
      case 1:  return reconHR2(u, t2.left[v]);
      case 2:  return reconHR2(u, t2.right[v]);
      default: return reconE(u, v);
    }
  }

  // ---- emit per-rooted-vertex bests: T1 block (ascending u), then T2 block (ascending v) -----
  std::vector<std::pair<W, std::vector<int>>> solve() {
    std::vector<std::pair<W, std::vector<int>>> out;
    for (int u = 1; u <= nv1; ++u) {
      W best = invalid_score(); int argv = 0;
      for (int v = 1; v <= nv2; ++v) { W e = E[idx(u, v)]; if (is_valid_score(e) && e > best) { best = e; argv = v; } }
      if (is_valid_score(best)) {
        std::vector<int> taxa = reconE(u, argv);
        if (taxa.size() > 1) out.emplace_back(best, std::move(taxa));
      }
    }
    for (int v = 1; v <= nv2; ++v) {
      W best = invalid_score(); int argu = 0;
      for (int u = 1; u <= nv1; ++u) { W e = E[idx(u, v)]; if (is_valid_score(e) && e > best) { best = e; argu = u; } }
      if (is_valid_score(best)) {
        std::vector<int> taxa = reconE(argu, v);
        if (taxa.size() > 1) out.emplace_back(best, std::move(taxa));
      }
    }
    return out;
  }
};

// quad_mast_solve(t1, t2): the drop-in replacement for centroid_mast_solve.
inline std::vector<std::pair<W, std::vector<int>>>
quad_mast_solve(const PhyloTree& t1, const PhyloTree& t2) {
  QuadDP dp(t1, t2);
  return dp.solve();
}

} // namespace maf::mast
