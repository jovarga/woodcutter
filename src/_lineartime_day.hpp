// _lineartime_day.hpp — linear-time reductions + cluster decomposition (NOT the exact track).
//
// Linear-time (O(t*n)) reduction + cluster-decomposition drivers.
//
// These replace the exact solver's O(n^2) cluster_sequence (minimal_clusters/compute_descendants)
// and its restrict_to_leaves-based reduce_exhaustive with Day's (1985) O(t*n) common-cluster scan
// and a mutable-forest worklist.  They only pay off on the LARGE instances that appear in the
// heuristic and lower-bound tracks.  They produce the SAME encodings the production reconstruction
// consumes (ClusterEntry / ClusterChoice / reconstruct_solution and SolverState / reconstruct), so
// they are drop-in.
//
// NAMING / TRUST: this is deliberately NOT named "heuristic_only".  The LOWER-BOUND track relies on
// these for its optimality guarantee, so their correctness matters and must be believed.  The EXACT
// track does NOT use them — it keeps its own trivial O(n^2) reductions/clusters (preprocessing.hpp /
// reduction_rule_stack.hpp), already trusted — so any doubt about the correctness of this file can
// NEVER cause an exact-track disqualification.  The heuristic track uses them purely for speed (no
// correctness requirement there).  Namespace: maf::lineartime.
#pragma once
#include "instance.hpp"
#include "preprocessing.hpp"          // ClusterEntry, ClusterChoice, restrict_to_leaves, add_pendant_leaf,
                                      // collapse_to_leaf, active_leaf_positions, solve_restricted_instance,
                                      // reconstruct_solution, ClusterSolver
#include "reduction_rule_stack.hpp"   // SolverState, reconstruct
#include <vector>
#include <set>
#include <algorithm>
#include <limits>
#include <tuple>
#include <cstdlib>

namespace maf { namespace lineartime {

// ---- Day's linear common-cluster scan (from _cluster_family) ---------------------------------

// _leaf_positions(tree, n): position of each leaf 1..n by an iterative DFS of `tree` (root = nv).
inline std::vector<int> day_leaf_positions(const DiGraph& tree, int n) {
  std::vector<int> pos(n + 1, 0);
  int cnt = 0;
  std::vector<int> stack{tree.nv()};
  while (!stack.empty()) {
    int v = stack.back(); stack.pop_back();
    if (v <= n) pos[v] = ++cnt;
    else for (int c : tree.out_neighbors(v)) stack.push_back(c);
  }
  return pos;
}

// _subtree_intervals(tree, n, pos): (lo,hi,size) of every node's leaf set under `pos`. One forward
// pass works because vertices are in reverse-topological order (every child id < its parent's).
inline void day_subtree_intervals(const DiGraph& tree, int n, const std::vector<int>& pos,
                                  std::vector<int>& lo, std::vector<int>& hi, std::vector<int>& sz) {
  int N = tree.nv();
  lo.assign(N + 1, 0); hi.assign(N + 1, 0); sz.assign(N + 1, 0);
  for (int v = 1; v <= N; ++v) {
    if (v <= n) { lo[v] = hi[v] = pos[v]; sz[v] = 1; }
    else {
      int l = std::numeric_limits<int>::max(), h = 0, s = 0;
      for (int c : tree.out_neighbors(v)) { l = std::min(l, lo[c]); h = std::max(h, hi[c]); s += sz[c]; }
      lo[v] = l; hi[v] = h; sz[v] = s;
    }
  }
}

inline long long day_interval_key(int lo, int hi, int n) { return (long long)(lo - 1) * (n + 1) + hi; }

// _cluster_family(maf): common clusters + laminar family in one pass. Fills:
//   iscommon[v]  — T1 clade v is common to all trees,
//   children     — node -> immediate family-children (common nodes whose nearest common ancestor is v),
//   lo1, sz1     — T1 subtree min-position / size,
//   leaf_at      — position -> leaf label,  root — T1 root vertex.
struct DayFamily {
  std::vector<char> iscommon;                 // over 1..nv(T1)
  std::vector<std::vector<int>> children;     // index 1..nv(T1)
  std::vector<int> lo1, sz1, leaf_at;
  int root = 0;
};

inline DayFamily day_cluster_family(const MAFP& maf) {
  int n = maf.n; int t = (int)maf.T.size();
  const DiGraph& T1 = maf.T[0];
  std::vector<int> pos = day_leaf_positions(T1, n);
  std::vector<int> lo1, hi1, sz1;
  day_subtree_intervals(T1, n, pos, lo1, hi1, sz1);

  std::vector<int> leaf_at(n + 1, 0);
  for (int leaf = 1; leaf <= n; ++leaf) leaf_at[pos[leaf]] = leaf;

  std::set<long long> common;
  for (int v = 1; v <= T1.nv(); ++v) common.insert(day_interval_key(lo1[v], hi1[v], n));
  for (int k = 1; k < t; ++k) {
    std::vector<int> lok, hik, szk;
    day_subtree_intervals(maf.T[k], n, pos, lok, hik, szk);
    std::set<long long> sk;
    for (int v = 1; v <= maf.T[k].nv(); ++v)
      if (szk[v] == hik[v] - lok[v] + 1) sk.insert(day_interval_key(lok[v], hik[v], n));
    std::set<long long> inter;
    for (long long key : common) if (sk.count(key)) inter.insert(key);
    common.swap(inter);
  }

  DayFamily fam;
  fam.iscommon.assign(T1.nv() + 1, 0);
  for (int v = 1; v <= T1.nv(); ++v)
    fam.iscommon[v] = common.count(day_interval_key(lo1[v], hi1[v], n)) ? 1 : 0;

  fam.root = T1.nv();
  fam.children.assign(T1.nv() + 1, {});
  std::vector<std::pair<int,int>> stack{{fam.root, 0}};   // (vertex, nearest common ancestor)
  while (!stack.empty()) {
    auto [v, cca] = stack.back(); stack.pop_back();
    int nc = cca;
    if (fam.iscommon[v]) { if (cca != 0) fam.children[cca].push_back(v); nc = v; }
    for (int c : T1.out_neighbors(v)) stack.push_back({c, nc});
  }
  fam.lo1 = std::move(lo1); fam.sz1 = std::move(sz1); fam.leaf_at = std::move(leaf_at);
  return fam;
}

// cluster_sequence_linear(maf) -> (seq, main_maf, main_label_map), same encoding as cluster_sequence:
// leaf_map[k] = positive original label, or -j (representative of seq[j], 1-based), or 0 (rho, last).
inline void cluster_sequence_linear(const MAFP& maf,
                                    std::vector<ClusterEntry>& seq, MAFP& main_maf,
                                    std::vector<int>& main_label_map) {
  int n = maf.n;
  DayFamily fam = day_cluster_family(maf);
  seq.clear();

  // Non-trivial common clusters (size 2..n-1), sorted ascending by size => children before parents.
  std::vector<int> cluster_nodes;
  for (int v = 1; v < (int)fam.iscommon.size(); ++v)
    if (fam.iscommon[v] && fam.sz1[v] >= 2 && fam.sz1[v] <= n - 1) cluster_nodes.push_back(v);
  std::stable_sort(cluster_nodes.begin(), cluster_nodes.end(),
                   [&](int a, int b) { return fam.sz1[a] < fam.sz1[b]; });
  std::vector<int> cluster_idx(fam.iscommon.size(), 0);   // cluster node -> 1-based seq index
  for (int i = 0; i < (int)cluster_nodes.size(); ++i) cluster_idx[cluster_nodes[i]] = i + 1;

  // family-child node -> (representative leaf, encoded label)
  auto child_rep_enc = [&](int u) -> std::pair<int,int> {
    if (u <= n) return {u, u};                            // leaf child: represents itself
    return {fam.leaf_at[fam.lo1[u]], -cluster_idx[u]};    // sub-cluster: leftmost leaf, encoded -idx
  };
  // Build the subinstance for a block from its family-children; with_rho adds the rho pendant leaf.
  auto build_block = [&](const std::vector<int>& kids, bool with_rho,
                         MAFP& out_maf, std::vector<int>& out_labels) {
    std::vector<std::pair<int,int>> reps_enc;
    for (int u : kids) reps_enc.push_back(child_rep_enc(u));
    std::sort(reps_enc.begin(), reps_enc.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
    std::vector<int> reps; reps.reserve(reps_enc.size());
    std::map<int,int> enc_of;
    for (auto& re : reps_enc) { reps.push_back(re.first); enc_of[re.first] = re.second; }
    auto [sub, posmap] = restrict_to_leaves(maf, reps);   // posmap[k] = reps[k]
    std::vector<int> labels; labels.reserve(sub.n);
    for (int k = 0; k < sub.n; ++k) labels.push_back(enc_of[posmap[k]]);
    if (with_rho) { out_maf = add_pendant_leaf(sub); labels.push_back(0); out_labels = std::move(labels); }
    else          { out_maf = std::move(sub);        out_labels = std::move(labels); }
  };

  for (int i = 0; i < (int)cluster_nodes.size(); ++i) {
    MAFP sub_rho; std::vector<int> lmap;
    build_block(fam.children[cluster_nodes[i]], /*with_rho=*/true, sub_rho, lmap);
    seq.push_back(ClusterEntry{std::move(sub_rho), std::move(lmap), -(i + 1)});
  }

  const std::vector<int>& rootkids = fam.children[fam.root];
  if (rootkids.empty()) {
    main_maf = maf;
    main_label_map.resize(n);
    for (int k = 0; k < n; ++k) main_label_map[k] = k + 1;
    return;
  }
  build_block(rootkids, /*with_rho=*/false, main_maf, main_label_map);
}

// solve_clusters_linear(maf, solver [, defer_second]): the LS11 minimum-weight-forest procedure driven
// by the linear cluster_sequence_linear. Mirrors preprocessing.hpp::solve_clusters (option A vs B
// independent), only the sequence source differs. Returns (partition, dual).
//
// Option A (keep ρ) is ALWAYS the locally-optimal cluster resolution: cost_A* <= cost_B*, because any
// B-forest plus ρ as a singleton is a valid A-forest (k_A* <= k_B*+1 => cost_A*=k_A*-1 <= k_B*=cost_B*).
// So we ALWAYS dedicate A first and take it as the incumbent choice; option B is preferred only for the
// DOWNSTREAM savings it unlocks (ρ singleton => aᵢ singleton in the main instance).  defer_mode selects
// whether/when to also dedicate option B:
//   * 0 = BOTH (default; exact / lower-bound): solve B for EVERY cluster — its dual is needed for the
//         sound min(dual_A,dual_B) composed lower bound, and its savings for best quality.
//   * 1 = A-ONLY: never solve B (heuristic FIRST pass — a fast, guaranteed-complete incumbent).
//   * 2 = ADAPTIVE WINDOW (heuristic refinement passes): solve B only while the remaining budget
//         comfortably covers finishing the rest of the pass.  The reserve is MACHINE-INDEPENDENT:
//             K = SAFETY * avg_time_per_cluster_measured_this_pass * (remaining_clusters + MAIN_WEIGHT)
//         so it auto-scales with instance SIZE (bigger blocks -> slower clusters -> larger K) and the
//         cluster-cascade COUNT (remaining_clusters), and transfers across machines (a raw seconds
//         constant would not).  When tight we keep the locally-optimal A and skip B, so the pass still
//         finishes and registers an improvement instead of being lost to check_budget().
// Bridge for PROGRESSIVE incumbent registration out of the per-block leaf solver.  ClusterSolver has a
// fixed signature (solver(sub_maf) -> (sol,dual)) and lives in the protected preprocessing.hpp, so the
// leaf cannot take a `consider` directly.  Instead the leaf (run_reduced_once) forwards every block-level
// improvement to THIS callback, and solve_clusters_linear swaps in a per-part closure that COMPOSES the
// block partition into a full (reduced-instance) forest and forwards it to the caller.  Empty => off.
// File-scope + single-threaded solve => no re-entrancy.
inline std::function<void(const std::vector<BitSet>&)>& cluster_progress_bridge() {
  static std::function<void(const std::vector<BitSet>&)> f;
  return f;
}

// Only clusters at least this large get a during-solve progress bridge (small clusters finish fast; a
// bridge on each would add O(#clusters^2) compose overhead for no SIGTERM benefit).  Env-overridable.
inline int cluster_progress_min_block() {
  static const int m = [](){ const char* e=std::getenv("MAF_PROG_MIN_BLOCK"); return (e&&*e)?std::atoi(e):256; }();
  return m;
}

// solve_clusters_linear(maf, solver, defer_mode, consider): as above, plus optional PROGRESSIVE
// registration.  `consider` (if set) receives full reduced-instance forests composed from every
// block-level improvement (during the whole-instance solve when there are no clusters, during the main
// core solve, and during any sufficiently large cluster solve).  All compositions are try/catch-guarded
// and the caller re-validates, so a malformed provisional forest can never crash or be published.
inline std::pair<std::vector<BitSet>, double>
solve_clusters_linear(const MAFP& maf, const ClusterSolver& solver, int defer_mode = 0,
                      const std::function<void(const std::vector<BitSet>&)>& consider = {}) {
  std::vector<ClusterEntry> seq; MAFP main_maf; std::vector<int> main_label_map;
  cluster_sequence_linear(maf, seq, main_maf, main_label_map);
  // Always drop the bridge on ANY exit (incl. an Interrupted throw mid-solve) so no stale closure with
  // dangling captures survives.  Clearing a std::function only destroys the closure (never dereferences
  // its captures), so this is safe even during unwinding.
  struct BridgeClearer { ~BridgeClearer(){ cluster_progress_bridge() = {}; } } _bridge_clearer;
  if (seq.empty()) {
    // No clusters: the block IS the whole (reduced) instance, so a block partition is already a full
    // forest — register each improvement directly (identity composition).
    if (consider) cluster_progress_bridge() = [&](const std::vector<BitSet>& b){ consider(b); };
    auto r = solver(maf);
    cluster_progress_bridge() = {};
    return r;
  }

  const int n = maf.n;
  auto clear_bridge = []{ cluster_progress_bridge() = {}; };
  // enc_map for a part: solve_restricted_instance uses restrict_to_leaves(maf, active) whose pos_map ==
  // the (sorted) active positions, so enc_map[k] == full_leaf_map[active[k]-1].  Matches exactly.
  auto enc_of = [](const std::vector<int>& full_map, const std::vector<int>& active) {
    std::vector<int> e; e.reserve(active.size());
    for (int p : active) e.push_back(full_map[p - 1]);
    return e;
  };
  auto trivial_choice = [](const std::vector<int>& lm) {
    std::vector<BitSet> s; s.reserve(lm.size());
    for (int k = 1; k <= (int)lm.size(); ++k) s.push_back(BitSet{k});
    return ClusterChoice{std::move(s), lm};
  };
  std::vector<BitSet> trivial_main;                       // all main leaves as singletons (provisional main)
  for (int k = 1; k <= main_maf.n; ++k) trivial_main.push_back(BitSet{k});

  // Dimensionless window tunables (env-overridable for tuning; baked defaults otherwise).
  // TODO(retune-300s): 1.3 / 15.0 were set from a 120 s proxy sweep (second-order effect); re-confirm at
  // the real 300 s budget that no large cluster-heavy instance fails to complete a pass.  See docs/TODO.md.
  static const double DEFER_SAFETY = [](){ const char* e=std::getenv("MAF_DEFER_SAFETY"); return (e&&*e)?std::atof(e):1.3; }();
  static const double DEFER_MAINW  = [](){ const char* e=std::getenv("MAF_DEFER_MAINW");  return (e&&*e)?std::atof(e):15.0; }();

  std::set<int> rho_singleton;
  std::vector<ClusterChoice> cluster_choices(seq.size());
  std::vector<double> cluster_duals(seq.size());
  std::vector<int> cluster_costs(seq.size());

  const double pass_start_rem = remaining_seconds();
  for (int i = 0; i < (int)seq.size(); ++i) {
    check_budget();
    const MAFP& sub_maf = seq[i].sub_maf;
    const std::vector<int>& leaf_map = seq[i].leaf_map;

    // Provisional-composition bridge for a large cluster's solve: compose its in-progress block sol with
    // the clusters resolved so far + all-singleton fallbacks for the rest + a singleton main.  Valid but
    // deliberately poor for the unsolved remainder; the caller's guard keeps only genuine improvements.
    const bool prog_here = (bool)consider && sub_maf.n >= cluster_progress_min_block();
    auto set_cluster_bridge = [&](const std::vector<int>& active) {
      if (!prog_here) return;
      std::vector<int> enc_i = enc_of(leaf_map, active);
      int ii = i;
      cluster_progress_bridge() = [&, ii, enc_i](const std::vector<BitSet>& blk) {
        try {
          std::vector<ClusterChoice> ctx(seq.size());
          for (int j = 0; j < (int)seq.size(); ++j) {
            if      (j < ii)  ctx[j] = cluster_choices[j];
            else if (j == ii) ctx[j] = ClusterChoice{blk, enc_i};
            else              ctx[j] = trivial_choice(seq[j].leaf_map);
          }
          consider(reconstruct_solution(ctx, trivial_main, main_label_map, n));
        } catch (...) {}
      };
    };

    // Dedicate option A (keep ρ) first — always locally optimal — and register it as the choice.
    std::vector<int> active_A = active_leaf_positions(leaf_map, rho_singleton, /*keep_rho=*/true);
    set_cluster_bridge(active_A);
    auto [sol_A, dual_A, map_A] = solve_restricted_instance(sub_maf, leaf_map, active_A, solver);
    clear_bridge();
    int cost_A = (int)sol_A.size() - 1;
    double dual_cost_A = std::max(std::ceil(dual_A - 1e-8) - 1.0, 0.0);
    cluster_choices[i] = ClusterChoice{sol_A, map_A};
    cluster_duals[i] = dual_cost_A; cluster_costs[i] = cost_A;

    // Decide whether to also dedicate option B (see defer_mode above).
    bool solve_b;
    if (defer_mode == 1)      solve_b = false;                // A-only pass
    else if (defer_mode == 2) {                               // adaptive machine-independent window
      double elapsed = pass_start_rem - remaining_seconds();
      double avg     = (i > 0) ? elapsed / (double)i : 0.0;  // measured time per cluster so far
      double K       = DEFER_SAFETY * avg * ((double)((int)seq.size() - i) + DEFER_MAINW);
      solve_b = (remaining_seconds() > K);
    } else                    solve_b = true;                 // 0 = BOTH (exact / lower-bound)

    if (solve_b) {
      std::vector<int> active_B = active_leaf_positions(leaf_map, rho_singleton, /*keep_rho=*/false);
      set_cluster_bridge(active_B);
      auto [sol_B, dual_B, map_B] = solve_restricted_instance(sub_maf, leaf_map, active_B, solver);
      clear_bridge();
      int cost_B = (int)sol_B.size();
      double dual_cost_B = std::max(std::ceil(dual_B - 1e-8), 0.0);
      cluster_duals[i] = std::min(dual_cost_A, dual_cost_B);   // valid composed global lower bound
      if (cost_B <= cost_A) {                                  // B wins (ties -> B: ρ_i singleton)
        cluster_choices[i] = ClusterChoice{sol_B, map_B};
        cluster_costs[i] = cost_B;
        rho_singleton.insert(i + 1);
      }
    }
  }

  std::vector<int> active_main = active_leaf_positions(main_label_map, rho_singleton, /*keep_rho=*/false);
  // Main-core bridge: cluster_choices are now COMPLETE, so every in-progress main partition composes into
  // a fully valid forest — this is the dominant progressive-registration case (the core is usually the
  // slowest part).  enc_of(main_label_map, active_main) == the main_map solve_restricted_instance returns.
  if (consider) {
    std::vector<int> enc_main = enc_of(main_label_map, active_main);
    cluster_progress_bridge() = [&, enc_main](const std::vector<BitSet>& blk) {
      try { consider(reconstruct_solution(cluster_choices, blk, enc_main, n)); } catch (...) {}
    };
  }
  auto [sol_main, dual_main, main_map] = solve_restricted_instance(main_maf, main_label_map, active_main, solver);
  clear_bridge();

  std::vector<BitSet> sol = reconstruct_solution(cluster_choices, sol_main, main_map, maf.n);
  double dual = 0.0; for (double d : cluster_duals) dual += d;
  dual += std::max(std::ceil(dual_main - 1e-8), 0.0);
  return {sol, std::max(dual, 0.0)};
}

// ---- Linear (O(t*n)) reduction rules (from reduce_linear) ------------------------------------

// Mutable binary forest: leaf removal + parent contraction is O(t), no O(n) restrict_to_leaves.
// Leaf ids 1..n are the ORIGINAL labels and persist (only leaves are removed).
struct MForest {
  std::vector<std::vector<int>> par, c1, c2;   // per tree, indexed by vertex id
  std::vector<int> root;
  int n = 0;
  std::vector<char> alive;                     // per leaf 1..n
  int t() const { return (int)par.size(); }
};

inline MForest make_mforest(const MAFP& maf) {
  int t = (int)maf.T.size(); int n = maf.n;
  MForest F; F.n = n; F.alive.assign(n + 1, 1);
  F.par.resize(t); F.c1.resize(t); F.c2.resize(t); F.root.assign(t, 0);
  for (int i = 0; i < t; ++i) {
    const DiGraph& T = maf.T[i]; int N = T.nv();
    F.par[i].assign(N + 1, 0); F.c1[i].assign(N + 1, 0); F.c2[i].assign(N + 1, 0);
    for (int v = 1; v <= N; ++v) {
      const auto& ins = T.in_neighbors(v);
      if (ins.empty()) { F.root[i] = v; F.par[i][v] = 0; } else F.par[i][v] = ins[0];
      const auto& outs = T.out_neighbors(v);
      if (!outs.empty()) { F.c1[i][v] = outs[0]; F.c2[i][v] = outs.size() >= 2 ? outs[1] : 0; }
    }
  }
  return F;
}

// leaf children of node u in tree i: sets (a,b,count).
inline void mf_leafkids(const MForest& F, int i, int u, int& a, int& b, int& cnt) {
  int x = F.c1[i][u], y = F.c2[i][u];
  bool lx = (x != 0 && x <= F.n), ly = (y != 0 && y <= F.n);
  if (lx && ly) { a = x; b = y; cnt = 2; }
  else if (lx)  { a = x; b = 0; cnt = 1; }
  else if (ly)  { a = y; b = 0; cnt = 1; }
  else          { a = 0; b = 0; cnt = 0; }
}
inline int mf_innerkid(const MForest& F, int i, int u) {
  int x = F.c1[i][u], y = F.c2[i][u];
  if (x != 0 && x > F.n) return x;
  if (y != 0 && y > F.n) return y;
  return 0;
}
inline int mf_sibleaf(const MForest& F, int i, int v) {   // sibling leaf of v, or 0
  int p = F.par[i][v]; if (p == 0) return 0;
  int s = F.c1[i][p] == v ? F.c2[i][p] : F.c1[i][p];
  return (s != 0 && s <= F.n) ? s : 0;
}
// Remove leaf b in every tree (contract its parent); fills per-tree (sibling, grandparent).
inline void mf_remove(MForest& F, int b, std::vector<std::pair<int,int>>& res) {
  res.assign(F.t(), {0, 0});
  for (int i = 0; i < F.t(); ++i) {
    int p = F.par[i][b];
    int s = F.c1[i][p] == b ? F.c2[i][p] : F.c1[i][p];
    int g = F.par[i][p];
    if (g == 0) { F.root[i] = s; F.par[i][s] = 0; }
    else { F.par[i][s] = g; if (F.c1[i][g] == p) F.c1[i][g] = s; else F.c2[i][g] = s; }
    res[i] = {s, g};
  }
  F.alive[b] = 0;
}

// common cherry: leaf a and its T1-sibling are siblings in EVERY tree -> (remove=max, keep=min).
inline bool mf_cherry(const MForest& F, int a, int& rem, int& keep) {
  int sib = mf_sibleaf(F, 0, a); if (sib == 0) return false;
  for (int i = 1; i < F.t(); ++i) if (mf_sibleaf(F, i, a) != sib) return false;
  rem = std::max(a, sib); keep = std::min(a, sib); return true;
}

// (4,3) common chain truncation rooted at chain-top x (t==2). Returns the 4th taxon y, or 0.
inline int mf_chain43(const MForest& F, int x) {
  int p1 = F.par[0][x], p2 = F.par[1][x];
  if (p1 == 0 || p2 == 0) return 0;
  int a, b, k1, k2;
  mf_leafkids(F, 0, p1, a, b, k1); mf_leafkids(F, 1, p2, a, b, k2);
  if (!(k1 == 1 && k2 == 1)) return 0;
  int a1 = mf_innerkid(F, 0, p1), a2 = mf_innerkid(F, 1, p2);
  if (a1 == 0 || a2 == 0) return 0;

  int l21, d21, k21, l22, d22, k22;
  mf_leafkids(F, 0, a1, l21, d21, k21); mf_leafkids(F, 1, a2, l22, d22, k22);
  if (!(k21 == 1 && k22 == 1 && l21 == l22)) return 0;
  a1 = mf_innerkid(F, 0, a1); a2 = mf_innerkid(F, 1, a2);
  if (a1 == 0 || a2 == 0) return 0;

  int e31, f31, k31, e32, f32, k32;
  mf_leafkids(F, 0, a1, e31, f31, k31); mf_leafkids(F, 1, a2, e32, f32, k32);
  if (k31 == 0 || k32 == 0) return 0;
  int next1 = 0, next2 = 0;
  if (k31 == 1 && k32 == 1) {
    if (e31 != e32) return 0;
    next1 = mf_innerkid(F, 0, a1); next2 = mf_innerkid(F, 1, a2);
    if (next1 == 0 || next2 == 0) return 0;
  } else if (k31 == 1 && k32 == 2) {
    if (!(e31 == e32 || e31 == f32)) return 0;
    next1 = mf_innerkid(F, 0, a1); next2 = a2; if (next1 == 0) return 0;
  } else if (k31 == 2 && k32 == 1) {
    if (!(e32 == e31 || e32 == f31)) return 0;
    next1 = a1; next2 = mf_innerkid(F, 1, a2); if (next2 == 0) return 0;
  } else return 0;

  int g21, g22, gk; mf_leafkids(F, 1, next2, g21, g22, gk);
  int c1, c2, ck; mf_leafkids(F, 0, next1, c1, c2, ck);
  for (int y : {c1, c2}) {
    if (y == 0) continue;
    if (gk >= 1 && y == g21) return y;
    if (gk == 2 && y == g22) return y;
  }
  return 0;
}

// 3-2 near-cherry around leaf p (t==2). Returns r (taxon deleted as a singleton), or 0.
inline int mf_threetwo(const MForest& F, int p) {
  int sib1 = mf_sibleaf(F, 0, p), sib2 = mf_sibleaf(F, 1, p);
  if (sib1 == 0 || sib2 == 0) return 0;
  if (sib1 == sib2) return 0;
  int ppar2 = F.par[1][p];
  if (ppar2 != 0 && F.par[1][ppar2] != 0 && F.par[1][ppar2] == F.par[1][sib1]) return sib2;
  int ppar1 = F.par[0][p];
  if (ppar1 != 0 && F.par[0][ppar1] != 0 && F.par[0][ppar1] == F.par[0][sib2]) return sib1;
  return 0;
}

// reduce_linear(maf) -> (reduced_maf, leaf_map, actions, k_adjustment).
// O(t*n) worklist of common-cherry / (4,3)-chain / 3-2 rules; actions + k_adjustment plug straight
// into the production SolverState/reconstruct (actions use original labels, which persist here).
struct LinearReduction {
  MAFP reduced_maf;
  std::vector<int> leaf_map;                           // leaf_map[k] = original label of reduced leaf k
  std::vector<std::tuple<bool,int,int>> actions;       // (is_singleton, orig_b, orig_a)
  int k_adjustment = 0;
};

inline LinearReduction reduce_linear(const MAFP& maf) {
  int n = maf.n; int t = (int)maf.T.size();
  MForest F = make_mforest(maf);
  LinearReduction R;
  std::vector<char> inq(n + 1, 1);
  std::vector<int> queue; queue.reserve(n);
  for (int v = 1; v <= n; ++v) queue.push_back(v);

  auto enq = [&](int v) { if (v >= 1 && v <= n && F.alive[v] && !inq[v]) { inq[v] = 1; queue.push_back(v); } };
  // neighbourhood re-enqueue: kept leaf + per-tree sibling and <=5 spine ancestors' leaf children
  auto reenqueue = [&](const std::vector<std::pair<int,int>>& res, int kept) {
    enq(kept);
    for (int i = 0; i < t; ++i) {
      int s = res[i].first, g = res[i].second;
      if (s <= n) enq(s);
      int u = g, steps = 0;
      while (u != 0 && steps < 5) {
        int a = F.c1[i][u], b = F.c2[i][u];
        if (a != 0 && a <= n) enq(a);
        if (b != 0 && b <= n) enq(b);
        u = F.par[i][u]; ++steps;
      }
    }
  };

  std::vector<std::pair<int,int>> res;
  while (!queue.empty()) {
    int v = queue.back(); queue.pop_back(); inq[v] = 0;
    if (!F.alive[v]) continue;

    int rem, keep;
    if (mf_cherry(F, v, rem, keep)) {
      mf_remove(F, rem, res);
      R.actions.emplace_back(false, rem, keep); reenqueue(res, keep); continue;
    }
    if (t == 2) {
      int y = mf_chain43(F, v);
      if (y != 0) { mf_remove(F, y, res); R.actions.emplace_back(false, y, v); reenqueue(res, v); continue; }
      int r = mf_threetwo(F, v);
      if (r != 0) { mf_remove(F, r, res); R.actions.emplace_back(true, r, 0); R.k_adjustment += 1; reenqueue(res, v); continue; }
    }
  }

  std::vector<int> alive_leaves;
  for (int v = 1; v <= n; ++v) if (F.alive[v]) alive_leaves.push_back(v);
  if ((int)alive_leaves.size() <= 1) {
    std::vector<int> keep = alive_leaves.empty() ? std::vector<int>{1} : alive_leaves;
    R.reduced_maf = restrict_to_leaves(maf, keep).first;
    R.leaf_map = keep;
    return R;
  }
  auto [rm, lmap] = restrict_to_leaves(maf, alive_leaves);
  R.reduced_maf = std::move(rm); R.leaf_map = std::move(lmap);
  return R;
}

// solve_with_reductions_linear(maf, solve_func): linear analogue of solve_with_reductions.
inline std::pair<std::vector<BitSet>, double>
solve_with_reductions_linear(const MAFP& maf, const ClusterSolver& solve_func) {
  LinearReduction R = reduce_linear(maf);

  std::vector<BitSet> sol_reduced; double dual = 0.0;
  if (R.reduced_maf.n == 0)      { sol_reduced = {}; dual = 0.0; }
  else if (R.reduced_maf.n == 1) { sol_reduced = {BitSet{1}}; dual = 0.0; }
  else                           { std::tie(sol_reduced, dual) = solve_func(R.reduced_maf); }

  // Build a SolverState the production reconstruct() understands: reduced maf + leaf_map + actions.
  SolverState state(R.reduced_maf);
  state.leaf_map = R.leaf_map;
  state.pending_actions = R.actions;
  state.k_adjustment = R.k_adjustment;

  std::vector<BitSet> sol = reconstruct(state, sol_reduced);
  return {sol, dual + (double)R.k_adjustment};
}

}} // namespace maf::lineartime
