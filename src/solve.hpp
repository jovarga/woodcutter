// The Branch-Cut-and-Price driver.
//   * BCPNode    — one B&B node: its master LP (MPModel), lpsol, clique cuts, and branch rows.
//   * BCPContext — the whole search: a vector of nodes (nullptr == pruned), the shared incumbent,
//                  the accumulated MIP column pool, and the mast cache.
// The production entry is now `solve_bcp = ctx_solve`, which processes the root (column generation +
// clique cuts + MIP heuristic) and then best-bound branch-and-price until every node is pruned.
//
// IMPLEMENTATION NOTES:
//   * MPModel::clone (getModel/passModel) preserves all row/column indices, so cut-row (cons) and
//     branch-row (branch_cons) indices survive a clone with no remapping.
//   * Strong branching (ctx_get_branch_vertex) uses LP OBJECTIVE values (obj0/obj1), which are
//     unique at the optimum even when the degenerate LP vertex is not — so scores match cross-API.
//   * v_saturation sums LP column values in canonical lex column order (deterministic FP summation),
//     matching the compute_cuts canonical-sort convention.
//   * The branch-constraint dual adjustment to the pricing internal weights is
//     iw(i,v) = -beta(i,v) + sum_{branch rows on (i,v)} dual(row)  (both cut and no-cut branches).
#pragma once
#include "core.hpp"
#include "instance.hpp"
#include "mpmodel.hpp"
#include "mast.hpp"
#include "mast/pricing.hpp"           // pricing_mast_solve (cut-pricing branch)
#include "columngeneration.hpp"
#include "cut_generation.hpp"         // compute_cuts (clique-cut separation)
#include "preprocessing.hpp"          // solve_clusters
#include "reduction_rule_stack.hpp"   // SolverState, reduce_exhaustive, reconstruct
#include <vector>
#include <utility>
#include <memory>
#include <limits>
#include <cmath>
#include <functional>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>

namespace maf {

// --- branching trace (env MAF_BTRACE=<file>): a compact per-decision validation log. Zero cost
// when the env var is unset. Lines:
//   R <primal> <dual> <ncuts>            after processing the root
//   B <node1based> <tree> <vertex>       each branch decision (node + strong-branch vertex)
//   N <node1based> <dual> <opt01>        each child node's dual after its column generation
//   P <node1based> <lpobj>               each node pruned
//   E <ncomp> <dual>                     final result
inline std::FILE*& btrace_file() { static std::FILE* f = nullptr; return f; }
inline void btrace_open() {
  const char* p = std::getenv("MAF_BTRACE");
  if (p && *p) btrace_file() = std::fopen(p, "w");
}
inline void btrace_close() { if (btrace_file()) { std::fclose(btrace_file()); btrace_file() = nullptr; } }
template <class... A> inline void BT(const char* fmt, A... a) {
  if (btrace_file()) { std::fprintf(btrace_file(), fmt, a...); std::fflush(btrace_file()); }
}

// A clique cut: an ordered list of (column, coefficient), kept insertion-ordered (initial clique
// columns ascending-by-index, then pricing-extension columns appended) so the member order fed to the
// pricing subproblem is deterministic. All coeffs are 1.0.
using Cut = std::vector<std::pair<BitSet, double>>;

// ---- BCPNode ---------------------------------------------------------------------------------
// One B&B node.  beta (inner-vertex prices) is not stored separately: node.beta[(i,v)] == mp.beta(i,v)
// after each LP solve, and every read of beta is preceded by an LP solve of this node's mp.
struct BCPNode {
  MPModel mp;
  std::vector<std::pair<BitSet, double>> lpsol;   // node.lpsol: column -> LP value (creation order)
  bool lpsol_isopt = false;
  bool lpsol_isvalid = true;

  std::vector<Cut> cuts;                            // node.cuts
  std::vector<double> rhss;                         // node.rhss
  std::vector<int> cons;                            // node.cons: cut index r (mp.cut_dual(r))

  // Branch rows.  branch_cons stores ABSOLUTE HiGHS row indices (stable across clone()).
  std::vector<std::pair<int,int>> branch_nodes;    // (tree_idx 1-based, inner vertex v)
  std::vector<double> branch_rhss;                 // {0,1} rhs placeholder
  std::vector<int> branch_cons;                    // absolute HiGHS row index of each branch row

  // copy_node(node): clone the master (fresh, unsolved) and copy all bookkeeping.
  BCPNode copy() const {
    BCPNode nn;
    nn.mp = mp.clone();
    nn.lpsol = lpsol;
    nn.lpsol_isopt = lpsol_isopt;
    nn.lpsol_isvalid = lpsol_isvalid;
    nn.cuts = cuts;
    nn.rhss = rhss;
    nn.cons = cons;
    nn.branch_nodes = branch_nodes;
    nn.branch_rhss = branch_rhss;
    nn.branch_cons = branch_cons;
    return nn;
  }
};

// ---- BCPContext -------------------------------------------------------------------------------
struct BCPContext {
  const MAFP* maf;
  std::vector<std::unique_ptr<BCPNode>> nodes;   // nullptr == pruned

  std::vector<BitSet> incumbent;                 // ctx.incumbent (a partition of the leaves)
  std::vector<BitSet> mip_vars;                  // ctx.mip_vars, kept insertion-ordered & unique
  std::unordered_set<BitSet> mip_vars_seen;

  explicit BCPContext(const MAFP& m) : maf(&m) {
    auto node = std::make_unique<BCPNode>();
    node->mp = MPModel(m, {}, /*binary=*/false);          // MPModel(maf, initvars=[])
    nodes.push_back(std::move(node));
    for (int x = 1; x <= m.n; ++x) incumbent.push_back(BitSet{x});   // [{x} for x in 1:n]
  }

  void mip_vars_add(const BitSet& s) {
    if (mip_vars_seen.insert(s).second) mip_vars.push_back(s);
  }
};

// Current incumbent size.
inline int ctx_primal(const BCPContext& ctx) { return (int)ctx.incumbent.size(); }

// Node LP objective. i is 0-based here.
inline double ctx_lpobj(const BCPContext& ctx, int i) {
  double s = 0.0;
  for (const auto& [var, val] : ctx.nodes[i]->lpsol) s += (double)(var.size() - 1) * val;
  return (double)ctx.maf->n - s;
}

// 0 unless every live node is LP-optimal; else min lpobj.
inline double ctx_dual(const BCPContext& ctx) {
  for (const auto& np : ctx.nodes) if (np && !np->lpsol_isopt) return 0.0;
  double m = std::numeric_limits<double>::infinity();
  for (int i = 0; i < (int)ctx.nodes.size(); ++i)
    if (ctx.nodes[i]) m = std::min(m, ctx_lpobj(ctx, i));
  return m;
}

// Solve the node master LP and refresh solution caches.
inline void ctx_rmp_lpsolve(BCPContext& ctx, int i) {
  BCPNode& node = *ctx.nodes[i];
  node.mp.solve();                              // mp_solve! — also refreshes mp.beta(i,v)

  node.lpsol.clear();
  node.lpsol.reserve(node.mp.col_sets.size());
  for (int r = 0; r < (int)node.mp.col_sets.size(); ++r)
    node.lpsol.emplace_back(node.mp.col_sets[r], node.mp.sol_value(r));

  // Add positive-valued columns to the shared MIP pool.
  for (int r = 0; r < (int)node.mp.col_sets.size(); ++r)
    if (node.mp.sol_value(r) >= 1e-8) ctx.mip_vars_add(node.mp.col_sets[r]);

  node.lpsol_isvalid = true;
}

// Whether a variable receives coefficient 1 in a clique cut.
inline bool covers_cut(const MAFP& maf, const Cut& cut, const BitSet& var) {
  std::vector<BitSet> var_inner(maf.T.size());
  for (size_t i = 0; i < maf.T.size(); ++i) {
    std::vector<int> iv = inner_vs(maf.T[i], var);
    var_inner[i] = BitSet(iv.begin(), iv.end());
  }
  for (const auto& [v, _w] : cut) {
    bool any_tree = false;
    for (size_t i = 0; i < maf.T.size(); ++i) {
      std::vector<int> iv = inner_vs(maf.T[i], v);
      BitSet vi(iv.begin(), iv.end());
      if (!isdisjoint(var_inner[i], vi)) { any_tree = true; break; }
    }
    if (!any_tree) return false;
  }
  return true;
}

// Build the pricing PhyloTrees for a node.  internal_weight(i,v) = -beta(i,v) + sum of branch-row
// duals on (i,v).  Both the cut branch (pricing internal_weight = beta = -node.beta + dual) and the
// no-cut branch (mast internal_weight = -beta_mast = -(node.beta - dual)) reduce to this formula.
inline std::vector<mast::PhyloTree> build_pricing_trees(const BCPContext& ctx, const BCPNode& node) {
  const int n = ctx.maf->n, T = (int)ctx.maf->T.size();
  // Per-tree branch-dual adjustment: adj[tree0][v] += dual(branch row) for each branch on (tree0+1,v).
  std::vector<std::unordered_map<int,double>> adj(T);
  for (size_t b = 0; b < node.branch_nodes.size(); ++b) {
    int tree0 = node.branch_nodes[b].first - 1;
    int v = node.branch_nodes[b].second;
    adj[tree0][v] += node.mp.row_dual_abs(node.branch_cons[b]);
  }
  std::vector<mast::PhyloTree> trees;
  trees.reserve(T);
  for (int i = 0; i < T; ++i) {
    std::vector<W> iw;
    iw.reserve(n - 1);
    for (int v = n + 1; v <= ctx.maf->T[i].nv(); ++v) {
      W w = -node.mp.beta(i, v);
      if (!adj[i].empty()) { auto it = adj[i].find(v); if (it != adj[i].end()) w += it->second; }
      iw.push_back(w);
    }
    trees.push_back(mast::build_phylo_tree(ctx.maf->T[i], iw, n));
  }
  return trees;
}

// Price new columns for node i.
inline void ctx_price(BCPContext& ctx, int i, bool use_clique_prices = true) {
  BCPNode& node = *ctx.nodes[i];
  const int n = ctx.maf->n;
  std::vector<int> added;                        // real indices of newly-added columns

  if (!node.cuts.empty() && use_clique_prices) {
    // Cut-pricing branch: price clique-cut extensions via pricing_mast.
    //
    // CORRECTNESS FIX (2026-07-04): pass the FULL clique member list to the pricer, NOT the
    // ee1233f "reduced_cliques" (members restricted to their BASIC columns).  The pricer activates a
    // cut's dual penalty when a candidate column conflicts with *all* passed members
    // (pricing.hpp:358, isdisjoint over required_variables).  The master, however, gives a new column
    // coefficient 1 in the cut iff it conflicts with *all* FULL members (covers_cut, below).  With
    // basic-only members (a subset), the penalty fires MORE eagerly than the true coefficient, so the
    // pricer UNDER-estimates reduced costs, wrongly certifies "no improving column", and column
    // generation stops at an INFLATED (invalid) LP bound.  On batch04_rgmul_0922 (t=4,n=16) that made
    // the exact solver prove 14 optimal when k*=13 (LP bound 13.2 > 13; missing columns {3,4},{3,5},
    // {4,5}).  Using the full member list makes the pricing "covers" test match covers_cut exactly, so
    // CG converges to the true LP relaxation. Cost: negligible (exact117 4.2s->5.1s; 0 k* changes
    // across the 150 exact2 instances).
    std::vector<std::vector<std::vector<int>>> cut_sets;   // active cuts -> member leaf-lists
    std::vector<W> cut_weights;
    for (size_t j = 0; j < node.cuts.size(); ++j) {
      double pen = node.mp.cut_dual(node.cons[j]);
      if (pen <= -1e-8) {                        // active_cuts
        std::vector<std::vector<int>> members;
        for (const auto& [col, _c] : node.cuts[j])
          members.push_back(col.to_vec());       // ALL members, ascending leaves
        cut_sets.push_back(std::move(members));
        cut_weights.push_back(pen);
      }
    }

    std::vector<mast::PhyloTree> trees = build_pricing_trees(ctx, node);
    auto result = mast::pricing_mast_solve(std::move(trees), cut_sets, cut_weights);
    std::vector<BitSet> new_vars_sets;
    for (const auto& s : result)
      if (s.objective >= 1 + 1e-8) new_vars_sets.emplace_back(s.taxa.begin(), s.taxa.end());

    added = node.mp.add_columns(new_vars_sets);

    // Extend each covered clique cut with the newly-priced columns. Add the member only if not
    // already present so a column re-priced across CG iterations can't append a duplicate member.
    // Inert for k*: a duplicate member has an identical conflict set.
    for (const BitSet& var : new_vars_sets)
      for (size_t j = 0; j < node.cuts.size(); ++j)
        if (covers_cut(*ctx.maf, node.cuts[j], var)) {
          node.mp.set_cut_coefficient(node.cons[j], var, 1.0);
          bool present = false;
          for (const auto& [col, _c] : node.cuts[j]) if (col == var) { present = true; break; }
          if (!present) node.cuts[j].emplace_back(var, 1.0);
        }

    // `isopt` requires both pricing and column insertion to produce no new variables.
    if (new_vars_sets.empty()) node.lpsol_isopt = true;
  } else {
    // No-cut branch: blacklist = all columns iff cuts nonempty (else empty).
    std::vector<std::vector<int>> blacklist;
    if (!node.cuts.empty())
      for (const BitSet& s : node.mp.col_sets) blacklist.push_back(s.to_vec());

    std::vector<mast::PhyloTree> trees = build_pricing_trees(ctx, node);
    std::vector<std::pair<W, BitSet>> priced = mast_solve(trees, n, blacklist);
    std::vector<W> best_ws;
    std::vector<BitSet> best_sets;
    best_ws.reserve(priced.size());
    best_sets.reserve(priced.size());
    for (auto& p : priced) { best_ws.push_back(p.first); best_sets.push_back(p.second); }

    std::vector<BitSet> new_vars = select_variables(best_ws, best_sets, node.mp.col_idx);
    added = node.mp.add_columns(new_vars);
    if (added.empty()) node.lpsol_isopt = true;
  }

  // Set the branch-row coefficient of each newly-priced column that covers a branched vertex.
  for (int r : added)
    for (size_t b = 0; b < node.branch_nodes.size(); ++b) {
      int tree0 = node.branch_nodes[b].first - 1;
      int v = node.branch_nodes[b].second;
      const std::vector<int>& iv = node.mp.col_inner[r][tree0];
      if (std::binary_search(iv.begin(), iv.end(), v))
        node.mp.set_row_coefficient_abs(node.branch_cons[b], r, 1.0);
    }
}

// Solve the RMP as a MIP, warm-started from the current incumbent.
inline void ctx_rmp_mipsolve(BCPContext& ctx) {
  auto [obj, sol] = solve_rmp(*ctx.maf, ctx.mip_vars, &ctx.incumbent);
  (void)obj;
  ctx.incumbent = sol;
}

// Column generation loop for node i.
inline void ctx_columngeneration(BCPContext& ctx, int i, bool use_clique_prices = true) {
  BCPNode& node = *ctx.nodes[i];
  while (!node.lpsol_isopt) {
    check_budget();
    ctx_rmp_lpsolve(ctx, i);
    ctx_price(ctx, i, use_clique_prices);
  }
}

// Separate clique cuts for node i.
inline void ctx_clique_cuts(BCPContext& ctx, int i, int max_cuts) {
  BCPNode& node = *ctx.nodes[i];
  auto [newcuts, newrhss] = compute_cuts(*ctx.maf, node.lpsol, max_cuts);
  for (size_t c = 0; c < newcuts.size(); ++c) {
    int cut_r = node.mp.add_cut_row(newcuts[c], newrhss[c]);
    Cut cut;
    for (const BitSet& col : newcuts[c]) cut.emplace_back(col, 1.0);   // ascending-index member order
    node.cuts.push_back(std::move(cut));
    node.rhss.push_back(newrhss[c]);
    node.cons.push_back(cut_r);
  }
  node.lpsol_isvalid = false;
  node.lpsol_isopt = false;
}

// LP-solve invalid nodes, return argmin lpobj (0-based).
inline int ctx_get_branch_node(BCPContext& ctx) {
  for (int i = 0; i < (int)ctx.nodes.size(); ++i)
    if (ctx.nodes[i] && !ctx.nodes[i]->lpsol_isvalid) {
      try { ctx_rmp_lpsolve(ctx, i); }
      catch (const NodeInfeasible&) { ctx.nodes[i].reset(); }   // infeasible branch -> prune
    }
  int best = -1;                                                // -1 = all nodes pruned
  double bestd = std::numeric_limits<double>::infinity();
  for (int i = 0; i < (int)ctx.nodes.size(); ++i) {
    if (!ctx.nodes[i]) continue;
    double d = ctx_lpobj(ctx, i);
    if (d < bestd) { bestd = d; best = i; }     // argmin: first minimal index
  }
  return best;
}

// cov[tree0][v] = sum of LP values of columns covering inner vertex v in tree tree0. Columns are
// summed in canonical lex order (deterministic FP).
inline std::vector<std::vector<double>> v_saturation(const BCPContext& ctx, const BCPNode& node) {
  const int T = (int)ctx.maf->T.size();
  const int nv = ctx.maf->T[0].nv();
  std::vector<std::vector<double>> cov(T, std::vector<double>(nv + 1, 0.0));
  std::vector<int> ordr;
  for (int r = 0; r < (int)node.lpsol.size(); ++r)
    if (node.lpsol[r].second > 1e-8) ordr.push_back(r);
  std::sort(ordr.begin(), ordr.end(),
            [&](int a, int b){ return node.lpsol[a].first < node.lpsol[b].first; });
  for (int r : ordr) {
    double val = node.lpsol[r].second;
    for (int i = 0; i < T; ++i)
      for (int v : node.mp.col_inner[r][i]) cov[i][v] += val;
  }
  return cov;
}

// Basis-aware branching score. For each (tree, inner vertex v):
//   cov           = v_saturation (LP mass of columns covering v)
//   row_basics    = # of BASIC columns covering v ;  nonrow_basics = (#basic cols) - row_basics
//   score         = 4*cov*(1-cov) * sqrt((1+row_basics)*(1+nonrow_basics))
// 4*cov*(1-cov) peaks at cov=0.5 (most fractional); the sqrt term favours vertices whose basic columns
// split more evenly across the two branch children (a stronger disjunction).  Leaves score 0 (cov==0).
// The basis (col_is_basic) is the one captured by node.mp's last LP solve, i.e. the CG optimum.
inline std::vector<std::vector<double>> row_branching_scores(const BCPContext& ctx, const BCPNode& node) {
  const int T = (int)ctx.maf->T.size();
  const int nv = ctx.maf->T[0].nv();
  std::vector<std::vector<double>> cov = v_saturation(ctx, node);   // [T][nv+1]
  int n_basic = 0;
  std::vector<std::vector<int>> row_basics(T, std::vector<int>(nv + 1, 0));
  for (int r = 0; r < (int)node.mp.col_sets.size(); ++r) {
    if (!node.mp.col_is_basic(r)) continue;
    ++n_basic;
    for (int i = 0; i < T; ++i)
      for (int v : node.mp.col_inner[r][i]) row_basics[i][v] += 1;
  }
  std::vector<std::vector<double>> score(T, std::vector<double>(nv + 1, 0.0));
  for (int i = 0; i < T; ++i)
    for (int v = 1; v <= nv; ++v) {
      double c = cov[i][v];
      double rb = (double)row_basics[i][v];
      double nrb = (double)n_basic - rb;
      score[i][v] = 4.0 * c * (1.0 - c) * std::sqrt((1.0 + rb) * (1.0 + nrb));
    }
  return score;
}

// Strong branching over the top `trials` vertices by row_branching_scores; returns the (tree_idx, v)
// maximizing (obj0-obj)*(obj1-obj).
inline std::pair<int,int> ctx_get_branch_vertex(BCPContext& ctx, int i, int trials = 15) {
  std::vector<std::vector<double>> score = row_branching_scores(ctx, *ctx.nodes[i]);
  const int T = (int)ctx.maf->T.size();
  const int nv = ctx.maf->T[0].nv();

  // Column-major order: tree index varies fastest.
  std::vector<std::pair<int,int>> ks;
  std::vector<double> key;
  ks.reserve((size_t)T * nv);
  key.reserve((size_t)T * nv);
  for (int v = 1; v <= nv; ++v)
    for (int t0 = 0; t0 < T; ++t0) {
      ks.emplace_back(t0 + 1, v);
      key.push_back(score[t0][v]);
    }
  // Stable descending order; ties keep ascending index order.
  std::vector<int> order(ks.size());
  for (int j = 0; j < (int)order.size(); ++j) order[j] = j;
  std::stable_sort(order.begin(), order.end(), [&](int a, int b){ return key[a] > key[b]; });

  ctx_rmp_lpsolve(ctx, i);   // refresh lpsol/basis
  BCPNode nn = ctx.nodes[i]->copy();
  MPModel& mp = nn.mp;
  // Strong-branching probes may hit an infeasible direction; treat its objective as +inf (that
  // direction prunes a whole subtree, so it is a strong branch).  Never fires when no probe is
  // infeasible (e.g. all exact2 instances), so their strong-branching choices are unchanged.
  auto probe = [&]() -> double { try { return mp.solve(); } catch (const NodeInfeasible&) { return std::numeric_limits<double>::infinity(); } };
  double obj = probe();

  int best_idx = order.empty() ? 0 : order[0];
  double best_score = -std::numeric_limits<double>::infinity();
  int lim = std::min(trials, (int)order.size());
  for (int oi = 0; oi < lim; ++oi) {
    int idx = order[oi];
    int tree_idx = ks[idx].first, v = ks[idx].second;
    std::vector<int> cover = mp.cols_covering(tree_idx - 1, v);
    // node value 0: node_expr <= 0
    int c0 = mp.add_branch_row(cover, /*ge1=*/false, /*track=*/false);
    double obj0 = probe();
    mp.set_row_bounds_abs(c0, -kHighsInf, 1.0);          // set_normalized_rhs(c0, 1)
    // node value 1: node_expr >= 1
    int c1 = mp.add_branch_row(cover, /*ge1=*/true, /*track=*/false);
    double obj1 = probe();
    mp.set_row_bounds_abs(c1, 0.0, kHighsInf);           // set_normalized_rhs(c1, 0)
    double score = (obj0 - obj) * (obj1 - obj);
    if (score > best_score) { best_score = score; best_idx = idx; }   // strict >: keep earlier
  }
  return ks[best_idx];
}

// Split node i on branch vertex k=(tree_idx,v):
//   node i        gets `sum_{cols covering (tree_idx,v)} a <= 0`   (0-branch, in place)
//   a fresh clone gets `... >= 1`                                  (1-branch, appended)
// Returns the pair (i, new_index) (0-based).
inline std::pair<int,int> ctx_branch(BCPContext& ctx, int i, std::pair<int,int> k) {
  BT("B %d %d %d\n", i + 1, k.first, k.second);
  BCPNode nn = ctx.nodes[i]->copy();
  int tree_idx = k.first, v = k.second;
  std::vector<int> cover_i  = ctx.nodes[i]->mp.cols_covering(tree_idx - 1, v);
  std::vector<int> cover_nn = nn.mp.cols_covering(tree_idx - 1, v);   // identical set (clone)

  int c0 = ctx.nodes[i]->mp.add_branch_row(cover_i, /*ge1=*/false, /*track=*/true);
  ctx.nodes[i]->branch_nodes.push_back(k);
  ctx.nodes[i]->branch_rhss.push_back(0.0);
  ctx.nodes[i]->branch_cons.push_back(c0);
  ctx.nodes[i]->lpsol_isvalid = false;
  ctx.nodes[i]->lpsol_isopt = false;

  int c1 = nn.mp.add_branch_row(cover_nn, /*ge1=*/true, /*track=*/true);
  nn.branch_nodes.push_back(k);
  nn.branch_rhss.push_back(0.0);
  nn.branch_cons.push_back(c1);
  nn.lpsol_isvalid = false;
  nn.lpsol_isopt = false;
  ctx.nodes.push_back(std::make_unique<BCPNode>(std::move(nn)));

  return {i, (int)ctx.nodes.size() - 1};
}

// Pick node + vertex, then branch.
inline std::pair<int,int> ctx_branch(BCPContext& ctx, int trials = 5) {
  int i = ctx_get_branch_node(ctx);
  if (i < 0) return {-1, -1};                  // all live nodes were pruned as infeasible
  std::pair<int,int> k = ctx_get_branch_vertex(ctx, i, trials);
  return ctx_branch(ctx, i, k);
}

// Root column generation + clique-cut separation loop.
inline void ctx_process_root(BCPContext& ctx, double progress_delta = 0.03) {
  ctx_columngeneration(ctx, 0);
  ctx_rmp_mipsolve(ctx);
  if (ctx_primal(ctx) - ctx_lpobj(ctx, 0) <= 1 - 1e-8) return;
  if ((int)ctx.maf->T.size() == 2) return;   // 2-tree: skip cuts

  const bool ucp = (int)ctx.maf->T.size() > 2;
  double olddual = ctx_lpobj(ctx, 0);          // dual BEFORE the initial cut round
  for (int r = 0; r < 1; ++r) {                // ONE initial cut round (was 3)
    ctx_clique_cuts(ctx, 0, /*max_cuts=*/3);
    ctx_columngeneration(ctx, 0, ucp);
  }
  ctx_rmp_mipsolve(ctx);
  while (ctx_primal(ctx) - ctx_lpobj(ctx, 0) > 1 - 1e-8 &&
         ctx_lpobj(ctx, 0) - olddual > progress_delta) {
    check_budget();
    olddual = ctx_lpobj(ctx, 0);
    int n_cuts = -1, it = 0;
    while (ctx_primal(ctx) - ctx_lpobj(ctx, 0) > 0.98 &&
           (int)ctx.nodes[0]->cuts.size() > n_cuts && it < 10) {   // inner cap 20 -> 10
      n_cuts = (int)ctx.nodes[0]->cuts.size();
      it += 1;
      ctx_clique_cuts(ctx, 0, /*max_cuts=*/4);
      ctx_rmp_lpsolve(ctx, 0);
    }
    ctx_columngeneration(ctx, 0, ucp);
    // A better MIP solution only helps if the dual crossed an integer AND the incumbent isn't optimal.
    if (std::ceil(ctx_lpobj(ctx, 0) + 1e-8) > std::ceil(olddual + 1e-8) &&
        ctx_primal(ctx) - ctx_lpobj(ctx, 0) > 1 - 1e-8) {
      ctx_rmp_mipsolve(ctx);
    }
  }
}

// Optional incumbent hook: if set, ctx_solve calls it each B&B iteration with the current incumbent
// forest and lower bound. Returning true asks the search to STOP (throws Interrupted). Default EMPTY
// => INERT: exact/heuristic never set it, so ctx_solve is byte-for-byte unchanged (only cost: one
// null-check per iteration). The lower-bound track's large-block solver arms it to stop once the
// incumbent is within the approximation factor.
inline std::function<bool(const std::vector<BitSet>&, double)> g_incumbent_hook{};

// Process root, then best-bound branch-and-price. The B&P loop includes three throughput
// optimizations:
//   * adaptive strong-branching effort: `trials = strong_branch_trials[min(it, end)]` — more probes
//     early (20), fewer later (down to 5);
//   * a PRUNE sweep BEFORE the primal (MIP) solve, so dominated nodes are dropped without a MIP;
//   * a CONDITIONAL MIP solve — skip it unless the dual just crossed an integer, or on a few early
//     iterations, or every mip_delay-th iteration (the primal rarely improves otherwise).
// `warm_cols` (default empty): seed the root master with a column pool so the BCP warm-starts from it
// rather than re-deriving every column from scratch.
inline std::pair<std::vector<BitSet>, double>
ctx_solve(const MAFP& maf,
          const std::vector<BitSet>& warm_cols = {},
          int mip_delay = 7,
          std::vector<int> strong_branch_trials = {20, 20, 15, 10, 10, 5}) {
  BCPContext ctx(maf);
  const bool ucp = (int)ctx.maf->T.size() > 2;   // use_clique_prices only for t > 2
  if (!warm_cols.empty()) {                        // warm-start the root master with the given pool
    std::vector<BitSet> add;
    for (const BitSet& c : warm_cols) if ((int)c.size() >= 2) add.push_back(c);
    if (!add.empty()) ctx.nodes[0]->mp.add_columns(add);
  }

  ctx_process_root(ctx);
  const double root_lb = ctx_lpobj(ctx, 0);      // valid running lower bound (used only by the hook)
  BT("R %d %.10g %d\n", ctx_primal(ctx), ctx_lpobj(ctx, 0), (int)ctx.nodes[0]->cuts.size());
  if (ctx_primal(ctx) - ctx_lpobj(ctx, 0) <= 1 - 1e-8) {
    BT("E %zu %.10g\n", ctx.incumbent.size(), ctx_lpobj(ctx, 0));
    return {ctx.incumbent, ctx_lpobj(ctx, 0)};
  }

  double dual_bound = std::numeric_limits<double>::infinity();
  auto all_nothing = [&]() { for (const auto& np : ctx.nodes) if (np) return false; return true; };
  // Column-generate a branch node; if its master is INFEASIBLE, that branch admits no agreement
  // forest, so PRUNE it (reset the node) rather than aborting the whole solve.  (Standard branch-and-
  // price.  Never triggers when the solve stays optimal — e.g. all exact2 instances — so those runs
  // are byte-for-byte unchanged.)
  auto cg_or_prune = [&](int idx) {
    try { ctx_columngeneration(ctx, idx, ucp); }
    catch (const NodeInfeasible&) { if (ctx.nodes[idx]) ctx.nodes[idx].reset(); }
  };
  // Prune every LP-optimal node whose bound can no longer beat the incumbent; record dual_bound.
  auto prune_dominated = [&]() {
    for (int idx = 0; idx < (int)ctx.nodes.size(); ++idx) {
      if (ctx.nodes[idx] && ctx.nodes[idx]->lpsol_isopt &&
          ctx_lpobj(ctx, idx) + 1 - 1e-8 > (double)ctx_primal(ctx)) {
        dual_bound = std::min(dual_bound, ctx_lpobj(ctx, idx));
        BT("P %d %.10g\n", idx + 1, ctx_lpobj(ctx, idx));
        ctx.nodes[idx].reset();
      }
    }
  };

  int it = 1;
  while (!all_nothing()) {
    check_budget();
    if (g_incumbent_hook &&
        g_incumbent_hook(ctx.incumbent, std::isfinite(dual_bound) ? std::max(dual_bound, root_lb) : root_lb))
      throw Interrupted{};
    double olddual = ctx_dual(ctx);
    int trials = strong_branch_trials[std::min(it, (int)strong_branch_trials.size()) - 1];
    auto [bi, bj] = ctx_branch(ctx, trials);
    if (bi < 0) continue;                       // all nodes pruned during selection -> loop exits
    cg_or_prune(bi);
    cg_or_prune(bj);
    // N lines log each child's bound after its CG.  Guard against a child pruned as infeasible
    // (null node) so the trace — and production — never dereferences a dead node.
    { double di = ctx.nodes[bi] ? ctx_lpobj(ctx, bi) : std::numeric_limits<double>::infinity();
      double dj = ctx.nodes[bj] ? ctx_lpobj(ctx, bj) : std::numeric_limits<double>::infinity();
      BT("N %d %.10g %d\n", bi + 1, di, ctx.nodes[bi] ? (int)ctx.nodes[bi]->lpsol_isopt : -1);
      BT("N %d %.10g %d\n", bj + 1, dj, ctx.nodes[bj] ? (int)ctx.nodes[bj]->lpsol_isopt : -1); }

    // Prune BEFORE the (possibly skipped) primal solve (3cb299e).
    prune_dominated();
    if (all_nothing()) break;

    // Conditional primal (MIP) solve: only when it can plausibly improve the incumbent.
    if (std::ceil(ctx_dual(ctx) + 1e-8) > std::ceil(olddual + 1e-8) ||
        it == 1 || it == 2 || it == 4 || it == 7 || it == 10 || (it % mip_delay == 0)) {
      ctx_rmp_mipsolve(ctx);
    }

    // Prune AFTER the primal solve.
    prune_dominated();
    it += 1;
  }
  BT("E %zu %.10g\n", ctx.incumbent.size(), dual_bound);
  return {ctx.incumbent, dual_bound};
}

// solve_bcp(maf) — the production sub-instance solver (driver entry).
// Returns (incumbent partition, dual bound).
inline std::pair<std::vector<BitSet>, double> solve_bcp(const MAFP& maf) {
  return ctx_solve(maf);
}

// Kernelize, solve the reduced instance via solve_func, lift back to original labels. The dual is
// adjusted by k_adjustment (one per ThreeTwo singleton).
inline std::pair<std::vector<BitSet>, double>
solve_with_reductions(const MAFP& maf, const ClusterSolver& solve_func,
    const std::vector<Rule>& sequence = {Rule::PendantSubtree, Rule::Chain43, Rule::ThreeTwo, Rule::PendantSubtree}) {
  SolverState state(maf);
  reduce_exhaustive(state, sequence);

  std::vector<BitSet> sol_reduced; double dual_reduced;
  if (state.maf.n == 0)      { sol_reduced = {};            dual_reduced = 0.0; }
  else if (state.maf.n == 1) { sol_reduced = {BitSet{1}};   dual_reduced = 0.0; }
  else                       { std::tie(sol_reduced, dual_reduced) = solve_func(state.maf); }

  std::vector<BitSet> sol = reconstruct(state, sol_reduced);
  return {sol, dual_reduced + (double)state.k_adjustment};
}

// solve_maf(maf): the production exact-track entry —
//   solve_with_reductions(maf; solve_func = m -> solve_clusters(m, solve_bcp)).
inline std::pair<std::vector<BitSet>, double> solve_maf(const MAFP& maf) {
  return solve_with_reductions(maf, [](const MAFP& m) { return solve_clusters(m, solve_bcp); });
}

} // namespace maf
