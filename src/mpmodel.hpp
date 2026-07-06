// The set-packing master problem over a persistent HiGHS model.
//
// CRITICAL: this must be an INCREMENTAL, WARM-STARTED wrapper: one persistent Highs object, columns
// and cut-rows added in deterministic order, the basis retained between solves (HiGHS warm-starts
// automatically on incremental edits), threads=1.
// Do NOT rebuild/reload the model per iteration.
//
// Master LP:
//   variables   a_s >= 0           one per column s (a taxa BitSet)
//   constraints disjoint[i,v]:  sum_{s : v in inner_vs(T[i],s)} a_s <= 1   (Steiner-vertex disjoint)
//   objective   Min  n - sum_s (|s|-1) a_s
// Column costs are -(|s|-1); the constant n is the objective offset.  All data is integer.
//
// BRANCHING: the B&B copies a node's master and appends branch rows
// (sum_{s covering (i,v)} a_s <= 0  or  >= 1). Cloning via getModel/passModel preserves row+column
// indices, so cut-row and branch-row indices stay valid across a clone. MPModel holds its Highs via
// unique_ptr so it is movable/clonable.
#pragma once
#include "core.hpp"
#include "instance.hpp"
#include "maf_util.hpp"
#include "budget.hpp"
#include "Highs.h"
#include <unordered_map>
#include <vector>
#include <utility>
#include <memory>
#include <algorithm>

namespace maf {

// Validation instrumentation: counts every LP + MIP solve. Reset to 0 by a driver before a solve to
// compare iteration counts. Has no effect on the solver itself.
inline long mp_solve_calls = 0;

struct MPModel {
  const MAFP* maf = nullptr;
  std::unique_ptr<Highs> hs;                 // persistent, warm-started (unique_ptr => MPModel movable)
  bool isbinary = false;

  // Columns. HiGHS column index of real column r is `real_base + r`; `real_base` is 1 when a dummy
  // column (fixed to 0) occupies index 0.
  int real_base = 0;
  std::vector<BitSet> col_sets;             // r -> taxa set, in creation order
  std::unordered_map<BitSet, int> col_idx;  // taxa set -> r

  // Basis bookkeeping, updated every solve.
  std::vector<int> nonbasic_counter;        // per real column r
  std::vector<char> wasbasic;               // per real column r

  // disjoint rows occupy HiGHS row indices 0 .. n_disjoint-1, in (i,v) order:
  // row_of(i,v) = i*n_inner + (v-(n+1)),  i in 0..t-1,  v in n+1..2n-1.  Cut rows are appended
  // (indices n_disjoint .. n_disjoint+n_cut_rows-1), then branch rows after those.
  int n_inner = 0;   // inner vertices per tree = n-1
  int n_disjoint = 0;// t*(n-1)
  int n_cut_rows = 0;
  int n_branch_rows = 0;

  // Last solve's raw HiGHS solution (for accessors below).
  double last_obj = 0.0;
  std::vector<double> last_col_value;
  std::vector<double> last_row_dual;
  std::vector<char> last_col_basic;         // per real column r: was BASIC in last solve's basis

  // inner_vs(T[i], s) for each column — cached so add/build don't recompute.
  // col_inner[r] is a vector over trees, each a list of inner vertices.
  std::vector<std::vector<std::vector<int>>> col_inner;

  // ---- construction -------------------------------------------------------------------------
  MPModel() = default;                       // empty shell used by clone()

  MPModel(const MAFP& maf_, const std::vector<BitSet>& vars, bool binary = false)
      : maf(&maf_), hs(std::make_unique<Highs>()), isbinary(binary) {
    hs->setOptionValue("output_flag", false);
    hs->setOptionValue("threads", 1);
    n_inner = maf->n - 1;
    n_disjoint = (int)maf->T.size() * n_inner;

    // Variables first (dummy iff empty, else one a_s per var), then constraints.
    if (vars.empty()) {
      // dummy == 0  (lb = ub = 0, cost 0) at column 0.
      HighsInt idx0 = 0; double val0 = 0.0;  // no nonzeros yet
      hs->addCol(0.0, 0.0, 0.0, 0, &idx0, &val0);
      real_base = 1;
    } else {
      real_base = 0;
    }
    hs->changeObjectiveOffset((double)maf->n);

    // Add the real columns (cost only; matrix coefficients are filled when rows are added below).
    for (const BitSet& s : vars) register_column(s, /*with_rows=*/false);

    // Add the disjoint rows in (i,v) order, with coefficient 1 for every column covering (i,v).
    // Build, per tree, the column lists touching each inner vertex.
    std::vector<std::vector<std::vector<int>>> rows_cols(maf->T.size(),
        std::vector<std::vector<int>>(n_inner));  // [i][v-(n+1)] -> HiGHS col indices
    for (int r = 0; r < (int)col_sets.size(); ++r)
      for (size_t i = 0; i < maf->T.size(); ++i)
        for (int v : col_inner[r][i]) rows_cols[i][v - (maf->n + 1)].push_back(real_base + r);

    for (size_t i = 0; i < maf->T.size(); ++i)
      for (int v = maf->n + 1; v <= 2 * maf->n - 1; ++v) {
        const std::vector<int>& cols = rows_cols[i][v - (maf->n + 1)];
        std::vector<double> ones(cols.size(), 1.0);
        std::vector<HighsInt> idx(cols.begin(), cols.end());
        hs->addRow(-kHighsInf, 1.0, (HighsInt)idx.size(), idx.data(), ones.data());
      }

    if (binary) set_integral_all();
  }

  // ---- clone: duplicate the whole master via getModel/passModel -----------------------------
  // The cloned model preserves EVERY column/row index (so cut-row and branch-row indices stay
  // valid), and starts WITHOUT a warm-start basis, so its first solve is cold. The last-solve caches
  // are intentionally left empty (the clone is unsolved until re-solved).
  MPModel clone() const {
    MPModel m;
    m.maf = maf;
    m.isbinary = isbinary;
    m.real_base = real_base;
    m.col_sets = col_sets;
    m.col_idx = col_idx;
    m.nonbasic_counter = nonbasic_counter;
    m.wasbasic = wasbasic;
    m.n_inner = n_inner;
    m.n_disjoint = n_disjoint;
    m.n_cut_rows = n_cut_rows;
    m.n_branch_rows = n_branch_rows;
    m.col_inner = col_inner;
    m.hs = std::make_unique<Highs>();
    m.hs->setOptionValue("output_flag", false);
    m.hs->setOptionValue("threads", 1);
    HighsModel model = hs->getModel();     // deep copy of the LP (bounds/costs/matrix/offset/integrality)
    m.hs->passModel(std::move(model));
    return m;
  }

  // Register a column's metadata + (optionally) its disjoint-row coefficients via addCol.
  // With with_rows=false the column is added with no matrix entries (used during construction,
  // where rows don't exist yet); with_rows=true it is the incremental mp_add_columns! path.
  void register_column(const BitSet& s, bool with_rows) {
    int r = (int)col_sets.size();
    col_sets.push_back(s);
    col_idx[s] = r;
    nonbasic_counter.push_back(0);
    wasbasic.push_back(0);
    // inner_vs(T[i], s) per tree.
    std::vector<std::vector<int>> inner(maf->T.size());
    for (size_t i = 0; i < maf->T.size(); ++i) inner[i] = inner_vs(maf->T[i], s);
    col_inner.push_back(inner);

    double cost = -(double)(s.size() - 1);
    if (with_rows) {
      std::vector<HighsInt> idx;
      for (size_t i = 0; i < maf->T.size(); ++i)
        for (int v : inner[i]) idx.push_back((HighsInt)row_of((int)i, v));
      std::vector<double> ones(idx.size(), 1.0);
      hs->addCol(cost, 0.0, kHighsInf, (HighsInt)idx.size(), idx.data(), ones.data());
    } else {
      HighsInt z = 0; double zv = 0.0;
      hs->addCol(cost, 0.0, kHighsInf, 0, &z, &zv);
    }
    if (isbinary) hs->changeColIntegrality(real_base + r, HighsVarType::kInteger);
  }

  int row_of(int i, int v) const { return i * n_inner + (v - (maf->n + 1)); }

  // mp_add_columns!: add new columns incrementally (disjoint-row coefficients included).
  // Returns the new columns' real indices (in creation order).
  std::vector<int> add_columns(const std::vector<BitSet>& new_vars) {
    std::vector<int> added;
    for (const BitSet& s : new_vars) {
      if (col_idx.count(s)) continue;  // already present
      added.push_back((int)col_sets.size());
      register_column(s, /*with_rows=*/true);
    }
    return added;
  }

  // ---- cut rows (clique cuts, appended after the n_disjoint disjoint rows) ------------------
  // add_subset_row_cut!: add a `<= rhs` row with coefficient 1.0 on each clique column.
  // Returns the cut index r (its HiGHS row index is n_disjoint + r).
  int add_cut_row(const std::vector<BitSet>& columns, double rhs) {
    std::vector<HighsInt> idx;
    for (const BitSet& col : columns) idx.push_back((HighsInt)(real_base + col_idx.at(col)));
    std::vector<double> ones(idx.size(), 1.0);
    hs->addRow(-kHighsInf, rhs, (HighsInt)idx.size(), idx.data(), ones.data());
    return n_cut_rows++;
  }

  // set_normalized_coefficient(cut row r, column, value): used when pricing extends a clique cut
  // with a newly-priced covering column.
  void set_cut_coefficient(int r, const BitSet& col, double value) {
    hs->changeCoeff((HighsInt)(n_disjoint + r), (HighsInt)(real_base + col_idx.at(col)), value);
  }

  // dual(ctx.cons[i]) — the LP dual of cut row r (the clique-cut penalty, <= 0 for a Min master).
  double cut_dual(int r) const { return last_row_dual[n_disjoint + r]; }

  // ---- branch rows (a3d68e4) ----------------------------------------------------------------
  // Real columns (r) covering inner vertex v in tree tree0 (0-based tree index), ascending r.
  std::vector<int> cols_covering(int tree0, int v) const {
    std::vector<int> r;
    for (int c = 0; c < (int)col_sets.size(); ++c) {
      const std::vector<int>& iv = col_inner[c][tree0];  // ascending
      if (std::binary_search(iv.begin(), iv.end(), v)) r.push_back(c);
    }
    return r;
  }

  // Add a branch row over the given real columns.  ge1 => `sum >= 1` (lower=1); else `sum <= 0`
  // (upper=0).  Returns the ABSOLUTE HiGHS row index (stable across clone).  n_branch_rows is
  // incremented only for real branch rows (track=true); strong-branching throwaway rows pass
  // track=false but still return their absolute index so they can be relaxed.
  int add_branch_row(const std::vector<int>& real_cols, bool ge1, bool track = true) {
    std::vector<HighsInt> idx; idx.reserve(real_cols.size());
    for (int c : real_cols) idx.push_back((HighsInt)(real_base + c));
    std::vector<double> ones(idx.size(), 1.0);
    double lower = ge1 ? 1.0 : -kHighsInf;
    double upper = ge1 ? kHighsInf : 0.0;
    hs->addRow(lower, upper, (HighsInt)idx.size(), idx.data(), ones.data());
    if (track) ++n_branch_rows;
    return (int)hs->getNumRow() - 1;   // absolute index of the row just added
  }

  // set_normalized_rhs analogue on an absolute row (used by strong branching to relax a trial row).
  void set_row_bounds_abs(int abs_row, double lower, double upper) {
    hs->changeRowBounds((HighsInt)abs_row, lower, upper);
  }
  // set_normalized_coefficient on an absolute row (branch-row coeff for a newly-priced column).
  void set_row_coefficient_abs(int abs_row, int real_col, double value) {
    hs->changeCoeff((HighsInt)abs_row, (HighsInt)(real_base + real_col), value);
  }
  double row_dual_abs(int abs_row) const { return last_row_dual[abs_row]; }

  void set_integral_all() {
    if (col_sets.empty()) return;
    hs->changeColsIntegrality(real_base, real_base + (HighsInt)col_sets.size() - 1, nullptr);
    std::vector<HighsVarType> t(col_sets.size(), HighsVarType::kInteger);
    hs->changeColsIntegrality(real_base, real_base + (HighsInt)col_sets.size() - 1, t.data());
  }

  // MIP warm-start. Seed HiGHS with a known feasible incumbent so the MIP solve can prune faster.
  // Every pooled column that is in the incumbent partition starts at 1.0, all other columns at 0.0.
  // The incumbent is a
  // valid partition => its columns are Steiner-disjoint => the start is packing-feasible.  PURE
  // SPEEDUP: the MIP OPTIMUM (= k*) is start-independent, so this never changes the objective or the
  // BCP control flow (which depends only on ctx_primal = k*); it can at most change WHICH of several
  // equally-optimal partitions HiGHS returns — an already-accepted, both-valid divergence.  Call after
  // the columns are registered (i.e. right after construction) and before solve().
  void set_mip_start(const std::vector<BitSet>& incumbent) {
    HighsSolution s;
    s.col_value.assign((size_t)hs->getNumCol(), 0.0);
    for (const BitSet& c : incumbent) {
      auto it = col_idx.find(c);
      if (it != col_idx.end()) s.col_value[real_base + it->second] = 1.0;
    }
    hs->setSolution(s);
  }

  // mp_solve!: optimize, capture solution, update basis bookkeeping.
  // Returns the objective; per-column values and row duals are stored for the accessors.
  double solve() {
    ++mp_solve_calls;
    // Cap this solve at the remaining time budget so we cannot blow past the deadline inside HiGHS.
    check_budget();
    hs->setOptionValue("time_limit", remaining_seconds());
    hs->run();
    // Exact-track safety: a non-optimal master must not be trusted.  Distinguish the two causes:
    //   * kInfeasible  -> this branch-and-price NODE is infeasible; throw NodeInfeasible so ctx_solve
    //     PRUNES it (standard B&P) instead of aborting.  (Cannot occur at the root — the packing
    //     master is always feasible at a=0 — so it never changes an instance that has no branching;
    //     in particular it is inert on all exact2 instances, which never create an infeasible node.)
    //   * anything else (time limit / error) -> genuine interruption; emit nothing.
    if (hs->getModelStatus() != HighsModelStatus::kOptimal) {
      if (hs->getModelStatus() == HighsModelStatus::kInfeasible && !budget_exceeded())
        throw NodeInfeasible{};
      throw Interrupted{};
    }
    const HighsInfo& info = hs->getInfo();
    const HighsSolution& sol = hs->getSolution();
    last_obj = info.objective_function_value;
    last_col_value = sol.col_value;
    last_row_dual = sol.row_dual;

    // nonbasic_counter++ for all; reset to 0 and mark wasbasic for basic columns.  last_col_basic
    // captures the current basis.
    const HighsBasis& basis = hs->getBasis();
    bool have_basis = basis.valid && (int)basis.col_status.size() >= real_base + (int)col_sets.size();
    last_col_basic.assign(col_sets.size(), 0);
    for (int r = 0; r < (int)col_sets.size(); ++r) {
      nonbasic_counter[r] += 1;
      if (have_basis && basis.col_status[real_base + r] == HighsBasisStatus::kBasic) {
        nonbasic_counter[r] = 0;
        wasbasic[r] = 1;
        last_col_basic[r] = 1;
      }
    }
    return last_obj;
  }

  // Whether real column r was BASIC in the last solve.
  bool col_is_basic(int r) const { return r >= 0 && r < (int)last_col_basic.size() && last_col_basic[r]; }

  // ---- accessors on the last solve ----------------------------------------------------------
  double sol_value(int r) const { return last_col_value[real_base + r]; }
  double sol_value(const BitSet& s) const { return last_col_value[real_base + col_idx.at(s)]; }
  // beta(i,v) = -dual(disjoint[i,v]) — the pricing penalty.
  double beta(int i, int v) const { return -last_row_dual[row_of(i, v)]; }
  double objective() const { return last_obj; }
};

} // namespace maf
