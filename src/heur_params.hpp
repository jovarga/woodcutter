// heur_params.hpp — runtime-tunable heuristic-track parameters (env-overridable).
//
// This header collects the hard-coded knobs of heuristics.hpp::run_heuristic / feasibility_pump into
// one struct read ONCE from the environment, with the current literals as defaults — so the DEFAULT
// behaviour is byte-for-byte unchanged and nothing needs a recompile to tune. A run_heuristic call
// reads `heur_params().<field>` at its use sites; set the corresponding MAF_HEUR_* env var to
// override.
//
// The lower-bound track (harness_lowerbound.hpp) calls feasibility_pump with its OWN explicit args and
// is intentionally left untouched here; it can adopt heur_params() if wanted.
#pragma once
#include <cstdlib>
#include <cstdint>

namespace maf { namespace heur {

struct HeurParams {
  // run_heuristic exact/heuristic transition: n <= this => mimic the exact solver on the (reduced or
  // cluster) block. The single biggest lever (per the heuristic review). Env: MAF_HEUR_MAX_TAXA_EXACT.
  int    max_taxa_exact;
  // feasibility_pump: no-improve stall limit (iterations), pump penalty per used column, base
  // LP size-weight and its restart-alternate value, base jitter and per-restart jitter increment.
  int    fp_stall;            // MAF_HEUR_FP_STALL           (25)
  double fp_pump_penalty;     // MAF_HEUR_FP_PUMP_PENALTY     (1.0)
  double fp_size_weight;      // MAF_HEUR_FP_SIZE_WEIGHT      (0.0)  base (single-pass + even restarts)
  double fp_size_weight_alt;  // MAF_HEUR_FP_SIZE_WEIGHT_ALT  (1.0)  odd-restart size-aware value
  double fp_jitter;           // MAF_HEUR_FP_JITTER           (0.4)  base tie-break magnitude
  double fp_jitter_step;      // MAF_HEUR_FP_JITTER_STEP      (0.2)  jitter += step*(restart%3)
  // fast-floor (mast_peel) column cap + wall-clock cap (fraction of remaining, ceilinged in seconds).
  int    floor_max_cols;      // MAF_HEUR_FLOOR_MAX_COLS      (8000)
  double floor_cap_frac;      // MAF_HEUR_FLOOR_CAP_FRAC      (0.3)
  double floor_cap_max;       // MAF_HEUR_FLOOR_CAP_MAX       (60.0)
  // fraction of the remaining budget the root column generation (ctx_process_root) may use.
  double root_cg_frac;        // MAF_HEUR_ROOT_CG_FRAC        (0.5)
  // Experimental pool enrichment, disabled by default.
  int    cg_all_cols;         // MAF_HEUR_CG_ALL_COLS         (0)   add every root-master column to FP pool
  int    cg_all_cap;          // MAF_HEUR_CG_ALL_CAP          (0)   0 = uncapped
  int    cg_all_min_size;     // MAF_HEUR_CG_ALL_MIN_SIZE     (2)
  int    cg_all_sort_size;    // MAF_HEUR_CG_ALL_SORT_SIZE    (0)   add larger columns first
  double near_rc_min_score;   // MAF_HEUR_NEAR_RC             (0.0) disabled; add priced cols with score >= this
  int    near_rc_cap;         // MAF_HEUR_NEAR_RC_CAP         (4000)
  int    dual_jitter_passes;  // MAF_HEUR_DUAL_JITTER         (0)   disabled
  double dual_jitter_mag;     // MAF_HEUR_DUAL_JITTER_MAG     (0.05)
  double dual_jitter_min_score; // MAF_HEUR_DUAL_JITTER_MIN   (1.0)
  int    dual_jitter_cap;     // MAF_HEUR_DUAL_JITTER_CAP     (4000)
  int    pool_log;            // MAF_HEUR_POOL_LOG            (0)   stderr pool stats
};

inline const HeurParams& heur_params() {
  static const HeurParams p = [] {
    auto ei = [](const char* k, int d)    { const char* e = std::getenv(k); return e && *e ? std::atoi(e) : d; };
    auto ed = [](const char* k, double d) { const char* e = std::getenv(k); return e && *e ? std::atof(e) : d; };
#ifdef MAX_TAXA_COUNT__EXACT_APPROACH
    const int def_max_taxa = MAX_TAXA_COUNT__EXACT_APPROACH;   // stay in sync with the #define default
#else
    const int def_max_taxa = 2000;
#endif
    HeurParams h;
    h.max_taxa_exact      = ei("MAF_HEUR_MAX_TAXA_EXACT",     def_max_taxa);
    h.fp_stall            = ei("MAF_HEUR_FP_STALL",           25);
    h.fp_pump_penalty     = ed("MAF_HEUR_FP_PUMP_PENALTY",    1.0);
    h.fp_size_weight      = ed("MAF_HEUR_FP_SIZE_WEIGHT",     0.0);
    h.fp_size_weight_alt  = ed("MAF_HEUR_FP_SIZE_WEIGHT_ALT", 1.0);
    h.fp_jitter           = ed("MAF_HEUR_FP_JITTER",          0.4);
    h.fp_jitter_step      = ed("MAF_HEUR_FP_JITTER_STEP",     0.2);
    h.floor_max_cols      = ei("MAF_HEUR_FLOOR_MAX_COLS",     8000);
    h.floor_cap_frac      = ed("MAF_HEUR_FLOOR_CAP_FRAC",     0.3);
    h.floor_cap_max       = ed("MAF_HEUR_FLOOR_CAP_MAX",      60.0);
    h.root_cg_frac        = ed("MAF_HEUR_ROOT_CG_FRAC",       0.5);
    h.cg_all_cols         = ei("MAF_HEUR_CG_ALL_COLS",        0);
    h.cg_all_cap          = ei("MAF_HEUR_CG_ALL_CAP",         0);
    h.cg_all_min_size     = ei("MAF_HEUR_CG_ALL_MIN_SIZE",    2);
    h.cg_all_sort_size    = ei("MAF_HEUR_CG_ALL_SORT_SIZE",   0);
    h.near_rc_min_score   = ed("MAF_HEUR_NEAR_RC",            0.0);
    h.near_rc_cap         = ei("MAF_HEUR_NEAR_RC_CAP",        4000);
    h.dual_jitter_passes  = ei("MAF_HEUR_DUAL_JITTER",        0);
    h.dual_jitter_mag     = ed("MAF_HEUR_DUAL_JITTER_MAG",    0.05);
    h.dual_jitter_min_score = ed("MAF_HEUR_DUAL_JITTER_MIN",  1.0);
    h.dual_jitter_cap     = ei("MAF_HEUR_DUAL_JITTER_CAP",    4000);
    h.pool_log            = ei("MAF_HEUR_POOL_LOG",           0);
    return h;
  }();
  return p;
}

}} // namespace maf::heur
