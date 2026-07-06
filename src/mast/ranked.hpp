// Ranked-solution machinery, blacklist handling, and score tie-breaking shared by MAST routines.
#pragma once
#include "trees.hpp"
#include <vector>
#include <unordered_set>
#include <algorithm>

namespace maf::mast {

// Strictly greater score wins; on ==, lexicographically smaller taxa wins; an invalid best is beaten
// by any valid candidate.
inline bool better_score(W cs, const std::vector<TaxonId>& ct, W bs, const std::vector<TaxonId>& bt) {
  if (!is_valid_score(bs)) return is_valid_score(cs);
  if (!is_valid_score(cs)) return false;
  return cs > bs || (cs == bs && lexless(ct, bt));
}

// ---- blacklist (SOL:238-257): a set keyed by the sorted taxa sequence -------------------------
struct VecHash {
  size_t operator()(const std::vector<int>& v) const {
    size_t h = 1469598103934665603ULL;
    for (int x : v) { h ^= (size_t)x; h *= 1099511628211ULL; }
    return h;
  }
};
using BlacklistSet = std::unordered_set<std::vector<int>, VecHash>;

inline BlacklistSet normalize_blacklist(const std::vector<std::vector<int>>& blacklist) {
  BlacklistSet s;
  for (std::vector<int> taxa : blacklist) { std::sort(taxa.begin(), taxa.end()); s.insert(std::move(taxa)); }
  return s;
}
inline bool is_blacklisted(const BlacklistSet& s, const std::vector<TaxonId>& taxa) {
  return s.find(taxa) != s.end();   // taxa is always sorted at call sites
}

// _insert_unique_solution_score_bounded! (SOL:328): keep `s` sorted by DESCENDING score, bounded to
// `limit`, dropping the lowest when full.  No taxa dedup (caller guarantees uniqueness of interest).
inline void insert_unique_solution_score_bounded(std::vector<MASTSolution>& s, MASTSolution sol, int limit) {
  if (limit < 1) return;
  if ((int)s.size() >= limit && sol.score <= s.back().score) return;
  int idx;
  if ((int)s.size() < limit) { s.push_back(MASTSolution{}); idx = (int)s.size() - 1; }
  else idx = (int)s.size() - 1;                      // full: overwrite the (now-dropped) last slot
  while (idx > 0 && sol.score > s[idx - 1].score) { s[idx] = s[idx - 1]; --idx; }
  s[idx] = std::move(sol);
}

// ---- _ScoreProductItem k-way merge heap (SOL:408-472) -----------------------------------------
// Array binary max-heap (1-based-style indexing emulated; here 0-based vectors with manual sift).
struct ScoreProductItem { W score; int left_idx; int right_idx; };
// precedes (SOL:414): higher score, then smaller left_idx, then smaller right_idx.
inline bool sp_precedes(const ScoreProductItem& a, const ScoreProductItem& b) {
  if (a.score != b.score) return a.score > b.score;
  if (a.left_idx != b.left_idx) return a.left_idx < b.left_idx;
  return a.right_idx < b.right_idx;
}
inline void sp_sift_up(std::vector<ScoreProductItem>& h, int idx) {   // idx is 0-based
  while (idx > 0) {
    int parent = (idx - 1) / 2;
    if (!sp_precedes(h[idx], h[parent])) break;
    std::swap(h[idx], h[parent]); idx = parent;
  }
}
inline void sp_sift_down(std::vector<ScoreProductItem>& h, int idx) {
  int n = (int)h.size();
  while (true) {
    int left = 2 * idx + 1;
    if (left >= n) break;
    int right = left + 1;
    int best = (right < n && sp_precedes(h[right], h[left])) ? right : left;
    if (!sp_precedes(h[best], h[idx])) break;
    std::swap(h[idx], h[best]); idx = best;
  }
}
inline void sp_push(std::vector<ScoreProductItem>& h, ScoreProductItem it) { h.push_back(it); sp_sift_up(h, (int)h.size() - 1); }
inline ScoreProductItem sp_pop(std::vector<ScoreProductItem>& h) {
  ScoreProductItem top = h[0];
  ScoreProductItem last = h.back(); h.pop_back();
  if (!h.empty()) { h[0] = last; sp_sift_down(h, 0); }
  return top;
}

// _combine_unique_solution_lists_score_bounded! (SOL:549): best-first merge of left x right (each
// sorted desc by score) into `dest`, taxa via sorted union, +score_delta, bounded to `limit`.
inline void combine_unique_solution_lists_score_bounded(
    std::vector<MASTSolution>& dest, const std::vector<MASTSolution>& left,
    const std::vector<MASTSolution>& right, W score_delta, int limit, const SolutionRef& ref) {
  if (left.empty() || right.empty() || limit < 1) return;
  std::vector<ScoreProductItem> heap;
  heap.reserve(left.size());
  for (int li = 0; li < (int)left.size(); ++li)
    sp_push(heap, {left[li].score + right[0].score + score_delta, li, 0});

  while (!heap.empty()) {
    ScoreProductItem item = sp_pop(heap);
    if ((int)dest.size() >= limit && item.score <= dest.back().score) break;
    const MASTSolution& ls = left[item.left_idx];
    const MASTSolution& rs = right[item.right_idx];
    std::vector<TaxonId> taxa = sorted_taxa_union(ls.taxa, rs.taxa);
    insert_unique_solution_score_bounded(dest, MASTSolution{item.score, std::move(taxa), ref}, limit);
    int nri = item.right_idx + 1;
    if (nri < (int)right.size())
      sp_push(heap, {ls.score + right[nri].score + score_delta, item.left_idx, nri});
  }
}

} // namespace maf::mast
