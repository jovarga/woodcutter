// core.hpp — foundational primitives shared by the solver.
// Provides the two pervasive data structures:
//   * maf::BitSet  (set of ints, ASCENDING iteration, word-backed)
//   * maf::DiGraph (directed graph with sorted adjacency)
#pragma once
#include <cstdint>
#include <vector>
#include <initializer_list>
#include <algorithm>
#include <functional>
#include <cassert>

namespace maf {

// The weight/score scalar type used throughout the solver (production fixes Float64).
// Also re-declared as maf::mast::W inside the MAST headers; both alias double, so unqualified `W`
// resolves consistently in either namespace.
using W = double;

// ---------------------------------------------------------------------------------------------
// BitSet: elements are nonnegative ints; iteration is in ASCENDING order, and downstream code relies
// on that for determinism.
// Backed by 64-bit words.  The word array is kept "trimmed" (no trailing zero words) so that
// operator== and std::hash are consistent for equal sets.
// ---------------------------------------------------------------------------------------------
class BitSet {
  std::vector<uint64_t> w_;
  static constexpr int WB = 64;

  void ensure_word(size_t nw) { if (w_.size() < nw) w_.resize(nw, 0); }
  void trim() { while (!w_.empty() && w_.back() == 0) w_.pop_back(); }

public:
  BitSet() = default;
  BitSet(std::initializer_list<int> xs) { for (int x : xs) insert(x); }
  template <class It> BitSet(It b, It e) { for (; b != e; ++b) insert((int)*b); }

  bool contains(int x) const {
    if (x < 0) return false;
    size_t i = (size_t)x / WB;
    return i < w_.size() && ((w_[i] >> (x % WB)) & 1ULL);
  }
  void insert(int x) { assert(x >= 0); ensure_word((size_t)x / WB + 1); w_[x / WB] |= (1ULL << (x % WB)); }
  void erase(int x) {
    if (x < 0) return;
    size_t i = (size_t)x / WB;
    if (i < w_.size()) { w_[i] &= ~(1ULL << (x % WB)); trim(); }
  }

  bool empty() const { return w_.empty(); }                 // trimmed => empty iff no words
  int size() const { int c = 0; for (uint64_t x : w_) c += __builtin_popcountll(x); return c; }

  // Ascending iteration helpers.
  template <class F> void for_each(F f) const {
    for (size_t i = 0; i < w_.size(); ++i) {
      uint64_t b = w_[i];
      while (b) { int t = __builtin_ctzll(b); f((int)(i * WB + t)); b &= b - 1; }
    }
  }
  std::vector<int> to_vec() const { std::vector<int> v; v.reserve(size()); for_each([&](int x){ v.push_back(x); }); return v; }
  int first() const {  // smallest element (caller ensures nonempty)
    for (size_t i = 0; i < w_.size(); ++i) if (w_[i]) return (int)(i * WB + __builtin_ctzll(w_[i]));
    return -1;
  }

  // In-place set operations.
  void union_with(const BitSet& o) {
    ensure_word(o.w_.size());
    for (size_t i = 0; i < o.w_.size(); ++i) w_[i] |= o.w_[i];
  }
  void intersect_with(const BitSet& o) {
    for (size_t i = 0; i < w_.size(); ++i) w_[i] &= (i < o.w_.size() ? o.w_[i] : 0ULL);
    trim();
  }
  void setdiff_with(const BitSet& o) {
    for (size_t i = 0; i < w_.size() && i < o.w_.size(); ++i) w_[i] &= ~o.w_[i];
    trim();
  }

  friend bool isdisjoint(const BitSet& a, const BitSet& b) {
    size_t n = std::min(a.w_.size(), b.w_.size());
    for (size_t i = 0; i < n; ++i) if (a.w_[i] & b.w_[i]) return false;
    return true;
  }

  bool operator==(const BitSet& o) const { return w_ == o.w_; }   // both trimmed => exact
  bool operator!=(const BitSet& o) const { return w_ != o.w_; }
  // Canonical total order by ascending element sequence (smallest differing element decides).
  bool operator<(const BitSet& o) const {
    size_t n = std::max(w_.size(), o.w_.size());
    for (size_t i = 0; i < n; ++i) {
      uint64_t a = i < w_.size() ? w_[i] : 0ULL, b = i < o.w_.size() ? o.w_[i] : 0ULL;
      if (a != b) { int t = __builtin_ctzll(a ^ b); return (b >> t) & 1ULL; }
    }
    return false;
  }

  const std::vector<uint64_t>& words() const { return w_; }

  // Ascending forward iterator so `for (int x : bitset)` works.
  // Elements are >= 0, so val_ == -1 uniquely marks end().
  class const_iterator {
    const std::vector<uint64_t>* w_; size_t wi_; uint64_t bits_; int val_;
    void advance() {
      while (bits_ == 0) {
        if (++wi_ >= w_->size()) { val_ = -1; return; }
        bits_ = (*w_)[wi_];
      }
      val_ = (int)(wi_ * WB + __builtin_ctzll(bits_));
      bits_ &= bits_ - 1;
    }
  public:
    const_iterator(const std::vector<uint64_t>* w, bool end) : w_(w), wi_(0), bits_(0), val_(-1) {
      if (!end && !w_->empty()) { bits_ = (*w_)[0]; advance(); }
    }
    int operator*() const { return val_; }
    const_iterator& operator++() { advance(); return *this; }
    bool operator!=(const const_iterator& o) const { return val_ != o.val_; }
    bool operator==(const const_iterator& o) const { return val_ == o.val_; }
  };
  const_iterator begin() const { return const_iterator(&w_, false); }
  const_iterator end()   const { return const_iterator(&w_, true); }
};

inline BitSet set_union(const BitSet& a, const BitSet& b)     { BitSet r = a; r.union_with(b);     return r; }
inline BitSet set_intersect(const BitSet& a, const BitSet& b) { BitSet r = a; r.intersect_with(b); return r; }
inline BitSet set_diff(const BitSet& a, const BitSet& b)      { BitSet r = a; r.setdiff_with(b);   return r; }

// ---------------------------------------------------------------------------------------------
// DiGraph: minimal directed graph used by the solver. Vertices are 1-based (1..nv). Adjacency lists
// are kept SORTED ascending after finalize() so traversal order is deterministic.
// ---------------------------------------------------------------------------------------------
class DiGraph {
  int nv_ = 0;
  std::vector<std::vector<int>> out_, in_;

public:
  DiGraph() = default;
  explicit DiGraph(int n) : nv_(n), out_(n + 1), in_(n + 1) {}

  int nv() const { return nv_; }
  void add_edge(int u, int v) { out_[u].push_back(v); in_[v].push_back(u); }
  void finalize() {  // keep adjacency sorted
    for (auto& a : out_) std::sort(a.begin(), a.end());
    for (auto& a : in_)  std::sort(a.begin(), a.end());
  }

  const std::vector<int>& out_neighbors(int v) const { return out_[v]; }
  const std::vector<int>& in_neighbors(int v) const { return in_[v]; }
  int out_degree(int v) const { return (int)out_[v].size(); }
  int in_degree(int v) const { return (int)in_[v].size(); }
  std::vector<int> all_neighbors(int v) const {
    std::vector<int> r = in_[v];
    r.insert(r.end(), out_[v].begin(), out_[v].end());
    return r;
  }
  // For a directed tree: parent (single in-neighbor) or 0 if root.
  int parent(int v) const { return in_[v].empty() ? 0 : in_[v][0]; }
};

} // namespace maf

template <> struct std::hash<maf::BitSet> {
  size_t operator()(const maf::BitSet& s) const {
    size_t h = 1469598103934665603ULL;
    for (uint64_t w : s.words()) { h ^= std::hash<uint64_t>()(w); h *= 1099511628211ULL; }
    return h;
  }
};
