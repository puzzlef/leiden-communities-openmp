#pragma once
#include <cstdint>
#include <utility>
#include <algorithm>
#include <vector>
#include "_ctypes.hxx"
#include "_queue.hxx"

using std::pair;
using std::vector;
using std::lower_bound;




// HELPER MACROS
// -------------
// Helps create bitsets.

#ifndef BITSET_TYPES
#define BITSET_TYPES(K, V, pairs) \
  using key_type   = K; \
  using value_type = V; \
  using entry_type = pair<K, V>; \
  using iterator       = decltype(pairs.begin()); \
  using const_iterator = decltype(pairs.cbegin());
#endif


#ifndef BITSET_DATA
#define BITSET_DATA(K, V, pairs) \
  inline       pair<K, V>* data()       noexcept { return pairs.data(); } \
  inline const pair<K, V>* data() const noexcept { return pairs.data(); }
#endif


#ifndef BITSET_ITERATOR
#define BITSET_ITERATOR(K, V, pairs) \
  inline auto begin() noexcept { return pairs.begin(); } \
  inline auto end()   noexcept { return pairs.end(); } \
  inline auto begin() const noexcept { return pairs.begin(); } \
  inline auto end()   const noexcept { return pairs.end(); }
#endif


#ifndef BITSET_SIZE
#define BITSET_SIZE(K, V, pairs)  \
  inline size_t size() const noexcept { return pairs.size(); } \
  inline bool empty()  const noexcept { return size() == 0; }
#endif


#ifndef BITSET_AT
#define BITSET_AT(K, V, pairs) \
  inline       auto& operator[](size_t i)       noexcept { return pairs[i]; } \
  inline const auto& operator[](size_t i) const noexcept { return pairs[i]; } \
  inline       auto&    entryAt(size_t i)       noexcept { return pairs[i]; } \
  inline const auto&    entryAt(size_t i) const noexcept { return pairs[i]; } \
  inline          K&      keyAt(size_t i)       noexcept { return pairs[i].first; } \
  inline const    K&      keyAt(size_t i) const noexcept { return pairs[i].first; } \
  inline          V&    valueAt(size_t i)       noexcept { return pairs[i].second; } \
  inline const    V&    valueAt(size_t i) const noexcept { return pairs[i].second; }
#endif


#ifndef BITSET_FCOMPARE
#define BITSET_FCOMPARE(K, V, fn, op) \
  auto fn = [](const auto& p, const auto& q) { return p.first op q.first; };

#define BITSET_FCOMPARES(K, V) \
  BITSET_FCOMPARE(K, V, fl, <) \
  BITSET_FCOMPARE(K, V, fe, ==)
#endif


#ifndef BITSET_FIND
#define BITSET_FIND(K, V) \
  inline auto find(K k)       noexcept { return locate_match(k); } \
  inline auto find(K k) const noexcept { return locate_match(k); }
#endif


#ifndef BITSET_ENTRIES
#define BITSET_ENTRIES(K, V, pairs) \
  /* dont change the keys! */ \
  inline auto values()  noexcept { return staticTransformIterable(pairs, PairSecond<K, V>()); } \
  inline auto entries() noexcept { return iterable(pairs); } \
  inline auto keys()    const noexcept { return staticTransformIterable(pairs, ConstPairFirst<K, V>());  } \
  inline auto values()  const noexcept { return staticTransformIterable(pairs, ConstPairSecond<K, V>()); } \
  inline auto entries() const noexcept { return iterable(pairs); }
#endif


#ifndef BITSET_FOREACH
#define BITSET_FOREACH(K, V, pairs) \
  /* dont change the keys! */ \
  template <class F> \
  inline void forEachValue(F fn) { for (auto& p : pairs) fn(p.second); } \
  template <class F> \
  inline void forEachEntry(F fn) { for (auto& p : pairs) fn(p); } \
  template <class F> \
  inline void forEach(F fn)      { for (auto& p : pairs) fn(p.first, p.second); } \
  template <class F> \
  inline void forEachKey(F fn)   const { for (const auto& p : pairs) fn(p.first); } \
  template <class F> \
  inline void forEachValue(F fn) const { for (const auto& p : pairs) fn(p.second); } \
  template <class F> \
  inline void forEachEntry(F fn) const { for (const auto& p : pairs) fn(p); } \
  template <class F> \
  inline void forEach(F fn)      const { for (const auto& p : pairs) fn(p.first, p.second); }
#endif


#ifndef BITSET_HAS
#define BITSET_HAS(K, V) \
  inline bool has(K k) const noexcept { \
    return locate_match(k) != end(); \
  }
#endif


#ifndef BITSET_GET
#define BITSET_GET(K, V) \
  inline V get(K k) const noexcept { \
    auto it = locate_match(k); \
    return it == end()? V() : (*it).second; \
  }
#endif


#ifndef BITSET_SET
#define BITSET_SET(K, V) \
  inline bool set(K k, V v) noexcept { \
    auto it = locate_match(k); \
    if (it == end()) return false; \
    (*it).second = v; \
    return true; \
  }
#endif




// LAZY BITSET
// -----------
// A lazy bitset that updates insertions and deletions upon calling update().
// It maintains keys in ascending order.

#ifndef LAZY_BITSET_LOCATE
#define LAZY_BITSET_LOCATE(K, V, f0, f1) \
  f0 auto locate_spot(K k) f1 { \
    auto fl = [](const auto& p, K k) { return p.first < k; }; \
    return lower_bound(begin(), end(), k, fl); \
  } \
  f0 auto locate_match(K k) f1 { \
    auto it = locate_spot(k); \
    return it == end() || (*it).first != k? end() : it; \
  }
#endif


template <class K=uint32_t, class V=NONE>
class LazyBitset {
  // Data.
  protected:
  vector<pair<K, V>> pairs;
  ssize_t unprocessed;

  // Basic operations.
  public:
  BITSET_TYPES(K, V, pairs)
  BITSET_ITERATOR(K, V, pairs)
  BITSET_SIZE(K, V, pairs)
  BITSET_AT(K, V, pairs)

  // Find operations.
  protected:
  LAZY_BITSET_LOCATE(K, V, inline, noexcept)
  LAZY_BITSET_LOCATE(K, V, inline, const noexcept)
  public:
  BITSET_FIND(K, V)

  // Access operations.
  public:
  BITSET_ENTRIES(K, V, pairs)
  BITSET_FOREACH(K, V, pairs)
  BITSET_HAS(K, V)
  BITSET_GET(K, V)
  BITSET_SET(K, V)

  // Update operations.
  protected:
  inline void updateRemove() {
    BITSET_FCOMPARES(K, V)
    size_t N = pairs.size(),  n  = N + unprocessed;
    auto  ib = pairs.begin(), ie = pairs.end(), im = ib + n;
    sort_values(im, ie, fl);
    auto  it = set_difference_inplace(ib, im, im, ie, fl, fe);
    pairs.resize(it - ib);
    unprocessed = 0;
  }

  inline void updateAdd(vector<pair<K, V>> *buf=nullptr) {
    BITSET_FCOMPARES(K, V)
    size_t N = pairs.size(), n = N - unprocessed, need = unprocessed + 4;
    if (!buf)  pairs.resize(N + need);
    else if (buf->size() < need) buf->resize(need);
    auto bb = buf? buf->begin() : pairs.begin() + N;
    auto be = buf? buf->end()   : pairs.end();
    auto ib = pairs.begin(), im = ib + n, ie = ib + N;
    sort_values(im, ie, fl);
    auto it = set_union_last_inplace(ib, im, im, ie, bb, be, fl, fe);
    pairs.resize(it - ib);
    unprocessed = 0;
  }

  public:
  inline void reserve(size_t n) {
    pairs.reserve(n);
  }

  inline void clear() noexcept {
    pairs.clear();
    unprocessed = 0;
  }

  inline void update(vector<pair<K, V>> *buf=nullptr) {
    if (unprocessed == 0) return;
    if (unprocessed  < 0) updateRemove();
    else updateAdd(buf);
  }

  inline void remove(K k, vector<pair<K, V>> *buf=nullptr) {
    if (unprocessed > 0) updateAdd(buf);
    pairs.push_back({k, V()});
    --unprocessed;
  }

  inline void add(K k, V v=V()) {
    if (unprocessed < 0) updateRemove();
    pairs.push_back({k, v});
    ++unprocessed;
  }
};

template <class K=uint32_t, class V=NONE>
inline auto lazyBitset(K _k=K(), V _v=V()) {
  return LazyBitset<K, V>();
}




// RETYPE
// ------

template <class K, class V, class KA=K, class VA=V>
constexpr auto retype(const LazyBitset<K, V>& x, KA _k=KA(), VA _v=VA()) {
  return LazyBitset<KA, VA>();
}




// WRITE
// -----

template <class B>
void writeBitset(ostream& a, const B& x) {
  a << "{\n";
  x.forEach([&](auto k, auto v) {
    a << "  " << k << ": " << v << "\n";
  });
  a << "}";
}

template <class K, class V>
inline void write(ostream& a, const LazyBitset<K, V>& x) { writeBitset(a, x); }

template <class K, class V>
inline ostream& operator<<(ostream& a, const LazyBitset<K, V>& x) { write(a, x); return a; }
