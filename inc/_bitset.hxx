#pragma once
#include <utility>
#include <vector>
#include <algorithm>
#include <cstdint>
#include "_algorithm.hxx"
#include "_ctypes.hxx"

using std::pair;
using std::vector;
using std::sort;
using std::lower_bound;




#pragma region CLASSES
/**
 * A Lazy bitset is a sparse integer key to value map that updates insertions
 * and deletions upon calling update(). It maintains keys in ascending order.
 * @tparam K key type
 * @tparam V value type
 */
template <class K=uint32_t, class V=NONE>
class LazyBitset {
  #pragma region TYPES
  public:
  /** The key type. */
  using key_type   = K;
  /** The value type. */
  using value_type = V;
  /** The entry type. */
  using entry_type = pair<K, V>;
  #pragma endregion


  #pragma region DATA
  protected:
  /** The pairs of keys and values. */
  vector<pair<K, V>> pairs;
  /** The number of unprocessed insertions and deletions (-ve). */
  ssize_t unprocessed;
  #pragma endregion


  #pragma region METHODS
  #pragma region SIZE
  public:
  /**
   * Get the number of entries in the bitset.
   * @returns |this|
   */
  inline size_t size() const noexcept {
    return pairs.size();  //  - abs(unprocessed)
  }

  /**
   * Check if the bitset is empty.
   * @returns |this| == 0
   */
  inline bool empty() const noexcept {
    return size() == 0;
  }
  #pragma endregion


  #pragma region AT
  public:
  /**
   * Get the entry at given index.
   * @param i index
   * @returns this[i]
   */
  inline pair<K, V> at(size_t i) const noexcept {
    return pairs[i];
  }

  /**
   * Get the key at given index.
   * @param i index
   * @returns this[i].first
   */
  inline K keyAt(size_t i) const noexcept {
    return pairs[i].first;
  }

  /**
   * Get the value at given index.
   * @param i index
   * @returns this[i].second
   */
  inline V valueAt(size_t i) const noexcept {
    return pairs[i].second;
  }
  #pragma endregion


  #pragma region FOREACH
  public:
  /**
   * Iterate over the entries in the bitset.
   * @param fp process function (key, value)
   */
  template <class F>
  inline void forEach(F fp) const noexcept {
    for (const auto& p : pairs)
      fp(p.first, p.second);
  }

  /**
   * Iterate over the keys in the bitset.
   * @param fp process function (key)
   */
  template <class F>
  inline void forEachKey(F fp) const noexcept {
    for (const auto& p : pairs)
      fp(p.first);
  }
  #pragma endregion


  #pragma region ITERATOR
  protected:
  /**
   * Get const iterator to the first entry.
   * @returns const& this[0]
   */
  inline auto begin() const noexcept {
    return pairs.begin();
  }

  /**
   * Get iterator to the first entry.
   * @returns & this[0]
   */
  inline auto begin() noexcept {
    return pairs.begin();
  }

  /**
   * Get const iterator to the end.
   * @returns const& this[|this|]
   */
  inline auto end() const noexcept {
    return pairs.end();  // pairs.begin() + size();
  }

  /**
   * Get iterator to the end.
   * @returns & this[|this|]
   */
  inline auto end() noexcept {
    return pairs.end();  // pairs.begin() + size();
  }
  #pragma endregion


  #pragma region FIND
  protected:
  /**
   * Find the entry with given key.
   * @param k key
   * @returns const& this[i] where this[i].first == k, or & this[|this|] if not found
   */
  inline auto find(K k) const noexcept {
    auto   fl = [](const auto& p, K k) { return p.first < k; };
    auto   it = lower_bound(begin(), end(), k, fl);
    return it == end() || (*it).first != k? end() : it;
  }

  /**
   * Find the entry with given key.
   * @param k key
   * @returns & this[i] where this[i].first == k, or & this[|this|] if not found
   */
  inline auto find(K k) noexcept {
    auto   fl = [](const auto& p, K k) { return p.first < k; };
    auto   it = lower_bound(begin(), end(), k, fl);
    return it == end() || (*it).first != k? end() : it;
  }
  #pragma endregion


  #pragma region ACCESS
  public:
  /**
   * Check if the bitset has an entry with given key.
   * @param k key
   * @returns has this[k]?
   */
  inline bool has(K k) const noexcept {
    auto   it = find(k);
    return it != end();
  }

  /**
   * Get the value of the entry with given key.
   * @param k key
   * @returns this[k]
   */
  inline V get(K k) const noexcept {
    auto   it = find(k);
    return it == end()? V() : (*it).second;
  }

  /**
   * Set the value of the entry with given key, if key is present.
   * @param k key
   * @param v value
   * @returns success?
   */
  inline bool set(K k, V v) noexcept {
    auto it = find(k);
    if  (it == end()) return false;
    (*it).second = v;
    return true;
  }
  #pragma endregion


  #pragma region UPDATE
  protected:
  /**
   * Update the bitset by sorting out all unprocessed deletions.
   * @note This is an expensive operation.
   */
  inline void updateRemove() {
    auto  fl = [](const auto& p, const auto& q) { return p.first <  q.first; };
    auto  fe = [](const auto& p, const auto& q) { return p.first == q.first; };
    size_t N = pairs.size();
    size_t n = N + unprocessed;
    auto  ib = pairs.begin();
    auto  ie = pairs.end();
    auto  im = ib + n;
    sort(im, ie, fl);
    auto  it = set_difference_inplace(ib, im, im, ie, fl, fe);
    pairs.resize(it - ib);
    unprocessed = 0;
  }

  /**
   * Update the bitset by sorting out all unprocessed insertions.
   * @note This is an expensive operation.
   */
  inline void updateAdd(vector<pair<K, V>> *buf=nullptr) {
    auto  fl = [](const auto& p, const auto& q) { return p.first <  q.first; };
    auto  fe = [](const auto& p, const auto& q) { return p.first == q.first; };
    size_t N = pairs.size();
    size_t n = N - unprocessed;
    size_t need  = unprocessed + 4;
    if (!buf)  pairs.resize(N + need);
    else if (buf->size() < need) buf->resize(need);
    auto bb = buf? buf->begin() : pairs.begin() + N;
    auto be = buf? buf->end()   : pairs.end();
    auto ib = pairs.begin();
    auto im = ib + n;
    auto ie = ib + N;
    sort(im, ie, fl);
    auto it = set_union_last_inplace(ib, im, im, ie, bb, be, fl, fe);
    pairs.resize(it - ib);
    unprocessed = 0;
  }

  public:
  /**
   * Reserve space for n entries.
   * @param n number of entries
   */
  inline void reserve(size_t n) {
    pairs.reserve(n);
  }

  /**
   * Clear the bitset.
   */
  inline void clear() noexcept {
    pairs.clear();
    unprocessed = 0;
  }

  /**
   * Update the bitset by sorting out all unprocessed insertions and deletions.
   * @note This is an expensive operation.
   */
  inline void update(vector<pair<K, V>> *buf=nullptr) {
    if (unprocessed == 0) return;
    if (unprocessed  < 0) updateRemove();
    else updateAdd(buf);
  }

  /**
   * Remove the entry with given key.
   * @param k key
   * @param buf buffer for unprocessed insertions
   * @note This operation is lazy.
   */
  inline void remove(K k, vector<pair<K, V>> *buf=nullptr) {
    if (unprocessed > 0) updateAdd(buf);
    pairs.push_back({k, V()});
    --unprocessed;
  }

  /**
   * Add the entry with given key.
   * @param k key
   * @param v value
   * @note This operation is lazy.
   */
  inline void add(K k, V v=V()) {
    if (unprocessed < 0) updateRemove();
    pairs.push_back({k, v});
    ++unprocessed;
  }
  #pragma endregion
  #pragma endregion
};
#pragma endregion




#pragma region METHODS
#pragma region WRITE
/**
 * Write a bitset to a stream.
 * @tparam B bitset type
 * @param a stream
 * @param x bitset
 */
template <class B>
void writeBitset(ostream& a, const B& x) {
  a << "{\n";
  x.forEach([&](auto k, auto v) {
    a << "  " << k << ": " << v << "\n";
  });
  a << "}";
}

/**
 * Write a lazy bitset to a stream.
 * @tparam K key type
 * @tparam V value type
 * @param a stream
 * @param x bitset
 */
template <class K, class V>
inline void write(ostream& a, const LazyBitset<K, V>& x) {
  writeBitset(a, x);
}

/**
 * Write a lazy bitset to a stream.
 * @tparam K key type
 * @tparam V value type
 * @param a stream
 * @param x bitset
 */
template <class K, class V>
inline ostream& operator<<(ostream& a, const LazyBitset<K, V>& x) {
  write(a, x);
  return a;
}
#pragma endregion
#pragma endregion
