#pragma once
#include <iterator>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include "_queue.hxx"

using std::iterator_traits;
using std::vector;
using std::unordered_map;
using std::all_of;
using std::count_if;
using std::equal;
using std::copy;
using std::sort;




// FIRST
// -----
// First position.

template <class I>
inline auto first_value(I ib, I ie) {
  using T = typename iterator_traits<I>::value_type;
  T a = ib != ie? *ib : T();
  return a;
}
template <class J>
inline auto firstValue(const J& x) {
  return first_value(x.begin(), x.end());
}




// ALL OF
// ------
// Is everything there?

template <class J, class F>
inline bool allOf(const J& x, F fn) {
  return all_of(x.begin(), x.end(), fn);
}




// ADJACENT FIND
// -------------

template <class I, class FE>
auto non_adjacent_find(I ib, I ie, FE fe) {
  if (ib==ie) return ie;
  // Compare adjacent elements, and if they
  // dont match, return the first one.
  I it = ib++;
  for (; ib!=ie; ++it, ++ib)
    if (!fe(*it, *ib)) return it;
  // The last element must be non-adjacent.
  return it;
}




// EQUAL
// -----
// Check if values match.

template <class IX, class IY>
inline bool equal_values(IX xb, IX xe, IY yb) {
  return equal(xb, xe, yb);
}
template <class IX, class IY>
inline bool equal_values(IX xb, IX xe, IY yb, IY ye) {
  return equal(xb, xe, yb, ye);
}
template <class IX, class IY, class FE>
inline bool equal_values(IX xb, IX xe, IY yb, FE fe) {
  return equal(xb, xe, yb, fe);
}
template <class IX, class IY, class FE>
inline bool equal_values(IX xb, IX xe, IY yb, IY ye, FE fe) {
  return equal(xb, xe, yb, ye, fe);
}
template <class JX, class JY>
inline bool equalValues(const JX& x, const JY& y) {
  return equal_values(x.begin(), x.end(), y.begin(), y.end());
}
template <class JX, class JY, class FE>
inline bool equalValues(const JX& x, const JY& y, FE fe) {
  return equal_values(x.begin(), x.end(), y.begin(), y.end(), fe);
}




// COUNT
// -----
// Count businesses in a sector.

template <class J, class F>
inline size_t countIf(const J& x, F fn) {
  return count_if(x.begin(), x.end(), fn);
}




// COPY
// ----

template <class I, class IA>
inline auto copy_values(I ib, I ie, IA ab) {
  return copy(ib, ie, ab);
}
template <class J, class JA>
inline size_t copyValues(const J& x, JA& a) {
  auto   it = copy_values(x.begin(), x.end(), a.begin());
  return it - a.begin();
}


template <class I, class T>
inline auto copy_append(I ib, I ie, vector<T>& a) {
  return a.insert(a.end(), ib, ie);
}
template <class J, class T>
inline size_t copyAppend(const J& x, vector<T>& a) {
  auto   it = copy_append(x.begin(), x.end(), a);
  return it - a.begin();
}


template <class I>
inline auto copy_vector(I ib, I ie) {
  using T = typename iterator_traits<I>::value_type; vector<T> a;
  copy_append(ib, ie, a);
  return a;
}
template <class J>
inline auto copyVector(const J& x) {
  return copy_vector(x.begin(), x.end());
}




// VALUE INDEX
// -----------
// Keep the address of each business (yellow pages).

template <class I, class M>
auto value_index(I ib, I ie, M& a) {
  size_t i = 0;
  for (; ib != ie; ++ib)
    a[*ib] = i++;
  return a;
}
template <class J, class M>
inline auto valueIndex(const J& x, M& a) {
  return value_index(x.begin(), x.end(), a);
}

template <class I>
inline auto value_index_unordered_map(I ib, I ie) {
  using K = typename iterator_traits<I>::value_type;
  unordered_map<K, size_t> a;
  return value_index(ib, ie, a);
}
template <class J>
inline auto valueIndexUnorderedMap(const J& x) {
  return value_index_unordered_map(x.begin(), x.end());
}




// SORT
// ----
// Arrange your portfolio by ROCE.

template <class I>
inline void sort_values(I ib, I ie) {
  sort(ib, ie);
}
template <class I, class FL>
inline void sort_values(I ib, I ie, FL fl) {
  sort(ib, ie, fl);
}
template <class J>
inline void sortValues(J& x) {
  sort_values(x.begin(), x.end());
}
template <class J, class FL>
inline void sortValues(J& x, FL fl) {
  sort_values(x.begin(), x.end(), fl);
}




// UNIQUE
// ------

template <class IX, class IA, class FE>
auto unique_last_copy(IX xb, IX xe, IA ab, FE fe) {
  if (xb==xe) return ab;
  // Compare adjacent elements, and
  // only copy non-matching ones.
  IX it = xb++;
  for (; xb!=xe; ++it, ++xb)
  { if (!fe(*it, *xb)) *(ab++) = *it; }
  // Copy the last unique element.
  *(ab++) = *it;
  return ab;
}
template <class IX, class IA>
inline auto unique_last_copy(IX xb, IX xe, IA ab) {
  auto fe = [](const auto& a, const auto& b) { return a == b; };
  return unique_last_copy(xb, xe, ab, fe);
}




// SET DIFFERENCE
// --------------

// Remove from `x`, elements given in `y`.
// Both `x` and `y` must be sorted.
template <class IX, class IY, class FL, class FE>
auto set_difference_inplace(IX xb, IX xe, IY yb, IY ye, FL fl, FE fe) {
  if (xb==xe || yb==ye) return xe;
  // Write-free loop when there is
  // nothing to remove.
  while (true) {
    while (fl(*xb, *yb))
    { if (++xb==xe) return xe; }
    if (fe(*xb, *(yb++))) break;
    if (yb==ye) return xe;
  }
  // There was a match, remove it.
  IX it = xb++;
  // Only one element needs removal.
  if (xb==xe) return it;
  // With-write loop when there are
  // more elements to remove.
   while (yb!=ye) {
    while (fl(*xb, *yb)) {
      *(it++) = *xb;
      if (++xb==xe) return it;
    }
    if (fe(*xb, *(yb++)))
    { if (++xb==xe) return it; }
  }
  // No more elements to remove.
  // Shift the remaining elements.
  return copy(xb, xe, it);
}
template <class IX, class IY>
inline auto set_difference_inplace(IX xb, IX xe, IY yb, IY ye) {
  auto fl = [](const auto& a, const auto& b) { return a <  b; };
  auto fe = [](const auto& a, const auto& b) { return a == b; };
  return set_difference_inplace(xb, xe, yb, ye, fl, fe);
}

template <class JX, class JY, class FL, class FE>
inline size_t setDifferenceInplace(JX& x, const JY& y, FL fl, FE fe) {
  auto   it = set_difference_inplace(x.begin(), x.end(), y.begin(), y.end(), fl, fe);
  return it - x.begin();
}
template <class JX, class JY>
inline size_t setDifferenceInplace(JX& x, const JY& y) {
  auto   it = set_difference_inplace(x.begin(), x.end(), y.begin(), y.end());
  return it - x.begin();
}




// SET UNION
// ---------

// Add elements from `y` into `x`, preferring the last in `y` among matching elements.
// Both `x` and `y` must be sorted. There must be sufficient space in `x` and `b` (buffer = |y|+2+1).
template <class IX, class IY, class IB, class FL, class FE>
auto set_union_last_inplace(IX xb, IX xe, IY yb, IY ye, IB bb, IB be, FL fl, FE fe) {
  if (yb==ye) return xe;
  if (xb==xe) return unique_last_copy(yb, ye, xb);
  // Deque-free loop when there
  // is nothing to insert.
  while (true) {
    while (fl(*xb, *yb))
      if (++xb==xe) return unique_last_copy(yb, ye, xb);
    if (!fe(*xb, *yb)) break;
    *xb = *yb;
    if (++yb==ye) return xe;
  }
  // Insert smaller element from `y`, after
  // saving the element in `x` to deque.
  auto q = unsized_deque_view(bb, be);
  IX it = xb;
  q.push_back(*(xb++));
  *it = *(yb++);
  // With-deque loop when elements
  // in `y` can be inserted into `x`.
  while (yb!=ye) {
    if (fe(*it, *yb)) *it = *(yb++);
    else {
      if (xb!=xe) q.push_back(*(xb++));
      *(++it) = !q.empty() && fl(q.front(), *yb)? q.pop_front() : *(yb++);
    }
  }
  // Continue until both `x` and
  // the deque are empty.
  while (true) {
    if (xb!=xe) q.push_back(*(xb++));
    if (q.empty()) break;
    *(++it) = q.pop_front();
  }
  return ++it;
}
template <class IX, class IY, class IB>
inline auto set_union_last_inplace(IX xb, IX xe, IY yb, IY ye, IB bb, IB be) {
  auto fl = [](const auto& a, const auto& b) { return a <  b; };
  auto fe = [](const auto& a, const auto& b) { return a == b; };
  return set_union_last_inplace(xb, xe, yb, ye, bb, be, fl, fe);
}

template <class JX, class JY, class JB, class FL, class FE>
inline size_t setUnionLastInplace(JX& x, const JY& y, JB& b, FL fl, FE fe) {
  auto   it = set_union_last_inplace(x.begin(), x.end(), y.begin(), y.end(), b.begin(), b.end(), fl, fe);
  return it - x.begin();
}
template <class JX, class JY, class JB>
inline size_t setUnionLastInplace(JX& x, const JY& y, JB& b) {
  auto   it = set_union_last_inplace(x.begin(), x.end(), y.begin(), y.end(), b.begin(), b.end());
  return it - x.begin();
}
