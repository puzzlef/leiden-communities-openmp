#pragma once
// Avoid compiler error/bug "error: structured binding refers to incomplete type"
#include <unordered_map>
#include <algorithm>
#include "_queue.hxx"

using std::copy;




#pragma region METHODS
/**
 * Find the first element that does not match its adjacent element.
 * @param ib begin iterator
 * @param ie end iterator
 * @param fe equality function (a, b)
 * @returns iterator to the first non-adjacent element
 */
template <class I, class FE>
inline I non_adjacent_find(I ib, I ie, FE fe) {
  if (ib==ie) return ie;
  // Compare adjacent elements, and if they
  // dont match, return the first one.
  I it = ib++;
  for (; ib!=ie; ++it, ++ib)
    if (!fe(*it, *ib)) return it;
  // The last element must be non-adjacent.
  return it;
}


/**
 * Find the first element that does not match its adjacent element.
 * @param ib begin iterator
 * @param ie end iterator
 * @returns iterator to the first non-adjacent element
 */
template <class I>
inline auto non_adjacent_find(I ib, I ie) {
  auto fe = [](const auto& a, const auto& b) { return a == b; };
  return non_adjacent_find(ib, ie, fe);
}




/**
 * Obtain the index of each element in a container.
 * @param ib begin iterator
 * @param ie end iterator
 * @param a map to store the index of each element (updated)
 */
template <class I, class M>
inline void value_index(I ib, I ie, M& a) {
  size_t i = 0;
  for (; ib != ie; ++ib)
    a[*ib] = i++;
  return a;
}




/**
 * Keep only the last unique element in a container.
 * @param xb begin iterator of input (updated)
 * @param xe end iterator of input (updated)
 * @param ab begin iterator of output
 * @param fe equality function (a, b)
 * @returns iterator to the end of output
 */
template <class IX, class IA, class FE>
inline IA unique_last_copy(IX xb, IX xe, IA ab, FE fe) {
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


/**
 * Keep only the last unique element in a container.
 * @param xb begin iterator of input (updated)
 * @param xe end iterator of input (updated)
 * @param ab begin iterator of output
 * @returns iterator to the end of output
 */
template <class IX, class IA>
inline IA unique_last_copy(IX xb, IX xe, IA ab) {
  auto fe = [](const auto& a, const auto& b) { return a == b; };
  return unique_last_copy(xb, xe, ab, fe);
}




/**
 * Keep elements in a container that are not in another container (in-place).
 * @param xb begin iterator of input (updated)
 * @param xe end iterator of input (updated)
 * @param yb begin iterator of elements to remove
 * @param ye end iterator of elements to remove
 * @param fl less-than function (a, b)
 * @param fe equality function (a, b)
 * @returns iterator to the end of updated input
 */
template <class IX, class IY, class FL, class FE>
inline IX set_difference_inplace(IX xb, IX xe, IY yb, IY ye, FL fl, FE fe) {
  // Remove from `x`, elements given in `y`.
  // Both `x` and `y` must be sorted.
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


/**
 * Keep elements in a container that are not in another container (in-place).
 * @param xb begin iterator of input (updated)
 * @param xe end iterator of input (updated)
 * @param yb begin iterator of elements to remove
 * @param ye end iterator of elements to remove
 * @returns iterator to the end of updated input
 */
template <class IX, class IY>
inline IX set_difference_inplace(IX xb, IX xe, IY yb, IY ye) {
  auto fl = [](const auto& a, const auto& b) { return a <  b; };
  auto fe = [](const auto& a, const auto& b) { return a == b; };
  return set_difference_inplace(xb, xe, yb, ye, fl, fe);
}




/**
 * Find the set union of two containers, keeping the last element among matching ones (in-place).
 * @param xb begin iterator of input (updated)
 * @param xe end iterator of input (updated)
 * @param yb begin iterator of elements to add
 * @param ye end iterator of elements to add
 * @param bb begin iterator of buffer, of size |y|+2+1
 * @param be end iterator of buffer
 * @param fl less-than function (a, b)
 * @param fe equality function (a, b)
 * @returns iterator to the end of updated input
 */
template <class IX, class IY, class IB, class FL, class FE>
inline IX set_union_last_inplace(IX xb, IX xe, IY yb, IY ye, IB bb, IB be, FL fl, FE fe) {
  // Add elements from `y` into `x`, preferring the last in `y` among matching elements.
  // Both `x` and `y` must be sorted. There must be sufficient space in `x` and `b` (buffer = |y|+2+1).
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


/**
 * Find the set union of two containers, keeping the last element among matching ones (in-place).
 * @param xb begin iterator of input (updated)
 * @param xe end iterator of input (updated)
 * @param yb begin iterator of elements to add
 * @param ye end iterator of elements to add
 * @param bb begin iterator of buffer
 * @param be end iterator of buffer
 * @returns iterator to the end of updated input
 */
template <class IX, class IY, class IB>
inline auto set_union_last_inplace(IX xb, IX xe, IY yb, IY ye, IB bb, IB be) {
  auto fl = [](const auto& a, const auto& b) { return a <  b; };
  auto fe = [](const auto& a, const auto& b) { return a == b; };
  return set_union_last_inplace(xb, xe, yb, ye, bb, be, fl, fe);
}
#pragma endregion
