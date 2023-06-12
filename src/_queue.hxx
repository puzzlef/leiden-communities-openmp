#pragma once
#include <iterator>

using std::iterator_traits;




// DEQUE VIEW
// ----------
// Bounded (circular) deque view capable of storing N elements.

template <class I>
class DequeView {
  // Data.
  protected:
  using T = typename iterator_traits<I>::value_type;
  const I xb, xe;
  I ib, ie;
  size_t n;

  // Basic operations.
  public:
  using value_type = T;
  inline size_t size()  const { return n; }
  inline bool   empty() const { return n==0; }

  // Read operations.
  public:
  inline auto back()  const { return ie==xb? *(xe-1) : *(ie-1); }
  inline auto front() const { return *ib; }

  // Write operations.
  public:
  inline void push_back(const T& v) {
    if (ib!=ie || n==0) ++n;
    *(ie++) = v;
    if (ie==xe) ie = xb;
  }

  inline void push_front(const T& v) {
    if (ib==xb) ib = xe;
    *(--ib) = v;
    if (ib!=ie || n==0) ++n;
  }

  inline auto pop_back() {
    if (ie==xb) ie = xe;
    if (n>0) --n;
    return *(--ie);
  }

  inline auto pop_front() {
    if (n>0) --n;
    auto v = *(ib++);
    if (ib==xe) ib = xb;
    return v;
  }

  // Lifetime operations.
  DequeView(I xb, I xe) :
  xb(xb), xe(xe), ib(xb), ie(xb), n(0) {}
};

template <class I>
inline auto deque_view(I xb, I xe) {
  return DequeView<I>(xb, xe);
}
template <class J>
inline auto dequeView(J& x) {
  return deque_view(x.begin(), x.end());
}




// UNSIZED DEQUE VIEW
// ------------------
// Bounded (circular) deque view capable of storing N-1 elements.

template <class I>
class UnsizedDequeView {
  // Data.
  protected:
  using T = typename iterator_traits<I>::value_type;
  const I xb, xe;
  I ib, ie;

  // Basic operations
  public:
  using value_type = T;
  inline bool empty() const { return ib==ie; }

  // Read operations.
  public:
  inline auto back()  const { return ie==xb? *(xe-1) : *(ie-1); }
  inline auto front() const { return *ib; }

  // Write operations.
  public:
  inline void push_back(const T& v) {
    *(ie++) = v;
    if (ie==xe) ie = xb;
  }

  inline void push_front(const T& v) {
    if (ib==xb) ib = xe;
    *(--ib) = v;
  }

  inline auto pop_back() {
    if (ie==xb) ie = xe;
    return *(--ie);
  }

  inline auto pop_front() {
    auto v = *(ib++);
    if (ib==xe) ib = xb;
    return v;
  }

  // Lifetime operations.
  UnsizedDequeView(I xb, I xe) :
  xb(xb), xe(xe), ib(xb), ie(xb) {}
};

template <class I>
inline auto unsized_deque_view(I xb, I xe) {
  return UnsizedDequeView<I>(xb, xe);
}
template <class J>
inline auto unsizedDequeView(J& x) {
  return unsized_deque_view(x.begin(), x.end());
}
