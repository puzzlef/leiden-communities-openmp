#pragma once
#include <cstddef>
#include <utility>
#include <iterator>
#include <vector>
#include <algorithm>
#include <type_traits>
#include "_utility.hxx"

using std::pair;
using std::ptrdiff_t;
using std::remove_reference_t;
using std::input_iterator_tag;
using std::output_iterator_tag;
using std::forward_iterator_tag;
using std::bidirectional_iterator_tag;
using std::random_access_iterator_tag;
using std::iterator_traits;
using std::vector;
using std::make_pair;
using std::distance;
using std::max;




// ITERATOR HELPER MACROS
// ----------------------
// Helps create iterators.

#ifndef PAREN_OPEN
#define PAREN_OPEN  (
#define PAREN_CLOSE )
#endif




// https://stackoverflow.com/questions/8054273/how-to-implement-an-stl-style-iterator-and-avoid-common-pitfalls
// class iterator {
//     iterator(const iterator&);
//     ~iterator();
//     iterator& operator=(const iterator&);
//     iterator& operator++(); //prefix increment
//     reference operator*() const;
//     friend void swap(iterator& lhs, iterator& rhs); //C++11 I think
// };
//
// class input_iterator : public virtual iterator {
//     iterator operator++(int); //postfix increment
//     value_type operator*() const;
//     pointer operator->() const;
//     friend bool operator==(const iterator&, const iterator&);
//     friend bool operator!=(const iterator&, const iterator&);
// };
// //once an input iterator has been dereferenced, it is
// //undefined to dereference one before that.
//
// class output_iterator : public virtual iterator {
//     reference operator*() const;
//     iterator operator++(int); //postfix increment
// };
// //dereferences may only be on the left side of an assignment
// //once an output iterator has been dereferenced, it is
// //undefined to dereference one before that.
//
// class forward_iterator : input_iterator, output_iterator {
//     forward_iterator();
// };
// //multiple passes allowed
//
// class bidirectional_iterator : forward_iterator {
//     iterator& operator--(); //prefix decrement
//     iterator operator--(int); //postfix decrement
// };
//
// class random_access_iterator : bidirectional_iterator {
//     friend bool operator<(const iterator&, const iterator&);
//     friend bool operator>(const iterator&, const iterator&);
//     friend bool operator<=(const iterator&, const iterator&);
//     friend bool operator>=(const iterator&, const iterator&);
//
//     iterator& operator+=(size_type);
//     friend iterator operator+(const iterator&, size_type);
//     friend iterator operator+(size_type, const iterator&);
//     iterator& operator-=(size_type);
//     friend iterator operator-(const iterator&, size_type);
//     friend difference_type operator-(iterator, iterator);
//
//     reference operator[](size_type) const;
// };




#ifndef ITERATOR_USING
#define ITERATOR_USING(cat, dif, val, ref, ptr) \
  using iterator_category = cat; \
  using difference_type   = dif; \
  using value_type        = val; \
  using reference         = ref; \
  using pointer           = ptr;

#define ITERATOR_USING_X(X) \
  using iterator_category = X::iterator_category; \
  using difference_type   = typename X::difference_type; \
  using value_type        = typename X::value_type; \
  using reference         = typename X::reference; \
  using pointer           = typename X::pointer;

#define ITERATOR_USING_XC(X, cat) \
  using iterator_category = cat; \
  using difference_type   = typename X::difference_type; \
  using value_type        = typename X::value_type; \
  using reference         = typename X::reference; \
  using pointer           = typename X::pointer;

#define ITERATOR_USING_XCD(X, cat, dif) \
  using iterator_category = cat; \
  using difference_type   = dif; \
  using value_type        = typename X::value_type; \
  using reference         = typename X::reference; \
  using pointer           = typename X::pointer;

#define ITERATOR_USING_XCVRP(X, cat, val, ref, ptr) \
  using iterator_category = cat; \
  using difference_type   = typename X::difference_type; \
  using value_type        = val; \
  using reference         = ref; \
  using pointer           = ptr;

#define ITERATOR_USING_XCRP(X, cat, ref, ptr) \
  using iterator_category = cat; \
  using difference_type   = typename X::difference_type; \
  using value_type        = remove_reference_t<ref>; \
  using reference         = ref; \
  using pointer           = ptr;
#endif


#ifndef ITERATOR_COPY
#define ITERATOR_COPY(f0, f1, cv, ce) \
  f0 iterator(const iterator& cv) f1 { ce; }
#endif

#ifndef ITERATOR_ASSIGN
#define ITERATOR_ASSIGN(f0, f1, av, ae) \
  f0 iterator& operator=(const iterator& av) f1 { ae; }
#endif

#ifndef ITERATOR_DESTROY
#define ITERATOR_DESTROY(f0, f1, de) \
  f0 ~iterator() f1 { de; }
#endif


#ifndef ITERATOR_DEREF
#define ITERATOR_DEREF_VALUE_DO(f0, f1, de) \
  f0 value_type operator*() f1 { de; }
#define ITERATOR_DEREF_VALUE(f0, f1, de) \
  ITERATOR_DEREF_VALUE_DO(f0, f1, return de)
#define ITERATOR_DEREF_DO(f0, f1, de) \
  f0 reference operator*() f1 { de; }
#define ITERATOR_DEREF(f0, f1, de) \
  ITERATOR_DEREF_DO(f0, f1, return de)
#endif

#ifndef ITERATOR_LOOKUP
#define ITERATOR_LOOKUP_DO(f0, f1, nv, le) \
  f0 reference operator[](difference_type nv) f1 { le; }
#define ITERATOR_LOOKUP(f0, f1, nv, le) \
  ITERATOR_LOOKUP_DO(f0, f1, nv, return le)
#endif

#ifndef ITERATOR_POINTER
#define ITERATOR_POINTER_DO(f0, f1, pe) \
  f0 pointer operator->() f1 { pe; }
#define ITERATOR_POINTER(f0, f1, pe) \
  ITERATOR_POINTER_DO(f0, f1, return pe)
#endif


#ifndef ITERATOR_INCREMENT
#define ITERATOR_INCREMENT_PRE(f0, f1, ie) \
  f0 iterator& operator++() f1 { ie; return *this; }
#define ITERATOR_INCREMENT_POST(f0, f1) \
  f0 iterator operator++(int) f1 { auto it = *this; ++(*this); return it; }
#define ITERATOR_INCREMENT(f0, f1, ie) \
  ITERATOR_INCREMENT_PRE(f0, f1, ie) \
  ITERATOR_INCREMENT_POST(f0, f1)
#endif

#ifndef ITERATOR_DECREMENT
#define ITERATOR_DECREMENT_PRE(f0, f1, de) \
  f0 iterator& operator--() f1 { de; return *this; }
#define ITERATOR_DECREMENT_POST(f0, f1)  \
  f0 iterator operator--(int) f1 { auto it = *this; --(*this); return it; }
#define ITERATOR_DECREMENT(f0, f1, de) \
  ITERATOR_DECREMENT_PRE(f0, f1, de) \
  ITERATOR_DECREMENT_POST(f0, f1)
#endif

#ifndef ITERATOR_NEXT
#define ITERATOR_NEXT_PRE(f0, f1, ie, de) \
  ITERATOR_INCREMENT_PRE(f0, f1, ie) \
  ITERATOR_DECREMENT_PRE(f0, f1, de)
#define ITERATOR_NEXT_POST(f0, f1) \
  ITERATOR_INCREMENT_POST(f0, f1) \
  ITERATOR_DECREMENT_POST(f0, f1)
#define ITERATOR_NEXT(f0, f1, ie, de) \
  ITERATOR_NEXT_PRE(f0, f1, ie, de) \
  ITERATOR_NEXT_POST(f0, f1)
#endif


#ifndef ITERATOR_ADVANCE
#define ITERATOR_ADVANCE_PLUS(f0, f1, nv, pe) \
  f0 iterator& operator+=(difference_type nv) f1 { pe; return *this; }
#define ITERATOR_ADVANCE_MINUS(f0, f1, nv, me) \
  f0 iterator& operator-=(difference_type nv) f1 { me; return *this; }
#define ITERATOR_ADVANCE(f0, f1, nv, pe, me) \
  ITERATOR_ADVANCE_PLUS(f0, f1, nv, pe) \
  ITERATOR_ADVANCE_MINUS(f0, f1, nv, me)
#endif


#ifndef ITERATOR_ARITHMETIC
#define ITERATOR_ARITHMETIC_PLUS_DO(f0, f1, iv, nv, pe) \
  f0 friend iterator operator+(const iterator& iv, difference_type nv) f1 { pe; } \
  f0 friend iterator operator+(difference_type nv, const iterator& iv) f1 { return iv+nv; }
#define ITERATOR_ARITHMETIC_PLUS(f0, f1, iv, nv, pe) \
  ITERATOR_ARITHMETIC_PLUS_DO(f0, f1, iv, nv, return pe)

#define ITERATOR_ARITHMETIC_MINUS_DO(f0, f1, iv, nv, me) \
  f0 friend iterator operator-(const iterator& iv, difference_type nv) f1 { me; }
#define ITERATOR_ARITHMETIC_MINUS(f0, f1, iv, nv, me) \
  ITERATOR_ARITHMETIC_MINUS_DO(f0, f1, iv, nv, return me)

#define ITERATOR_ARITHMETIC_DO(f0, f1, iv, nv, pe, me)  \
  ITERATOR_ARITHMETIC_PLUS_DO(f0, f1, iv, nv, pe) \
  ITERATOR_ARITHMETIC_MINUS_DO(f0, f1, iv, nv, me)
#define ITERATOR_ARITHMETIC(f0, f1, iv, nv, pe, me)  \
  ITERATOR_ARITHMETIC_PLUS(f0, f1, iv, nv, pe) \
  ITERATOR_ARITHMETIC_MINUS(f0, f1, iv, nv, me)
#endif


#ifndef ITERATOR_DIFFERENCE
#define ITERATOR_DIFFERENCE_DO(f0, f1, lv, rv, de) \
  f0 friend difference_type operator-(iterator lv, iterator rv) f1 { de; }
#define ITERATOR_DIFFERENCE(f0, f1, lv, rv, de) \
  ITERATOR_DIFFERENCE_DO(f0, f1, lv, rv, return de)
#endif


#ifndef ITERATOR_SWAP
#define ITERATOR_SWAP(f0, f1, lv, rv, se) \
  f0 friend void swap(iterator& lv, iterator& rv) f1 { se; }
#endif


#ifndef ITERATOR_COMPARE
#define ITERATOR_COMPARE_EQ_DO(f0, f1, lv, rv, ce) \
  f0 friend bool operator==(const iterator& lv, const iterator& rv) f1 { ce; }
#define ITERATOR_COMPARE_NE_DO(f0, f1, lv, rv, ce) \
  f0 friend bool operator!=(const iterator& lv, const iterator& rv) f1 { ce; }
#define ITERATOR_COMPARE_GT_DO(f0, f1, lv, rv, ce) \
  f0 friend bool operator>(const iterator& lv, const iterator& rv) f1 { ce; }
#define ITERATOR_COMPARE_LT_DO(f0, f1, lv, rv, ce) \
  f0 friend bool operator<(const iterator& lv, const iterator& rv) f1 { ce; }
#define ITERATOR_COMPARE_GE_DO(f0, f1, lv, rv, ce) \
  f0 friend bool operator>=(const iterator& lv, const iterator& rv) f1 { ce; }
#define ITERATOR_COMPARE_LE_DO(f0, f1, lv, rv, ce) \
  f0 friend bool operator<=(const iterator& lv, const iterator& rv) f1 { ce; }

#define ITERATOR_COMPARE_EQ(f0, f1, lv, rv, ce) \
  ITERATOR_COMPARE_EQ_DO(f0, f1, lv, rv, return ce)
#define ITERATOR_COMPARE_NE(f0, f1, lv, rv, ce) \
  ITERATOR_COMPARE_NE_DO(f0, f1, lv, rv, return ce)
#define ITERATOR_COMPARE_GT(f0, f1, lv, rv, ce) \
  ITERATOR_COMPARE_GT_DO(f0, f1, lv, rv, return ce)
#define ITERATOR_COMPARE_LT(f0, f1, lv, rv, ce) \
  ITERATOR_COMPARE_LT_DO(f0, f1, lv, rv, return ce)
#define ITERATOR_COMPARE_GE(f0, f1, lv, rv, ce) \
  ITERATOR_COMPARE_GE_DO(f0, f1, lv, rv, return ce)
#define ITERATOR_COMPARE_LE(f0, f1, lv, rv, ce) \
  ITERATOR_COMPARE_LE_DO(f0, f1, lv, rv, return ce)

#define ITERATOR_COMPARE_EQNE_DO(f0, f1, lv, rv, eq, ne) \
  ITERATOR_COMPARE_EQ_DO(f0, f1, lv, rv, eq) \
  ITERATOR_COMPARE_NE_DO(f0, f1, lv, rv, ne)
#define ITERATOR_COMPARE_EQNE(f0, f1, lv, rv, eq, ne) \
  ITERATOR_COMPARE_EQ(f0, f1, lv, rv, eq) \
  ITERATOR_COMPARE_NE(f0, f1, lv, rv, ne)
#define ITERATOR_COMPARE_EQNE_INT(f0, f1, lv, rv, ce) \
  ITERATOR_COMPARE_EQ(f0, f1, lv, rv, (ce)==0) \
  ITERATOR_COMPARE_NE(f0, f1, lv, rv, (ce)!=0)

#define ITERATOR_COMPARE_GTLT_DO(f0, f1, lv, rv, gt, lt) \
  ITERATOR_COMPARE_GT_DO(f0, f1, lv, rv, gt) \
  ITERATOR_COMPARE_LT_DO(f0, f1, lv, rv, lt)
#define ITERATOR_COMPARE_GTLT(f0, f1, lv, rv, gt, lt) \
  ITERATOR_COMPARE_GT(f0, f1, lv, rv, gt) \
  ITERATOR_COMPARE_LT(f0, f1, lv, rv, lt)
#define ITERATOR_COMPARE_GTLT_INT(f0, f1, lv, rv, ce) \
  ITERATOR_COMPARE_GT(f0, f1, lv, rv, (ce)>0) \
  ITERATOR_COMPARE_LT(f0, f1, lv, rv, (ce)<0)

#define ITERATOR_COMPARE_GELE_DO(f0, f1, lv, rv, ge, le) \
  ITERATOR_COMPARE_GE_DO(f0, f1, lv, rv, ge) \
  ITERATOR_COMPARE_LE_DO(f0, f1, lv, rv, le)
#define ITERATOR_COMPARE_GELE(f0, f1, lv, rv, ge, le) \
  ITERATOR_COMPARE_GE(f0, f1, lv, rv, ge) \
  ITERATOR_COMPARE_LE(f0, f1, lv, rv, le)
#define ITERATOR_COMPARE_GELE_INT(f0, f1, lv, rv, ce) \
  ITERATOR_COMPARE_GE(f0, f1, lv, rv, (ce)>0) \
  ITERATOR_COMPARE_LE(f0, f1, lv, rv, (ce)<0)

#define ITERATOR_COMPARE_DO(f0, f1, lv, rv, eq, ne, gt, lt, ge, le) \
  ITERATOR_COMPARE_EQ_DO(f0, f1, lv, rv, eq) \
  ITERATOR_COMPARE_NE_DO(f0, f1, lv, rv, ne) \
  ITERATOR_COMPARE_GT_DO(f0, f1, lv, rv, gt) \
  ITERATOR_COMPARE_LT_DO(f0, f1, lv, rv, lt) \
  ITERATOR_COMPARE_GE_DO(f0, f1, lv, rv, ge) \
  ITERATOR_COMPARE_LE_DO(f0, f1, lv, rv, le)
#define ITERATOR_COMPARE(f0, f1, lv, rv, eq, ne, gt, lt, ge, le) \
  ITERATOR_COMPARE_EQ(f0, f1, lv, rv, eq) \
  ITERATOR_COMPARE_NE(f0, f1, lv, rv, ne) \
  ITERATOR_COMPARE_GT(f0, f1, lv, rv, gt) \
  ITERATOR_COMPARE_LT(f0, f1, lv, rv, lt) \
  ITERATOR_COMPARE_GE(f0, f1, lv, rv, ge) \
  ITERATOR_COMPARE_LE(f0, f1, lv, rv, le)
#define ITERATOR_COMPARE_INT(f0, f1, lv, rv, ce) \
  ITERATOR_COMPARE_EQ(f0, f1, lv, rv, (ce)==0) \
  ITERATOR_COMPARE_NE(f0, f1, lv, rv, (ce)!=0) \
  ITERATOR_COMPARE_GT(f0, f1, lv, rv, (ce)>0) \
  ITERATOR_COMPARE_LT(f0, f1, lv, rv, (ce)<0) \
  ITERATOR_COMPARE_GE(f0, f1, lv, rv, (ce)>=0) \
  ITERATOR_COMPARE_LE(f0, f1, lv, rv, (ce)<=0)
#endif




// ITERABLE HELPER MACROS
// ----------------------
// Helps create iterables.

#ifndef ITERABLE_ITERATOR
#define ITERABLE_NAME_DO(f0, f1, fn, fe) \
  f0 auto fn() f1 { fe; }
#define ITERABLE_BEGIN_DO(f0, f1, be) \
  f0 auto begin() f1 { be; }
#define ITERABLE_END_DO(f0, f1, ee) \
  f0 auto end()   f1 { ee; }

#define ITERABLE_CNAME_DO(f0, f1, fn, fe) \
  f0 auto fn()    const f1 { fe; } \
  f0 auto c##fn() const f1 { fe; }
#define ITERABLE_CBEGIN_DO(f0, f1, be) \
  f0 auto begin()  const f1 { be; } \
  f0 auto cbegin() const f1 { be; }
#define ITERABLE_CEND_DO(f0, f1, ee) \
  f0 auto end()  const f1 { ee; } \
  f0 auto cend() const f1 { ee; }

#define ITERABLE_NAMES_DO(f0, f1, fn, fe) \
  ITERABLE_NAME_DO(f0, f1, fn, fe) \
  ITERABLE_CNAME_DO(f0, f1, fn, fe)
#define ITERABLE_BEGINS_DO(f0, f1, be) \
  ITERABLE_BEGIN_DO(f0, f1, be) \
  ITERABLE_CBEGIN_DO(f0, f1, be)
#define ITERABLE_ENDS_DO(f0, f1, ee) \
  ITERABLE_END_DO(f0, f1, ee) \
  ITERABLE_CEND_DO(f0, f1, ee)

#define ITERABLE_NAME(f0, f1, fn, fe) \
  ITERABLE_NAME_DO(f0, f1, fn, return fe)
#define ITERABLE_BEGIN(f0, f1, be) \
  ITERABLE_BEGIN_DO(f0, f1, return be)
#define ITERABLE_END(f0, f1, ee) \
  ITERABLE_END_DO(f0, f1, return ee)

#define ITERABLE_CNAME(f0, f1, fn, fe) \
  ITERABLE_CNAME_DO(f0, f1, fn, return fe)
#define ITERABLE_CBEGIN(f0, f1, be) \
  ITERABLE_CBEGIN_DO(f0, f1, return be)
#define ITERABLE_CEND(f0, f1, ee) \
  ITERABLE_CEND_DO(f0, f1, return ee)

#define ITERABLE_NAMES(f0, f1, fn, fe) \
  ITERABLE_NAMES_DO(f0, f1, fn, return fe)
#define ITERABLE_BEGINS(f0, f1, be) \
  ITERABLE_BEGINS_DO(f0, f1, return be)
#define ITERABLE_ENDS(f0, f1, ee) \
  ITERABLE_ENDS_DO(f0, f1, return ee)

#define ITERABLE_ITERATOR_DO(f0, f1, be, ee) \
  ITERABLE_BEGIN_DO(f0, f1, be) \
  ITERABLE_END_DO(f0, f1, ee)
#define ITERABLE_ITERATOR(f0, f1, be, ee) \
  ITERABLE_BEGIN(f0, f1, be) \
  ITERABLE_END(f0, f1, ee)

#define ITERABLE_CITERATOR_DO(f0, f1, be, ee) \
  ITERABLE_CBEGIN_DO(f0, f1, be) \
  ITERABLE_CEND_DO(f0, f1, ee)
#define ITERABLE_CITERATOR(f0, f1, be, ee) \
  ITERABLE_CBEGIN(f0, f1, be) \
  ITERABLE_CEND(f0, f1, ee)

#define ITERABLE_ITERATORS_DO(f0, f1, be, ee) \
  ITERABLE_BEGINS_DO(f0, f1, be) \
  ITERABLE_ENDS_DO(f0, f1, ee)
#define ITERABLE_ITERATORS(f0, f1, be, ee) \
  ITERABLE_BEGINS(f0, f1, be) \
  ITERABLE_ENDS(f0, f1, ee)

#define ITERABLE_ITERATOR_DEFAULT(ib, ie) \
  ITERABLE_ITERATOR(inline, const, ib, ie)
#define ITERABLE_CITERATOR_DEFAULT(ib, ie) \
  ITERABLE_CITERATOR(inline, const, ib, ie)
#define ITERABLE_ITERATORS_DEFAULT(ib, ie) \
  ITERABLE_ITERATORS(inline, const, ib, ie)
#endif


#ifndef ITERABLE_SIZES
#define ITERABLE_SIZE_DO(f0, f1, se) \
  f0 size_t size() f1 { se; }
#define ITERABLE_SIZE(f0, f1, se) \
  ITERABLE_SIZE_DO(f0, f1, return se)

#define ITERABLE_EMPTY_DO(f0, f1, ee) \
  f0 bool empty() f1 { ee; }
#define ITERABLE_EMPTY(f0, f1, ee) \
  ITERABLE_EMPTY_DO(f0, f1, return ee)

#define ITERABLE_SIZES_DO(f0, f1, se, ee) \
  ITERABLE_SIZE_DO(f0, f1, se) \
  ITERABLE_EMPTY_DO(f0, f1, ee)
#define ITERABLE_SIZES(f0, f1, se, ee) \
  ITERABLE_SIZE(f0, f1, se) \
  ITERABLE_EMPTY(f0, f1, ee)
#define ITERABLE_SIZES_DEFAULT(ib, ie) \
  ITERABLE_SIZES(inline, const, distance(ib, ie), ib == ie)
#endif




// ITERABLE
// --------
// Holds begin and end iterators.

template <class I>
class Iterable {
  const I ib, ie;
  public:
  using iterator = I;
  Iterable(I ib, I ie) noexcept : ib(ib), ie(ie) {}
  ITERABLE_ITERATOR_DEFAULT(ib, ie)
  ITERABLE_SIZES_DEFAULT(ib, ie)
};

template <class I>
inline auto iterable(I ib, I ie) noexcept {
  return Iterable<I>(ib, ie);
}
template <class J>
inline auto iterable(const J& x) {
  return iterable(x.begin(), x.end());
}




// SIZED ITERABLE
// --------------

template <class I>
class SizedIterable {
  const I ib, ie;
  const size_t N;
  public:
  SizedIterable(I ib, I ie, size_t N) noexcept : ib(ib), ie(ie), N(N) {}
  SizedIterable(I ib, I ie) : ib(ib), ie(ie), N(distance(ib, ie)) {}
  public:
  ITERABLE_ITERATOR_DEFAULT(ib, ie)
  ITERABLE_SIZES(inline, const, N, N == 0)
};

template <class I>
inline auto sized_iterable(I ib, I ie, size_t N) noexcept {
  return SizedIterable<I>(ib, ie, N);
}
template <class I>
inline auto sized_iterable(I ib, I ie) {
  return SizedIterable<I>(ib, ie);
}
template <class J>
inline auto sizedIterable(const J& x, size_t N) {
  return sized_iterable(x.begin(), x.end(), N);
}
template <class J>
inline auto sizedIterable(const J& x) {
  return sized_iterable(x.begin(), x.end());
}




// PAIRED ITERABLE
// ---------------

#define PAIRED_ITERABLE_SHORT_TYPES(I) \
  using K = decltype((*ib).first); \
  using V = decltype((*ib).second);


template <class I>
class PairedIterable {
  const I ib, ie;
  PAIRED_ITERABLE_SHORT_TYPES(I)
  public:
  PairedIterable(I ib, I ie) noexcept : ib(ib), ie(ie) {}
  inline auto firsts()  noexcept { return staticTransformIterable(*this, PairFirst<K, V>()); }
  inline auto seconds() noexcept { return staticTransformIterable(*this, PairSecond<K, V>()); }
  inline auto keys()    noexcept { return firsts(); }
  inline auto values()  noexcept { return seconds(); }
  inline auto pairs()   noexcept { return *this; }
};

template <class I>
class ConstPairedIterable {
  const I ib, ie;
  PAIRED_ITERABLE_SHORT_TYPES(I)
  public:
  ConstPairedIterable(I ib, I ie) noexcept : ib(ib), ie(ie) {}
  inline auto firsts()  noexcept { return staticTransformIterable(*this, ConstPairFirst<K, V>()); }
  inline auto seconds() noexcept { return staticTransformIterable(*this, ConstPairSecond<K, V>()); }
  inline auto keys()    noexcept { return firsts(); }
  inline auto values()  noexcept { return seconds(); }
  inline auto pairs()   noexcept { return *this; }
};

template <class I>
class PairedValueIterable {
  const I ib, ie;
  PAIRED_ITERABLE_SHORT_TYPES(I)
  public:
  PairedValueIterable(I ib, I ie) noexcept : ib(ib), ie(ie) {}
  inline auto firsts()  noexcept { return staticTransformIterable(*this, PairFirstValue<K, V>()); }
  inline auto seconds() noexcept { return staticTransformIterable(*this, PairSecondValue<K, V>()); }
  inline auto keys()    noexcept { return firsts(); }
  inline auto values()  noexcept { return seconds(); }
  inline auto pairs()   noexcept { return *this; }
};


template <class I>
inline auto paired_iterable(I ib, I ie) noexcept {
  return PairedIterable<I>(ib, ie);
}
template <class I>
inline auto const_paired_iterable(I ib, I ie) noexcept {
  return ConstPairedIterable<I>(ib, ie);
}
template <class I>
inline auto paired_value_iterable(I ib, I ie) noexcept {
  return PairedValueIterable<I>(ib, ie);
}

template <class J>
inline auto pairedIterable(const J& x) {
  return paired_iterable(x.begin(), x.end());
}
template <class J>
inline auto constPairedIterable(const J& x) {
  return const_paired_iterable(x.begin(), x.end());
}
template <class J>
inline auto pairedValueIterable(const J& x) {
  return paired_value_iterable(x.begin(), x.end());
}




// SIZE
// ----

template <class I>
inline size_t size(const SizedIterable<I>& x) {
  return x.size();
}




// FAST SIZE
// ---------
// Fast size in O(1).

template <class T>
inline size_t fastSize(const vector<T>& x) noexcept {
  return x.size();
}
template <class I>
inline size_t fastSize(const SizedIterable<I>& x) noexcept {
  return x.size();
}
template <class J>
inline size_t fastSize(const J& x) noexcept {
  return size_t(-1);
}




// SLICE
// -----

template <class J>
inline auto sliceIterable(const J& x, size_t i) {
  auto b = x.begin(), e = x.end();
  return iterable(b+i<e? b+i:e, e);
}
template <class J>
inline auto sliceIterable(const J& x, size_t i, size_t I) {
  auto b = x.begin(), e = x.end();
  return iterable(b+i<e? b+i:e, b+I<e? b+I:e);
}




// DEFAULT ITERATOR
// ----------------
// Return default value of type, always.

template <class T>
class DefaultIterator {
  using iterator = DefaultIterator;
  const T x;
  public:
  ITERATOR_USING(random_access_iterator_tag, ptrdiff_t, T, const T&, const T*)
  DefaultIterator() noexcept : x() {}
  ITERATOR_DEREF(inline, const noexcept, x)
  ITERATOR_POINTER(inline, const noexcept, &x)
  ITERATOR_LOOKUP(inline, const noexcept, i, x)
  ITERATOR_NEXT(inline, noexcept,,)
  ITERATOR_ADVANCE(inline, noexcept, i,,)
  ITERATOR_ARITHMETIC(inline, noexcept, l, i, l, l)
  ITERATOR_DIFFERENCE(inline, noexcept, l, r, 0)
  ITERATOR_COMPARE_INT(inline, noexcept, l, r, 0)
};

template <class T>
class DefaultValueIterator {
  using iterator = DefaultValueIterator;
  public:
  ITERATOR_USING(random_access_iterator_tag, ptrdiff_t, T, T, const T*)
  DefaultValueIterator() noexcept {}
  ITERATOR_DEREF_VALUE(inline, const noexcept, T())
  ITERATOR_LOOKUP(inline, const noexcept, i, T())
  ITERATOR_NEXT(inline, noexcept,,)
  ITERATOR_ADVANCE(inline, noexcept, i,,)
  ITERATOR_ARITHMETIC(inline, noexcept, l, i, l, l)
  ITERATOR_DIFFERENCE(inline, noexcept, l, r, 0)
  ITERATOR_COMPARE_INT(inline, noexcept, l, r, 0)
};

template <class T>
inline auto default_iterator(const T& _) noexcept {
  return DefaultIterator<T>();
}
template <class T>
inline auto default_value_iterator(const T& _) noexcept {
  return DefaultValueIterator<T>();
}




// RANGE ITERATOR
// --------------

template <class T>
class RangeIterator {
  using iterator = RangeIterator;
  T v;
  public:
  ITERATOR_USING(random_access_iterator_tag, T, T, T, const T*)
  RangeIterator(T v) noexcept : v(v) {}
  ITERATOR_DEREF(inline, const noexcept, v)
  ITERATOR_POINTER(inline, const noexcept, &v)
  ITERATOR_LOOKUP(inline, const noexcept, n, v+n)
  ITERATOR_NEXT(inline, noexcept, ++v, --v)
  ITERATOR_ADVANCE(inline, noexcept, n, v+=n, v-=n)
  ITERATOR_ARITHMETIC(inline, noexcept, l, n, iterator(l.v+n), iterator(l.v-n))
  ITERATOR_DIFFERENCE(inline, noexcept, l, r, l.v-r.v)
  ITERATOR_COMPARE_INT(inline, noexcept, l, r, l.v-r.v)
};

template <class T>
class RangeFactor {
  const T v, DV;
  public:
  RangeFactor(T v, T DV) noexcept : v(v), DV(DV) {}
  inline T operator()(size_t n) const noexcept { return v+DV*n; }
};


template <class T>
inline size_t rangeSize(T v, T V, T DV=1) {
  return max(size_t(0), size_t(ceilDiv(V-v, DV)));
}
template <class T>
inline size_t rangeLast(T v, T V, T DV=1) {
  return v + DV*(rangeSize(v, V, DV) - 1);
}

template <class T>
inline auto rangeIterator(T v) noexcept { return RangeIterator<T>(v); }

template <class T>
inline auto rangeIterable(T V) noexcept {
  RangeIterator<T> b(0), e(V);
  return iterable(b, e);
}
template <class T>
inline auto rangeIterable(T v, T V, T DV=1) {
  auto x = rangeIterable(rangeSize(v, V, DV));
  return transformIterable(x, [=](size_t n) { return v+DV*n; });
}

template <class T>
inline auto rangeVector(T V) noexcept {
  auto x = rangeIterable(V);
  return vector<T>(x.begin(), x.end());
}
template <class T>
inline auto rangeVector(T v, T V, T DV=1) {
  auto x = rangeIterable(v, V, DV);
  return vector<T>(x.begin(), x.end());
}




// CIRCULAR ITERATOR
// -----------------

#ifndef CIRCULAR_ITERATOR_ACCESS
#define CIRCULAR_ITERATOR_ACCESS(I, ib, ie, it) \
  inline I valueBegin() const noexcept { return ib; } \
  inline I valueEnd()   const noexcept { return ie; } \
  inline I value()      const noexcept { return it; }
#endif

#ifndef CIRCULAR_ITERABLE_ACCESS
#define CIRCULAR_ITERABLE_ACCESS(I, xb, xe) \
  inline auto values() const noexcept { return Iterable<I>(xb, xe); }
#endif


template <class I>
class InputCircularIterator {
  using iterator = InputCircularIterator;
  const I xb, xe; I it;
  public:
  ITERATOR_USING_XC(I, input_iterator_tag)
  InputCircularIterator(I xb, I xe, I it) noexcept : xb(xb), xe(xe), it(it) {}
  ITERATOR_DEREF_VALUE(inline, const, *it)
  ITERATOR_INCREMENT(inline,, ++it; if (it == xe) it = xb)
  ITERATOR_COMPARE_EQNE(inline,, l, r, l.it==r.it, l.it!=r.it)
  CIRCULAR_ITERATOR_ACCESS(I, xb, xe, it)
};


template <class I>
class InputCircularIterable {
  const I xb, xe;
  const I ib, ie;
  public:
  InputCircularIterable(I xb, I xe, I ib, I ie) noexcept : xb(xb), xe(xe), ib(ib), ie(ie) {}
  inline auto begin() const { return InputCircularIterator<I>(xb, xe, ib); }
  inline auto end()   const { return InputCircularIterator<I>(xb, xe, ie); }
  CIRCULAR_ITERABLE_ACCESS(I, xb, xe)
  ITERABLE_SIZES_DEFAULT(ib, ie)
};

template <class I>
inline auto input_circular_iterator(I ib, I ie, I it) noexcept {
  return InputCircularIterator<I>(ib, ie, it);
}

template <class I>
inline auto input_circular_iterable(I xb, I xe, I ib, I ie) noexcept {
  return InputCircularIterable<I>(xb, xe, ib, ie);
}
template <class J>
inline auto inputCircularIterable(const J& x, size_t i, size_t I) noexcept {
  return input_circular_iterable(x.begin(), x.end(), x.begin()+i, x.begin()+I);
}




// PAIR ITERATOR
// -------------

#ifndef PAIR_ITERATOR_ACCESS
#define PAIR_ITERATOR_ACCESS(I0, I1, i0, i1) \
  inline I0 first()  const noexcept { return i0; } \
  inline I1 second() const noexcept { return i1; } \
  inline I0 key()    const noexcept { return i0; } \
  inline I1 value()  const noexcept { return i1; }
#endif


template <class I0, class I1>
class InputPairIterator {
  using iterator = InputPairIterator;
  using V0 = typename I0::value_type;
  using V1 = typename I1::value_type;
  using VP = pair<V0, V1>;
  I0 i0; I1 i1;
  public:
  ITERATOR_USING_XCVRP(I0, input_iterator_tag, VP, VP, const VP*)
  InputPairIterator(I0 i0, I1 i1) noexcept : i0(i0), i1(i1) {}
  ITERATOR_DEREF_VALUE(inline, const, make_pair(*i0, *i1))
  ITERATOR_INCREMENT(inline,, ++i0; ++i1)
  ITERATOR_COMPARE_EQNE(inline,, l, r, l.i0==r.i0, l.i0!=r.i0)
  PAIR_ITERATOR_ACCESS(I0, I1, i0, i1)
};

template <class I0, class I1>
class OutputPairIterator {
  using iterator = OutputPairIterator;
  using V0 = typename I0::value_type;
  using V1 = typename I1::value_type;
  using VP = pair<V0, V1>;
  I0 i0; I1 i1;
  public:
  ITERATOR_USING_XCVRP(I0, output_iterator_tag, VP, VP, const VP*)
  OutputPairIterator(I0 i0, I1 i1) noexcept : i0(i0), i1(i1) {}
  ITERATOR_DEREF(inline, const, make_pair(*i0, *i1))
  ITERATOR_INCREMENT(inline,, ++i0; ++i1)
  PAIR_ITERATOR_ACCESS(I0, I1, i0, i1)
};

template <class I0, class I1>
class ForwardPairIterator {
  using iterator = ForwardPairIterator;
  using V0 = typename I0::value_type;
  using V1 = typename I1::value_type;
  using VP = pair<V0, V1>;
  I0 i0; I1 i1;
  public:
  ITERATOR_USING_XCVRP(I0, forward_iterator_tag, VP, VP, const VP*)
  ForwardPairIterator(I0 i0, I1 i1) noexcept : i0(i0), i1(i1) {}
  ITERATOR_DEREF(inline, const, make_pair(*i0, *i1))
  ITERATOR_INCREMENT(inline,, ++i0; ++i1)
  ITERATOR_COMPARE_EQNE(inline,, l, r, l.i0==r.i0, l.i0!=r.i0)
  PAIR_ITERATOR_ACCESS(I0, I1, i0, i1)
};

template <class I0, class I1>
class BidirectionalPairIterator {
  using iterator = BidirectionalPairIterator;
  using V0 = typename I0::value_type;
  using V1 = typename I1::value_type;
  using VP = pair<V0, V1>;
  I0 i0; I1 i1;
  public:
  ITERATOR_USING_XCVRP(I0, bidirectional_iterator_tag, VP, VP, const VP*)
  BidirectionalPairIterator(I0 i0, I1 i1) noexcept : i0(i0), i1(i1) {}
  ITERATOR_DEREF(inline, const, make_pair(*i0, *i1))
  ITERATOR_NEXT(inline,, ++i0; ++i1, --i0; --i1)
  ITERATOR_COMPARE_EQNE(inline,, l, r, l.i0==r.i0, l.i0!=r.i0)
  PAIR_ITERATOR_ACCESS(I0, I1, i0, i1)
};

template <class I0, class I1>
class RandomAccessPairIterator {
  using iterator = RandomAccessPairIterator;
  using V0 = typename I0::value_type;
  using V1 = typename I1::value_type;
  using VP = pair<V0, V1>;
  I0 i0; I1 i1;
  public:
  ITERATOR_USING_XCVRP(I0, random_access_iterator_tag, VP, VP, const VP*)
  RandomAccessPairIterator(I0 i0, I1 i1) noexcept : i0(i0), i1(i1) {}
  ITERATOR_DEREF(inline, const, make_pair(*i0, *i1))
  ITERATOR_LOOKUP(inline, const, i, make_pair(i0[i], i1[i]))
  ITERATOR_NEXT(inline,, ++i0; ++i1, --i0; --i1)
  ITERATOR_ADVANCE(inline,, n, i0+=n; i1+=n, i0-=n; i1-=n)
  ITERATOR_ARITHMETIC(inline,, l, n, iterator(l.i0+n, l.i1+n), iterator(l.i0-n, l.i1-n))
  ITERATOR_DIFFERENCE(inline,, l, r, l.i0-r.i0)
  ITERATOR_COMPARE_EQNE(inline,, l, r, l.i0==r.i0, l.i0!=r.i0)
  ITERATOR_COMPARE_GTLT(inline,, l, r, l.i0>r.i0, l.i0<r.i0)
  ITERATOR_COMPARE_GELE(inline,, l, r, l.i0>=r.i0, l.i0<=r.i0)
  PAIR_ITERATOR_ACCESS(I0, I1, i0, i1)
};

template <class I0, class I1>
inline auto input_pair_iterator(I0 i0, I1 i1) noexcept {
  return InputPairIterator<I0, I1>(i0, i1);
}
template <class I0, class I1>
inline auto output_pair_iterator(I0 i0, I1 i1) noexcept {
  return OutputPairIterator<I0, I1>(i0, i1);
}
template <class I0, class I1>
inline auto forward_pair_iterator(I0 i0, I1 i1) noexcept {
  return ForwardPairIterator<I0, I1>(i0, i1);
}
template <class I0, class I1>
inline auto bidirectional_pair_iterator(I0 i0, I1 i1) noexcept {
  return BidirectionalPairIterator<I0, I1>(i0, i1);
}
template <class I0, class I1>
inline auto random_access_pair_iterator(I0 i0, I1 i1) noexcept {
  return RandomAccessPairIterator<I0, I1>(i0, i1);
}

template <class I0, class I1>
inline auto input_pair_iterable(I0 i0b, I0 i0e, I1 i1b, I1 i1e) noexcept {
  InputPairIterator<I0, I1> b(i0b, i1b), e(i0e, i1e);
  return iterable(b, e);
}
template <class I0, class I1>
inline auto output_pair_iterable(I0 i0b, I0 i0e, I1 i1b, I1 i1e) noexcept {
  OutputPairIterator<I0, I1> b(i0b, i1b), e(i0e, i1e);
  return iterable(b, e);
}
template <class I0, class I1>
inline auto forward_pair_iterable(I0 i0b, I0 i0e, I1 i1b, I1 i1e) noexcept {
  ForwardPairIterator<I0, I1> b(i0b, i1b), e(i0e, i1e);
  return iterable(b, e);
}
template <class I0, class I1>
inline auto bidirectional_pair_iterable(I0 i0b, I0 i0e, I1 i1b, I1 i1e) noexcept {
  BidirectionalPairIterator<I0, I1> b(i0b, i1b), e(i0e, i1e);
  return iterable(b, e);
}
template <class I0, class I1>
inline auto random_access_pair_iterable(I0 i0b, I0 i0e, I1 i1b, I1 i1e) noexcept {
  RandomAccessPairIterator<I0, I1> b(i0b, i1b), e(i0e, i1e);
  return iterable(b, e);
}

template <class J0, class J1>
inline auto inputPairIterable(const J0& x0, const J1& x1) {
  return input_pair_iterable(x0.begin(), x0.end(), x1.begin(), x1.end());
}
template <class J0, class J1>
inline auto outputPairIterable(const J0& x0, const J1& x1) {
  return output_pair_iterable(x0.begin(), x0.end(), x1.begin(), x1.end());
}
template <class J0, class J1>
inline auto forwardPairIterable(const J0& x0, const J1& x1) {
  return forward_pair_iterable(x0.begin(), x0.end(), x1.begin(), x1.end());
}
template <class J0, class J1>
inline auto bidirectionalPairIterable(const J0& x0, const J1& x1) {
  return bidirectional_pair_iterable(x0.begin(), x0.end(), x1.begin(), x1.end());
}
template <class J0, class J1>
inline auto randomAccessPairIterable(const J0& x0, const J1& x1) {
  return random_access_pair_iterable(x0.begin(), x0.end(), x1.begin(), x1.end());
}

template <class I0, class I1>
inline auto pair_iterator(I0 i0, I1 i1) noexcept {
  return random_access_pair_iterator(i0, i1);
}

template <class I0, class I1>
inline auto pair_iterable(I0 i0b, I0 i0e, I1 i1b, I1 i1e) noexcept {
  return random_access_pair_iterable(i0b, i0e, i1b, i1e);
}
template <class J0, class J1>
inline auto pairIterable(const J0& x0, const J1& x1) {
  return pair_iterable(x0.begin(), x0.end(), x1.begin(), x1.end());
}




// FILTER ITERATOR
// ---------------

#ifndef FILTER_ITERATOR_ACCESS
#define FILTER_ITERATOR_ACCESS(I, F, it, ie, fn) \
  inline I value()     const noexcept { return it; } \
  inline I valueEnd()  const noexcept { return ie; } \
  inline F predicate() const noexcept { return fn; }
#endif

#ifndef FILTER_ITERABLE_ACCESS
#define FILTER_ITERABLE_ACCESS(I, F, ib, ie, fn) \
  inline auto values() const noexcept { return Iterable<I>(ib, ie); } \
  inline F predicate() const noexcept { return fn; }
#endif


template <class I, class F>
class InputFilterIterator {
  using iterator = InputFilterIterator;
  I it; const I ie; const F fn;
  inline void next() { while (it!=ie && !fn(*it)) ++it; }
  public:
  ITERATOR_USING_XC(I, input_iterator_tag);
  InputFilterIterator(I it, I ie, F fn) : it(it), ie(ie), fn(fn) { next(); }
  ITERATOR_DEREF_VALUE(inline, const, *it)
  ITERATOR_POINTER(inline, const noexcept, it.I::operator->())
  ITERATOR_INCREMENT(inline,, ++it; next())
  ITERATOR_COMPARE_EQNE(inline,, l, r, l.it==r.it, l.it!=r.it)
  FILTER_ITERATOR_ACCESS(I, F, it, ie, fn)
};

template <class I, class F>
class ForwardFilterIterator {
  using iterator = ForwardFilterIterator;
  I it; const I ie; const F fn;
  inline void next() { while (it!=ie && !fn(*it)) ++it; }
  public:
  ITERATOR_USING_XC(I, forward_iterator_tag);
  ForwardFilterIterator(I it, I ie, F fn) : it(it), ie(ie), fn(fn) { next(); }
  ITERATOR_DEREF(inline, const, *it)
  ITERATOR_POINTER(inline, const noexcept, it.I::operator->())
  ITERATOR_INCREMENT(inline,, ++it; next())
  ITERATOR_COMPARE_EQNE(inline,, l, r, l.it==r.it, l.it!=r.it)
  FILTER_ITERATOR_ACCESS(I, F, it, ie, fn)
};

template <class I, class F>
class InputFilterIterable {
  const I ib, ie;
  const F fn;
  public:
  InputFilterIterable(I ib, I ie, F fn) noexcept : ib(ib), ie(ie), fn(fn) {}
  inline auto begin() const { return InputFilterIterator<I, F>(ib, ie, fn); }
  inline auto end()   const { return InputFilterIterator<I, F>(ie, ie, fn); }
  FILTER_ITERABLE_ACCESS(I, F, ib, ie, fn)
  ITERABLE_SIZES_DEFAULT(ib, ie)
};

template <class I, class F>
class ForwardFilterIterable {
  const I ib, ie;
  const F fn;
  public:
  ForwardFilterIterable(I ib, I ie, F fn) noexcept : ib(ib), ie(ie), fn(fn) {}
  inline auto begin() const { return ForwardFilterIterator<I, F>(ib, ie, fn); }
  inline auto end()   const { return ForwardFilterIterator<I, F>(ie, ie, fn); }
  FILTER_ITERABLE_ACCESS(I, F, ib, ie, fn)
  ITERABLE_SIZES_DEFAULT(ib, ie)
};

template <class I, class F>
inline auto input_filter_iterator(I it, I ie, F fn) {
  return InputFilterIterator<I, F>(it, ie, fn);
}
template <class I, class F>
inline auto forward_filter_iterator(I it, I ie, F fn) {
  return ForwardFilterIterator<I, F>(it, ie, fn);
}

template <class I, class F>
inline auto input_filter_iterable(I ib, I ie, F fn) {
  return InputFilterIterable<I, F>(ib, ie, fn);
}
template <class I, class F>
inline auto forward_filter_iterable(I ib, I ie, F fn) {
  return ForwardFilterIterable<I, F>(ib, ie, fn);
}

template <class J, class F>
inline auto inputFilterIterable(const J& x, F fn) {
  return input_filter_iterable(x.begin(), x.end(), fn);
}
template <class J, class F>
inline auto forwardFilterIterable(const J& x, F fn) {
  return forward_filter_iterable(x.begin(), x.end(), fn);
}

template <class I, class F>
inline auto filter_iterator(I it, I ie, F fn) {
  return forward_filter_iterator(it, ie, fn);
}

template <class I, class F>
inline auto filter_iterable(I ib, I ie, F fn) {
  return forward_filter_iterable(ib, ie, fn);
}
template <class J, class F>
inline auto filterIterable(const J& x, F fn) {
  return filter_iterable(x.begin(), x.end(), fn);
}




// CONDITIONAL ITERATOR
// --------------------

#ifndef CONDITIONAL_ITERATOR_ACCESS
#define CONDITIONAL_ITERATOR_ACCESS(I, IC, it, ie, ic) \
  inline I  value()     const noexcept { return it; } \
  inline I  valueEnd()  const noexcept { return ie; } \
  inline IC condition() const noexcept { return ic; }
#endif

#ifndef CONDITIONAL_ITERABLE_ACCESS
#define CONDITIONAL_ITERABLE_ACCESS(I, IC, ib, ie, ic) \
  inline auto values()     const noexcept { return Iterable<I>(ib, ie); } \
  inline auto conditions() const noexcept { return Iterable<I>(ic, ic+distance(ib, ie)); }
#endif


template <class I, class IC>
class InputConditionalIterator {
  using iterator = InputConditionalIterator;
  I it; const I ie; IC ic;
  inline void next() { while (it!=ie && !(*ic)) { ++it; ++ic; } }
  public:
  ITERATOR_USING_XC(I, input_iterator_tag);
  InputConditionalIterator(I it, I ie, IC ic) : it(it), ie(ie), ic(ic) { next(); }
  ITERATOR_DEREF_VALUE(inline, const, *it)
  ITERATOR_POINTER(inline, const noexcept, it.I::operator->())
  ITERATOR_INCREMENT(inline,, ++it; ++ic; next())
  ITERATOR_COMPARE_EQNE(inline,, l, r, l.it==r.it, l.it!=r.it)
  CONDITIONAL_ITERATOR_ACCESS(I, IC, it, it, ic)
};

template <class I, class IC>
class ForwardConditionalIterator {
  using iterator = ForwardConditionalIterator;
  I it; const I ie; IC ic;
  inline void next() { while (it!=ie && !(*ic)) { ++it; ++ic; } }
  public:
  ITERATOR_USING_XC(I, forward_iterator_tag);
  ForwardConditionalIterator(I it, I ie, IC ic) : it(it), ie(ie), ic(ic) { next(); }
  ITERATOR_DEREF(inline, const, *it)
  ITERATOR_POINTER(inline, const noexcept, it.I::operator->())
  ITERATOR_INCREMENT(inline,, ++it; ++ic; next())
  ITERATOR_COMPARE_EQNE(inline,, l, r, l.it==r.it, l.it!=r.it)
  CONDITIONAL_ITERATOR_ACCESS(I, IC, it, it, ic)
};

template <class I, class IC>
class InputConditionalIterable {
  const I ib, ie;
  const IC ic;
  public:
  InputConditionalIterable(I ib, I ie, IC ic) noexcept : ib(ib), ie(ie), ic(ic) {}
  inline auto begin() const { return InputConditionalIterator<I, IC>(ib, ie, ic); }
  inline auto end()   const { return InputConditionalIterator<I, IC>(ie, ie, ic); }
  CONDITIONAL_ITERABLE_ACCESS(I, IC, ib, ie, ic)
  ITERABLE_SIZES_DEFAULT(ib, ie)
};

template <class I, class IC>
class ForwardConditionalIterable {
  const I ib, ie;
  const IC ic;
  public:
  ForwardConditionalIterable(I ib, I ie, IC ic) noexcept : ib(ib), ie(ie), ic(ic) {}
  inline auto begin() const { return ForwardConditionalIterator<I, IC>(ib, ie, ic); }
  inline auto end()   const { return ForwardConditionalIterator<I, IC>(ie, ie, ic); }
  CONDITIONAL_ITERABLE_ACCESS(I, IC, ib, ie, ic)
  ITERABLE_SIZES_DEFAULT(ib, ie)
};

template <class I, class IC>
inline auto input_conditional_iterator(I it, I ie, IC ic) {
  return InputConditionalIterator<I, IC>(it, ie, ic);
}
template <class I, class IC>
inline auto forward_conditional_iterator(I it, I ie, IC ic) {
  return ForwardConditionalIterator<I, IC>(it, ie, ic);
}

template <class I, class IC>
inline auto input_conditional_iterable(I ib, I ie, IC ic) {
  return InputConditionalIterable<I, IC>(ib, ie, ic);
}
template <class I, class IC>
inline auto forward_conditional_iterable(I ib, I ie, IC ic) {
  return ForwardConditionalIterable<I, IC>(ib, ie, ic);
}

template <class J, class JC>
inline auto inputConditionalIterable(const J& x, const JC& c) {
  return input_conditional_iterable(x.begin(), x.end(), c.begin());
}
template <class J, class JC>
inline auto forwardConditionalIterable(const J& x, const JC& c) {
  return forward_conditional_iterable(x.begin(), x.end(), c.begin());
}

template <class I, class IC>
inline auto conditional_iterator(I it, I ie, IC ic) {
  return forward_conditional_iterator(it, ie, ic);
}

template <class I, class IC>
inline auto conditional_iterable(I ib, I ie, IC ic) {
  return forward_conditional_iterable(ib, ie, ic);
}
template <class J, class JC>
inline auto conditionalIterable(const J& x, const JC& c) {
  return conditional_iterable(x.begin(), x.end(), c.begin());
}




// TRANSFORM ITERATOR
// ------------------

#ifndef TRANSFORM_ITERATOR_ACCESS
#define TRANSFORM_ITERATOR_ACCESS(I, F, it, fn) \
  inline I value()     const noexcept { return it; } \
  inline F operation() const noexcept { return fn; }
#endif

#ifndef TRANSFORM_ITERABLE_ACCESS
#define TRANSFORM_ITERABLE_ACCESS(I, F, ib, ie, fn) \
  inline auto values() const noexcept { return Iterable<I>(ib, ie); } \
  inline F operation() const noexcept { return fn; }
#endif


template <class I, class F>
class InputTransformIterator {
  using iterator = InputTransformIterator;
  I it; const F fn;
  public:
  ITERATOR_USING_XCRP(I, input_iterator_tag, decltype(fn(*it)), const value_type*)
  InputTransformIterator(I it, F fn) noexcept : it(it), fn(fn) {}
  ITERATOR_DEREF_VALUE(inline, const, fn(*it))
  ITERATOR_INCREMENT(inline,, ++it)
  ITERATOR_COMPARE_EQNE(inline,, l, r, l.it==r.it, l.it!=r.it)
  TRANSFORM_ITERATOR_ACCESS(I, F, it, fn)
};

template <class I, class F>
class OutputTransformIterator {
  using iterator = OutputTransformIterator;
  I it; const F fn;
  public:
  ITERATOR_USING_XCRP(I, output_iterator_tag, decltype(fn(*it)), const value_type*)
  OutputTransformIterator(I it, F fn) noexcept : it(it), fn(fn) {}
  ITERATOR_DEREF(inline, const, fn(*it))
  ITERATOR_INCREMENT(inline,, ++it)
  TRANSFORM_ITERATOR_ACCESS(I, F, it, fn)
};

template <class I, class F>
class ForwardTransformIterator {
  using iterator = ForwardTransformIterator;
  I it; const F fn;
  public:
  ITERATOR_USING_XCRP(I, forward_iterator_tag, decltype(fn(*it)), const value_type*)
  ForwardTransformIterator(I it, F fn) noexcept : it(it), fn(fn) {}
  ITERATOR_DEREF(inline, const, fn(*it))
  ITERATOR_INCREMENT(inline,, ++it)
  ITERATOR_COMPARE_EQNE(inline,, l, r, l.it==r.it, l.it!=r.it)
  TRANSFORM_ITERATOR_ACCESS(I, F, it, fn)
};

template <class I, class F>
class BidirectionalTransformIterator {
  using iterator = BidirectionalTransformIterator;
  I it; const F fn;
  public:
  ITERATOR_USING_XCRP(I, bidirectional_iterator_tag, decltype(fn(*it)), const value_type*)
  BidirectionalTransformIterator(I it, F fn) noexcept : it(it), fn(fn) {}
  ITERATOR_DEREF(inline, const, fn(*it))
  ITERATOR_NEXT(inline,, ++it, --it)
  ITERATOR_COMPARE_EQNE(inline,, l, r, l.it==r.it, l.it!=r.it)
  TRANSFORM_ITERATOR_ACCESS(I, F, it, fn)
};

template <class I, class F>
class RandomAccessTransformIterator {
  using iterator = RandomAccessTransformIterator;
  I it; const F fn;
  public:
  ITERATOR_USING_XCRP(I, random_access_iterator_tag, decltype(fn(*it)), const value_type*)
  RandomAccessTransformIterator(I it, F fn) noexcept : it(it), fn(fn) {}
  ITERATOR_DEREF(inline, const, fn(*it))
  ITERATOR_LOOKUP(inline, const, i, fn(it[i]))
  ITERATOR_NEXT(inline,, ++it, --it)
  ITERATOR_ADVANCE(inline,, n, it+=n, it-=n)
  ITERATOR_ARITHMETIC(inline,, l, n, iterator(l.it+n, l.fn), iterator(l.it-n, l.fn))
  ITERATOR_DIFFERENCE(inline,, l, r, l.it-r.it)
  ITERATOR_COMPARE_EQNE(inline,, l, r, l.it==r.it, l.it!=r.it)
  ITERATOR_COMPARE_GTLT(inline,, l, r, l.it>r.it, l.it<r.it)
  ITERATOR_COMPARE_GELE(inline,, l, r, l.it>=r.it, l.it<=r.it)
  TRANSFORM_ITERATOR_ACCESS(I, F, it, fn)
};

template <class I, class F>
class InputTransformIterable {
  const I ib, ie;
  const F fn;
  public:
  InputTransformIterable(I ib, I ie, F fn) noexcept : ib(ib), ie(ie), fn(fn) {}
  inline auto begin() const { return InputTransformIterator<I, F>(ib, fn); }
  inline auto end()   const { return InputTransformIterator<I, F>(ie, fn); }
  TRANSFORM_ITERABLE_ACCESS(I, F, ib, ie, fn)
  ITERABLE_SIZES_DEFAULT(ib, ie)
};

template <class I, class F>
class OutputTransformIterable {
  const I ib, ie;
  const F fn;
  public:
  OutputTransformIterable(I ib, I ie, F fn) noexcept : ib(ib), ie(ie), fn(fn) {}
  inline auto begin() const { return OutputTransformIterator<I, F>(ib, fn); }
  inline auto end()   const { return OutputTransformIterator<I, F>(ie, fn); }
  TRANSFORM_ITERABLE_ACCESS(I, F, ib, ie, fn)
  ITERABLE_SIZES_DEFAULT(ib, ie)
};

template <class I, class F>
class ForwardTransformIterable {
  const I ib, ie;
  const F fn;
  public:
  ForwardTransformIterable(I ib, I ie, F fn) noexcept : ib(ib), ie(ie), fn(fn) {}
  inline auto begin() const { return ForwardTransformIterator<I, F>(ib, fn); }
  inline auto end()   const { return ForwardTransformIterator<I, F>(ie, fn); }
  TRANSFORM_ITERABLE_ACCESS(I, F, ib, ie, fn)
  ITERABLE_SIZES_DEFAULT(ib, ie)
};

template <class I, class F>
class BidirectionalTransformIterable {
  const I ib, ie;
  const F fn;
  public:
  BidirectionalTransformIterable(I ib, I ie, F fn) noexcept : ib(ib), ie(ie), fn(fn) {}
  inline auto begin() const { return BidirectionalTransformIterator<I, F>(ib, fn); }
  inline auto end()   const { return BidirectionalTransformIterator<I, F>(ie, fn); }
  TRANSFORM_ITERABLE_ACCESS(I, F, ib, ie, fn)
  ITERABLE_SIZES_DEFAULT(ib, ie)
};

template <class I, class F>
class RandomAccessTransformIterable {
  const I ib, ie;
  const F fn;
  public:
  RandomAccessTransformIterable(I ib, I ie, F fn) noexcept : ib(ib), ie(ie), fn(fn) {}
  inline auto begin() const { return RandomAccessTransformIterator<I, F>(ib, fn); }
  inline auto end()   const { return RandomAccessTransformIterator<I, F>(ie, fn); }
  TRANSFORM_ITERABLE_ACCESS(I, F, ib, ie, fn)
  ITERABLE_SIZES_DEFAULT(ib, ie)
};

template <class I, class F>
inline auto input_transform_iterator(I it, F fn) noexcept {
  return InputTransformIterator<I, F>(it, fn);
}
template <class I, class F>
inline auto output_transform_iterator(I it, F fn) noexcept {
  return OutputTransformIterator<I, F>(it, fn);
}
template <class I, class F>
inline auto forward_transform_iterator(I it, F fn) noexcept {
  return ForwardTransformIterator<I, F>(it, fn);
}
template <class I, class F>
inline auto bidirectional_transform_iterator(I it, F fn) noexcept {
  return BidirectionalTransformIterator<I, F>(it, fn);
}
template <class I, class F>
inline auto random_access_transform_iterator(I it, F fn) noexcept {
  return RandomAccessTransformIterator<I, F>(it, fn);
}

template <class I, class F>
inline auto input_transform_iterable(I ib, I ie, F fn) noexcept {
  return InputTransformIterable<I, F>(ib, ie, fn);
}
template <class I, class F>
inline auto output_transform_iterable(I ib, I ie, F fn) noexcept {
  return OutputTransformIterable<I, F>(ib, ie, fn);
}
template <class I, class F>
inline auto forward_transform_iterable(I ib, I ie, F fn) noexcept {
  return ForwardTransformIterable<I, F>(ib, ie, fn);
}
template <class I, class F>
inline auto bidirectional_transform_iterable(I ib, I ie, F fn) noexcept {
  return BidirectionalTransformIterable<I, F>(ib, ie, fn);
}
template <class I, class F>
inline auto random_access_transform_iterable(I ib, I ie, F fn) noexcept {
  return RandomAccessTransformIterable<I, F>(ib, ie, fn);
}

template <class J, class F>
inline auto inputTransformIterable(const J& x, F fn) {
  return input_transform_iterable(x.begin(), x.end(), fn);
}
template <class J, class F>
inline auto outputTransformIterable(const J& x, F fn) {
  return output_transform_iterable(x.begin(), x.end(), fn);
}
template <class J, class F>
inline auto forwardTransformIterable(const J& x, F fn) {
  return forward_transform_iterable(x.begin(), x.end(), fn);
}
template <class J, class F>
inline auto bidirectionalTransformIterable(const J& x, F fn) {
  return bidirectional_transform_iterable(x.begin(), x.end(), fn);
}
template <class J, class F>
inline auto randomAccessTransformIterable(const J& x, F fn) {
  return random_access_transform_iterable(x.begin(), x.end(), fn);
}

template <class I, class F>
inline auto transform_iterator(I it, F fn) noexcept {
  return random_access_transform_iterator(it, fn);
}

template <class I, class F>
inline auto transform_iterable(I ib, I ie, F fn) noexcept {
  return random_access_transform_iterable(ib, ie, fn);
}
template <class J, class F>
inline auto transformIterable(const J& x, F fn) {
  return transform_iterable(x.begin(), x.end(), fn);
}




// STATIC TRANSFORM ITERATOR
// -------------------------

#ifndef STATIC_TRANSFORM_ITERATOR
#define STATIC_TRANSFORM_ITERATOR(I, F, it) \
  inline I value()     const noexcept { return it; } \
  inline F operation() const noexcept { return F(); }
#endif


template <class I, class F>
class InputStaticTransformIterator {
  using iterator = InputStaticTransformIterator;
  I it;
  public:
  ITERATOR_USING_XCRP(I, input_iterator_tag, decltype(F()(*it)), const value_type*)
  InputStaticTransformIterator(I it) noexcept : it(it) {}
  ITERATOR_DEREF_VALUE(inline, const, F()(*it))
  ITERATOR_INCREMENT(inline,, ++it)
  ITERATOR_COMPARE_EQNE(inline,, l, r, l.it==r.it, l.it!=r.it)
  STATIC_TRANSFORM_ITERATOR(I, F, it)
};

template <class I, class F>
class OutputStaticTransformIterator {
  using iterator = OutputStaticTransformIterator;
  I it;
  public:
  ITERATOR_USING_XCRP(I, output_iterator_tag, decltype(F()(*it)), const value_type*)
  OutputStaticTransformIterator(I it) noexcept : it(it) {}
  ITERATOR_DEREF(inline, const, F()(*it))
  ITERATOR_INCREMENT(inline,, ++it)
  STATIC_TRANSFORM_ITERATOR(I, F, it)
};

template <class I, class F>
class ForwardStaticTransformIterator {
  using iterator = ForwardStaticTransformIterator;
  I it;
  public:
  ITERATOR_USING_XCRP(I, forward_iterator_tag, decltype(F()(*it)), const value_type*)
  ForwardStaticTransformIterator(I it) noexcept : it(it) {}
  ITERATOR_DEREF(inline, const, F()(*it))
  ITERATOR_INCREMENT(inline,, ++it)
  ITERATOR_COMPARE_EQNE(inline,, l, r, l.it==r.it, l.it!=r.it)
  STATIC_TRANSFORM_ITERATOR(I, F, it)
};

template <class I, class F>
class BidirectionalStaticTransformIterator {
  using iterator = BidirectionalStaticTransformIterator;
  I it;
  public:
  ITERATOR_USING_XCRP(I, bidirectional_iterator_tag, decltype(F()(*it)), const value_type*)
  BidirectionalStaticTransformIterator(I it) noexcept : it(it) {}
  ITERATOR_DEREF(inline, const, F()(*it))
  ITERATOR_NEXT(inline,, ++it, --it)
  ITERATOR_COMPARE_EQNE(inline,, l, r, l.it==r.it, l.it!=r.it)
  STATIC_TRANSFORM_ITERATOR(I, F, it)
};

template <class I, class F>
class RandomAccessStaticTransformIterator {
  using iterator = RandomAccessStaticTransformIterator;
  I it;
  public:
  ITERATOR_USING_XCRP(I, random_access_iterator_tag, decltype(F()(*it)), const value_type*)
  RandomAccessStaticTransformIterator(I it) noexcept : it(it) {}
  ITERATOR_DEREF(inline, const, F()(*it))
  ITERATOR_LOOKUP(inline, const, i, F()(it[i]))
  ITERATOR_NEXT(inline,, ++it, --it)
  ITERATOR_ADVANCE(inline,, n, it+=n, it-=n)
  ITERATOR_ARITHMETIC(inline,, l, n, iterator(l.it+n), iterator(l.it-n))
  ITERATOR_DIFFERENCE(inline,, l, r, l.it-r.it)
  ITERATOR_COMPARE_EQNE(inline,, l, r, l.it==r.it, l.it!=r.it)
  ITERATOR_COMPARE_GTLT(inline,, l, r, l.it>r.it, l.it<r.it)
  ITERATOR_COMPARE_GELE(inline,, l, r, l.it>=r.it, l.it<=r.it)
  STATIC_TRANSFORM_ITERATOR(I, F, it)
};

template <class I, class F>
inline auto input_static_transform_iterator(I it, F _) noexcept {
  return InputStaticTransformIterator<I, F>(it);
}
template <class I, class F>
inline auto output_static_transform_iterator(I it, F _) noexcept {
  return OutputStaticTransformIterator<I, F>(it);
}
template <class I, class F>
inline auto forward_static_transform_iterator(I it, F _) noexcept {
  return ForwardStaticTransformIterator<I, F>(it);
}
template <class I, class F>
inline auto bidirectional_static_transform_iterator(I it, F _) noexcept {
  return BidirectionalStaticTransformIterator<I, F>(it);
}
template <class I, class F>
inline auto random_access_static_transform_iterator(I it, F _) noexcept {
  return RandomAccessStaticTransformIterator<I, F>(it);
}

template <class I, class F>
inline auto input_static_transform_iterable(I ib, I ie, F _) noexcept {
  InputStaticTransformIterator<I, F> b(ib), e(ie);
  return iterable(b, e);
}
template <class I, class F>
inline auto output_static_transform_iterable(I ib, I ie, F _) noexcept {
  OutputStaticTransformIterator<I, F> b(ib), e(ie);
  return iterable(b, e);
}
template <class I, class F>
inline auto forward_static_transform_iterable(I ib, I ie, F _) noexcept {
  ForwardStaticTransformIterator<I, F> b(ib), e(ie);
  return iterable(b, e);
}
template <class I, class F>
inline auto bidirectional_static_transform_iterable(I ib, I ie, F _) noexcept {
  BidirectionalStaticTransformIterator<I, F> b(ib), e(ie);
  return iterable(b, e);
}
template <class I, class F>
inline auto random_access_static_transform_iterable(I ib, I ie, F _) noexcept {
  RandomAccessStaticTransformIterator<I, F> b(ib), e(ie);
  return iterable(b, e);
}

template <class J, class F>
inline auto inputStaticTransformIterable(const J& x, F _) {
  return input_static_transform_iterable(x.begin(), x.end(), _);
}
template <class J, class F>
inline auto outputStaticTransformIterable(const J& x, F _) {
  return output_static_transform_iterable(x.begin(), x.end(), _);
}
template <class J, class F>
inline auto forwardStaticTransformIterable(const J& x, F _) {
  return forward_static_transform_iterable(x.begin(), x.end(), _);
}
template <class J, class F>
inline auto bidirectionalStaticTransformIterable(const J& x, F _) {
  return bidirectional_static_transform_iterable(x.begin(), x.end(), _);
}
template <class J, class F>
inline auto randomAccessStaticTransformIterable(const J& x, F _) {
  return random_access_static_transform_iterable(x.begin(), x.end(), _);
}

template <class I, class F>
inline auto static_transform_iterator(I it, F _) noexcept {
  return random_access_static_transform_iterator(it, _);
}

template <class I, class F>
inline auto static_transform_iterable(I ib, I ie, F _) noexcept {
  return random_access_static_transform_iterable(ib, ie, _);
}
template <class J, class F>
inline auto staticTransformIterable(const J& x, F _) {
  return static_transform_iterable(x.begin(), x.end(), _);
}




// TERNARY ITERATOR
// ----------------
// Select iterator by boolean.

#ifndef TERNARY_ITERATOR_ACCESS
#define TERNARY_ITERATOR_ACCESS(I1, I0, sel, i1, i0) \
  inline bool select() const noexcept { return sel; } \
  inline I1 truthy()   const noexcept { return i1; } \
  inline I0 falsy()    const noexcept { return i0; }
#endif


template <class I1, class I0>
class InputTernaryIterator {
  using iterator = InputTernaryIterator;
  using T = typename I1::value_type;
  const bool sel; I1 i1; I0 i0;
  public:
  ITERATOR_USING_XC(I1, input_iterator_tag)
  InputTernaryIterator(bool sel, I1 i1, I0 i0) noexcept : sel(sel), i1(i1), i0(i0) {}
  ITERATOR_DEREF_VALUE(inline, const, sel? *i1 : *i0)
  ITERATOR_POINTER(inline, const, sel? i1.I1::operator->() : i0.I0::operator->())
  ITERATOR_INCREMENT(inline,, if (sel) ++i1; else ++i0)
  ITERATOR_COMPARE_EQ(inline,, l, r, l.sel==r.sel && (l.sel? l.i1==r.i1 : l.i0==r.i0))
  ITERATOR_COMPARE_NE(inline,, l, r, l.sel==r.sel && (l.sel? l.i1!=r.i1 : l.i0!=r.i0))
  TERNARY_ITERATOR_ACCESS(I1, I0, sel, i1, i0)
};

template <class I1, class I0>
class OutputTernaryIterator {
  using iterator = OutputTernaryIterator;
  using T = typename I1::value_type;
  const bool sel; I1 i1; I0 i0;
  public:
  ITERATOR_USING_XC(I1, output_iterator_tag)
  OutputTernaryIterator(bool sel, I1 i1, I0 i0) noexcept : sel(sel), i1(i1), i0(i0) {}
  ITERATOR_DEREF(inline, const, sel? *i1 : *i0)
  ITERATOR_INCREMENT(inline,, if (sel) ++i1; else ++i0)
  TERNARY_ITERATOR_ACCESS(I1, I0, sel, i1, i0)
};

template <class I1, class I0>
class ForwardTernaryIterator {
  using iterator = ForwardTernaryIterator;
  using T = typename I1::value_type;
  const bool sel; I1 i1; I0 i0;
  public:
  ITERATOR_USING_XC(I1, forward_iterator_tag)
  ForwardTernaryIterator(bool sel, I1 i1, I0 i0) noexcept : sel(sel), i1(i1), i0(i0) {}
  ITERATOR_DEREF(inline, const, sel? *i1 : *i0)
  ITERATOR_POINTER(inline, const, sel? i1.I1::operator->() : i0.I0::operator->())
  ITERATOR_INCREMENT(inline,, if (sel) ++i1; else ++i0)
  ITERATOR_COMPARE_EQ(inline,, l, r, l.sel==r.sel && (l.sel? l.i1==r.i1 : l.i0==r.i0))
  ITERATOR_COMPARE_NE(inline,, l, r, l.sel==r.sel && (l.sel? l.i1!=r.i1 : l.i0!=r.i0))
  TERNARY_ITERATOR_ACCESS(I1, I0, sel, i1, i0)
};

template <class I1, class I0>
class BidirectionalTernaryIterator {
  using iterator = BidirectionalTernaryIterator;
  using T = typename I1::value_type;
  const bool sel; I1 i1; I0 i0;
  public:
  ITERATOR_USING_XC(I1, bidirectional_iterator_tag)
  BidirectionalTernaryIterator(bool sel, I1 i1, I0 i0) noexcept : sel(sel), i1(i1), i0(i0) {}
  ITERATOR_DEREF(inline, const, sel? *i1 : *i0)
  ITERATOR_POINTER(inline, const, sel? i1.I1::operator->() : i0.I0::operator->())
  ITERATOR_INCREMENT(inline,, if (sel) ++i1; else ++i0)
  ITERATOR_DECREMENT(inline,, if (sel) --i1; else --i0)
  ITERATOR_COMPARE_EQ(inline,, l, r, l.sel==r.sel && (l.sel? l.i1==r.i1 : l.i0==r.i0))
  ITERATOR_COMPARE_NE(inline,, l, r, l.sel==r.sel && (l.sel? l.i1!=r.i1 : l.i0!=r.i0))
  TERNARY_ITERATOR_ACCESS(I1, I0, sel, i1, i0)
};

template <class I1, class I0>
class RandomAccessTernaryIterator {
  using iterator = RandomAccessTernaryIterator;
  using T = typename I1::value_type;
  const bool sel; I1 i1; I0 i0;
  public:
  ITERATOR_USING_XC(I1, random_access_iterator_tag)
  RandomAccessTernaryIterator(bool sel, I1 i1, I0 i0) noexcept : sel(sel), i1(i1), i0(i0) {}
  ITERATOR_DEREF(inline, const, sel? *i1 : *i0)
  ITERATOR_POINTER(inline, const, sel? i1.I1::operator->() : i0.I0::operator->())
  ITERATOR_LOOKUP(inline, const, i, sel? i1[i] : i0[i])
  ITERATOR_NEXT(inline,, if (sel) ++i1; else ++i0, if (sel) --i1; else --i0)
  ITERATOR_ADVANCE(inline,, n, if (sel) i1+=n; else i0+=n, if (sel) i1-=n; else i0-=n)
  ITERATOR_ARITHMETIC(inline,, l, n, iterator(l.sel, l.sel? l.i1+n:l.i1, l.sel? l.i0:l.i0+n), iterator(l.sel, l.sel? l.i1-n:l.i1, l.sel? l.i0:l.i0-n))
  ITERATOR_DIFFERENCE(inline,, l, r, l.sel==r.sel? (l.sel? l.i1-r.i1 : l.i0-r.i0) : 0)
  ITERATOR_COMPARE_EQ(inline,, l, r, l.sel==r.sel && (l.sel? l.i1==r.i1 : l.i0==r.i0))
  ITERATOR_COMPARE_NE(inline,, l, r, l.sel==r.sel && (l.sel? l.i1!=r.i1 : l.i0!=r.i0))
  ITERATOR_COMPARE_GT(inline,, l, r, l.sel==r.sel && (l.sel? l.i1>r.i1 : l.i0>r.i0))
  ITERATOR_COMPARE_LT(inline,, l, r, l.sel==r.sel && (l.sel? l.i1<r.i1 : l.i0<r.i0))
  ITERATOR_COMPARE_GE(inline,, l, r, l.sel==r.sel && (l.sel? l.i1>=r.i1 : l.i0>=r.i0))
  ITERATOR_COMPARE_LE(inline,, l, r, l.sel==r.sel && (l.sel? l.i1<=r.i1 : l.i0<=r.i0))
  TERNARY_ITERATOR_ACCESS(I1, I0, sel, i1, i0)
};

template <class I1, class I0>
inline auto input_ternary_iterator(bool sel, I1 i1, I0 i0) noexcept {
  return InputTernaryIterator<I1, I0>(sel, i1, i0);
}
template <class I1, class I0>
inline auto output_ternary_iterator(bool sel, I1 i1, I0 i0) noexcept {
  return OutputTernaryIterator<I1, I0>(sel, i1, i0);
}
template <class I1, class I0>
inline auto forward_ternary_iterator(bool sel, I1 i1, I0 i0) noexcept {
  return ForwardTernaryIterator<I1, I0>(sel, i1, i0);
}
template <class I1, class I0>
inline auto bidirectional_ternary_iterator(bool sel, I1 i1, I0 i0) noexcept {
  return BidirectionalTernaryIterator<I1, I0>(sel, i1, i0);
}
template <class I1, class I0>
inline auto random_access_ternary_iterator(bool sel, I1 i1, I0 i0) noexcept {
  return RandomAccessTernaryIterator<I1, I0>(sel, i1, i0);
}

template <class I1, class I0>
inline auto input_ternary_iterable(bool sel, I1 i1b, I1 i1e, I0 i0b, I0 i0e) noexcept {
  auto b = InputTernaryIterator<I1, I0>(sel, i1b, i0b);
  auto e = InputTernaryIterator<I1, I0>(sel, i1e, i0e);
  return iterable(b, e);
}
template <class I1, class I0>
inline auto output_ternary_iterable(bool sel, I1 i1b, I1 i1e, I0 i0b, I0 i0e) noexcept {
  auto b = OutputTernaryIterator<I1, I0>(sel, i1b, i0b);
  auto e = OutputTernaryIterator<I1, I0>(sel, i1e, i0e);
  return iterable(b, e);
}
template <class I1, class I0>
inline auto forward_ternary_iterable(bool sel, I1 i1b, I1 i1e, I0 i0b, I0 i0e) noexcept {
  auto b = ForwardTernaryIterator<I1, I0>(sel, i1b, i0b);
  auto e = ForwardTernaryIterator<I1, I0>(sel, i1e, i0e);
  return iterable(b, e);
}
template <class I1, class I0>
inline auto bidirectional_ternary_iterable(bool sel, I1 i1b, I1 i1e, I0 i0b, I0 i0e) noexcept {
  auto b = BidirectionalTernaryIterator<I1, I0>(sel, i1b, i0b);
  auto e = BidirectionalTernaryIterator<I1, I0>(sel, i1e, i0e);
  return iterable(b, e);
}
template <class I1, class I0>
inline auto random_access_ternary_iterable(bool sel, I1 i1b, I1 i1e, I0 i0b, I0 i0e) noexcept {
  auto b = RandomAccessTernaryIterator<I1, I0>(sel, i1b, i0b);
  auto e = RandomAccessTernaryIterator<I1, I0>(sel, i1e, i0e);
  return iterable(b, e);
}

template <class J1, class J0>
inline auto inputTernaryIterable(bool sel, const J1& x1, const J0& x0) {
  return input_ternary_iterable(sel, x1.begin(), x1.end(), x0.begin(), x0.end());
}
template <class J1, class J0>
inline auto outputTernaryIterable(bool sel, const J1& x1, const J0& x0) {
  return output_ternary_iterable(sel, x1.begin(), x1.end(), x0.begin(), x0.end());
}
template <class J1, class J0>
inline auto forwardTernaryIterable(bool sel, const J1& x1, const J0& x0) {
  return forward_ternary_iterable(sel, x1.begin(), x1.end(), x0.begin(), x0.end());
}
template <class J1, class J0>
inline auto bidirectionalTernaryIterable(bool sel, const J1& x1, const J0& x0) {
  return bidirectional_ternary_iterable(sel, x1.begin(), x1.end(), x0.begin(), x0.end());
}
template <class J1, class J0>
inline auto randomAccessTernaryIterable(bool sel, const J1& x1, const J0& x0) {
  return random_access_ternary_iterable(sel, x1.begin(), x1.end(), x0.begin(), x0.end());
}

template <class I1, class I0>
inline auto ternary_iterator(bool sel, I1 i1, I0 i0) noexcept {
  return random_access_transform_iterator(sel, i1, i0);
}

template <class I1, class I0>
inline auto ternary_iterable(bool sel, I1 i1b, I1 i1e, I0 i0b, I0 i0e) noexcept {
  return random_access_transform_iterable(sel, i1b, i1e, i0b, i0e);
}
template <class J1, class J0>
inline auto ternaryIterable(bool sel, const J1& x1, const J0& x0) {
  return ternary_iterable(sel, x1.begin(), x1.end(), x0.begin(), x0.end());
}




// SELECT ITERATOR
// ---------------
// Select iterator by index.
// Can be done using tuples.
