#pragma once
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <vector>
#include "_debug.hxx"

#ifdef OPENMP
#include <omp.h>
#endif

using std::vector;
using std::abs;
using std::min;
using std::max;
using std::fill;




// VECTOR 2D
// ---------

template <class T>
using vector2d = vector<vector<T>>;




// GATHER VALUES
// -------------

template <class TA, class TX, class IS, class FM>
inline void gatherValuesW(TA *a, const TX *x, const IS& is, FM fm) {
  ASSERT(a && x);
  size_t j = 0;
  for (auto i : is)
    a[j++] = TA(fm(x[i]));
}
template <class TA, class TX, class IS>
inline void gatherValuesW(TA *a, const TX *x, const IS& is) {
  auto fm = [](const auto& v) { return v; };
  gatherValuesW(a, x, is, fm);
}

template <class TA, class TX, class IS, class FM>
inline void gatherValuesW(vector<TA>& a, const vector<TX>& x, const IS& is, FM fm) {
  gatherValuesW(a.data(), x.data(), is, fm);
}
template <class TA, class TX, class IS>
inline void gatherValuesW(vector<TA>& a, const vector<TX>& x, const IS& is) {
  gatherValuesW(a.data(), x.data(), is);
}

template <class IS>
inline void gatherValuesW(vector<bool>& a, const vector<bool>& x, const IS& is) {
  size_t j = 0;
  for (auto i : is)
    a[j++] = x[i];
}


#ifdef OPENMP
template <class TA, class TX, class IS, class FM>
inline void gatherValuesOmpW(TA *a, const TX *x, const IS& is, FM fm) {
  ASSERT(a && x);
  size_t N = is.size();
  #pragma omp parallel for schedule(auto)
  for (size_t j=0; j<N; ++j)
    a[j] = TA(fm(x[is[j]]));
}
template <class TA, class TX, class IS>
inline void gatherValuesOmpW(TA *a, const TX *x, const IS& is) {
  auto fm = [](const auto& v) { return v; };
  gatherValuesOmpW(a, x, is, fm);
}

template <class TA, class TX, class IS, class FM>
inline void gatherValuesOmpW(vector<TA>& a, const vector<TX>& x, const IS& is, FM fm) {
  gatherValuesOmpW(a.data(), x.data(), is, fm);
}
template <class TA, class TX, class IS>
inline void gatherValuesOmpW(vector<TA>& a, const vector<TX>& x, const IS& is) {
  gatherValuesOmpW(a.data(), x.data(), is);
}
#endif




// SCATTER VALUES
// --------------

template <class TA, class TX, class IS, class FM>
inline void scatterValuesW(TA *a, const TX *x, const IS& is, FM fm) {
  ASSERT(a && x);
  size_t j = 0;
  for (auto i : is)
    a[i] = TA(fm(x[j++]));
}
template <class TA, class TX, class IS>
inline void scatterValuesW(TA *a, const TX *x, const IS& is) {
  auto fm = [](const auto& v) { return v; };
  scatterValuesW(a, x, is, fm);
}

template <class TA, class TX, class IS, class FM>
inline void scatterValuesW(vector<TA>& a, const vector<TX>& x, const IS& is, FM fm) {
  scatterValuesW(a.data(), x.data(), is, fm);
}
template <class TA, class TX, class IS>
inline void scatterValuesW(vector<TA>& a, const vector<TX>& x, const IS& is) {
  scatterValuesW(a.data(), x.data(), is);
}


#ifdef OPENMP
template <class TA, class TX, class IS, class FM>
inline void scatterValuesOmpW(TA *a, const TX *x, const IS& is, FM fm) {
  ASSERT(a && x);
  size_t N = is.size();
  #pragma omp parallel for schedule(auto)
  for (size_t j=0; j<N; ++j)
    a[is[j]] = TA(fm(x[j]));
}
template <class TA, class TX, class IS>
inline void scatterValuesOmpW(TA *a, const TX *x, const IS& is) {
  auto fm = [](const auto& v) { return v; };
  scatterValuesOmpW(a, x, is, fm);
}

template <class TA, class TX, class IS, class FM>
inline void scatterValuesOmpW(vector<TA>& a, const vector<TX>& x, const IS& is, FM fm) {
  scatterValuesOmpW(a.data(), x.data(), is, fm);
}
template <class TA, class TX, class IS>
inline void scatterValuesOmpW(vector<TA>& a, const vector<TX>& x, const IS& is) {
  scatterValuesOmpW(a.data(), x.data(), is);
}
#endif




// SCATTER OR
// ----------

template <class TA, class TX, class IS>
inline void scatterOrW(TA *a, const TX *x, const IS& is) {
  ASSERT(a && x);
  size_t j = 0;
  for (auto i : is)
    a[i] |= TA(x[j++]);
}
template <class TA, class TX, class IS>
inline void scatterOrW(vector<TA>& a, const vector<TX>& x, const IS& is) {
  scatterOrW(a.data(), x.data(), is);
}


#ifdef OPENMP
template <class TA, class TX, class IS>
inline void scatterOrOmpW(TA *a, const TX *x, const IS& is) {
  ASSERT(a && x);
  size_t N = is.size();
  #pragma omp parallel for schedule(auto)
  for (size_t j=0; j<N; ++j)
    a[is[j]] |= TA(x[j]);
}
template <class TA, class TX, class IS>
inline void scatterOrOmpW(vector<TA>& a, const vector<TX>& x, const IS& is) {
  scatterOrOmpW(a.data(), x.data(), is);
}
#endif




// VALUE INDICES
// -------------

template <class TA, class TX, class FM>
inline void valueIndicesW(vector2d<TA>& a, const vector<TX>& x, FM fm) {
  size_t N = x.size();
  for (size_t i=0; i<N; ++i) {
    TX  v = fm(x[i]);
    if (v <  TX()) continue;
    if (v >= a.size()) a.resize(v+1);
    a[v].push_back(TA(i));
  }
}
template <class TA, class TX, class FM>
inline vector2d<TA> valueIndicesAs(const vector<TX>& x, FM fm) {
  vector2d<TA> a; valueIndicesW(a, x, fm);
  return a;
}




// FILL VALUE
// ----------

template <class T>
inline void fillValueU(T *a, size_t N, const T& v) {
  ASSERT(a);
  fill(a, a+N, v);
}
template <class T>
inline void fillValueU(vector<T>& a, const T& v) {
  fill(a.begin(), a.end(), v);
}


#ifdef OPENMP
template <class T>
inline void fillValueOmpU(T *a, size_t N, const T& v) {
  ASSERT(a);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a[i] = v;
}
template <class T>
inline void fillValueOmpU(vector<T>& a, const T& v) {
  fillValueOmpU(a.data(), a.size(), v);
}
inline void fillValueOmpU(vector<bool>& a, const bool& v) {
  fill(a.begin(), a.end(), v);
}
#endif




// ADD VALUE
// ---------

template <class T>
inline void addValueU(T *a, size_t N, const T& v) {
  ASSERT(a);
  for (size_t i=0; i<N; ++i)
    a[i] += v;
}
template <class T>
inline void addValueU(vector<T>& a, const T& v) {
  addValueU(a.data(), a.size(), v);
}


#ifdef OPENMP
template <class T>
inline void addValueOmpU(T *a, size_t N, const T& v) {
  ASSERT(a);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a[i] += v;
}
template <class T>
inline void addValueOmpU(vector<T>& a, const T& v) {
  addValueOmpU(a.data(), a.size(), v);
}
#endif




// COPY VALUES
// -----------

template <class TA, class TX>
inline void copyValuesW(TA *a, const TX *x, size_t N) {
  ASSERT(a && x);
  for (size_t i=0; i<N; ++i)
    a[i] = x[i];
}
template <class TA, class TX>
inline void copyValuesW(vector<TA>& a, const vector<TX>& x) {
  return copyValuesW(a.data(), x.data(), x.size());
}


#ifdef OPENMP
template <class TA, class TX>
inline void copyValuesOmpW(TA *a, const TX *x, size_t N) {
  ASSERT(a && x);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a[i] = x[i];
}
template <class TA, class TX>
inline void copyValuesOmpW(vector<TA>& a, const vector<TX>& x) {
  return copyValuesOmpW(a.data(), x.data(), x.size());
}
#endif




// MULTIPLY VALUE
// --------------

template <class TA, class TX, class TV>
inline void multiplyValueW(TA *a, const TX *x, TV v, size_t N) {
  ASSERT(a && x);
  for (size_t i=0; i<N; ++i)
    a[i] = TA(x[i] * v);
}
template <class TA, class TX, class TV>
inline void multiplyValueW(vector<TA>& a, const vector<TX>& x, TV v) {
  multiplyValueW(a.data(), x.data(), v, x.size());
}


#ifdef OPENMP
template <class TA, class TX, class TV>
inline void multiplyValueOmpW(TA *a, const TX *x, TV v, size_t N) {
  ASSERT(a && x);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a[i] = TA(x[i] * v);
}
template <class TA, class TX, class TV>
inline void multiplyValueOmpW(vector<TA>& a, const vector<TX>& x, TV v) {
  multiplyValueOmpW(a.data(), x.data(), v, x.size());
}
#endif




// MULTIPLY VALUES
// ---------------

template <class TA, class TX, class TY>
inline void multiplyValuesW(TA *a, const TX *x, const TY *y, size_t N) {
  ASSERT(a && x && y);
  for (size_t i=0; i<N; ++i)
    a[i] = TA(x[i] * y[i]);
}
template <class TA, class TX, class TY>
inline void multiplyValuesW(vector<TA>& a, const vector<TX>& x, const vector<TY>& y) {
  multiplyValuesW(a.data(), x.data(), y.data(), x.size());
}


#ifdef OPENMP
template <class TA, class TX, class TY>
inline void multiplyValuesOmpW(TA *a, const TX *x, const TY *y, size_t N) {
  ASSERT(a && x && y);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a[i] = TA(x[i] * y[i]);
}
template <class TA, class TX, class TY>
inline void multiplyValuesOmpW(vector<TA>& a, const vector<TX>& x, const vector<TY>& y) {
  multiplyValuesOmpW(a.data(), x.data(), y.data(), x.size());
}
#endif




// SUM VALUES
// ----------


template <class TX, class TA=TX>
inline TA sumValues(const TX *x, size_t N, TA a=TA()) {
  ASSERT(x);
  for (size_t i=0; i<N; ++i)
    a += TA(x[i]);
  return a;
}
template <class TX, class TA=TX>
inline TA sumValues(const vector<TX>& x, TA a=TA()) {
  return sumValues(x.data(), x.size(), a);
}


#ifdef OPENMP
template <class TX, class TA=TX>
inline TA sumValuesOmp(const TX *x, size_t N, TA a=TA()) {
  ASSERT(x);
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (size_t i=0; i<N; ++i)
    a += TA(x[i]);
  return a;
}
template <class TX, class TA=TX>
inline TA sumValuesOmp(const vector<TX>& x, TA a=TA()) {
  return sumValuesOmp(x.data(), x.size(), a);
}
#endif




// L1-NORM
// -------

template <class TX, class TA=TX>
inline TA l1Norm(const TX *x, size_t N, TA a=TA()) {
  ASSERT(x);
  for (size_t i=0; i<N; ++i)
    a += TA(abs(x[i]));
  return a;
}
template <class TX, class TA=TX>
inline TA l1Norm(const vector<TX>& x, TA a=TA()) {
  return l1Norm(x.data(), x.size(), a);
}


#ifdef OPENMP
template <class TX, class TA=TX>
inline TA l1NormOmp(const TX *x, size_t N, TA a=TA()) {
  ASSERT(x);
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (size_t i=0; i<N; ++i)
    a += TA(abs(x[i]));
  return a;
}
template <class TX, class TA=TX>
inline TA l1NormOmp(const vector<TX>& x, TA a=TA()) {
  return l1NormOmp(x.data(), x.size(), a);
}
#endif




// L1-NORM DELTA
// -------------

template <class TX, class TY, class TA=TX>
inline TA l1NormDelta(const TX *x, const TY *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  for (size_t i=0; i<N; ++i)
    a += TA(abs(x[i] - y[i]));
  return a;
}
template <class TX, class TY, class TA=TX>
inline TA l1NormDelta(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return l1NormDelta(x.data(), y.data(), x.size(), a);
}


#ifdef OPENMP
template <class TX, class TY, class TA=TX>
inline TA l1NormOmpDelta(const TX *x, const TY *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (size_t i=0; i<N; ++i)
    a += TA(abs(x[i] - y[i]));
  return a;
}
template <class TX, class TY, class TA=TX>
inline TA l1NormOmpDelta(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return l1NormOmpDelta(x.data(), y.data(), x.size(), a);
}
#endif




// L1-NORM DELTA AT
// ----------------

template <class TX, class TY, class TI, class TA=TX>
inline TA l1NormDeltaAt(const TX *x, const TY *y, const TI *is, size_t IS, TA a=TA()) {
  ASSERT(x && y && is);
  for (size_t l=0; l<IS; ++l) {
    TI i = is[l];
    a += TA(abs(x[i] - y[i]));
  }
  return a;
}
template <class TX, class TY, class TI, class TA=TX>
inline TA l1NormDeltaAt(const vector<TX>& x, const vector<TY>& y, const vector<TI>& is, TA a=TA()) {
  return l1NormDeltaAt(x.data(), y.data(), is.data(), is.size(), a);
}


#ifdef OPENMP
template <class TX, class TY, class TI, class TA=TX>
inline TA l1NormDeltaAtOmp(const TX *x, const TY *y, const TI *is, size_t IS, TA a=TA()) {
  ASSERT(x && y && is);
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (size_t l=0; l<IS; ++l) {
    TI i = is[l];
    a += TA(abs(x[i] - y[i]));
  }
  return a;
}
template <class TX, class TY, class TI, class TA=TX>
inline TA l1NormDeltaAtOmp(const vector<TX>& x, const vector<TY>& y, const vector<TI>& is, TA a=TA()) {
  return l1NormDeltaAtOmp(x.data(), y.data(), is.data(), is.size(), a);
}
#endif




// L2-NORM
// -------

template <class TX, class TA=TX>
inline TA l2Norm(const TX *x, size_t N, TA a=TA()) {
  ASSERT(x);
  for (size_t i=0; i<N; ++i)
    a += TA(x[i]) * TA(x[i]);
  return a;
}
template <class TX, class TA=TX>
inline TA l2Norm(const vector<TX>& x, TA a=TA()) {
  return l2Norm(x.data(), x.size(), a);
}


#ifdef OPENMP
template <class TX, class TA=TX>
inline TA l2NormOmp(const TX *x, size_t N, TA a=TA()) {
  ASSERT(x);
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (size_t i=0; i<N; ++i)
    a += TA(x[i]) * TA(x[i]);
  return a;
}
template <class TX, class TA=TX>
inline TA l2NormOmp(const vector<TX>& x, TA a=TA()) {
  return l2NormOmp(x.data(), x.size(), a);
}
#endif




// L2-NORM DELTA
// -------------

template <class TX, class TY, class TA=TX>
inline TA l2NormDelta(const TX *x, const TY *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  for (size_t i=0; i<N; ++i)
    a += TA(x[i] - y[i]) * TA(x[i] - y[i]);
  return a;
}
template <class TX, class TY, class TA=TX>
inline TA l2NormDelta(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return l2NormDelta(x.data(), y.data(), x.size(), a);
}


#ifdef OPENMP
template <class TX, class TY, class TA=TX>
inline TA l2NormDeltaOmp(const TX *x, const TY *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (size_t i=0; i<N; ++i)
    a += TA(x[i] - y[i]) * TA(x[i] - y[i]);
  return a;
}
template <class TX, class TY, class TA=TX>
inline TA l2NormDeltaOmp(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return l2NormDeltaOmp(x.data(), y.data(), x.size(), a);
}
#endif




// L2-NORM DELTA AT
// ----------------

template <class TX, class TY, class TI, class TA=TX>
inline TA l2NormDeltaAt(const TX *x, const TY *y, const TI *is, size_t IS, TA a=TA()) {
  ASSERT(x && y && is);
  for (size_t l=0; l<IS; ++l) {
    TI i = is[l];
    a += TA(x[i] - y[i]) * TA(x[i] - y[i]);
  }
  return a;
}
template <class TX, class TY, class TI, class TA=TX>
inline TA l2NormDeltaAt(const vector<TX>& x, const vector<TY>& y, const vector<TI>& is, TA a=TA()) {
  return l2NormDeltaAt(x.data(), y.data(), is.data(), is.size(), a);
}


#ifdef OPENMP
template <class TX, class TY, class TI, class TA=TX>
inline TA l2NormOmp(const TX *x, const TY *y, const TI *is, size_t IS, TA a=TA()) {
  ASSERT(x && y && is);
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (size_t l=0; l<IS; ++l) {
    TI i = is[l];
    a += TA(x[i] - y[i]) * TA(x[i] - y[i]);
  }
  return a;
}
template <class TX, class TY, class TI, class TA=TX>
inline TA l2NormOmp(const vector<TX>& x, const vector<TY>& y, const vector<TI>& is, TA a=TA()) {
  return l2NormOmp(x.data(), y.data(), is.data(), is.size(), a);
}
#endif




// LI-NORM
// -------

template <class TX, class TA=TX>
inline TA liNorm(const TX *x, size_t N, TA a=TA()) {
  ASSERT(x);
  for (size_t i=0; i<N; ++i)
    a = max(a, TA(abs(x[i])));
  return a;
}
template <class TX, class TA=TX>
inline TA liNorm(const vector<TX>& x, TA a=TA()) {
  return liNorm(x.data(), x.size(), a);
}


#ifdef OPENMP
template <class TX, class TA=TX>
inline TA liNormOmp(const TX *x, size_t N, TA a=TA()) {
  ASSERT(x);
  #pragma omp parallel for schedule(auto) reduction(max:a)
  for (size_t i=0; i<N; ++i)
    a = max(a, TA(abs(x[i])));
  return a;
}
template <class TX, class TA=TX>
inline TA liNormOmp(const vector<TX>& x, TA a=TA()) {
  return liNormOmp(x.data(), x.size(), a);
}
#endif




// LI-NORM DELTA
// -------------

template <class TX, class TY, class TA=TX>
inline TA liNormDelta(const TX *x, const TY *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  for (size_t i=0; i<N; ++i)
    a = max(a, TA(abs(x[i] - y[i])));
  return a;
}
template <class TX, class TY, class TA=TX>
inline TA liNormDelta(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return liNormDelta(x.data(), y.data(), x.size(), a);
}


#ifdef OPENMP
template <class TX, class TY, class TA=TX>
inline TA liNormDeltaOmp(const TX *x, const TY *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  #pragma omp parallel for schedule(auto) reduction(max:a)
  for (size_t i=0; i<N; ++i)
    a = max(a, TA(abs(x[i] - y[i])));
  return a;
}
template <class TX, class TY, class TA=TX>
inline TA liNormDeltaOmp(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return liNormDeltaOmp(x.data(), y.data(), x.size(), a);
}
#endif




// LI-NORM DELTA AT
// ----------------

template <class TX, class TY, class TI, class TA=TX>
inline TA liNormDeltaAt(const TX *x, const TY *y, const TI *is, size_t IS, TA a=TA()) {
  ASSERT(x && y && is);
  for (size_t l=0; l<IS; ++l) {
    TI i = is[l];
    a = max(a, TA(abs(x[i] - y[i])));
  }
  return a;
}
template <class TX, class TY, class TI, class TA=TX>
inline TA liNormDeltaAt(const vector<TX>& x, const vector<TY>& y, const vector<TI>& is, TA a=TA()) {
  return liNormDeltaAt(x.data(), y.data(), is.data(), is.size(), a);
}


#ifdef OPENMP
template <class TX, class TY, class TI, class TA=TX>
inline TA liNormOmp(const TX *x, const TY *y, const TI *is, size_t IS, TA a=TA()) {
  ASSERT(x && y && is);
  #pragma omp parallel for schedule(auto) reduction(max:a)
  for (size_t l=0; l<IS; ++l) {
    TI i = is[l];
    a = max(a, TA(abs(x[i] - y[i])));
  }
  return a;
}
template <class TX, class TY, class TI, class TA=TX>
inline TA liNormOmp(const vector<TX>& x, const vector<TY>& y, const vector<TI>& is, TA a=TA()) {
  return liNormOmp(x.data(), y.data(), is.data(), is.size(), a);
}
#endif




// INCLUSIVE SCAN
// --------------

template <class TA, class TX>
inline TA inclusiveScanW(TA *a, const TX *x, size_t N, TA acc=TA()) {
  ASSERT(a && x);
  for (size_t i=0; i<N; ++i) {
    acc += x[i];
    a[i] = acc;
  }
  return acc;
}
template <class TA, class TX>
inline TA inclusiveScanW(vector<TA>& a, const vector<TX>& x, TA acc=TA()) {
  return inclusiveScanW(a.data(), x.data(), x.size(), acc);
}


#ifdef OPENMP
template <class TA, class TX>
inline TA inclusiveScanOmpW(TA *a, TA *buf, const TX *x, size_t N, TA acc=TA()) {
  ASSERT(a && x);
  // Each thread computes a local scan of its chunk of the input array,
  // and then the local scans are combined into a global scan.
  int H = omp_get_max_threads();
  fillValueU(buf, H, TA());
  #pragma omp parallel
  {
    int T = omp_get_num_threads();
    int t = omp_get_thread_num();
    size_t chunkSize = (N + T - 1) / T;
    size_t i = t * chunkSize;
    size_t I = min(i + chunkSize, N);
    buf[t]   = inclusiveScanW(a+i, x+i, I-i);
  }
  // The global scan is computed on the local scans.
  inclusiveScanW(buf, buf, H);
  // The global scan is added to the local scans.
  #pragma omp parallel
  {
    int T = omp_get_num_threads();
    int t = omp_get_thread_num();
    size_t chunkSize = (N + T - 1) / T;
    size_t i = t * chunkSize;
    size_t I = min(i + chunkSize, N);
    addValueU(a+i, I-i, t==0? acc : buf[t-1] + acc);
  }
  return buf[H-1] + acc;
}
template <class TA, class TX>
inline TA inclusiveScanOmpW(vector<TA>& a, vector<TA>& buf, const vector<TX>& x, TA acc=TA()) {
  return inclusiveScanOmpW(a.data(), buf.data(), x.data(), x.size(), acc);
}
#endif




// EXCLUSIVE SCAN
// --------------

template <class TA, class TX>
inline TA exclusiveScanW(TA *a, const TX *x, size_t N, TA acc=TA()) {
  ASSERT(a && x);
  for (size_t i=0; i<N; ++i) {
    TA t = x[i];
    a[i] = acc;
    acc += t;
  }
  return acc;
}
template <class TA, class TX>
inline TA exclusiveScanW(vector<TA>& a, const vector<TX>& x, TA acc=TA()) {
  return exclusiveScanW(a.data(), x.data(), x.size(), acc);
}


#ifdef OPENMP
template <class TA, class TX>
inline TA exclusiveScanOmpW(TA *a, TA *buf, const TX *x, size_t N, TA acc=TA()) {
  ASSERT(a && x);
  // Each thread computes a local scan of its chunk of the input array,
  // and then the local scans are combined into a global scan.
  int H = omp_get_max_threads();
  fillValueU(buf, H, TA());
  #pragma omp parallel
  {
    int T = omp_get_num_threads();
    int t = omp_get_thread_num();
    size_t chunkSize = (N + T - 1) / T;
    size_t i = t * chunkSize;
    size_t I = min(i + chunkSize, N);
    buf[t]   = exclusiveScanW(a+i, x+i, I-i);
  }
  // The global scan is computed on the local scans.
  inclusiveScanW(buf, buf, H);
  // The global scan is added to the local scans.
  #pragma omp parallel
  {
    int T = omp_get_num_threads();
    int t = omp_get_thread_num();
    size_t chunkSize = (N + T - 1) / T;
    size_t i = t * chunkSize;
    size_t I = min(i + chunkSize, N);
    addValueU(a+i, I-i, t==0? acc : buf[t-1] + acc);
  }
  return buf[H-1] + acc;
}
template <class TA, class TX>
inline TA exclusiveScanOmpW(vector<TA>& a, vector<TA>& buf, const vector<TX>& x, TA acc=TA()) {
  return exclusiveScanOmpW(a.data(), buf.data(), x.data(), x.size(), acc);
}
#endif
