#pragma once
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include "_debug.hxx"
#ifdef OPENMP
#include <omp.h>
#endif

using std::vector;
using std::fill;
using std::abs;
using std::min;
using std::max;




#pragma region TYPES
/**
 * Represents a 2D vector.
 * @tparam T element type
 */
template <class T>
using vector2d = vector<vector<T>>;
#pragma endregion




#pragma region METHODS
#pragma region GATHER VALUES
/**
 * Gather values at specified indices of an array into another array.
 * @param a output array (updated)
 * @param x input array
 * @param is indices
 * @param fm mapping function (value)
 */
template <class TA, class TX, class IS, class FM>
inline void gatherValuesW(TA *a, const TX *x, const IS& is, FM fm) {
  ASSERT(a && x);
  size_t j = 0;
  for (auto i : is)
    a[j++] = TA(fm(x[i]));
}

/**
 * Gather values at specified indices of an array into another array.
 * @param a output array (updated)
 * @param x input array
 * @param is indices
 */
template <class TA, class TX, class IS>
inline void gatherValuesW(TA *a, const TX *x, const IS& is) {
  auto fm = [](const auto& v) { return v; };
  gatherValuesW(a, x, is, fm);
}

/**
 * Gather values at specified indices of a vector into another vector.
 * @param a output vector (updated)
 * @param x input vector
 * @param is indices
 * @param fm mapping function (value)
 */
template <class TA, class TX, class IS, class FM>
inline void gatherValuesW(vector<TA>& a, const vector<TX>& x, const IS& is, FM fm) {
  gatherValuesW(a.data(), x.data(), is, fm);
}

/**
 * Gather values at specified indices of a vector into another vector.
 * @param a output vector (updated)
 * @param x input vector
 * @param is indices
 */
template <class TA, class TX, class IS>
inline void gatherValuesW(vector<TA>& a, const vector<TX>& x, const IS& is) {
  gatherValuesW(a.data(), x.data(), is);
}

/**
 * Gather values at specified indices of a boolean vector into another boolean vector.
 * @param a output boolean vector (updated)
 * @param x input boolean vector
 * @param is indices
 */
template <class IS>
inline void gatherValuesW(vector<bool>& a, const vector<bool>& x, const IS& is) {
  size_t j = 0;
  for (auto i : is)
    a[j++] = x[i];
}


#ifdef OPENMP
/**
 * Gather values at specified indices of an array into another array, in parallel.
 * @param a output array (updated)
 * @param x input array
 * @param is indices
 * @param fm mapping function (value)
 */
template <class TA, class TX, class IS, class FM>
inline void gatherValuesOmpW(TA *a, const TX *x, const IS& is, FM fm) {
  ASSERT(a && x);
  size_t N = is.size();
  #pragma omp parallel for schedule(auto)
  for (size_t j=0; j<N; ++j)
    a[j] = TA(fm(x[is[j]]));
}

/**
 * Gather values at specified indices of an array into another array, in parallel.
 * @param a output array (updated)
 * @param x input array
 * @param is indices
 */
template <class TA, class TX, class IS>
inline void gatherValuesOmpW(TA *a, const TX *x, const IS& is) {
  auto fm = [](const auto& v) { return v; };
  gatherValuesOmpW(a, x, is, fm);
}

/**
 * Gather values at specified indices of a vector into another vector, in parallel.
 * @param a output vector (updated)
 * @param x input vector
 * @param is indices
 * @param fm mapping function (value)
 */
template <class TA, class TX, class IS, class FM>
inline void gatherValuesOmpW(vector<TA>& a, const vector<TX>& x, const IS& is, FM fm) {
  gatherValuesOmpW(a.data(), x.data(), is, fm);
}

/**
 * Gather values at specified indices of a vector into another vector, in parallel.
 * @param a output vector (updated)
 * @param x input vector
 * @param is indices
 */
template <class TA, class TX, class IS>
inline void gatherValuesOmpW(vector<TA>& a, const vector<TX>& x, const IS& is) {
  gatherValuesOmpW(a.data(), x.data(), is);
}
#endif
#pragma endregion




#pragma region SCATTER VALUES
/**
 * Scatter values of an array into another array at specified indices.
 * @param a output array (updated)
 * @param x input array
 * @param is indices
 * @param fm mapping function (value)
 */
template <class TA, class TX, class IS, class FM>
inline void scatterValuesW(TA *a, const TX *x, const IS& is, FM fm) {
  ASSERT(a && x);
  size_t j = 0;
  for (auto i : is)
    a[i] = TA(fm(x[j++]));
}

/**
 * Scatter values of an array into another array at specified indices.
 * @param a output array (updated)
 * @param x input array
 * @param is indices
 */
template <class TA, class TX, class IS>
inline void scatterValuesW(TA *a, const TX *x, const IS& is) {
  auto fm = [](const auto& v) { return v; };
  scatterValuesW(a, x, is, fm);
}

/**
 * Scatter values of a vector into another vector at specified indices.
 * @param a output vector (updated)
 * @param x input vector
 * @param is indices
 * @param fm mapping function (value)
 */
template <class TA, class TX, class IS, class FM>
inline void scatterValuesW(vector<TA>& a, const vector<TX>& x, const IS& is, FM fm) {
  scatterValuesW(a.data(), x.data(), is, fm);
}

/**
 * Scatter values of a vector into another vector at specified indices.
 * @param a output vector (updated)
 * @param x input vector
 * @param is indices
 */
template <class TA, class TX, class IS>
inline void scatterValuesW(vector<TA>& a, const vector<TX>& x, const IS& is) {
  scatterValuesW(a.data(), x.data(), is);
}

/**
 * Scatter values of a boolean vector into another boolean vector at specified indices.
 * @param a output boolean vector (updated)
 * @param x input boolean vector
 * @param is indices
 */
template <class IS>
inline void scatterValuesW(vector<bool>& a, const vector<bool>& x, const IS& is) {
  size_t j = 0;
  for (auto i : is)
    a[i] = x[j++];
}


#ifdef OPENMP
/**
 * Scatter values of an array into another array at specified indices, in parallel.
 * @param a output array (updated)
 * @param x input array
 * @param is indices
 * @param fm mapping function (value)
 */
template <class TA, class TX, class IS, class FM>
inline void scatterValuesOmpW(TA *a, const TX *x, const IS& is, FM fm) {
  ASSERT(a && x);
  size_t N = is.size();
  #pragma omp parallel for schedule(auto)
  for (size_t j=0; j<N; ++j)
    a[is[j]] = TA(fm(x[j]));
}

/**
 * Scatter values of an array into another array at specified indices, in parallel.
 * @param a output array (updated)
 * @param x input array
 * @param is indices
 */
template <class TA, class TX, class IS>
inline void scatterValuesOmpW(TA *a, const TX *x, const IS& is) {
  auto fm = [](const auto& v) { return v; };
  scatterValuesOmpW(a, x, is, fm);
}

/**
 * Scatter values of a vector into another vector at specified indices, in parallel.
 * @param a output vector (updated)
 * @param x input vector
 * @param is indices
 * @param fm mapping function (value)
 */
template <class TA, class TX, class IS, class FM>
inline void scatterValuesOmpW(vector<TA>& a, const vector<TX>& x, const IS& is, FM fm) {
  scatterValuesOmpW(a.data(), x.data(), is, fm);
}

/**
 * Scatter values of a vector into another vector at specified indices, in parallel.
 * @param a output vector (updated)
 * @param x input vector
 * @param is indices
 */
template <class TA, class TX, class IS>
inline void scatterValuesOmpW(vector<TA>& a, const vector<TX>& x, const IS& is) {
  scatterValuesOmpW(a.data(), x.data(), is);
}
#endif
#pragma endregion




#pragma region SCATTER OR
/**
 * Scatter values of an array into another array at specified indices with OR operation.
 * @param a output array (updated)
 * @param x input array
 * @param is indices
 */
template <class TA, class TX, class IS>
inline void scatterOrW(TA *a, const TX *x, const IS& is) {
  ASSERT(a && x);
  size_t j = 0;
  for (auto i : is)
    a[i] |= TA(x[j++]);
}

/**
 * Scatter values of a vector into another vector at specified indices with OR operation.
 * @param a output vector (updated)
 * @param x input vector
 * @param is indices
 */
template <class TA, class TX, class IS>
inline void scatterOrW(vector<TA>& a, const vector<TX>& x, const IS& is) {
  scatterOrW(a.data(), x.data(), is);
}


#ifdef OPENMP
/**
 * Scatter values of an array into another array at specified indices with OR operation, in parallel.
 * @param a output array (updated)
 * @param x input array
 * @param is indices
 */
template <class TA, class TX, class IS>
inline void scatterOrOmpW(TA *a, const TX *x, const IS& is) {
  ASSERT(a && x);
  size_t N = is.size();
  #pragma omp parallel for schedule(auto)
  for (size_t j=0; j<N; ++j)
    a[is[j]] |= TA(x[j]);
}

/**
 * Scatter values of a vector into another vector at specified indices with OR operation, in parallel.
 * @param a output vector (updated)
 * @param x input vector
 * @param is indices
 */
template <class TA, class TX, class IS>
inline void scatterOrOmpW(vector<TA>& a, const vector<TX>& x, const IS& is) {
  scatterOrOmpW(a.data(), x.data(), is);
}
#endif
#pragma endregion




#pragma region VALUE INDICES
/**
 * Get indices of values in an array.
 * @param a output array of indices for each value (updated)
 * @param x input array
 * @param fm mapping function (value)
 */
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
#pragma endregion




#pragma region FILL VALUE
/**
 * Fill an array with a value.
 * @param a output array (a[i] = v, updated)
 * @param N size of array
 * @param v value to fill
 */
template <class T>
inline void fillValueU(T *a, size_t N, const T& v) {
  ASSERT(a);
  fill(a, a+N, v);
}

/**
 * Fill a vector with a value.
 * @param a output vector (a[i] = v, updated)
 * @param v value to fill
 */
template <class T>
inline void fillValueU(vector<T>& a, const T& v) {
  fill(a.begin(), a.end(), v);
}


#ifdef OPENMP
/**
 * Fill an array with a value in parallel.
 * @param a output array (a[i] = v, updated)
 * @param N size of array
 * @param v value to fill
 */
template <class T>
inline void fillValueOmpU(T *a, size_t N, const T& v) {
  ASSERT(a);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a[i] = v;
}

/**
 * Fill a vector with a value in parallel.
 * @param a output vector (a[i] = v, updated)
 * @param v value to fill
 */
template <class T>
inline void fillValueOmpU(vector<T>& a, const T& v) {
  fillValueOmpU(a.data(), a.size(), v);
}
inline void fillValueOmpU(vector<bool>& a, const bool& v) {
  fill(a.begin(), a.end(), v);
}
#endif
#pragma endregion




#pragma region ADD VALUE
/**
 * Add a value to each element of an array.
 * @param a output array (a[i] += v, updated)
 * @param N size of array
 * @param v value to add
 */
template <class T>
inline void addValueU(T *a, size_t N, const T& v) {
  ASSERT(a);
  for (size_t i=0; i<N; ++i)
    a[i] += v;
}

/**
 * Add a value to each element of a vector.
 * @param a output vector (a[i] += v, updated)
 * @param v value to add
 */
template <class T>
inline void addValueU(vector<T>& a, const T& v) {
  addValueU(a.data(), a.size(), v);
}


#ifdef OPENMP
/**
 * Add a value to each element of an array in parallel.
 * @param a output array (a[i] += v, updated)
 * @param N size of array
 * @param v value to add
 */
template <class T>
inline void addValueOmpU(T *a, size_t N, const T& v) {
  ASSERT(a);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a[i] += v;
}

/**
 * Add a value to each element of a vector in parallel.
 * @param a output vector (a[i] += v, updated)
 * @param v value to add
 */
template <class T>
inline void addValueOmpU(vector<T>& a, const T& v) {
  addValueOmpU(a.data(), a.size(), v);
}
#endif
#pragma endregion




#pragma region COPY VALUES
/**
 * Copy values from an array to another output array.
 * @param a output array (a[i] = x[i], updated)
 * @param x input array
 * @param N size of arrays
 */
template <class TA, class TX>
inline void copyValuesW(TA *a, const TX *x, size_t N) {
  ASSERT(a && x);
  for (size_t i=0; i<N; ++i)
    a[i] = x[i];
}

/**
 * Copy values from a vector to another output vector.
 * @param a output vector (a[i] = x[i], updated)
 * @param x input vector
 */
template <class TA, class TX>
inline void copyValuesW(vector<TA>& a, const vector<TX>& x) {
  return copyValuesW(a.data(), x.data(), x.size());
}


#ifdef OPENMP
/**
 * Copy values from an array to another output array in parallel.
 * @param a output array (a[i] = x[i], updated)
 * @param x input array
 * @param N size of arrays
 */
template <class TA, class TX>
inline void copyValuesOmpW(TA *a, const TX *x, size_t N) {
  ASSERT(a && x);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a[i] = x[i];
}

/**
 * Copy values from a vector to another output vector in parallel.
 * @param a output vector (a[i] = x[i], updated)
 * @param x input vector
 */
template <class TA, class TX>
inline void copyValuesOmpW(vector<TA>& a, const vector<TX>& x) {
  return copyValuesOmpW(a.data(), x.data(), x.size());
}
#endif
#pragma endregion




#pragma region MULTIPLY VALUE
/**
 * Multiply a value to each element of an array and store the result in an output array.
 * @param a output array (a[i] = x[i] * v, updated)
 * @param x input array
 * @param v value to multiply
 * @param N size of arrays
 */
template <class TA, class TX, class TV>
inline void multiplyValueW(TA *a, const TX *x, TV v, size_t N) {
  ASSERT(a && x);
  for (size_t i=0; i<N; ++i)
    a[i] = TA(x[i] * v);
}

/**
 * Multiply a value to each element of a vector and store the result in an output vector.
 * @param a output vector (a[i] = x[i] * v, updated)
 * @param x input vector
 * @param v value to multiply
 */
template <class TA, class TX, class TV>
inline void multiplyValueW(vector<TA>& a, const vector<TX>& x, TV v) {
  multiplyValueW(a.data(), x.data(), v, x.size());
}


#ifdef OPENMP
/**
 * Multiply a value to each element of an array and store the result in an output array, in parallel.
 * @param a output array (a[i] = x[i] * v, updated)
 * @param x input array
 * @param v value to multiply
 * @param N size of arrays
 */
template <class TA, class TX, class TV>
inline void multiplyValueOmpW(TA *a, const TX *x, TV v, size_t N) {
  ASSERT(a && x);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a[i] = TA(x[i] * v);
}

/**
 * Multiply a value to each element of a vector and store the result in an output vector, in parallel.
 * @param a output vector (a[i] = x[i] * v, updated)
 * @param x input vector
 * @param v value to multiply
 */
template <class TA, class TX, class TV>
inline void multiplyValueOmpW(vector<TA>& a, const vector<TX>& x, TV v) {
  multiplyValueOmpW(a.data(), x.data(), v, x.size());
}
#endif
#pragma endregion




#pragma region MULTIPLY VALUES
/**
 * Multiply two arrays element-wise and store the result in an output array.
 * @param a output array (a[i] = x[i] * y[i], updated)
 * @param x first array
 * @param y second array
 * @param N size of arrays
 */
template <class TA, class TX, class TY>
inline void multiplyValuesW(TA *a, const TX *x, const TY *y, size_t N) {
  ASSERT(a && x && y);
  for (size_t i=0; i<N; ++i)
    a[i] = TA(x[i] * y[i]);
}

/**
 * Multiply two vectors element-wise and store the result in an output vector.
 * @param a output vector (a[i] = x[i] * y[i], updated)
 * @param x first vector
 * @param y second vector
 */
template <class TA, class TX, class TY>
inline void multiplyValuesW(vector<TA>& a, const vector<TX>& x, const vector<TY>& y) {
  multiplyValuesW(a.data(), x.data(), y.data(), x.size());
}


#ifdef OPENMP
/**
 * Multiply two arrays element-wise and store the result in an output array, in parallel.
 * @param a output array (a[i] = x[i] * y[i], updated)
 * @param x first array
 * @param y second array
 * @param N size of arrays
 */
template <class TA, class TX, class TY>
inline void multiplyValuesOmpW(TA *a, const TX *x, const TY *y, size_t N) {
  ASSERT(a && x && y);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a[i] = TA(x[i] * y[i]);
}

/**
 * Multiply two vectors element-wise and store the result in an output vector, in parallel.
 * @param a output vector (a[i] = x[i] * y[i], updated)
 * @param x first vector
 * @param y second vector
 */
template <class TA, class TX, class TY>
inline void multiplyValuesOmpW(vector<TA>& a, const vector<TX>& x, const vector<TY>& y) {
  multiplyValuesOmpW(a.data(), x.data(), y.data(), x.size());
}
#endif
#pragma endregion




#pragma region SUM VALUES
/**
 * Compute the sum of values in an array.
 * @param x an array
 * @param N size of array
 * @param a initial value
 * @returns sum of values
 */
template <class TX, class TA=TX>
inline TA sumValues(const TX *x, size_t N, TA a=TA()) {
  ASSERT(x);
  for (size_t i=0; i<N; ++i)
    a += TA(x[i]);
  return a;
}

/**
 * Compute the sum of values in a vector.
 * @param x a vector
 * @param a initial value
 * @returns sum of values
 */
template <class TX, class TA=TX>
inline TA sumValues(const vector<TX>& x, TA a=TA()) {
  return sumValues(x.data(), x.size(), a);
}


#ifdef OPENMP
/**
 * Compute the sum of values in an array in parallel.
 * @param x an array
 * @param N size of array
 * @param a initial value
 * @returns sum of values
 */
template <class TX, class TA=TX>
inline TA sumValuesOmp(const TX *x, size_t N, TA a=TA()) {
  ASSERT(x);
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (size_t i=0; i<N; ++i)
    a += TA(x[i]);
  return a;
}

/**
 * Compute the sum of values in a vector in parallel.
 * @param x a vector
 * @param a initial value
 * @returns sum of values
 */
template <class TX, class TA=TX>
inline TA sumValuesOmp(const vector<TX>& x, TA a=TA()) {
  return sumValuesOmp(x.data(), x.size(), a);
}
#endif
#pragma endregion




#pragma region COUNT VALUE
/**
 * Count the number of times a value appears in an array.
 * @param x an array
 * @param N size of array
 * @param v value to count
 * @returns number of times value appears
 */
template <class TX>
inline size_t countValue(const TX *x, size_t N, const TX& v) {
  ASSERT(x);
  size_t a = 0;
  for (size_t i=0; i<N; ++i)
    if (x[i]==v) ++a;
  return a;
}

/**
 * Count the number of times a value appears in an array.
 * @param x a vector
 * @param v value to count
 * @returns number of times value appears
 */
template <class TX>
inline size_t countValue(const vector<TX>& x, const TX& v) {
  return countValue(x.data(), x.size(), v);
}


#ifdef OPENMP
/**
 * Count the number of times a value appears in an array.
 * @param x an array
 * @param N size of array
 * @param v value to count
 * @returns number of times value appears
 */
template <class TX>
inline size_t countValueOmp(const TX *x, size_t N, const TX& v) {
  ASSERT(x);
  size_t a = 0;
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (size_t i=0; i<N; ++i)
    if (x[i]==v) ++a;
  return a;
}

/**
 * Count the number of times a value appears in an array.
 * @param x a vector
 * @param v value to count
 * @returns number of times value appears
 */
template <class TX>
inline size_t countValueOmp(const vector<TX>& x, const TX& v) {
  return countValueOmp(x.data(), x.size(), v);
}
#endif
#pragma endregion




#pragma region L1-NORM
/**
 * Compute the L1-norm of an array.
 * @param x an array
 * @param N size of array
 * @param a initial value
 * @returns ||x||_1
 */
template <class TX, class TA=TX>
inline TA l1Norm(const TX *x, size_t N, TA a=TA()) {
  ASSERT(x);
  for (size_t i=0; i<N; ++i)
    a += TA(abs(x[i]));
  return a;
}

/**
 * Compute the L1-norm of a vector.
 * @param x a vector
 * @param a initial value
 * @returns ||x||_1
 */
template <class TX, class TA=TX>
inline TA l1Norm(const vector<TX>& x, TA a=TA()) {
  return l1Norm(x.data(), x.size(), a);
}


#ifdef OPENMP
/**
 * Compute the L1-norm of an array in parallel.
 * @param x an array
 * @param N size of array
 * @param a initial value
 * @returns ||x||_1
 */
template <class TX, class TA=TX>
inline TA l1NormOmp(const TX *x, size_t N, TA a=TA()) {
  ASSERT(x);
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (size_t i=0; i<N; ++i)
    a += TA(abs(x[i]));
  return a;
}

/**
 * Compute the L1-norm of a vector in parallel.
 * @param x a vector
 * @param a initial value
 * @returns ||x||_1
 */
template <class TX, class TA=TX>
inline TA l1NormOmp(const vector<TX>& x, TA a=TA()) {
  return l1NormOmp(x.data(), x.size(), a);
}
#endif
#pragma endregion




#pragma region L1-NORM DELTA
/**
 * Compute the L1-norm of the difference of two arrays.
 * @param x an array
 * @param y another array
 * @param N size of arrays
 * @param a initial value
 * @returns ||x-y||_1
 */
template <class TX, class TY, class TA=TX>
inline TA l1NormDelta(const TX *x, const TY *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  for (size_t i=0; i<N; ++i)
    a += TA(abs(x[i] - y[i]));
  return a;
}

/**
 * Compute the L1-norm of the difference of two vectors.
 * @param x a vector
 * @param y another vector
 * @param a initial value
 * @returns ||x-y||_1
 */
template <class TX, class TY, class TA=TX>
inline TA l1NormDelta(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return l1NormDelta(x.data(), y.data(), x.size(), a);
}


#ifdef OPENMP
/**
 * Compute the L1-norm of the difference of two arrays in parallel.
 * @param x an array
 * @param y another array
 * @param N size of arrays
 * @param a initial value
 * @returns ||x-y||_1
 */
template <class TX, class TY, class TA=TX>
inline TA l1NormDeltaOmp(const TX *x, const TY *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (size_t i=0; i<N; ++i)
    a += TA(abs(x[i] - y[i]));
  return a;
}

/**
 * Compute the L1-norm of the difference of two vectors in parallel.
 * @param x a vector
 * @param y another vector
 * @param a initial value
 * @returns ||x-y||_1
 */
template <class TX, class TY, class TA=TX>
inline TA l1NormDeltaOmp(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return l1NormDeltaOmp(x.data(), y.data(), x.size(), a);
}
#endif
#pragma endregion




#pragma region L1-NORM DELTA AT
/**
 * Compute the L1-norm of the difference of two arrays at given indices.
 * @param x an array
 * @param y another array
 * @param is indices
 * @param IS number of indices
 * @param a initial value
 * @returns ||x-y||_1
 */
template <class TX, class TY, class TI, class TA=TX>
inline TA l1NormDeltaAt(const TX *x, const TY *y, const TI *is, size_t IS, TA a=TA()) {
  ASSERT(x && y && is);
  for (size_t l=0; l<IS; ++l) {
    TI i = is[l];
    a += TA(abs(x[i] - y[i]));
  }
  return a;
}

/**
 * Compute the L1-norm of the difference of two vectors at given indices.
 * @param x a vector
 * @param y another vector
 * @param is indices
 * @param a initial value
 * @returns ||x-y||_1
 */
template <class TX, class TY, class TI, class TA=TX>
inline TA l1NormDeltaAt(const vector<TX>& x, const vector<TY>& y, const vector<TI>& is, TA a=TA()) {
  return l1NormDeltaAt(x.data(), y.data(), is.data(), is.size(), a);
}


#ifdef OPENMP
/**
 * Compute the L1-norm of the difference of two arrays at given indices in parallel.
 * @param x an array
 * @param y another array
 * @param is indices
 * @param IS number of indices
 * @param a initial value
 * @returns ||x-y||_1
 */
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

/**
 * Compute the L1-norm of the difference of two vectors at given indices in parallel.
 * @param x a vector
 * @param y another vector
 * @param is indices
 * @param a initial value
 * @returns ||x-y||_1
 */
template <class TX, class TY, class TI, class TA=TX>
inline TA l1NormDeltaAtOmp(const vector<TX>& x, const vector<TY>& y, const vector<TI>& is, TA a=TA()) {
  return l1NormDeltaAtOmp(x.data(), y.data(), is.data(), is.size(), a);
}
#endif
#pragma endregion




#pragma region L2-NORM
/**
 * Compute the L2-norm of an array.
 * @param x an array
 * @param N size of array
 * @param a initial value
 * @returns ||x||_2
 */
template <class TX, class TA=TX>
inline TA l2Norm(const TX *x, size_t N, TA a=TA()) {
  ASSERT(x);
  for (size_t i=0; i<N; ++i)
    a += TA(x[i]) * TA(x[i]);
  return a;
}

/**
 * Compute the L2-norm of a vector.
 * @param x a vector
 * @param a initial value
 * @returns ||x||_2
 */
template <class TX, class TA=TX>
inline TA l2Norm(const vector<TX>& x, TA a=TA()) {
  return l2Norm(x.data(), x.size(), a);
}


#ifdef OPENMP
/**
 * Compute the L2-norm of an array in parallel.
 * @param x an array
 * @param N size of array
 * @param a initial value
 * @returns ||x||_2
 */
template <class TX, class TA=TX>
inline TA l2NormOmp(const TX *x, size_t N, TA a=TA()) {
  ASSERT(x);
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (size_t i=0; i<N; ++i)
    a += TA(x[i]) * TA(x[i]);
  return a;
}

/**
 * Compute the L2-norm of a vector in parallel.
 * @param x a vector
 * @param a initial value
 * @returns ||x||_2
 */
template <class TX, class TA=TX>
inline TA l2NormOmp(const vector<TX>& x, TA a=TA()) {
  return l2NormOmp(x.data(), x.size(), a);
}
#endif
#pragma endregion




#pragma region L2-NORM DELTA
/**
 * Compute the L2-norm of the difference of two arrays.
 * @param x an array
 * @param y another array
 * @param N size of arrays
 * @param a initial value
 * @returns ||x-y||_2
 */
template <class TX, class TY, class TA=TX>
inline TA l2NormDelta(const TX *x, const TY *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  for (size_t i=0; i<N; ++i)
    a += TA(x[i] - y[i]) * TA(x[i] - y[i]);
  return a;
}

/**
 * Compute the L2-norm of the difference of two vectors.
 * @param x a vector
 * @param y another vector
 * @param a initial value
 * @returns ||x-y||_2
 */
template <class TX, class TY, class TA=TX>
inline TA l2NormDelta(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return l2NormDelta(x.data(), y.data(), x.size(), a);
}


#ifdef OPENMP
/**
 * Compute the L2-norm of the difference of two arrays in parallel.
 * @param x an array
 * @param y another array
 * @param N size of arrays
 * @param a initial value
 * @returns ||x-y||_2
 */
template <class TX, class TY, class TA=TX>
inline TA l2NormDeltaOmp(const TX *x, const TY *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (size_t i=0; i<N; ++i)
    a += TA(x[i] - y[i]) * TA(x[i] - y[i]);
  return a;
}

/**
 * Compute the L2-norm of the difference of two vectors in parallel.
 * @param x a vector
 * @param y another vector
 * @param a initial value
 * @returns ||x-y||_2
 */
template <class TX, class TY, class TA=TX>
inline TA l2NormDeltaOmp(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return l2NormDeltaOmp(x.data(), y.data(), x.size(), a);
}
#endif
#pragma endregion




#pragma region L2-NORM DELTA AT
/**
 * Compute the L2-norm of the difference of two arrays at given indices.
 * @param x an array
 * @param y another array
 * @param is indices
 * @param IS number of indices
 * @param a initial value
 * @returns ||x-y||_2
 */
template <class TX, class TY, class TI, class TA=TX>
inline TA l2NormDeltaAt(const TX *x, const TY *y, const TI *is, size_t IS, TA a=TA()) {
  ASSERT(x && y && is);
  for (size_t l=0; l<IS; ++l) {
    TI i = is[l];
    a += TA(x[i] - y[i]) * TA(x[i] - y[i]);
  }
  return a;
}

/**
 * Compute the L2-norm of the difference of two vectors at given indices.
 * @param x a vector
 * @param y another vector
 * @param is indices
 * @param a initial value
 * @returns ||x-y||_2
 */
template <class TX, class TY, class TI, class TA=TX>
inline TA l2NormDeltaAt(const vector<TX>& x, const vector<TY>& y, const vector<TI>& is, TA a=TA()) {
  return l2NormDeltaAt(x.data(), y.data(), is.data(), is.size(), a);
}


#ifdef OPENMP
/**
 * Compute the L2-norm of the difference of two arrays at given indices in parallel.
 * @param x an array
 * @param y another array
 * @param is indices
 * @param IS number of indices
 * @param a initial value
 * @returns ||x-y||_2
 */
template <class TX, class TY, class TI, class TA=TX>
inline TA l2NormDeltaAtOmp(const TX *x, const TY *y, const TI *is, size_t IS, TA a=TA()) {
  ASSERT(x && y && is);
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (size_t l=0; l<IS; ++l) {
    TI i = is[l];
    a += TA(x[i] - y[i]) * TA(x[i] - y[i]);
  }
  return a;
}

/**
 * Compute the L2-norm of the difference of two vectors at given indices in parallel.
 * @param x a vector
 * @param y another vector
 * @param is indices
 * @param a initial value
 * @returns ||x-y||_2
 */
template <class TX, class TY, class TI, class TA=TX>
inline TA l2NormDeltaAtOmp(const vector<TX>& x, const vector<TY>& y, const vector<TI>& is, TA a=TA()) {
  return l2NormDeltaAtOmp(x.data(), y.data(), is.data(), is.size(), a);
}
#endif
#pragma endregion




#pragma region LI-NORM
/**
 * Compute the L∞-norm of an array.
 * @param x an array
 * @param N size of array
 * @param a initial value
 * @returns ||x||_inf
 */
template <class TX, class TA=TX>
inline TA liNorm(const TX *x, size_t N, TA a=TA()) {
  ASSERT(x);
  for (size_t i=0; i<N; ++i)
    a = max(a, TA(abs(x[i])));
  return a;
}

/**
 * Compute the L∞-norm of a vector.
 * @param x a vector
 * @param a initial value
 * @returns ||x||_inf
 */
template <class TX, class TA=TX>
inline TA liNorm(const vector<TX>& x, TA a=TA()) {
  return liNorm(x.data(), x.size(), a);
}


#ifdef OPENMP
/**
 * Compute the L∞-norm of an array in parallel.
 * @param x an array
 * @param N size of array
 * @param a initial value
 * @returns ||x||_inf
 */
template <class TX, class TA=TX>
inline TA liNormOmp(const TX *x, size_t N, TA a=TA()) {
  ASSERT(x);
  #pragma omp parallel for schedule(auto) reduction(max:a)
  for (size_t i=0; i<N; ++i)
    a = max(a, TA(abs(x[i])));
  return a;
}

/**
 * Compute the L∞-norm of a vector in parallel.
 * @param x a vector
 * @param a initial value
 * @returns ||x||_inf
 */
template <class TX, class TA=TX>
inline TA liNormOmp(const vector<TX>& x, TA a=TA()) {
  return liNormOmp(x.data(), x.size(), a);
}
#endif
#pragma endregion




#pragma region LI-NORM DELTA
/**
 * Compute the L∞-norm of the difference of two arrays.
 * @param x an array
 * @param y another array
 * @param N size of arrays
 * @param a initial value
 * @returns ||x-y||_inf
 */
template <class TX, class TY, class TA=TX>
inline TA liNormDelta(const TX *x, const TY *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  for (size_t i=0; i<N; ++i)
    a = max(a, TA(abs(x[i] - y[i])));
  return a;
}

/**
 * Compute the L∞-norm of the difference of two vectors.
 * @param x a vector
 * @param y another vector
 * @param a initial value
 * @returns ||x-y||_inf
 */
template <class TX, class TY, class TA=TX>
inline TA liNormDelta(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return liNormDelta(x.data(), y.data(), x.size(), a);
}


#ifdef OPENMP
/**
 * Compute the L∞-norm of the difference of two arrays in parallel.
 * @param x an array
 * @param y another array
 * @param N size of arrays
 * @param a initial value
 * @returns ||x-y||_inf
 */
template <class TX, class TY, class TA=TX>
inline TA liNormDeltaOmp(const TX *x, const TY *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  #pragma omp parallel for schedule(auto) reduction(max:a)
  for (size_t i=0; i<N; ++i)
    a = max(a, TA(abs(x[i] - y[i])));
  return a;
}

/**
 * Compute the L∞-norm of the difference of two vectors in parallel.
 * @param x a vector
 * @param y another vector
 * @param a initial value
 * @returns ||x-y||_inf
 */
template <class TX, class TY, class TA=TX>
inline TA liNormDeltaOmp(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return liNormDeltaOmp(x.data(), y.data(), x.size(), a);
}
#endif
#pragma endregion




#pragma region LI-NORM DELTA AT
/**
 * Compute the L∞-norm of the difference of two arrays at given indices.
 * @param x an array
 * @param y another array
 * @param is indices
 * @param IS number of indices
 * @param a initial value
 * @returns ||x[..is] - y[..is]||_inf
 */
template <class TX, class TY, class TI, class TA=TX>
inline TA liNormDeltaAt(const TX *x, const TY *y, const TI *is, size_t IS, TA a=TA()) {
  ASSERT(x && y && is);
  for (size_t l=0; l<IS; ++l) {
    TI i = is[l];
    a = max(a, TA(abs(x[i] - y[i])));
  }
  return a;
}

/**
 * Compute the L∞-norm of the difference of two vectors at given indices.
 * @param x a vector
 * @param y another vector
 * @param is indices
 * @param a initial value
 * @returns ||x[..is] - y[..is]||_inf
 */
template <class TX, class TY, class TI, class TA=TX>
inline TA liNormDeltaAt(const vector<TX>& x, const vector<TY>& y, const vector<TI>& is, TA a=TA()) {
  return liNormDeltaAt(x.data(), y.data(), is.data(), is.size(), a);
}


#ifdef OPENMP
/**
 * Compute the L∞-norm of the difference of two arrays at given indices in parallel.
 * @param x an array
 * @param y another array
 * @param is indices
 * @param IS number of indices
 * @param a initial value
 * @returns ||x[..is] - y[..is]||_inf
 */
template <class TX, class TY, class TI, class TA=TX>
inline TA liNormDeltaAtOmp(const TX *x, const TY *y, const TI *is, size_t IS, TA a=TA()) {
  ASSERT(x && y && is);
  #pragma omp parallel for schedule(auto) reduction(max:a)
  for (size_t l=0; l<IS; ++l) {
    TI i = is[l];
    a = max(a, TA(abs(x[i] - y[i])));
  }
  return a;
}

/**
 * Compute the L∞-norm of the difference of two vectors at given indices in parallel.
 * @param x a vector
 * @param y another vector
 * @param is indices
 * @param a initial value
 * @returns ||x[..is] - y[..is]||_inf
 */
template <class TX, class TY, class TI, class TA=TX>
inline TA liNormDeltaAtOmp(const vector<TX>& x, const vector<TY>& y, const vector<TI>& is, TA a=TA()) {
  return liNormDeltaAtOmp(x.data(), y.data(), is.data(), is.size(), a);
}
#endif
#pragma endregion




#pragma region INCLUSIVE SCAN
/**
 * Perform inclusive scan of an array into another array.
 * @param a output array (updated)
 * @param x input array
 * @param N size of arrays
 * @param acc initial value
 * @returns final value
 */
template <class TA, class TX>
inline TA inclusiveScanW(TA *a, const TX *x, size_t N, TA acc=TA()) {
  ASSERT(a && x);
  for (size_t i=0; i<N; ++i) {
    acc += x[i];
    a[i] = acc;
  }
  return acc;
}

/**
 * Perform inclusive scan of a vector into another vector.
 * @param a output vector (updated)
 * @param x input vector
 * @param acc initial value
 * @returns final value
 */
template <class TA, class TX>
inline TA inclusiveScanW(vector<TA>& a, const vector<TX>& x, TA acc=TA()) {
  return inclusiveScanW(a.data(), x.data(), x.size(), acc);
}


#ifdef OPENMP
/**
 * Perform inclusive scan of an array into another array in parallel.
 * @param a output array (updated)
 * @param x input array
 * @param N size of arrays
 * @param acc initial value
 * @returns final value
 */
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
    size_t i = min(t * chunkSize, N);
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
    size_t i = min(t * chunkSize, N);
    size_t I = min(i + chunkSize, N);
    addValueU(a+i, I-i, t==0? acc : buf[t-1] + acc);
  }
  return buf[H-1] + acc;
}

/**
 * Perform inclusive scan of a vector into another vector in parallel.
 * @param a output vector (updated)
 * @param x input vector
 * @param acc initial value
 * @returns final value
 */
template <class TA, class TX>
inline TA inclusiveScanOmpW(vector<TA>& a, vector<TA>& buf, const vector<TX>& x, TA acc=TA()) {
  return inclusiveScanOmpW(a.data(), buf.data(), x.data(), x.size(), acc);
}
#endif
#pragma endregion




#pragma region EXCLUSIVE SCAN
/**
 * Perform exclusive scan of an array into another array.
 * @param a output array (updated)
 * @param x input array
 * @param N size of arrays
 * @param acc initial value
 * @returns final value
 */
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

/**
 * Perform exclusive scan of a vector into another vector.
 * @param a output vector (updated)
 * @param x input vector
 * @param acc initial value
 * @returns final value
 */
template <class TA, class TX>
inline TA exclusiveScanW(vector<TA>& a, const vector<TX>& x, TA acc=TA()) {
  return exclusiveScanW(a.data(), x.data(), x.size(), acc);
}


#ifdef OPENMP
/**
 * Perform exclusive scan of an array into another array in parallel.
 * @param a output array (updated)
 * @param x input array
 * @param N size of arrays
 * @param acc initial value
 * @returns final value
 */
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
    size_t i = min(t * chunkSize, N);
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
    size_t i = min(t * chunkSize, N);
    size_t I = min(i + chunkSize, N);
    addValueU(a+i, I-i, t==0? acc : buf[t-1] + acc);
  }
  return buf[H-1] + acc;
}

/**
 * Perform exclusive scan of a vector into another vector in parallel.
 * @param a output vector (updated)
 * @param x input vector
 * @param acc initial value
 * @returns final value
 */
template <class TA, class TX>
inline TA exclusiveScanOmpW(vector<TA>& a, vector<TA>& buf, const vector<TX>& x, TA acc=TA()) {
  return exclusiveScanOmpW(a.data(), buf.data(), x.data(), x.size(), acc);
}
#endif
#pragma endregion
#pragma endregion
