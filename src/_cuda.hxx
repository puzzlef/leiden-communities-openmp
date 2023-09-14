#pragma once
#include <limits>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>
#include "_debug.hxx"
#include "_cmath.hxx"

using std::numeric_limits;
using std::vector;
using std::min;
using std::max;
using std::abs;
using std::ceil;
using std::fprintf;
using std::exit;




#pragma region TYPES
/** 64-bit signed integer (CUDA specific). */
typedef long long int          int64_cu;
/** 64-bit unsigned integer (CUDA specific). */
typedef unsigned long long int uint64_cu;
// - https://stackoverflow.com/a/32862733/1413259
#pragma endregion




#pragma region KEYWORDS
#ifndef __global__
/** CUDA kernel function. */
#define __global__

/** CUDA host function. */
#define __host__

/** CUDA device function. */
#define __device__

/** CUDA shared memory. */
#define __shared__

/**
 * Synchronize all threads in a block.
 */
#define __syncthreads()
#endif
#pragma endregion




#pragma region LAUNCH CONFIG
#ifndef BLOCK_LIMIT_CUDA
/** Maximum number of threads per block. */
#define BLOCK_LIMIT_CUDA         1024
/** Maximum number of threads per block, when using a map-like kernel. */
#define BLOCK_LIMIT_MAP_CUDA     256
/** Maximum number of threads per block, when using a reduce-like kernel. */
#define BLOCK_LIMIT_REDUCE_CUDA  256
#endif


#ifndef GRID_LIMIT_CUDA
/** Maximum number of blocks per grid. */
#define GRID_LIMIT_CUDA          2147483647  // 2^31 - 1
/** Maximum number of blocks per grid, when using a map-like kernel. */
#define GRID_LIMIT_MAP_CUDA      65535
/** Maximum number of blocks per grid, when using a reduce-like kernel. */
#define GRID_LIMIT_REDUCE_CUDA   1024
#endif


/**
 * Get the block size for kernel launch, based on number of elements to process.
 * @tparam COARSE an element-per-block kernel?
 * @param N number of elements to process
 * @param BLIM block limit
 * @returns block size
 */
template <bool COARSE=false>
inline int blockSizeCu(size_t N, int BLIM=BLOCK_LIMIT_CUDA) noexcept {
  return COARSE? BLIM : int(min(N, size_t(BLIM)));
}


/**
 * Get the grid size for kernel launch, based on number of elements to process.
 * @tparam COARSE an element-per-block kernel?
 * @param N number of elements to process
 * @param B block size
 * @param GLIM grid limit
 * @returns grid size
 */
template <bool COARSE=false>
inline int gridSizeCu(size_t N, int B, int GLIM=GRID_LIMIT_CUDA) noexcept {
  return COARSE? int(min(N, size_t(GLIM))) : int(min(ceilDiv(N, size_t(B)), size_t(GLIM)));
}


/**
 * Get the number of elements produced by a reduce-like kernel.
 * @tparam COARSE an element-per-block kernel?
 * @param N number of elements to process
 * @returns number of reduced elements
 */
template <bool COARSE=false>
inline int reduceSizeCu(size_t N) noexcept {
  const int B = blockSizeCu<COARSE>(N,   BLOCK_LIMIT_REDUCE_CUDA);
  const int G = gridSizeCu <COARSE>(N, B, GRID_LIMIT_REDUCE_CUDA);
  return G;
}
#pragma endregion




#pragma region TRY
#ifndef TRY_CUDA
/**
 * Log error on CUDA function call failure.
 * @param err error code
 * @param exp expression string
 * @param func current function name
 * @param line current line number
 * @param file current file name
 */
void tryFailedCuda(cudaError err, const char* exp, const char* func, int line, const char* file) {
  if (err == cudaSuccess) return;
  fprintf(stderr,
    "ERROR: %s: %s\n"
    "  in expression %s\n"
    "  at %s:%d in %s\n",
    cudaGetErrorName(err), cudaGetErrorString(err), exp, func, line, file);
  exit(err);
}

/**
 * Try to execute a CUDA function call.
 * @param exp expression to execute
 */
#define TRY_CUDA(exp)  do { cudaError err = exp; if (err != cudaSuccess) tryFailedCuda(err, #exp, __func__, __LINE__, __FILE__); } while (0)

/**
 * Try to execute a CUDA function call only if build mode is error or higher.
 * @param exp expression to execute
 **/
#define TRY_CUDAE(exp)  PERFORME(TRY_CUDA(exp))

/**
 * Try to execute a CUDA function call only if build mode is warning or higher.
 * @param exp expression to execute
 **/
#define TRY_CUDAW(exp)  PERFORMW(TRY_CUDA(exp))

/**
 * Try to execute a CUDA function call only if build mode is info or higher.
 * @param exp expression to execute
 **/
#define TRY_CUDAI(exp)  PERFORMI(TRY_CUDA(exp))

/**
 * Try to execute a CUDA function call only if build mode is debug or higher.
 * @param exp expression to execute
 **/
#define TRY_CUDAD(exp)  PERFORMD(TRY_CUDA(exp))

/**
 * Try to execute a CUDA function call only if build mode is trace.
 * @param exp expression to execute
 **/
#define TRY_CUDAT(exp)  PERFORMT(TRY_CUDA(exp))
#endif
#pragma endregion




#pragma region UNUSED
/**
 * Mark CUDA variable as unused.
 */
template <class T>
inline void __device__ unusedCuda(T&&) {}


#ifndef UNUSED_CUDA
/**
 * Mark CUDA variable as unused.
 * @param x variable to mark as unused
 */
#define UNUSED_CUDA(x)  unusedCuda(x)
#endif
#pragma endregion




#pragma region DEFINE
#ifndef DEFINE_CUDA
/**
 * Define thread, block variables for CUDA.
 * @param t thread index
 * @param b block index
 * @param B block size
 * @param G grid size
 */
#define DEFINE_CUDA(t, b, B, G) \
  const int t = threadIdx.x; \
  const int b = blockIdx.x; \
  const int B = blockDim.x; \
  const int G = gridDim.x; \
  UNUSED_CUDA(t); \
  UNUSED_CUDA(b); \
  UNUSED_CUDA(B); \
  UNUSED_CUDA(G)


/**
 * Define 2D thread, block variables for CUDA.
 * @param tx thread x index
 * @param ty thread y index
 * @param bx block x index
 * @param by block y index
 * @param BX block x size
 * @param BY block y size
 * @param GX grid x size
 * @param GY grid y size
 */
#define DEFINE2D_CUDA(tx, ty, bx, by, BX, BY, GX, GY) \
  const int tx = threadIdx.x; \
  const int ty = threadIdx.y; \
  const int bx = blockIdx.x; \
  const int by = blockIdx.y; \
  const int BX = blockDim.x; \
  const int BY = blockDim.y; \
  const int GX = gridDim.x;  \
  const int GY = gridDim.y; \
  UNUSED_CUDA(tx); \
  UNUSED_CUDA(ty); \
  UNUSED_CUDA(bx); \
  UNUSED_CUDA(by); \
  UNUSED_CUDA(BX); \
  UNUSED_CUDA(BY); \
  UNUSED_CUDA(GX); \
  UNUSED_CUDA(GY)
#endif
#pragma endregion




#pragma region METHODS
#pragma region READ
/**
 * Read a value from global memory.
 * @param v address of value
 * @returns value
 */
template <class T>
inline T readValueCu(const T *v) {
  ASSERT(v);
  T vH;
  TRY_CUDA( cudaMemcpy(&vH, v, sizeof(T), cudaMemcpyDeviceToHost) );
  return vH;
}


/**
 * Read values from global memory.
 * @param v address of values
 * @param N number of values
 * @returns values as vector
 */
template <class T>
inline vector<T> readValuesCu(const T *v, size_t N) {
  ASSERT(v);
  vector<T> vH(N);
  TRY_CUDA( cudaMemcpy(vH.data(), v, N * sizeof(T), cudaMemcpyDeviceToHost) );
  return vH;
}
#pragma endregion




#pragma region SWAP
/**
 * Swap two values in device memory [device function].
 * @param a first value (updated)
 * @param b second value (updated)
 */
template <class T>
inline void __device__ swapCudU(T& a, T& b) {
  const T t = a;
  a = b;
  b = t;
}
#pragma endregion




#pragma region CEIL DIV
/**
 * Get the ceiling of a division [device function].
 * @param x the dividend
 * @param y the divisor
 * @returns ceil(x/y)
 */
template <class T>
inline T __device__ ceilDivCud(T x, T y) {
  ASSERT(y);
  return (x + y-1) / y;
}

/**
 * Get the ceiling of a division [device function].
 * @param x the dividend
 * @param y the divisor
 * @returns ceil(x/y)
 */
template <>
inline float __device__ ceilDivCud<float>(float x, float y) {
  return ceil(x/y);
}

/**
 * Get the ceiling of a division [device function].
 * @param x the dividend
 * @param y the divisor
 * @returns ceil(x/y)
 */
template <>
inline double __device__ ceilDivCud<double>(double x, double y) {
  return ceil(x/y);
}
#pragma endregion




#pragma region POW2
/**
 * Get the next power of 2 of a value [device function].
 * @param x the value
 * @returns next power of 2 of x
 */
template <class T>
inline T __device__ nextPow2Cud(T x) {
  return T(1) << (__clz(T()) - __clz(x));
}
#pragma endregion




#pragma region COPY
/**
 * Copy values from one array to another [device function].
 * @param a destination array (output)
 * @param x source array
 * @param N size of source array
 * @param i start index
 * @param DI index stride
 */
template <class T>
inline void __device__ copyValuesCudW(T *a, const T *x, size_t N, size_t i, size_t DI) {
  ASSERT(a && x && DI);
  for (; i<N; i+=DI)
    a[i] = x[i];
}


/**
 * Copy values from one array to another [kernel].
 * @param a destination array (output)
 * @param x source array
 * @param N size of source array
 */
template <class T>
void __global__ copyValuesCukW(T *a, const T *x, size_t N) {
  ASSERT(a && x);
  DEFINE_CUDA(t, b, B, G);
  copyValuesCudW(a, x, N, B*b+t, G*B);
}


/**
 * Copy values from one array to another.
 * @param a destination array (output)
 * @param x source array
 * @param N size of source array
 */
template <class T>
inline void copyValuesCuW(T *a, const T *x, size_t N) {
  ASSERT(a && x);
  const int B = blockSizeCu(N,   BLOCK_LIMIT_MAP_CUDA);
  const int G = gridSizeCu (N, B, GRID_LIMIT_MAP_CUDA);
  copyValuesCukW<<<G, B>>>(a, x, N);
}
#pragma endregion




#pragma region FILL
/**
 * Fill array with a value [device function].
 * @param a array to fill (output)
 * @param N size of array
 * @param v value to fill with
 * @param i start index
 * @param DI index stride
 */
template <class T>
inline void __device__ fillValueCudW(T *a, size_t N, T v, size_t i, size_t DI) {
  ASSERT(a && DI);
  for (; i<N; i+=DI)
    a[i] = v;
}


/**
 * Fill array with a value [kernel].
 * @param a array to fill (output)
 * @param N size of array
 * @param v value to fill with
 */
template <class T>
void __global__ fillValueCukW(T *a, size_t N, T v) {
  ASSERT(a);
  DEFINE_CUDA(t, b, B, G);
  fillValueCudW(a, N, v, B*b+t, G*B);
}


/**
 * Fill array with a value.
 * @param a array to fill (output)
 * @param N size of array
 * @param v value to fill with
 */
template <class T>
inline void fillValueCuW(T *a, size_t N, T v) {
  ASSERT(a);
  const int B = blockSizeCu(N,   BLOCK_LIMIT_MAP_CUDA);
  const int G = gridSizeCu (N, B, GRID_LIMIT_MAP_CUDA);
  fillValueCukW<<<G, B>>>(a, N, v);
}
#pragma endregion




#pragma region SUM
/**
 * Compute the sum of values in an array, from a thread [device function].
 * @param x array to sum
 * @param N size of array
 * @param i start index
 * @param DI index stride
 * @returns sum of values in x[i..DI..N]
 */
template <class T>
inline T __device__ sumValuesThreadCud(const T *x, size_t N, size_t i, size_t DI) {
  ASSERT(x && DI);
  T a = T();
  for (; i<N; i+=DI)
    a += x[i];
  return a;
}


/**
 * Compute the sum of values in an array, within a block [device function].
 * @param a array to sum (updated, a[0] is the result)
 * @param N size of array
 * @param i thread index
 */
template <class T>
inline void __device__ sumValuesBlockReduceCudU(T *a, size_t N, size_t i) {
  ASSERT(a);
  // Reduce values in a to a[0] in reverse binary tree fashion.
  for (; N>1;) {
    size_t DN = (N+1)/2;
    if (i<N/2) a[i] += a[DN+i];
    __syncthreads();
    N = DN;
  }
}


/**
 * Compute the sum of values in an array [kernel].
 * @tparam CACHE size of shared memory cache
 * @param a partial result array (output)
 * @param x array to sum
 * @param N size of array to sum
 */
template <int CACHE=BLOCK_LIMIT_REDUCE_CUDA, class T>
void __global__ sumValuesCukW(T *a, const T *x, size_t N) {
  ASSERT(a && x);
  DEFINE_CUDA(t, b, B, G);
  __shared__ T cache[CACHE];
  // Store per-thread sum in shared cache (for further reduction).
  cache[t] = sumValuesThreadCud(x, N, B*b+t, G*B);
  // Wait for all threads within the block to finish.
  __syncthreads();
  // Reduce the sum in the cache to a single value in reverse binary tree fashion.
  sumValuesBlockReduceCudU(cache, B, t);
  // Store this per-block sum into a partial result array.
  if (t==0) a[b] = cache[0];
}


/**
 * Compute the sum of values in an array, using memcpy approach.
 * @param a partial result array (output)
 * @param x array to sum
 * @param N size of array to sum
 */
template <class T>
inline void sumValuesMemcpyCuW(T *a, const T *x, size_t N) {
  ASSERT(a && x);
  const int B = blockSizeCu(N,   BLOCK_LIMIT_REDUCE_CUDA);
  const int G = gridSizeCu (N, B, GRID_LIMIT_REDUCE_CUDA);
  sumValuesCukW<<<G, B>>>(a, x, N);
}


/**
 * Compute the sum of values in an array, using inplace approach.
 * @param a result array (output, a[0] is the result)
 * @param x array to sum
 * @param N size of array to sum
 */
template <class T>
inline void sumValuesInplaceCuW(T *a, const T *x, size_t N) {
  ASSERT(a && x);
  const int B = blockSizeCu(N ,  BLOCK_LIMIT_REDUCE_CUDA);
  const int G = gridSizeCu (N, B, GRID_LIMIT_REDUCE_CUDA);
  sumValuesCukW<<<G, B>>>(a, x, N);
  TRY_CUDA( cudaDeviceSynchronize() );
  sumValuesCukW<GRID_LIMIT_REDUCE_CUDA><<<1, G>>>(a, a, G);
}
#pragma endregion




#pragma region LI-NORM
/**
 * Compute the L∞-norm of an array, from a thread [device function].
 * @param x array to compute on
 * @param N size of array
 * @param i start index
 * @param DI index stride
 * @returns ||x[i..DI..N]||_∞
 */
template <class T>
inline T __device__ liNormThreadCud(const T *x, size_t N, size_t i, size_t DI) {
  ASSERT(x && DI);
  T a = T();  // TODO: use numeric_limits<T>::min()?
  for (; i<N; i+=DI)
    a = max(a, x[i]);
  return a;
}


/**
 * Compute the L∞-norm of an array, within a block [device function].
 * @param a array to compute on (updated, a[0] is the result)
 * @param N size of array
 * @param i thread index
 */
template <class T>
inline void __device__ liNormBlockReduceCudU(T *a, size_t N, size_t i) {
  ASSERT(a);
  // Reduce values in a to a[0] in reverse binary tree fashion.
  for (; N>1;) {
    size_t DN = (N+1)/2;
    if (i<N/2) a[i] = max(a[i], a[DN+i]);
    __syncthreads();
    N = DN;
  }
}


/**
 * Compute the L∞-norm of an array [kernel].
 * @tparam CACHE size of shared memory cache
 * @param a partial result array (output)
 * @param x array to compute on
 * @param N size of array to compute on
 */
template <int CACHE=BLOCK_LIMIT_REDUCE_CUDA, class T>
void __global__ liNormCukW(T *a, const T *x, size_t N) {
  ASSERT(a && x);
  DEFINE_CUDA(t, b, B, G);
  __shared__ T cache[CACHE];
  // Store per-thread L∞-norm in shared cache (for further reduction).
  cache[t] = liNormThreadCud(x, N, B*b+t, G*B);
  // Wait for all threads within the block to finish.
  __syncthreads();
  // Reduce the L∞-norms in cache to a single value in reverse binary tree fashion.
  liNormBlockReduceCudU(cache, B, t);
  // Store this per-block L∞-norm into a partial result array.
  if (t==0) a[b] = cache[0];
}


/**
 * Compute the L∞-norm of an array, using memcpy approach.
 * @param a partial result array (output)
 * @param x array to compute on
 * @param N size of array to compute on
 */
template <class T>
inline void liNormMemcpyCuW(T *a, const T *x, size_t N) {
  ASSERT(a && x);
  const int B = blockSizeCu(N,   BLOCK_LIMIT_REDUCE_CUDA);
  const int G = gridSizeCu (N, B, GRID_LIMIT_REDUCE_CUDA);
  liNormCukW<<<G, B>>>(a, x, N);
}


/**
 * Compute the L∞-norm of an array, using inplace approach.
 * @param a result array (output, a[0] is the result)
 * @param x array to compute on
 * @param N size of array to compute on
 */
template <class T>
inline void liNormInplaceCuW(T *a, const T *x, size_t N) {
  ASSERT(a && x);
  const int B = blockSizeCu(N,   BLOCK_LIMIT_REDUCE_CUDA);
  const int G = gridSizeCu (N, B, GRID_LIMIT_REDUCE_CUDA);
  liNormCukW<<<G, B>>>(a, x, N);
  TRY_CUDA( cudaDeviceSynchronize() );
  liNormCukW<GRID_LIMIT_REDUCE_CUDA><<<1, G>>>(a, a, G);
}
#pragma endregion




#pragma region LI-NORM DELTA
/**
 * Compute L∞-norm of the difference between two arrays, from a thread [device function].
 * @param x first array
 * @param y second array
 * @param N size of each array
 * @param i start index
 * @param DI index stride
 * @returns ||x[i..DI..N] - y[i..DI..N]||_∞
 */
template <class T>
inline T __device__ liNormDeltaThreadCud(const T *x, const T *y, size_t N, size_t i, size_t DI) {
  ASSERT(x && y && DI);
  T a = T();  // TODO: use numeric_limits<T>::min()?
  for (; i<N; i+=DI)
    a = max(a, abs(x[i] - y[i]));
  return a;
}


/**
 * Compute L∞-norm of the difference between two arrays [kernel].
 * @tparam CACHE size of shared memory cache
 * @param a partial result array (output)
 * @param x first array
 * @param y second array
 * @param N size of each array
 */
template <int CACHE=BLOCK_LIMIT_REDUCE_CUDA, class T>
void __global__ liNormDeltaCukW(T *a, const T *x, const T *y, size_t N) {
  ASSERT(a && x);
  DEFINE_CUDA(t, b, B, G);
  __shared__ T cache[CACHE];
  // Store per-thread delta L∞-norm in shared cache (for further reduction).
  cache[t] = liNormDeltaThreadCud(x, y, N, B*b+t, G*B);
  // Wait for all threads within the block to finish.
  __syncthreads();
  // Reduce the delta L∞-norms in cache to a single value in reverse binary tree fashion.
  liNormBlockReduceCudU(cache, B, t);
  // Store this per-block delta L∞-norm into a partial result array.
  if (t==0) a[b] = cache[0];
}


/**
 * Compute L∞-norm of the difference between two arrays, using memcpy approach.
 * @param a partial result array (output)
 * @param x first array
 * @param y second array
 * @param N size of each array
 */
template <class T>
inline void liNormDeltaMemcpyCuW(T *a, const T *x, const T *y, size_t N) {
  ASSERT(a && x && y);
  const int B = blockSizeCu(N,   BLOCK_LIMIT_REDUCE_CUDA);
  const int G = gridSizeCu (N, B, GRID_LIMIT_REDUCE_CUDA);
  liNormDeltaCukW<<<G, B>>>(a, x, y, N);
}


/**
 * Compute L∞-norm of the difference between two arrays, using inplace approach.
 * @param a result array (output, a[0] is the result)
 * @param x first array
 * @param y second array
 * @param N size of each array
 */
template <class T>
inline void liNormDeltaInplaceCuW(T *a, const T *x, const T *y, size_t N) {
  ASSERT(a && x && y);
  const int B = blockSizeCu(N,   BLOCK_LIMIT_REDUCE_CUDA);
  const int G = gridSizeCu (N, B, GRID_LIMIT_REDUCE_CUDA);
  liNormDeltaCukW<<<G, B>>>(a, x, y, N);
  TRY_CUDA( cudaDeviceSynchronize() );
  liNormCukW<GRID_LIMIT_REDUCE_CUDA><<<1, G>>>(a, a, G);
}
#pragma endregion
#pragma endregion
