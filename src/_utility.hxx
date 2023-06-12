#pragma once
#include <utility>
#include <chrono>
#ifdef MPI
#include "_mpi.hxx"
#endif

using std::pair;
using std::chrono::microseconds;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;




// PAIR
// ----

template <class K, class V>
struct PairFirst  { inline K& operator()(pair<K, V>& x) noexcept { return x.first; } };
template <class K, class V>
struct PairSecond { inline V& operator()(pair<K, V>& x) noexcept { return x.second; } };
template <class K, class V>
struct ConstPairFirst  { inline const K& operator()(const pair<K, V>& x) noexcept { return x.first; } };
template <class K, class V>
struct ConstPairSecond { inline const V& operator()(const pair<K, V>& x) noexcept { return x.second; } };
template <class K, class V>
struct PairFirstValue  { inline K operator()(const pair<K, V>& x) noexcept { return x.first; } };
template <class K, class V>
struct PairSecondValue { inline V operator()(const pair<K, V>& x) noexcept { return x.second; } };




// MEASURE DURATION
// ----------------

/** Get current time. */
inline auto timeNow() {
  return high_resolution_clock::now();
}

/** Get time duration in milliseconds. */
template <class T>
inline float duration(const T& start, const T& stop) {
  auto a = duration_cast<microseconds>(stop - start);
  return a.count()/1000.0f;
}

/** Get time duration in milliseconds. */
template <class T>
inline float duration(const T& start) {
  auto stop = timeNow();
  return duration(start, stop);
}


template <class F>
inline float measureDuration(F fn, int N=1) {
  auto start = timeNow();
  for (int i=0; i<N; ++i)
    fn();
  auto stop  = timeNow();
  return duration(start, stop)/N;
}

#ifdef MPI
template <class F>
inline float measureDurationMpi(F fn, int N=1) {
  double total = 0;
  for (int i=0; i<N; ++i) {
    // Match up with other processes before start.
    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();
    fn();
    // Let all processes complete together.
    MPI_Barrier(MPI_COMM_WORLD);
    double stop  = MPI_Wtime();
    total += stop - start;
  }
  // Report in milliseconds.
  return float(total*1000/N);
}
#endif


template <class F>
inline float measureDurationMarked(F fn, int N=1) {
  float total = 0;
  for (int i=0; i<N; ++i)
    fn([&](auto fm) { total += measureDuration(fm); });
  return total/N;
}

#ifdef MPI
template <class F>
inline float measureDurationMarkedMpi(F fn, int N=1) {
  float total = 0;
  for (int i=0; i<N; ++i)
    fn([&](auto fm) { total += measureDurationMpi(fm); });
  return total/N;
}
#endif




// RETRY
// -----

template <class F>
inline bool retry(F fn, int N=2) {
  for (int i=0; i<N; ++i)
    if (fn()) return true;
  return false;
}




// MOVE
// ----
// Conditional move, otherwise value.

#define CMOVE(c, t, f) \
  ((c)? move(t) : (f))

#define CMOVE_VECTOR(t, f) \
  CMOVE(!(t).empty(), t, f)

#define CMOVE_GRAPH(t, f) \
  CMOVE((t).order()>0, t, f)
