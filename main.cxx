#include <utility>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include "inc/main.hxx"

using namespace std;




#pragma region CONFIGURATION
#ifndef TYPE
/** Type of edge weights. */
#define TYPE float
#endif
#ifndef MAX_THREADS
/** Maximum number of threads to use. */
#define MAX_THREADS 64
#endif
#ifndef REPEAT_METHOD
/** Number of times to repeat each method. */
#define REPEAT_METHOD 5
#endif
#pragma endregion




// HELPERS
// -------

template <class G, class R>
inline double getModularity(const G& x, const R& a, double M) {
  auto fc = [&](auto u) { return a.membership[u]; };
  return modularityByOmp(x, fc, M, 1.0);
}


template <class K, class W>
inline float refinementTime(const LouvainResult<K, W>& a) {
  return 0;
}
template <class K, class W>
inline float refinementTime(const LeidenResult<K, W>& a) {
  return a.refinementTime;
}




// PERFORM EXPERIMENT
// ------------------

template <class G>
void runExperiment(const G& x) {
  int repeat = REPEAT_METHOD;
  double   M = edgeWeightOmp(x)/2;
  // Follow a specific result logging format, which can be easily parsed later.
  auto flog = [&](const auto& ans, const char *technique) {
    printf(
      "{%09.1fms, %09.1fms mark, %09.1fms init, %09.1fms firstpass, %09.1fms locmove, %09.1fms refine, %09.1fms aggr, %.3e aff, %04d iters, %03d passes, %01.9f modularity, %zu/%zu disconnected} %s\n",
      ans.time, ans.markingTime, ans.initializationTime, ans.firstPassTime, ans.localMoveTime, refinementTime(ans), ans.aggregationTime,
      double(ans.affectedVertices), ans.iterations, ans.passes, getModularity(x, ans, M),
      countValue(communitiesDisconnectedOmp(x, ans.membership), char(1)),
      communities(x, ans.membership).size(), technique
    );
  };
  // Get community memberships on original graph (static).
  {
    auto a0 = louvainStaticOmp(x, {repeat});
    flog(a0, "louvainStaticOmp");
    auto b0 = leidenStaticOmp(x, {repeat});
    flog(b0, "leidenStaticOmp");
  }
}


int main(int argc, char **argv) {
  using K = uint32_t;
  using V = TYPE;
  install_sigsegv();
  char *file     = argv[1];
  bool symmetric = argc>2? stoi(argv[2]) : false;
  bool weighted  = argc>3? stoi(argv[3]) : false;
  omp_set_num_threads(MAX_THREADS);
  LOG("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  LOG("Loading graph %s ...\n", file);
  DiGraph<K, None, V> x;
  readMtxOmpW(x, file, weighted); LOG(""); println(x);
  if (!symmetric) { x = symmetricizeOmp(x); LOG(""); print(x); printf(" (symmetricize)\n"); }
  runExperiment(x);
  printf("\n");
  return 0;
}
