#include <utility>
#include <random>
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
  using K = typename G::key_type;
  using V = typename G::edge_value_type;
  random_device dev;
  default_random_engine rnd(dev());
  int repeat = REPEAT_METHOD;
  double   M = edgeWeightOmp(x)/2;
  // Follow a specific result logging format, which can be easily parsed later.
  auto flog = [&](const auto& ans, const char *technique) {
    printf(
      "{%09.1fms, %09.1fms mark, %09.1fms init, %09.1fms firstpass, %09.1fms locmove, %09.1fms refine, %09.1fms aggr, %09.1fms split, %.3e aff, %04d iters, %03d passes, %01.9f modularity, %zu/%zu disconnected} %s\n",
      ans.time, ans.markingTime, ans.initializationTime, ans.firstPassTime, ans.localMoveTime, refinementTime(ans), ans.aggregationTime, ans.splittingTime,
      double(ans.affectedVertices), ans.iterations, ans.passes, getModularity(x, ans, M),
      countValue(communitiesDisconnectedOmp(x, ans.membership), char(1)),
      communities(x, ans.membership).size(), technique
    );
  };
  // Get community memberships on original graph (static).
  {
    auto a0 = louvainStaticOmp(x, {repeat});
    flog(a0, "louvainStaticOmp");
  }
  {
    auto a1 = louvainSplitLastStaticOmp<1>(x, {repeat});
    flog(a1, "louvainSplitLastStaticOmp1");
    auto a2 = louvainSplitLastStaticOmp<2>(x, {repeat});
    flog(a2, "louvainSplitLastStaticOmp2");
    // auto a3 = louvainSplitLastStaticOmp<3>(x, {repeat});
    // flog(a3, "louvainSplitLastStaticOmp3");
    auto a4 = louvainSplitLastStaticOmp<4>(x, {repeat});
    flog(a4, "louvainSplitLastStaticOmp4");
  }
  {
    auto a5 = louvainSplitIterationStaticOmp<1>(x, {repeat});
    flog(a5, "louvainSplitIterationStaticOmp1");
    auto a6 = louvainSplitIterationStaticOmp<2>(x, {repeat});
    flog(a6, "louvainSplitIterationStaticOmp2");
    // auto a7 = louvainSplitIterationStaticOmp<3>(x, {repeat});
    // flog(a7, "louvainSplitIterationStaticOmp3");
    auto a8 = louvainSplitIterationStaticOmp<4>(x, {repeat});
    flog(a8, "louvainSplitIterationStaticOmp4");
  }
  {
    auto b0 = leidenStaticOmp<false>(rnd, x, {repeat});
    flog(b0, "leidenStaticOmpGreedy");
  }
  {
    auto b1 = leidenSplitLastStaticOmp<1, false>(rnd, x, {repeat});
    flog(b1, "leidenSplitLastStaticOmpGreedy1");
    auto b2 = leidenSplitLastStaticOmp<2, false>(rnd, x, {repeat});
    flog(b2, "leidenSplitLastStaticOmpGreedy2");
    // auto b3 = leidenSplitLastStaticOmp<3, false>(rnd, x, {repeat});
    // flog(b3, "leidenSplitLastStaticOmpGreedy3");
    auto b4 = leidenSplitLastStaticOmp<4, false>(rnd, x, {repeat});
    flog(b4, "leidenSplitLastStaticOmpGreedy4");
  }
  {
    auto b5 = leidenSplitIterationStaticOmp<1, false>(rnd, x, {repeat});
    flog(b5, "leidenSplitIterationStaticOmpGreedy1");
    auto b6 = leidenSplitIterationStaticOmp<2, false>(rnd, x, {repeat});
    flog(b6, "leidenSplitIterationStaticOmpGreedy2");
    // auto b7 = leidenSplitIterationStaticOmp<3, false>(rnd, x, {repeat});
    // flog(b7, "leidenSplitIterationStaticOmpGreedy3");
    auto b8 = leidenSplitIterationStaticOmp<4, false>(rnd, x, {repeat});
    flog(b8, "leidenSplitIterationStaticOmpGreedy4");
  }
  {
    auto c0 = leidenStaticOmp<false>(rnd, x, {repeat, 1.0, 1e-06, 1.0, 10.0, 100, 100});
    flog(c0, "leidenStaticOmpGreedyMedium");
  }
  {
    auto c1 = leidenSplitLastStaticOmp<1, false>(rnd, x, {repeat, 1.0, 1e-06, 1.0, 10.0, 100, 100});
    flog(c1, "leidenSplitLastStaticOmpGreedyMedium1");
    auto c2 = leidenSplitLastStaticOmp<2, false>(rnd, x, {repeat, 1.0, 1e-06, 1.0, 10.0, 100, 100});
    flog(c2, "leidenSplitLastStaticOmpGreedyMedium2");
    // auto c3 = leidenSplitLastStaticOmp<3, false>(rnd, x, {repeat, 1.0, 1e-06, 1.0, 10.0, 100, 100});
    // flog(c3, "leidenSplitLastStaticOmpGreedyMedium3");
    auto c4 = leidenSplitLastStaticOmp<4, false>(rnd, x, {repeat, 1.0, 1e-06, 1.0, 10.0, 100, 100});
    flog(c4, "leidenSplitLastStaticOmpGreedyMedium4");
  }
  {
    auto c5 = leidenSplitIterationStaticOmp<1, false>(rnd, x, {repeat, 1.0, 1e-06, 1.0, 10.0, 100, 100});
    flog(c5, "leidenSplitIterationStaticOmpGreedyMedium1");
    auto c6 = leidenSplitIterationStaticOmp<2, false>(rnd, x, {repeat, 1.0, 1e-06, 1.0, 10.0, 100, 100});
    flog(c6, "leidenSplitIterationStaticOmpGreedyMedium2");
    // auto c7 = leidenSplitIterationStaticOmp<3, false>(rnd, x, {repeat, 1.0, 1e-06, 1.0, 10.0, 100, 100});
    // flog(c7, "leidenSplitIterationStaticOmpGreedyMedium3");
    auto c8 = leidenSplitIterationStaticOmp<4, false>(rnd, x, {repeat, 1.0, 1e-06, 1.0, 10.0, 100, 100});
    flog(c8, "leidenSplitIterationStaticOmpGreedyMedium4");
  }
  {
    auto d0 = leidenStaticOmp<false>(rnd, x, {repeat, 1.0, 1e-10, 1.0, 1.00, 100, 100});
    flog(d0, "leidenStaticOmpGreedyHeavy");
  }
  {
    auto d1 = leidenSplitLastStaticOmp<1, false>(rnd, x, {repeat, 1.0, 1e-10, 1.0, 1.00, 100, 100});
    flog(d1, "leidenSplitLastStaticOmpGreedyHeavy1");
    auto d2 = leidenSplitLastStaticOmp<2, false>(rnd, x, {repeat, 1.0, 1e-10, 1.0, 1.00, 100, 100});
    flog(d2, "leidenSplitLastStaticOmpGreedyHeavy2");
    // auto d3 = leidenSplitLastStaticOmp<3, false>(rnd, x, {repeat, 1.0, 1e-10, 1.0, 1.00, 100, 100});
    // flog(d3, "leidenSplitLastStaticOmpGreedyHeavy3");
    auto d4 = leidenSplitLastStaticOmp<4, false>(rnd, x, {repeat, 1.0, 1e-10, 1.0, 1.00, 100, 100});
    flog(d4, "leidenSplitLastStaticOmpGreedyHeavy4");
  }
  {
    auto d5 = leidenSplitIterationStaticOmp<1, false>(rnd, x, {repeat, 1.0, 1e-10, 1.0, 1.00, 100, 100});
    flog(d5, "leidenSplitIterationStaticOmpGreedyHeavy1");
    auto d6 = leidenSplitIterationStaticOmp<2, false>(rnd, x, {repeat, 1.0, 1e-10, 1.0, 1.00, 100, 100});
    flog(d6, "leidenSplitIterationStaticOmpGreedyHeavy2");
    // auto d7 = leidenSplitIterationStaticOmp<3, false>(rnd, x, {repeat, 1.0, 1e-10, 1.0, 1.00, 100, 100});
    // flog(d7, "leidenSplitIterationStaticOmpGreedyHeavy3");
    auto d8 = leidenSplitIterationStaticOmp<4, false>(rnd, x, {repeat, 1.0, 1e-10, 1.0, 1.00, 100, 100});
    flog(d8, "leidenSplitIterationStaticOmpGreedyHeavy4");
  }
  {
    auto e0 = leidenStaticOmp<true> (rnd, x, {repeat});
    flog(e0, "leidenStaticOmpRandom");
  }
  {
    auto e1 = leidenSplitLastStaticOmp<1, true>(rnd, x, {repeat});
    flog(e1, "leidenSplitLastStaticOmpRandom1");
    auto e2 = leidenSplitLastStaticOmp<2, true>(rnd, x, {repeat});
    flog(e2, "leidenSplitLastStaticOmpRandom2");
    // auto e3 = leidenSplitLastStaticOmp<3, true>(rnd, x, {repeat});
    // flog(e3, "leidenSplitLastStaticOmpRandom3");
    auto e4 = leidenSplitLastStaticOmp<4, true>(rnd, x, {repeat});
    flog(e4, "leidenSplitLastStaticOmpRandom4");
  }
  {
    auto e5 = leidenSplitIterationStaticOmp<1, true>(rnd, x, {repeat});
    flog(e5, "leidenSplitIterationStaticOmpRandom1");
    auto e6 = leidenSplitIterationStaticOmp<2, true>(rnd, x, {repeat});
    flog(e6, "leidenSplitIterationStaticOmpRandom2");
    // auto e7 = leidenSplitIterationStaticOmp<3, true>(rnd, x, {repeat});
    // flog(e7, "leidenSplitIterationStaticOmpRandom3");
    auto e8 = leidenSplitIterationStaticOmp<4, true>(rnd, x, {repeat});
    flog(e8, "leidenSplitIterationStaticOmpRandom4");
  }
  {
    auto f0 = leidenStaticOmp<true> (rnd, x, {repeat, 1.0, 1e-06, 1.0, 10.0, 100, 100});
    flog(f0, "leidenStaticOmpRandomMedium");
  }
  {
    auto f1 = leidenSplitLastStaticOmp<1, true>(rnd, x, {repeat, 1.0, 1e-06, 1.0, 10.0, 100, 100});
    flog(f1, "leidenSplitLastStaticOmpRandomMedium1");
    auto f2 = leidenSplitLastStaticOmp<2, true>(rnd, x, {repeat, 1.0, 1e-06, 1.0, 10.0, 100, 100});
    flog(f2, "leidenSplitLastStaticOmpRandomMedium2");
    // auto f3 = leidenSplitLastStaticOmp<3, true>(rnd, x, {repeat, 1.0, 1e-06, 1.0, 10.0, 100, 100});
    // flog(f3, "leidenSplitLastStaticOmpRandomMedium3");
    auto f4 = leidenSplitLastStaticOmp<4, true>(rnd, x, {repeat, 1.0, 1e-06, 1.0, 10.0, 100, 100});
    flog(f4, "leidenSplitLastStaticOmpRandomMedium4");
  }
  {
    auto f5 = leidenSplitIterationStaticOmp<1, true>(rnd, x, {repeat, 1.0, 1e-06, 1.0, 10.0, 100, 100});
    flog(f5, "leidenSplitIterationStaticOmpRandomMedium1");
    auto f6 = leidenSplitIterationStaticOmp<2, true>(rnd, x, {repeat, 1.0, 1e-06, 1.0, 10.0, 100, 100});
    flog(f6, "leidenSplitIterationStaticOmpRandomMedium2");
    // auto f7 = leidenSplitIterationStaticOmp<3, true>(rnd, x, {repeat, 1.0, 1e-06, 1.0, 10.0, 100, 100});
    // flog(f7, "leidenSplitIterationStaticOmpRandomMedium3");
    auto f8 = leidenSplitIterationStaticOmp<4, true>(rnd, x, {repeat, 1.0, 1e-06, 1.0, 10.0, 100, 100});
    flog(f8, "leidenSplitIterationStaticOmpRandomMedium4");
  }
  {
    auto g0 = leidenStaticOmp<true> (rnd, x, {repeat, 1.0, 1e-10, 1.0, 1.00, 100, 100});
    flog(g0, "leidenStaticOmpRandomHeavy");
  }
  {
    auto g1 = leidenSplitLastStaticOmp<1, true>(rnd, x, {repeat, 1.0, 1e-10, 1.0, 1.00, 100, 100});
    flog(g1, "leidenSplitLastStaticOmpRandomHeavy1");
    auto g2 = leidenSplitLastStaticOmp<2, true>(rnd, x, {repeat, 1.0, 1e-10, 1.0, 1.00, 100, 100});
    flog(g2, "leidenSplitLastStaticOmpRandomHeavy2");
    // auto g3 = leidenSplitLastStaticOmp<3, true>(rnd, x, {repeat, 1.0, 1e-10, 1.0, 1.00, 100, 100});
    // flog(g3, "leidenSplitLastStaticOmpRandomHeavy3");
    auto g4 = leidenSplitLastStaticOmp<4, true>(rnd, x, {repeat, 1.0, 1e-10, 1.0, 1.00, 100, 100});
    flog(g4, "leidenSplitLastStaticOmpRandomHeavy4");
  }
  {
    auto g5 = leidenSplitIterationStaticOmp<1, true>(rnd, x, {repeat, 1.0, 1e-10, 1.0, 1.00, 100, 100});
    flog(g5, "leidenSplitIterationStaticOmpRandomHeavy1");
    auto g6 = leidenSplitIterationStaticOmp<2, true>(rnd, x, {repeat, 1.0, 1e-10, 1.0, 1.00, 100, 100});
    flog(g6, "leidenSplitIterationStaticOmpRandomHeavy2");
    // auto g7 = leidenSplitIterationStaticOmp<3, true>(rnd, x, {repeat, 1.0, 1e-10, 1.0, 1.00, 100, 100});
    // flog(g7, "leidenSplitIterationStaticOmpRandomHeavy3");
    auto g8 = leidenSplitIterationStaticOmp<4, true>(rnd, x, {repeat, 1.0, 1e-10, 1.0, 1.00, 100, 100});
    flog(g8, "leidenSplitIterationStaticOmpRandomHeavy4");
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
