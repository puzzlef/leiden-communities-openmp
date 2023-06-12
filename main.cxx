#include <cstdint>
#include <cstdio>
#include <utility>
#include <random>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include "src/main.hxx"

using namespace std;




// Fixed config
#ifndef TYPE
#define TYPE float
#endif
#ifndef MAX_THREADS
#define MAX_THREADS 64
#endif
#ifndef REPEAT_BATCH
#define REPEAT_BATCH 5
#endif
#ifndef REPEAT_METHOD
#define REPEAT_METHOD 1
#endif




// HELPERS
// -------

template <class G, class K>
inline double getModularity(const G& x, const LouvainResult<K>& a, double M) {
  auto fc = [&](auto u) { return a.membership[u]; };
  return modularityByOmp(x, fc, M, 1.0);
}




// GENERATE BATCH
// --------------

template <class G, class R>
inline auto addRandomEdges(G& a, R& rnd, size_t batchSize, size_t i, size_t n) {
  using K = typename G::key_type;
  using V = typename G::edge_value_type;
  int retries = 5;
  vector<tuple<K, K, V>> insertions;
  auto fe = [&](auto u, auto v, auto w) {
    a.addEdge(u, v, w);
    a.addEdge(v, u, w);
    insertions.push_back(make_tuple(u, v, w));
    insertions.push_back(make_tuple(v, u, w));
    return true;
  };
  for (size_t l=0; l<batchSize; ++l)
    retry([&]() { return addRandomEdge(a, rnd, i, n, V(1), fe); }, retries);
  updateOmpU(a);
  return insertions;
}


template <class G, class R>
auto removeRandomEdges(G& a, R& rnd, size_t batchSize, size_t i, size_t n) {
  using K = typename G::key_type;
  int retries = 5;
  vector<tuple<K, K>> deletions;
  auto fe = [&](auto u, auto v) {
    a.removeEdge(u, v);
    a.removeEdge(v, u);
    deletions.push_back(make_tuple(u, v));
    deletions.push_back(make_tuple(v, u));
    return true;
  };
  for (size_t l=0; l<batchSize; ++l)
    retry([&]() { return removeRandomEdge(a, rnd, i, n, fe); }, retries);
  updateOmpU(a);
  return deletions;
}




// PERFORM EXPERIMENT
// ------------------

template <class G, class R, class F>
inline void runAbsoluteBatches(const G& x, R& rnd, F fn) {
  size_t d = BATCH_DELETIONS_BEGIN;
  size_t i = BATCH_INSERTIONS_BEGIN;
  for (int epoch=0;; ++epoch) {
    for (int r=0; r<REPEAT_BATCH; ++r) {
      auto y  = duplicate(x);
      for (int sequence=0; sequence<BATCH_LENGTH; ++sequence) {
      auto deletions  = removeRandomEdges(y, rnd, d, 1, x.span()-1);
      auto insertions = addRandomEdges   (y, rnd, i, 1, x.span()-1);
        fn(y, deletions, insertions, sequence, epoch);
      }
    }
    if (d>=BATCH_DELETIONS_END && i>=BATCH_INSERTIONS_END) break;
    d BATCH_DELETIONS_STEP;
    i BATCH_INSERTIONS_STEP;
    d = min(d, size_t(BATCH_DELETIONS_END));
    i = min(i, size_t(BATCH_INSERTIONS_END));
  }
}


template <class G, class R, class F>
inline void runRelativeBatches(const G& x, R& rnd, F fn) {
  double d = BATCH_DELETIONS_BEGIN;
  double i = BATCH_INSERTIONS_BEGIN;
  for (int epoch=0;; ++epoch) {
    for (int r=0; r<REPEAT_BATCH; ++r) {
      auto y  = duplicate(x);
      for (int sequence=0; sequence<BATCH_LENGTH; ++sequence) {
      auto deletions  = removeRandomEdges(y, rnd, size_t(d * x.size()/2), 1, x.span()-1);
      auto insertions = addRandomEdges   (y, rnd, size_t(i * x.size()/2), 1, x.span()-1);
        fn(y, deletions, insertions, sequence, epoch);
      }
    }
    if (d>=BATCH_DELETIONS_END && i>=BATCH_INSERTIONS_END) break;
    d BATCH_DELETIONS_STEP;
    i BATCH_INSERTIONS_STEP;
    d = min(d, double(BATCH_DELETIONS_END));
    i = min(i, double(BATCH_INSERTIONS_END));
  }
}


template <class G, class R, class F>
inline void runBatches(const G& x, R& rnd, F fn) {
  if (BATCH_UNIT=="%") runRelativeBatches(x, rnd, fn);
  else runAbsoluteBatches(x, rnd, fn);
}


template <class F>
inline void runThreadsWithBatch(int epoch, F fn) {
  int t = NUM_THREADS_BEGIN;
  for (int l=0; l<epoch && t<=NUM_THREADS_END; ++l)
    t NUM_THREADS_STEP;
  omp_set_num_threads(t);
  fn(t);
  omp_set_num_threads(MAX_THREADS);
}


template <class F>
inline void runThreadsAll(F fn) {
  for (int t=NUM_THREADS_BEGIN; t<=NUM_THREADS_END; t NUM_THREADS_STEP) {
    omp_set_num_threads(t);
    fn(t);
    omp_set_num_threads(MAX_THREADS);
  }
}


template <class F>
inline void runThreads(int epoch, F fn) {
  if (NUM_THREADS_MODE=="with-batch") runThreadsWithBatch(epoch, fn);
  else runThreadsAll(fn);
}


template <class G>
void runExperiment(const G& x) {
  using K = typename G::key_type;
  using V = typename G::edge_value_type;
  random_device dev;
  default_random_engine rnd(dev());
  int repeat  = REPEAT_METHOD;
  int retries = 5;
  vector<K> *init = nullptr;
  double M = edgeWeightOmp(x)/2;
  // Get community memberships on original graph (static).
  auto b0 = louvainStaticOmp(x, init);
  #if BATCH_LENGTH>1
  vector<K> B2, B3, B4;
  vector<K> C2, C3, C4;
  vector<K> D2, D3, D4;
  vector<K> E2, E3, E4;
  #else
  const auto& B2 = b0.membership;
  const auto& B3 = b0.membership;
  const auto& B4 = b0.membership;
  const auto& C2 = b0.membership;
  const auto& C3 = b0.membership;
  const auto& C4 = b0.membership;
  const auto& D2 = b0.membership;
  const auto& D3 = b0.membership;
  const auto& D4 = b0.membership;
  const auto& E2 = b0.membership;
  const auto& E3 = b0.membership;
  const auto& E4 = b0.membership;
  #endif
  // Get community memberships on updated graph (dynamic).
  runBatches(x, rnd, [&](const auto& y, const auto& deletions, const auto& insertions, int sequence, int epoch) {
    double M = edgeWeightOmp(y)/2;
    // Follow a specific result logging format, which can be easily parsed later.
    auto glog = [&](const auto& ans, const char *technique, int numThreads) {
      printf(
        "{-%.3e/+%.3e [%04d] batch, %03d threads} -> "
        "{%09.1f/%09.1fms, %04d iters, %03d passes, %01.9f modularity} %s\n",
        double(deletions.size()), double(insertions.size()), sequence, numThreads,
        ans.preprocessingTime, ans.time, ans.iterations, ans.passes, getModularity(y, ans, M), technique
      );
    };
    #if BATCH_LENGTH>1
    if (sequence==0) {
      B2 = b0.membership;
      B3 = b0.membership;
      B4 = b0.membership;
      C2 = b0.membership;
      C3 = b0.membership;
      C4 = b0.membership;
      D2 = b0.membership;
      D3 = b0.membership;
      D4 = b0.membership;
      E2 = b0.membership;
      E3 = b0.membership;
      E4 = b0.membership;
    }
    #endif
    // Adjust number of threads.
    runThreads(epoch, [&](int numThreads) {
      auto flog = [&](const auto& ans, const char *technique) {
        glog(ans, technique, numThreads);
      };
      // Find static Louvain.
      auto b1 = louvainStaticOmp<false, false>(y, init, {repeat});
      flog(b1, "louvainStaticOmp");
      auto c1 = louvainStaticOmp<true, false> (y, init, {repeat});
      flog(c1, "louvainStaticOmp<JUMP>");
      auto d1 = louvainStaticOmp<false, true> (y, init, {repeat});
      flog(d1, "louvainStaticOmp<REFINE>");
      auto e1 = louvainStaticOmp<true, true>  (y, init, {repeat});
      flog(e1, "louvainStaticOmp<JUMP, REFINE>");
      // Find naive-dynamic Louvain.
      auto b2 = louvainStaticOmp<false, false>(y, &B2, {repeat});
      flog(b2, "louvainNaiveDynamicOmp");
      auto c2 = louvainStaticOmp<true, false> (y, &C2, {repeat});
      flog(c2, "louvainNaiveDynamicOmp<JUMP>");
      auto d2 = louvainStaticOmp<false, true> (y, &D2, {repeat});
      flog(d2, "louvainNaiveDynamicOmp<REFINE>");
      auto e2 = louvainStaticOmp<true, true>  (y, &E2, {repeat});
      flog(e2, "louvainNaiveDynamicOmp<JUMP, REFINE>");
      // Find frontier based dynamic Louvain.
      auto b4 = louvainDynamicFrontierOmp<false, false>(y, deletions, insertions, &B4, {repeat});
      flog(b4, "louvainDynamicFrontierOmp");
      auto c4 = louvainDynamicFrontierOmp<true, false> (y, deletions, insertions, &C4, {repeat});
      flog(c4, "louvainDynamicFrontierOmp<JUMP>");
      auto d4 = louvainDynamicFrontierOmp<false, true> (y, deletions, insertions, &D4, {repeat});
      flog(d4, "louvainDynamicFrontierOmp<REFINE>");
      auto e4 = louvainDynamicFrontierOmp<true, true>  (y, deletions, insertions, &E4, {repeat});
      flog(e4, "louvainDynamicFrontierOmp<JUMP, REFINE>");
      // Find delta-screening based dynamic Louvain.
      auto b3 = louvainDynamicDeltaScreeningOmp<false, false>(y, deletions, insertions, &B3, {repeat});
      flog(b3, "louvainDynamicDeltaScreeningOmp");
      auto c3 = louvainDynamicDeltaScreeningOmp<true, false> (y, deletions, insertions, &C3, {repeat});
      flog(c3, "louvainDynamicDeltaScreeningOmp<JUMP>");
      auto d3 = louvainDynamicDeltaScreeningOmp<false, true> (y, deletions, insertions, &D3, {repeat});
      flog(d3, "louvainDynamicDeltaScreeningOmp<REFINE>");
      auto e3 = louvainDynamicDeltaScreeningOmp<true, true>  (y, deletions, insertions, &E3, {repeat});
      flog(e3, "louvainDynamicDeltaScreeningOmp<JUMP,REFINE>");
      #if BATCH_LENGTH>1
      B2 = b2.membership;
      B3 = b3.membership;
      B4 = b4.membership;
      C2 = c2.membership;
      C3 = c3.membership;
      C4 = c4.membership;
      D2 = d2.membership;
      D3 = d3.membership;
      D4 = d4.membership;
      E2 = e2.membership;
      E3 = e3.membership;
      E4 = e4.membership;
      #endif
    });
  });
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
  OutDiGraph<K, None, V> x;
  readMtxOmpW(x, file, weighted); LOG(""); println(x);
  if (!symmetric) { x = symmetricizeOmp(x); LOG(""); print(x); printf(" (symmetricize)\n"); }
  runExperiment(x);
  printf("\n");
  return 0;
}
