#pragma once
#include <utility>
#include <tuple>
#include <vector>
#include <random>
#include <algorithm>
#include "_main.hxx"
#include "Graph.hxx"
#include "properties.hxx"
#include "csr.hxx"
#ifdef OPENMP
#include <omp.h>
#endif

using std::tuple;
using std::vector;
using std::uniform_int_distribution;
using std::make_pair;
using std::move;
using std::swap;
using std::get;
using std::min;
using std::max;




#pragma region TYPES
/**
 * Options for Leiden algorithm.
 */
struct LeidenOptions {
  #pragma region DATA
  /** Number of times to repeat the algorithm [1]. */
  int repeat;
  /** Resolution parameter for modularity [1]. */
  double resolution;
  /** Tolerance for convergence [1e-2]. */
  double tolerance;
  /** Tolerance for aggregation [0.8]. */
  double aggregationTolerance;
  /** Tolerance drop factor after each pass [10]. */
  double toleranceDrop;
  /** Maximum number of iterations per pass [20]. */
  int maxIterations;
  /** Maximum number of passes [10]. */
  int maxPasses;
  #pragma endregion


  #pragma region CONSTRUCTORS
  /**
   * Define options for Leiden algorithm.
   * @param repeat number of times to repeat the algorithm [1]
   * @param resolution resolution parameter for modularity [1]
   * @param tolerance tolerance for convergence [1e-2]
   * @param aggregationTolerance tolerance for aggregation [0.8]
   * @param toleranceDrop tolerance drop factor after each pass [10]
   * @param maxIterations maximum number of iterations per pass [20]
   * @param maxPasses maximum number of passes [10]
   */
  LeidenOptions(int repeat=1, double resolution=1, double tolerance=1e-2, double aggregationTolerance=0.8, double toleranceDrop=10, int maxIterations=20, int maxPasses=10) :
  repeat(repeat), resolution(resolution), tolerance(tolerance), aggregationTolerance(aggregationTolerance), toleranceDrop(toleranceDrop), maxIterations(maxIterations), maxPasses(maxPasses) {}
  #pragma endregion
};


/** Weight to be used in hashtable. */
#define LEIDEN_WEIGHT_TYPE double




/**
 * Result of Leiden algorithm.
 * @tparam K key type (vertex-id)
 */
template <class K>
struct LeidenResult {
  #pragma region DATA
  /** Community membership each vertex belongs to. */
  vector<K> membership;
  /** Number of iterations performed. */
  int iterations;
  /** Number of passes performed. */
  int passes;
  /** Time spent in milliseconds. */
  float time;
  /** Time spent in milliseconds for preprocessing. */
  float preprocessingTime;
  /** Time spent in milliseconds in first pass. */
  float firstPassTime;
  /** Time spent in milliseconds in local-moving phase. */
  float localMoveTime;
  /** Time spent in milliseconds in aggregation phase. */
  float aggregationTime;
  #pragma endregion


  #pragma region CONSTRUCTORS
  /**
   * Result of Leiden algorithm.
   * @param membership community membership each vertex belongs to
   * @param iterations number of iterations performed
   * @param passes number of passes performed
   * @param time time spent in milliseconds
   * @param preprocessingTime time spent in milliseconds for preprocessing
   * @param firstPassTime time spent in milliseconds in first pass
   * @param localMoveTime time spent in milliseconds in local-moving phase
   * @param aggregationTime time spent in milliseconds in aggregation phase
   */
  LeidenResult(vector<K>&& membership, int iterations=0, int passes=0, float time=0, float preprocessingTime=0, float firstPassTime=0, float localMoveTime=0, float aggregationTime=0) :
  membership(membership), iterations(iterations), passes(passes), time(time), preprocessingTime(preprocessingTime), firstPassTime(firstPassTime), localMoveTime(localMoveTime), aggregationTime(aggregationTime) {}


  /**
   * Result of Leiden algorithm.
   * @param membership community membership each vertex belongs to (moved)
   * @param iterations number of iterations performed
   * @param passes number of passes performed
   * @param time time spent in milliseconds
   * @param preprocessingTime time spent in milliseconds for preprocessing
   * @param firstPassTime time spent in milliseconds in first pass
   * @param localMoveTime time spent in milliseconds in local-moving phase
   * @param aggregationTime time spent in milliseconds in aggregation phase
   */
  LeidenResult(vector<K>& membership, int iterations=0, int passes=0, float time=0, float preprocessingTime=0, float firstPassTime=0, float localMoveTime=0, float aggregationTime=0) :
  membership(move(membership)), iterations(iterations), passes(passes), time(time), preprocessingTime(preprocessingTime), firstPassTime(firstPassTime), localMoveTime(localMoveTime), aggregationTime(aggregationTime) {}
  #pragma endregion
};
#pragma endregion




#pragma region METHODS
#pragma region HASHTABLES
/**
 * Allocate a number of hashtables.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param S size of each hashtable
 */
template <class K, class W>
inline void leidenAllocateHashtablesW(vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, size_t S) {
  size_t N = vcs.size();
  for (size_t i=0; i<N; ++i) {
    vcs[i]   = new vector<K>();
    vcout[i] = new vector<W>(S);
  }
}


/**
 * Free a number of hashtables.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 */
template <class K, class W>
inline void leidenFreeHashtablesW(vector<vector<K>*>& vcs, vector<vector<W>*>& vcout) {
  size_t N = vcs.size();
  for (size_t i=0; i<N; ++i) {
    delete vcs[i];
    delete vcout[i];
  }
}
#pragma endregion




#pragma region RANDOM NUMBER GENERATORS
/**
 * Allocate a number of random number generators.
 * @param rng per-thread random number generators (updated)
 * @param rnd random number generator for seeding
 */
template <class RND>
inline void leidenAllocateRngsW(vector<xorshift32_engine*>& rng, RND& rnd) {
  uniform_int_distribution<uint32_t> dis(0, UINT32_MAX);
  size_t N = rng.size();
  for (size_t i=0; i<N; ++i)
    rng[i] = new xorshift32_engine(dis(rnd));
}


/**
 * Free a number of random number generators.
 * @param rng per-thread random number generators (updated)
 */
inline void leidenFreeRngsW(vector<xorshift32_engine*>& rng) {
  size_t N = rng.size();
  for (size_t i=0; i<N; ++i)
    delete rng[i];
}
#pragma endregion




#pragma region INITIALIZE
/**
 * Find the total edge weight of each vertex.
 * @param vtot total edge weight of each vertex (updated, must be initialized)
 * @param x original graph
 */
template <class G, class W>
inline void leidenVertexWeightsW(vector<W>& vtot, const G& x) {
  x.forEachVertexKey([&](auto u) {
    x.forEachEdge(u, [&](auto v, auto w) {
      vtot[u] += w;
    });
  });
}


#ifdef OPENMP
/**
 * Find the total edge weight of each vertex.
 * @param vtot total edge weight of each vertex (updated, must be initialized)
 * @param x original graph
 */
template <class G, class W>
inline void leidenVertexWeightsOmpW(vector<W>& vtot, const G& x) {
  using  K = typename G::key_type;
  size_t S = x.span();
  #pragma omp parallel for schedule(dynamic, 2048)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    x.forEachEdge(u, [&](auto v, auto w) { vtot[u] += w; });
  }
}
#endif


/**
 * Find the total edge weight of each community.
 * @param ctot total edge weight of each community (updated, must be initialized)
 * @param x original graph
 * @param vcom community each vertex belongs to
 * @param vtot total edge weight of each vertex
 */
template <class G, class K, class W>
inline void leidenCommunityWeightsW(vector<W>& ctot, const G& x, const vector<K>& vcom, const vector<W>& vtot) {
  x.forEachVertexKey([&](auto u) {
    K c = vcom[u];
    ctot[c] += vtot[u];
  });
}


#ifdef OPENMP
/**
 * Find the total edge weight of each community.
 * @param ctot total edge weight of each community (updated, must be initialized)
 * @param x original graph
 * @param vcom community each vertex belongs to
 * @param vtot total edge weight of each vertex
 */
template <class G, class K, class W>
inline void leidenCommunityWeightsOmpW(vector<W>& ctot, const G& x, const vector<K>& vcom, const vector<W>& vtot) {
  size_t S = x.span();
  #pragma omp parallel for schedule(static, 2048)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    K c = vcom[u];
    #pragma omp atomic
    ctot[c] += vtot[u];
  }
}
#endif


/**
 * Initialize community weights such that each vertex is its own community.
 * @param ctot total edge weight of each community (updated)
 * @param x original graph
 * @param vtot total edge weight of each vertex
 */
template <class G, class W>
inline void leidenInitializeCommunityWeightsW(vector<W>& ctot, const G& x, const vector<W>& vtot) {
  x.forEachVertexKey([&](auto u) {
    ctot[u] = vtot[u];
  });
}


#ifdef OPENMP
/**
 * Initialize community weights such that each vertex is its own community.
 * @param ctot total edge weight of each community (updated)
 * @param x original graph
 * @param vtot total edge weight of each vertex
 */
template <class G, class W>
inline void leidenInitializeCommunityWeightsOmpW(vector<W>& ctot, const G& x, const vector<W>& vtot) {
  using  K = typename G::key_type;
  size_t S = x.span();
  #pragma omp parallel for schedule(auto)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    ctot[u] = vtot[u];
  }
}
#endif


/**
 * Initialize communities such that each vertex is its own community.
 * @param vcom community each vertex belongs to (updated, must be initialized)
 * @param vcob community bound each vertex belongs to (updated, must be initialized)
 * @param ctot total edge weight of each community (updated, must be initialized)
 * @param x original graph
 * @param vtot total edge weight of each vertex
 */
template <class G, class K, class W>
inline void leidenInitializeW(vector<K>& vcom, vector<K>& vcob, vector<W>& ctot, const G& x, const vector<W>& vtot) {
  x.forEachVertexKey([&](auto u) {
    vcom[u] = u;
    vcob[u] = u;
    ctot[u] = vtot[u];
  });
}


#ifdef OPENMP
/**
 * Initialize communities such that each vertex is its own community.
 * @param vcom community each vertex belongs to (updated, must be initialized)
 * @param vcob community bound each vertex belongs to (updated, must be initialized)
 * @param ctot total edge weight of each community (updated, must be initialized)
 * @param x original graph
 * @param vtot total edge weight of each vertex
 */
template <class G, class K, class W>
inline void leidenInitializeOmpW(vector<K>& vcom, vector<K>& vcob, vector<W>& ctot, const G& x, const vector<W>& vtot) {
  size_t S = x.span();
  #pragma omp parallel for schedule(static, 2048)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    vcom[u] = u;
    vcob[u] = u;
    ctot[u] = vtot[u];
  }
}
#endif


/**
 * Initialize communities from given initial communities.
 * @param vcom community each vertex belongs to (updated, must be initialized)
 * @param vcob community bound each vertex belongs to (updated, must be initialized)
 * @param ctot total edge weight of each community (updated, must be initialized)
 * @param x original graph
 * @param vtot total edge weight of each vertex
 * @param q initial community each vertex belongs to
 */
template <class G, class K, class W>
inline void leidenInitializeFromW(vector<K>& vcom, vector<K>& vcob, vector<W>& ctot, const G& x, const vector<W>& vtot, const vector<K>& q) {
  x.forEachVertexKey([&](auto u) {
    K c = q[u];
    vcom[u]  = u;
    vcob[u]  = c;
    ctot[c] += vtot[u];
  });
}


#ifdef OPENMP
/**
 * Initialize communities from given initial communities.
 * @param vcom community each vertex belongs to (updated, must be initialized)
 * @param vcob community bound each vertex belongs to (updated, must be initialized)
 * @param ctot total edge weight of each community (updated, must be initialized)
 * @param x original graph
 * @param vtot total edge weight of each vertex
 * @param q initial community each vertex belongs to
 */
template <class G, class K, class W>
inline void leidenInitializeFromOmpW(vector<K>& vcom, vector<K>& vcob, vector<W>& ctot, const G& x, const vector<W>& vtot, const vector<K>& q) {
  size_t S = x.span();
  #pragma omp parallel for schedule(static, 2048)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    K c = q[u];
    vcom[u]  = u;
    vcob[u]  = c;
    #pragma omp atomic
    ctot[c] += vtot[u];
  }
}
#endif
#pragma endregion




#pragma region CHANGE COMMUNITY
/**
 * Scan an edge community connected to a vertex.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param u given vertex
 * @param v outgoing edge vertex
 * @param w outgoing edge weight
 * @param vcom community each vertex belongs to
 * @param vcob community bound each vertex belongs to
 */
template <bool SELF=false, bool REFINE=false, class K, class V, class W>
inline void leidenScanCommunityW(vector<K>& vcs, vector<W>& vcout, K u, K v, V w, const vector<K>& vcom, const vector<K>& vcob) {
  if (!SELF && u==v) return;
  if (REFINE && vcob[u]!=vcob[v]) return;
  K c = vcom[v];
  if (!vcout[c]) vcs.push_back(c);
  vcout[c] += w;
}


/**
 * Scan an edge community connected to a vertex.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param u given vertex
 * @param v outgoing edge vertex
 * @param w outgoing edge weight
 * @param vcom community each vertex belongs to
 */
template <bool SELF=false, class K, class V, class W>
inline void leidenScanCommunityW(vector<K>& vcs, vector<W>& vcout, K u, K v, V w, const vector<K>& vcom) {
  leidenScanCommunityW<SELF, false>(vcs, vcout, u, v, w, vcom, vcom);
}


/**
 * Scan communities connected to a vertex.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param x original graph
 * @param u given vertex
 * @param vcom community each vertex belongs to
 * @param vcob community bound each vertex belongs to
 */
template <bool SELF=false, bool REFINE=false, class G, class K, class W>
inline void leidenScanCommunitiesW(vector<K>& vcs, vector<W>& vcout, const G& x, K u, const vector<K>& vcom, const vector<K>& vcob) {
  x.forEachEdge(u, [&](auto v, auto w) { leidenScanCommunityW<SELF, REFINE>(vcs, vcout, u, v, w, vcom, vcob); });
}


/**
 * Scan communities connected to a vertex.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param x original graph
 * @param u given vertex
 * @param vcom community each vertex belongs to
 */
template <bool SELF=false, class G, class K, class W>
inline void leidenScanCommunitiesW(vector<K>& vcs, vector<W>& vcout, const G& x, K u, const vector<K>& vcom) {
  leidenScanCommunitiesW<SELF, false>(vcs, vcout, x, u, vcom, vcom);
}


/**
 * Clear communities scan data.
 * @param vcs total edge weight from vertex u to community C (updated)
 * @param vcout communities vertex u is linked to (updated)
 */
template <class K, class W>
inline void leidenClearScanW(vector<K>& vcs, vector<W>& vcout) {
  for (K c : vcs)
    vcout[c] = W();
  vcs.clear();
}


/**
 * Choose connected community with best delta modularity.
 * @param x original graph
 * @param u given vertex
 * @param vcom community each vertex belongs to
 * @param vtot total edge weight of each vertex
 * @param ctot total edge weight of each community
 * @param vcs communities vertex u is linked to
 * @param vcout total edge weight from vertex u to community C
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns [best community, delta modularity]
 */
template <bool SELF=false, class G, class K, class W>
inline auto leidenChooseCommunityGreedy(const G& x, K u, const vector<K>& vcom, const vector<W>& vtot, const vector<W>& ctot, const vector<K>& vcs, const vector<W>& vcout, double M, double R) {
  K cmax = K(), d = vcom[u];
  W emax = W();
  for (K c : vcs) {
    if (!SELF && c==d) continue;
    W e = deltaModularity(vcout[c], vcout[d], vtot[u], ctot[c], ctot[d], M, R);
    if (e>emax) { emax = e; cmax = c; }
  }
  return make_pair(cmax, emax);
}


/**
 * Choose connected community with best delta modularity.
 * @param rng random number generator
 * @param x original graph
 * @param u given vertex
 * @param vcom community each vertex belongs to
 * @param vtot total edge weight of each vertex
 * @param ctot total edge weight of each community
 * @param vcs communities vertex u is linked to
 * @param vcout total edge weight from vertex u to community C
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns [best community, delta modularity]
 */
template <bool SELF=false, class G, class K, class W>
inline auto leidenChooseCommunityRandom(xorshift32_engine& rng, const G& x, K u, const vector<K>& vcom, const vector<W>& vtot, const vector<W>& ctot, const vector<K>& vcs, const vector<W>& vcout, double M, double R) {
  K cmax = K(), d = vcom[u];
  W emax = W(), esum = W(), etil = W();
  for (K c : vcs) {
    if (!SELF && c==d) continue;
    W e = deltaModularity(vcout[c], vcout[d], vtot[u], ctot[c], ctot[d], M, R);
    if (e>0) esum += e;
  }
  W esel = ((rng() & 0xFFFF) / W(65536.0)) * esum;
  for (K c : vcs) {
    if (!SELF && c==d) continue;
    W e = deltaModularity(vcout[c], vcout[d], vtot[u], ctot[c], ctot[d], M, R);
    if (e>0) { etil += e; cmax = c; emax = e; }
    if (esel > etil) break;
  }
  return make_pair(cmax, emax);
}


/**
 * Choose connected community with best delta modularity.
 * @param rng random number generator
 * @param x original graph
 * @param u given vertex
 * @param vcom community each vertex belongs to
 * @param vtot total edge weight of each vertex
 * @param ctot total edge weight of each community
 * @param vcs communities vertex u is linked to
 * @param vcout total edge weight from vertex u to community C
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns [best community, delta modularity]
 */
template <bool SELF=false, bool RANDOM=false, class G, class K, class W>
inline auto leidenChooseCommunity(xorshift32_engine& rng, const G& x, K u, const vector<K>& vcom, const vector<W>& vtot, const vector<W>& ctot, const vector<K>& vcs, const vector<W>& vcout, double M, double R) {
  if (RANDOM) return leidenChooseCommunityRandom<SELF>(rng, x, u, vcom, vtot, ctot, vcs, vcout, M, R);
  else        return leidenChooseCommunityGreedy<SELF>     (x, u, vcom, vtot, ctot, vcs, vcout, M, R);
}


/**
 * Move vertex to another community C.
 * @param vcom community each vertex belongs to (updated)
 * @param ctot total edge weight of each community (updated)
 * @param x original graph
 * @param u given vertex
 * @param c community to move to
 * @param vtot total edge weight of each vertex
 */
template <class G, class K, class W>
inline void leidenChangeCommunityW(vector<K>& vcom, vector<W>& ctot, const G& x, K u, K c, const vector<W>& vtot) {
  K d = vcom[u];
  ctot[d] -= vtot[u];
  ctot[c] += vtot[u];
  vcom[u] = c;
}


#ifdef OPENMP
/**
 * Move vertex to another community C.
 * @param vcom community each vertex belongs to (updated)
 * @param ctot total edge weight of each community (updated)
 * @param x original graph
 * @param u given vertex
 * @param c community to move to
 * @param vtot total edge weight of each vertex
 */
template <class G, class K, class W>
inline void leidenChangeCommunityOmpW(vector<K>& vcom, vector<W>& ctot, const G& x, K u, K c, const vector<W>& vtot) {
  K d = vcom[u];
  #pragma omp atomic
  ctot[d] -= vtot[u];
  #pragma omp atomic
  ctot[c] += vtot[u];
  vcom[u] = c;
}
#endif
#pragma endregion




#pragma region LOCAL-MOVING PHASE
/**
 * Leiden algorithm's local moving phase.
 * @param vcom community each vertex belongs to (initial, updated)
 * @param ctot total edge weight of each community (precalculated, updated)
 * @param vaff is vertex affected flag (updated)
 * @param vcs communities vertex u is linked to (temporary buffer, updated)
 * @param vcout total edge weight from vertex u to community C (temporary buffer, updated)
 * @param rng random number generator
 * @param x original graph
 * @param vcob community bound each vertex belongs to
 * @param vtot total edge weight of each vertex
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @param L max iterations
 * @param fc has local moving phase converged?
 * @returns iterations performed (0 if converged already)
 */
template <bool REFINE=false, bool RANDOM=false, class G, class K, class W, class B, class FC>
inline int leidenMoveW(vector<K>& vcom, vector<W>& ctot, vector<B>& vaff, vector<K>& vcs, vector<W>& vcout, xorshift32_engine& rng, const G& x, const vector<K>& vcob, const vector<W>& vtot, double M, double R, int L, FC fc) {
  int l = 0;
  W  el = W();
  for (; l<L;) {
    el = W();
    x.forEachVertexKey([&](auto u) {
      if (!vaff[u]) return;
      leidenClearScanW(vcs, vcout);
      leidenScanCommunitiesW<false, REFINE>(vcs, vcout, x, u, vcom, vcob);
      auto [c, e] = leidenChooseCommunity<false, RANDOM>(rng, x, u, vcom, vtot, ctot, vcs, vcout, M, R);
      if (c)      { leidenChangeCommunityW(vcom, ctot, x, u, c, vtot); x.forEachEdgeKey(u, [&](auto v) { vaff[v] = B(1); }); }
      vaff[u] = B();
      el += e;  // l1-norm
    });
    if (fc(el, l++)) break;
  }
  return l>1 || el? l : 0;
}


#ifdef OPENMP
/**
 * Leiden algorithm's local moving phase.
 * @param vcom community each vertex belongs to (initial, updated)
 * @param ctot total edge weight of each community (precalculated, updated)
 * @param vaff is vertex affected flag (updated)
 * @param vcs communities vertex u is linked to (temporary buffer, updated)
 * @param vcout total edge weight from vertex u to community C (temporary buffer, updated)
 * @param rng random number generators
 * @param x original graph
 * @param vcob community bound each vertex belongs to
 * @param vtot total edge weight of each vertex
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @param L max iterations
 * @param fc has local moving phase converged?
 * @returns iterations performed (0 if converged already)
 */
template <bool REFINE=false, bool RANDOM=false, class G, class K, class W, class B, class FC>
inline int leidenMoveOmpW(vector<K>& vcom, vector<W>& ctot, vector<B>& vaff, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, vector<xorshift32_engine*>& rng, const G& x, const vector<K>& vcob, const vector<W>& vtot, double M, double R, int L, FC fc) {
  size_t S = x.span();
  int l = 0;
  W  el = W();
  for (; l<L;) {
    el = W();
    #pragma omp parallel for schedule(dynamic, 2048) reduction(+:el)
    for (K u=0; u<S; ++u) {
      int t = omp_get_thread_num();
      if (!x.hasVertex(u)) continue;
      if (!vaff[u]) continue;
      leidenClearScanW(*vcs[t], *vcout[t]);
      leidenScanCommunitiesW<false, REFINE>(*vcs[t], *vcout[t], x, u, vcom, vcob);
      auto [c, e] = leidenChooseCommunity<false, RANDOM>(*rng[t], x, u, vcom, vtot, ctot, *vcs[t], *vcout[t], M, R);
      if (c)      { leidenChangeCommunityOmpW(vcom, ctot, x, u, c, vtot); x.forEachEdgeKey(u, [&](auto v) { vaff[v] = B(1); }); }
      vaff[u] = B();
      el += e;  // l1-norm
    }
    if (fc(el, l++)) break;
  }
  return l>1 || el? l : 0;
}
#endif
#pragma endregion




#pragma region COMMUNITY PROPERTIES
/**
 * Examine if each community exists.
 * @param a does each community exist (updated)
 * @param x original graph
 * @param vcom community each vertex belongs to
 * @returns number of communities
 */
template <class G, class K, class A>
inline size_t leidenCommunityExistsW(vector<A>& a, const G& x, const vector<K>& vcom) {
  size_t C = 0;
  fillValueU(a, A());
  x.forEachVertexKey([&](auto u) {
    K c = vcom[u];
    if (!a[c]) ++C;
    a[c] = A(1);
  });
  return C;
}


#ifdef OPENMP
/**
 * Examine if each community exists.
 * @param a does each community exist (updated)
 * @param x original graph
 * @param vcom community each vertex belongs to
 * @returns number of communities
 */
template <class G, class K, class A>
inline size_t leidenCommunityExistsOmpW(vector<A>& a, const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  size_t C = 0;
  fillValueOmpU(a, A());
  #pragma omp parallel for schedule(static, 2048) reduction(+:C)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    K c = vcom[u];
    A m = A();
    #pragma omp atomic capture
    { m = a[c]; a[c] = A(1); }
    if (!m) ++C;
  }
  return C;
}
#endif




/**
 * Find the total degree of each community.
 * @param a total degree of each community (updated)
 * @param x original graph
 * @param vcom community each vertex belongs to
 */
template <class G, class K, class A>
inline void leidenCommunityTotalDegreeW(vector<A>& a, const G& x, const vector<K>& vcom) {
  fillValueU(a, A());
  x.forEachVertexKey([&](auto u) {
    K c   = vcom[u];
    a[c] += x.degree(u);
  });
}


#ifdef OPENMP
/**
 * Find the total degree of each community.
 * @param a total degree of each community (updated)
 * @param x original graph
 * @param vcom community each vertex belongs to
 */
template <class G, class K, class A>
inline void leidenCommunityTotalDegreeOmpW(vector<A>& a, const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  fillValueOmpU(a, A());
  #pragma omp parallel for schedule(static, 2048)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    K c = vcom[u];
    #pragma omp atomic
    a[c] += x.degree(u);
  }
}
#endif




/**
 * Find the number of vertices in each community.
 * @param a number of vertices belonging to each community (updated)
 * @param x original graph
 * @param vcom community each vertex belongs to
 */
template <class G, class K, class A>
inline void leidenCountCommunityVerticesW(vector<A>& a, const G& x, const vector<K>& vcom) {
  fillValueU(a, A());
  x.forEachVertexKey([&](auto u) {
    K c = vcom[u];
    ++a[c];
  });
}


#ifdef OPENMP
/**
 * Find the number of vertices in each community.
 * @param a number of vertices belonging to each community (updated)
 * @param x original graph
 * @param vcom community each vertex belongs to
 */
template <class G, class K, class A>
inline void leidenCountCommunityVerticesOmpW(vector<A>& a, const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  fillValueOmpU(a, A());
  #pragma omp parallel for schedule(static, 2048)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    K c = vcom[u];
    #pragma omp atomic
    ++a[c];
  }
}
#endif




/**
 * Find the vertices in each community.
 * @param coff csr offsets for vertices belonging to each community (updated)
 * @param cdeg number of vertices in each community (updated)
 * @param cedg vertices belonging to each community (updated)
 * @param x original graph
 * @param vcom community each vertex belongs to
 */
template <class G, class K>
inline void leidenCommunityVerticesW(vector<K>& coff, vector<K>& cdeg, vector<K>& cedg, const G& x, const vector<K>& vcom) {
  size_t C = coff.size() - 1;
  leidenCountCommunityVerticesW(coff, x, vcom);
  coff[C] = exclusiveScanW(coff.data(), coff.data(), C);
  fillValueU(cdeg, K());
  x.forEachVertexKey([&](auto u) {
    K c = vcom[u];
    csrAddEdgeU(cdeg, cedg, coff, c, u);
  });
}


#ifdef OPENMP
/**
 * Find the vertices in each community.
 * @param coff csr offsets for vertices belonging to each community (updated)
 * @param cdeg number of vertices in each community (updated)
 * @param cedg vertices belonging to each community (updated)
 * @param bufk buffer for exclusive scan of size |threads| (scratch)
 * @param x original graph
 * @param vcom community each vertex belongs to
 */
template <class G, class K>
inline void leidenCommunityVerticesOmpW(vector<K>& coff, vector<K>& cdeg, vector<K>& cedg, vector<K>& bufk, const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  size_t C = coff.size() - 1;
  leidenCountCommunityVerticesOmpW(coff, x, vcom);
  coff[C] = exclusiveScanOmpW(coff.data(), bufk.data(), coff.data(), C);
  fillValueOmpU(cdeg, K());
  #pragma omp parallel for schedule(static, 2048)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    K c = vcom[u];
    csrAddEdgeOmpU(cdeg, cedg, coff, c, u);
  }
}
#endif
#pragma endregion




#pragma region LOOKUP COMMUNITIES
/**
 * Update community membership in a tree-like fashion (to handle aggregation).
 * @param a output community each vertex belongs to (updated)
 * @param vcom community each vertex belongs to (at this aggregation level)
 */
template <class K>
inline void leidenLookupCommunitiesU(vector<K>& a, const vector<K>& vcom) {
  for (auto& v : a)
    v = vcom[v];
}


#ifdef OPENMP
/**
 * Update community membership in a tree-like fashion (to handle aggregation).
 * @param a output community each vertex belongs to (updated)
 * @param vcom community each vertex belongs to (at this aggregation level)
 */
template <class K>
inline void leidenLookupCommunitiesOmpU(vector<K>& a, const vector<K>& vcom) {
  size_t S = a.size();
  #pragma omp parallel for schedule(static, 2048)
  for (size_t u=0; u<S; ++u)
    a[u] = vcom[a[u]];
}
#endif
#pragma endregion




#pragma region AGGREGATION PHASE
/**
 * Aggregate outgoing edges of each community.
 * @param ydeg degree of each community (updated)
 * @param yedg vertex ids of outgoing edges of each community (updated)
 * @param ywei weights of outgoing edges of each community (updated)
 * @param vcs communities vertex u is linked to (temporary buffer, updated)
 * @param vcout total edge weight from vertex u to community C (temporary buffer, updated)
 * @param x original graph
 * @param vcom community each vertex belongs to
 * @param coff offsets for vertices belonging to each community
 * @param cedg vertices belonging to each community
 * @param yoff offsets for vertices belonging to each community
 */
template <class G, class K, class W>
inline void leidenAggregateEdgesW(vector<K>& ydeg, vector<K>& yedg, vector<W>& ywei, vector<K>& vcs, vector<W>& vcout, const G& x, const vector<K>& vcom, const vector<K>& coff, const vector<K>& cedg, const vector<size_t>& yoff) {
  size_t C = coff.size() - 1;
  fillValueU(ydeg, K());
  for (K c=0; c<C; ++c) {
    K n = csrDegree(coff, c);
    if (n==0) continue;
    leidenClearScanW(vcs, vcout);
    csrForEachEdgeKey(coff, cedg, c, [&](auto u) {
      leidenScanCommunitiesW<true>(vcs, vcout, x, u, vcom);
    });
    for (auto d : vcs)
      csrAddEdgeU(ydeg, yedg, ywei, yoff, c, d, vcout[d]);
  }
}


#ifdef OPENMP
/**
 * Aggregate outgoing edges of each community.
 * @param ydeg degree of each community (updated)
 * @param yedg vertex ids of outgoing edges of each community (updated)
 * @param ywei weights of outgoing edges of each community (updated)
 * @param vcs communities vertex u is linked to (temporary buffer, updated)
 * @param vcout total edge weight from vertex u to community C (temporary buffer, updated)
 * @param x original graph
 * @param vcom community each vertex belongs to
 * @param coff offsets for vertices belonging to each community
 * @param cedg vertices belonging to each community
 * @param yoff offsets for vertices belonging to each community
 */
template <class G, class K, class W>
inline void leidenAggregateEdgesOmpW(vector<K>& ydeg, vector<K>& yedg, vector<W>& ywei, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, const vector<K>& vcom, const vector<K>& coff, const vector<K>& cedg, const vector<size_t>& yoff) {
  size_t C = coff.size() - 1;
  fillValueOmpU(ydeg, K());
  #pragma omp parallel for schedule(dynamic, 2048)
  for (K c=0; c<C; ++c) {
    int t = omp_get_thread_num();
    K   n = csrDegree(coff, c);
    if (n==0) continue;
    leidenClearScanW(*vcs[t], *vcout[t]);
    csrForEachEdgeKey(coff, cedg, c, [&](auto u) {
      leidenScanCommunitiesW<true>(*vcs[t], *vcout[t], x, u, vcom);
    });
    for (auto d : *vcs[t])
      csrAddEdgeU(ydeg, yedg, ywei, yoff, c, d, (*vcout[t])[d]);
  }
}
#endif


/**
 * Re-number communities such that they are numbered 0, 1, 2, ...
 * @param vcom community each vertex belongs to (updated)
 * @param cext does each community exist (updated)
 * @param x original graph
 * @returns number of communities
 */
template <class G, class K>
inline size_t leidenRenumberCommunitiesW(vector<K>& vcom, vector<K>& cext, const G& x) {
  size_t C = exclusiveScanW(cext, cext);
  leidenLookupCommunitiesU(vcom, cext);
  return C;
}


#ifdef OPENMP
/**
 * Re-number communities such that they are numbered 0, 1, 2, ...
 * @param vcom community each vertex belongs to (updated)
 * @param cext does each community exist (updated)
 * @param bufk buffer for exclusive scan of size |threads| (scratch)
 * @param x original graph
 * @returns number of communities
 */
template <class G, class K>
inline size_t leidenRenumberCommunitiesOmpW(vector<K>& vcom, vector<K>& cext, vector<K>& bufk, const G& x) {
  size_t C = exclusiveScanOmpW(cext, bufk, cext);
  leidenLookupCommunitiesOmpU(vcom, cext);
  return C;
}
#endif


/**
 * Leiden algorithm's community aggregation phase.
 * @param yoff offsets for vertices belonging to each community (updated)
 * @param ydeg degree of each community (updated)
 * @param yedg vertex ids of outgoing edges of each community (updated)
 * @param ywei weights of outgoing edges of each community (updated)
 * @param vcs communities vertex u is linked to (temporary buffer, updated)
 * @param vcout total edge weight from vertex u to community C (temporary buffer, updated)
 * @param x original graph
 * @param vcom community each vertex belongs to
 * @param coff offsets for vertices belonging to each community
 * @param cedg vertices belonging to each community
 */
template <class G, class K, class W>
inline void leidenAggregateW(vector<size_t>& yoff, vector<K>& ydeg, vector<K>& yedg, vector<W>& ywei, vector<K>& vcs, vector<W>& vcout, const G& x, const vector<K>& vcom, vector<K>& coff, vector<K>& cedg) {
  size_t C = coff.size() - 1;
  leidenCommunityTotalDegreeW(yoff, x, vcom);
  yoff[C] = exclusiveScanW(yoff.data(), yoff.data(), C);
  leidenAggregateEdgesW(ydeg, yedg, ywei, vcs, vcout, x, vcom, coff, cedg, yoff);
}


#ifdef OPENMP
/**
 * Leiden algorithm's community aggregation phase.
 * @param yoff offsets for vertices belonging to each community (updated)
 * @param ydeg degree of each community (updated)
 * @param yedg vertex ids of outgoing edges of each community (updated)
 * @param ywei weights of outgoing edges of each community (updated)
 * @param bufs buffer for exclusive scan of size |threads| (scratch)
 * @param vcs communities vertex u is linked to (temporary buffer, updated)
 * @param vcout total edge weight from vertex u to community C (temporary buffer, updated)
 * @param x original graph
 * @param vcom community each vertex belongs to
 * @param coff offsets for vertices belonging to each community
 * @param cedg vertices belonging to each community
 */
template <class G, class K, class W>
inline void leidenAggregateOmpW(vector<size_t>& yoff, vector<K>& ydeg, vector<K>& yedg, vector<W>& ywei, vector<size_t>& bufs, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, const vector<K>& vcom, vector<K>& coff, vector<K>& cedg) {
  size_t C = coff.size() - 1;
  leidenCommunityTotalDegreeOmpW(yoff, x, vcom);
  yoff[C] = exclusiveScanOmpW(yoff.data(), bufs.data(), yoff.data(), C);
  leidenAggregateEdgesOmpW(ydeg, yedg, ywei, vcs, vcout, x, vcom, coff, cedg, yoff);
}
#endif
#pragma endregion




#pragma region ENVIRONMENT SETUP
/**
 * Setup and perform the Leiden algorithm.
 * @param rnd random number generator
 * @param x original graph
 * @param q initial community each vertex belongs to
 * @param o leiden options
 * @param fm marking affected vertices / preprocessing to be performed (vaff)
 * @returns leiden result
 */
template <bool RANDOM=false, class FLAG=char, class RND, class G, class K, class FM>
inline auto leidenInvoke(RND& rnd, const G& x, const vector<K> *q, const LeidenOptions& o, FM fm) {
  using  W = LEIDEN_WEIGHT_TYPE;
  using  B = FLAG;
  double R = o.resolution;
  int    L = o.maxIterations, l = 0;
  int    P = o.maxPasses, p = 0;
  size_t X = x.size();
  size_t S = x.span();
  double M = edgeWeight(x)/2;
  vector<B> vaff(S);
  vector<K> vcob(S), vcom(S), a(S);
  vector<W> vtot(S), ctot(S);
  vector<K> vcs;
  vector<W> vcout(S);
  vector<xorshift32_engine*> rng(1);
  leidenAllocateRngsW(rng, rnd);
  size_t Z = max(size_t(o.aggregationTolerance * X), X);
  size_t Y = max(size_t(o.aggregationTolerance * Z), Z);
  DiGraphCsr<K, None, None, K> cv(S, S);
  DiGraphCsr<K, None, W> y(S, Y);  // y(S, X)
  DiGraphCsr<K, None, W> z(S, Z);  // z(S, X)
  float tm = 0, tp = 0, tl = 0, ta = 0;
  float t  = measureDurationMarked([&](auto mark) {
    double E  = o.tolerance;
    auto   fc = [&](double el, int l) { return el<=E; };
    fillValueU(vcob, K());
    fillValueU(vcom, K());
    fillValueU(vtot, W());
    fillValueU(ctot, W());
    fillValueU(a, K());
    cv.respan(S);
    y .respan(S);
    z .respan(S);
    mark([&]() {
      tm += measureDuration([&]() { fm(vaff); });
      auto t0 = timeNow(), t1 = t0;
      leidenVertexWeightsW(vtot, x);
      if (q) leidenInitializeFromW(vcom, vcob, ctot, x, vtot, *q);
      else   leidenInitializeW    (vcom, vcob, ctot, x, vtot);
      for (l=0, p=0; M>0 && P>0;) {
        if (p==1) t1 = timeNow();
        bool isFirst = p==0;
        int m = 0;
        tl += measureDuration([&]() {
          if (isFirst) m += leidenMoveW<false, RANDOM>(vcob, ctot, vaff, vcs, vcout, *rng[0], x, vcob, vtot, M, R, L, fc);
          else         m += leidenMoveW<false, RANDOM>(vcob, ctot, vaff, vcs, vcout, *rng[0], y, vcob, vtot, M, R, L, fc);
          if (isFirst) fillValueOmpU(vaff.data(), x.order(), B(1));
          else         fillValueOmpU(vaff.data(), y.order(), B(1));
          if (isFirst) leidenInitializeCommunityWeightsW(ctot, x, vtot);
          else         leidenInitializeCommunityWeightsW(ctot, y, vtot);
          if (isFirst) m += leidenMoveW<true,  RANDOM>(vcom, ctot, vaff, vcs, vcout, *rng[0], x, vcob, vtot, M, R, L, fc);
          else         m += leidenMoveW<true,  RANDOM>(vcom, ctot, vaff, vcs, vcout, *rng[0], y, vcob, vtot, M, R, L, fc);
        });
        l += max(m, 1); ++p;
        if (m<=1 || p>=P) break;
        size_t GN = isFirst? x.order() : y.order();
        size_t GS = isFirst? x.span()  : y.span();
        size_t CN = 0;
        if (isFirst) CN = leidenCommunityExistsW(cv.degrees, x, vcom);
        else         CN = leidenCommunityExistsW(cv.degrees, y, vcom);
        if (double(CN)/GN >= o.aggregationTolerance) break;
        if (isFirst) leidenRenumberCommunitiesW(vcom, cv.degrees, x);
        else         leidenRenumberCommunitiesW(vcom, cv.degrees, y);
        if (isFirst) copyValuesW(a, vcom);
        else         leidenLookupCommunitiesU(a, vcom);
        cv.respan(CN); z.respan(CN);
        if (isFirst) leidenCommunityVerticesW(cv.offsets, cv.degrees, cv.edgeKeys, x, vcom);
        else         leidenCommunityVerticesW(cv.offsets, cv.degrees, cv.edgeKeys, y, vcom);
        ta += measureDuration([&]() {
          if (isFirst) leidenAggregateW(z.offsets, z.degrees, z.edgeKeys, z.edgeValues, vcs, vcout, x, vcom, cv.offsets, cv.edgeKeys);
          else         leidenAggregateW(z.offsets, z.degrees, z.edgeKeys, z.edgeValues, vcs, vcout, y, vcom, cv.offsets, cv.edgeKeys);
        });
        swap(y, z);
        // fillValueU(vcob.data(), CN, K());
        // fillValueU(vcom.data(), CN, K());
        // fillValueU(ctot.data(), CN, W());
        fillValueU(vtot.data(), CN, W());
        fillValueU(vaff.data(), CN, B(1));
        leidenVertexWeightsW(vtot, y);
        leidenInitializeW(vcom, vcob, ctot, y, vtot);
        E /= o.toleranceDrop;
      }
      if (p<=1) copyValuesW(a, vcom);
      else      leidenLookupCommunitiesU(a, vcom);
      if (p<=1) t1 = timeNow();
      tp += duration(t0, t1);
    });
  }, o.repeat);
  leidenFreeRngsW(rng);
  return LeidenResult<K>(a, l, p, t, tm/o.repeat, tp/o.repeat, tl/o.repeat, ta/o.repeat);
}


#ifdef OPENMP
/**
 * Setup and perform the Leiden algorithm.
 * @param rnd random number generator
 * @param x original graph
 * @param q initial community each vertex belongs to
 * @param o leiden options
 * @param fm marking affected vertices / preprocessing to be performed (vaff)
 * @returns leiden result
 */
template <bool RANDOM=false, class FLAG=char, class RND, class G, class K, class FM>
inline auto leidenInvokeOmp(RND& rnd, const G& x, const vector<K> *q, const LeidenOptions& o, FM fm) {
  using  W = LEIDEN_WEIGHT_TYPE;
  using  B = FLAG;
  double R = o.resolution;
  int    L = o.maxIterations, l = 0;
  int    P = o.maxPasses, p = 0;
  size_t X = x.size();
  size_t S = x.span();
  double M = edgeWeightOmp(x)/2;
  int    T = omp_get_max_threads();
  vector<B> vaff(S);
  vector<K> vcob(S), vcom(S), a(S);
  vector<W> vtot(S), ctot(S);
  vector<K> bufk(T);
  vector<size_t> bufs(T);
  vector<vector<K>*> vcs(T);
  vector<vector<W>*> vcout(T);
  vector<xorshift32_engine*> rng(T);
  leidenAllocateHashtablesW(vcs, vcout, S);
  leidenAllocateRngsW(rng, rnd);
  size_t Z = max(size_t(o.aggregationTolerance * X), X);
  size_t Y = max(size_t(o.aggregationTolerance * Z), Z);
  DiGraphCsr<K, None, None, K> cv(S, S);
  DiGraphCsr<K, None, W> y(S, Y);  // y(S, X)
  DiGraphCsr<K, None, W> z(S, Z);  // z(S, X)
  float tm = 0, tp = 0, tl = 0, ta = 0;
  float t  = measureDurationMarked([&](auto mark) {
    double E  = o.tolerance;
    auto   fc = [&](double el, int l) { return el<=E; };
    fillValueOmpU(vcob, K());
    fillValueOmpU(vcom, K());
    fillValueOmpU(vtot, W());
    fillValueOmpU(ctot, W());
    fillValueOmpU(a, K());
    cv.respan(S);
    y .respan(S);
    z .respan(S);
    mark([&]() {
      tm += measureDuration([&]() { fm(vaff); });
      auto t0 = timeNow(), t1 = t0;
      leidenVertexWeightsOmpW(vtot, x);
      if (q) leidenInitializeFromOmpW(vcom, vcob, ctot, x, vtot, *q);
      else   leidenInitializeOmpW    (vcom, vcob, ctot, x, vtot);
      for (l=0, p=0; M>0 && P>0;) {
        if (p==1) t1 = timeNow();
        bool isFirst = p==0;
        int m = 0;
        tl += measureDuration([&]() {
          if (isFirst) m += leidenMoveOmpW<false, RANDOM>(vcob, ctot, vaff, vcs, vcout, rng, x, vcob, vtot, M, R, L, fc);
          else         m += leidenMoveOmpW<false, RANDOM>(vcob, ctot, vaff, vcs, vcout, rng, y, vcob, vtot, M, R, L, fc);
          if (isFirst) fillValueOmpU(vaff.data(), x.order(), B(1));
          else         fillValueOmpU(vaff.data(), y.order(), B(1));
          if (isFirst) leidenInitializeCommunityWeightsOmpW(ctot, x, vtot);
          else         leidenInitializeCommunityWeightsOmpW(ctot, y, vtot);
          if (isFirst) m += leidenMoveOmpW<true,  RANDOM>(vcom, ctot, vaff, vcs, vcout, rng, x, vcob, vtot, M, R, L, fc);
          else         m += leidenMoveOmpW<true,  RANDOM>(vcom, ctot, vaff, vcs, vcout, rng, y, vcob, vtot, M, R, L, fc);
        });
        l += max(m, 1); ++p;
        if (m<=1 || p>=P) break;
        size_t GN = isFirst? x.order() : y.order();
        size_t GS = isFirst? x.span()  : y.span();
        size_t CN = 0;
        if (isFirst) CN = leidenCommunityExistsOmpW(cv.degrees, x, vcom);
        else         CN = leidenCommunityExistsOmpW(cv.degrees, y, vcom);
        if (double(CN)/GN >= o.aggregationTolerance) break;
        if (isFirst) leidenRenumberCommunitiesOmpW(vcom, cv.degrees, bufk, x);
        else         leidenRenumberCommunitiesOmpW(vcom, cv.degrees, bufk, y);
        if (isFirst) copyValuesOmpW(a, vcom);
        else         leidenLookupCommunitiesOmpU(a, vcom);
        cv.respan(CN); z.respan(CN);
        if (isFirst) leidenCommunityVerticesOmpW(cv.offsets, cv.degrees, cv.edgeKeys, bufk, x, vcom);
        else         leidenCommunityVerticesOmpW(cv.offsets, cv.degrees, cv.edgeKeys, bufk, y, vcom);
        ta += measureDuration([&]() {
          if (isFirst) leidenAggregateOmpW(z.offsets, z.degrees, z.edgeKeys, z.edgeValues, bufs, vcs, vcout, x, vcom, cv.offsets, cv.edgeKeys);
          else         leidenAggregateOmpW(z.offsets, z.degrees, z.edgeKeys, z.edgeValues, bufs, vcs, vcout, y, vcom, cv.offsets, cv.edgeKeys);
        });
        swap(y, z);
        // fillValueOmpU(vcob.data(), CN, K());
        // fillValueOmpU(vcom.data(), CN, K());
        // fillValueOmpU(ctot.data(), CN, W());
        fillValueOmpU(vtot.data(), CN, W());
        fillValueOmpU(vaff.data(), CN, B(1));
        leidenVertexWeightsOmpW(vtot, y);
        leidenInitializeOmpW(vcom, vcob, ctot, y, vtot);
        E /= o.toleranceDrop;
      }
      if (p<=1) copyValuesOmpW(a, vcom);
      else      leidenLookupCommunitiesOmpU(a, vcom);
      if (p<=1) t1 = timeNow();
      tp += duration(t0, t1);
    });
  }, o.repeat);
  leidenFreeHashtablesW(vcs, vcout);
  leidenFreeRngsW(rng);
  return LeidenResult<K>(a, l, p, t, tm/o.repeat, tp/o.repeat, tl/o.repeat, ta/o.repeat);
}
#endif
#pragma endregion




#pragma region STATIC/NAIVE-DYNAMIC APPROACH
/**
 * Obtain the community membership of each vertex with Static/Naive-dynamic Leiden.
 * @param rnd random number generator
 * @param x original graph
 * @param q initial community each vertex belongs to
 * @param o leiden options
 * @returns leiden result
 */
template <bool RANDOM=false, class FLAG=char, class RND, class G, class K>
inline auto leidenStatic(RND& rnd, const G& x, const vector<K>* q=nullptr, const LeidenOptions& o={}) {
  auto fm = [](auto& vaff) { fillValueU(vaff, FLAG(1)); };
  return leidenInvoke<RANDOM, FLAG>(rnd, x, q, o, fm);
}


#ifdef OPENMP
/**
 * Obtain the community membership of each vertex with Static/Naive-dynamic Leiden.
 * @param rnd random number generator
 * @param x original graph
 * @param q initial community each vertex belongs to
 * @param o leiden options
 * @returns leiden result
 */
template <bool RANDOM=false, class FLAG=char, class RND, class G, class K>
inline auto leidenStaticOmp(RND& rnd, const G& x, const vector<K>* q=nullptr, const LeidenOptions& o={}) {
  auto fm = [](auto& vaff) { fillValueOmpU(vaff, FLAG(1)); };
  return leidenInvokeOmp<RANDOM, FLAG>(rnd, x, q, o, fm);
}
#endif
#pragma endregion




#pragma region DYNAMIC FRONTIER APPROACH
/**
 * Find the vertices which should be processed upon a batch of edge insertions and deletions.
 * @param vertices vertex affected flags (output)
 * @param y updated graph
 * @param deletions edge deletions for this batch update (undirected)
 * @param insertions edge insertions for this batch update (undirected)
 * @param vcom community each vertex belongs to
 */
template <class B, class G, class K, class V>
inline void leidenAffectedVerticesFrontierW(vector<B>& vertices, const G& y, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom) {
  fillValueU(vertices, B());
  for (const auto& [u, v] : deletions) {
    if (vcom[u] != vcom[v]) continue;
    vertices[u]  = 1;
  }
  for (const auto& [u, v, w] : insertions) {
    if (vcom[u] == vcom[v]) continue;
    vertices[u]  = 1;
  }
}


#ifdef OPENMP
/**
 * Find the vertices which should be processed upon a batch of edge insertions and deletions.
 * @param vertices vertex affected flags (output)
 * @param y updated graph
 * @param deletions edge deletions for this batch update (undirected)
 * @param insertions edge insertions for this batch update (undirected)
 * @param vcom community each vertex belongs to
 */
template <class B, class G, class K, class V>
inline void leidenAffectedVerticesFrontierOmpW(vector<B>& vertices, const G& y, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom) {
  fillValueOmpU(vertices, B());
  size_t D = deletions.size();
  size_t I = insertions.size();
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<D; ++i) {
    K u = get<0>(deletions[i]);
    K v = get<1>(deletions[i]);
    if (vcom[u] != vcom[v]) continue;
    vertices[u]  = 1;
  }
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<I; ++i) {
    K u = get<0>(insertions[i]);
    K v = get<1>(insertions[i]);
    if (vcom[u] == vcom[v]) continue;
    vertices[u]  = 1;
  }
}
#endif




/**
 * Obtain the community membership of each vertex with Dynamic Frontier Leiden.
 * @param y updated graph
 * @param deletions edge deletions in batch update
 * @param insertions edge insertions in batch update
 * @param q initial communities
 * @param o leiden options
 * @returns leiden result
 */
template <class FLAG=char, class G, class K, class V>
inline auto leidenDynamicFrontier(const G& y, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const LeidenOptions& o={}) {
  const vector<K>& vcom = *q;
  auto fm = [&](auto& vaff) { leidenAffectedVerticesFrontierW(vaff, y, deletions, insertions, vcom); };
  return leidenInvoke<FLAG>(y, q, o, fm);
}


#ifdef OPENMP
/**
 * Obtain the community membership of each vertex with Dynamic Frontier Leiden.
 * @param y updated graph
 * @param deletions edge deletions in batch update
 * @param insertions edge insertions in batch update
 * @param q initial communities
 * @param o leiden options
 * @returns leiden result
 */
template <class FLAG=char, class G, class K, class V>
inline auto leidenDynamicFrontierOmp(const G& y, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const LeidenOptions& o={}) {
  const vector<K>& vcom = *q;
  auto fm = [&](auto& vaff) { leidenAffectedVerticesFrontierOmpW(vaff, y, deletions, insertions, vcom); };
  return leidenInvokeOmp<FLAG>(y, q, o, fm);
}
#endif
#pragma endregion
#pragma endregion
