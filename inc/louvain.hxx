#pragma once
#include <utility>
#include <tuple>
#include <vector>
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
using std::make_pair;
using std::move;
using std::swap;
using std::get;
using std::min;
using std::max;




#pragma region TYPES
/**
 * Options for Louvain algorithm.
 */
struct LouvainOptions {
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
   * Define options for Louvain algorithm.
   * @param repeat number of times to repeat the algorithm [1]
   * @param resolution resolution parameter for modularity [1]
   * @param tolerance tolerance for convergence [1e-2]
   * @param aggregationTolerance tolerance for aggregation [0.8]
   * @param toleranceDrop tolerance drop factor after each pass [10]
   * @param maxIterations maximum number of iterations per pass [20]
   * @param maxPasses maximum number of passes [10]
   */
  LouvainOptions(int repeat=1, double resolution=1, double tolerance=1e-2, double aggregationTolerance=0.8, double toleranceDrop=10, int maxIterations=20, int maxPasses=10) :
  repeat(repeat), resolution(resolution), tolerance(tolerance), aggregationTolerance(aggregationTolerance), toleranceDrop(toleranceDrop), maxIterations(maxIterations), maxPasses(maxPasses) {}
  #pragma endregion
};


/** Weight to be used in hashtable. */
#define LOUVAIN_WEIGHT_TYPE double




/**
 * Result of Louvain algorithm.
 * @tparam K key type (vertex-id)
 * @tparam W weight type
 */
template <class K, class W=LOUVAIN_WEIGHT_TYPE>
struct LouvainResult {
  #pragma region DATA
  /** Community membership each vertex belongs to. */
  vector<K> membership;
  /** Total edge weight of each vertex. */
  vector<W> vertexWeight;
  /** Total edge weight of each community. */
  vector<W> communityWeight;
  /** Number of iterations performed. */
  int iterations;
  /** Number of passes performed. */
  int passes;
  /** Time spent in milliseconds. */
  float time;
  /** Time spent in milliseconds for initial marking of affected vertices. */
  float markingTime;
  /** Time spent in initializing community memberships and total vertex/community weights. */
  float initializationTime;
  /** Time spent in milliseconds in first pass. */
  float firstPassTime;
  /** Time spent in milliseconds in local-moving phase. */
  float localMoveTime;
  /** Time spent in milliseconds in aggregation phase. */
  float aggregationTime;
  /** Number of vertices initially marked as affected. */
  size_t affectedVertices;
  #pragma endregion


  #pragma region CONSTRUCTORS
  /**
   * Result of Louvain algorithm.
   * @param membership community membership each vertex belongs to
   * @param vertexWeight total edge weight of each vertex
   * @param communityWeight total edge weight of each community
   * @param iterations number of iterations performed
   * @param passes number of passes performed
   * @param time time spent in milliseconds
   * @param markingTime time spent in milliseconds for initial marking of affected vertices
   * @param initializationTime time spent in initializing community memberships and total vertex/community weights
   * @param firstPassTime time spent in milliseconds in first pass
   * @param localMoveTime time spent in milliseconds in local-moving phase
   * @param aggregationTime time spent in milliseconds in aggregation phase
   * @param affectedVertices number of vertices initially marked as affected
   */
  LouvainResult(vector<K>&& membership, vector<W>&& vertexWeight, vector<W>&& communityWeight, int iterations=0, int passes=0, float time=0, float markingTime=0, float initializationTime=0, float firstPassTime=0, float localMoveTime=0, float aggregationTime=0, size_t affectedVertices=0) :
  membership(membership), vertexWeight(vertexWeight), communityWeight(communityWeight), iterations(iterations), passes(passes), time(time), markingTime(markingTime), initializationTime(initializationTime), firstPassTime(firstPassTime), localMoveTime(localMoveTime), aggregationTime(aggregationTime), affectedVertices(affectedVertices) {}


  /**
   * Result of Louvain algorithm.
   * @param membership community membership each vertex belongs to (moved)
   * @param vertexWeight total edge weight of each vertex (moved)
   * @param communityWeight total edge weight of each community (moved)
   * @param iterations number of iterations performed
   * @param passes number of passes performed
   * @param time time spent in milliseconds
   * @param markingTime time spent in milliseconds for initial marking of affected vertices
   * @param initializationTime time spent in initializing community memberships and total vertex/community weights
   * @param firstPassTime time spent in milliseconds in first pass
   * @param localMoveTime time spent in milliseconds in local-moving phase
   * @param aggregationTime time spent in milliseconds in aggregation phase
   * @param affectedVertices number of vertices initially marked as affected
   */
  LouvainResult(vector<K>& membership, vector<W>& vertexWeight, vector<W>& communityWeight, int iterations=0, int passes=0, float time=0, float markingTime=0, float initializationTime=0, float firstPassTime=0, float localMoveTime=0, float aggregationTime=0, size_t affectedVertices=0) :
  membership(move(membership)), vertexWeight(move(vertexWeight)), communityWeight(move(communityWeight)), iterations(iterations), passes(passes), time(time), markingTime(markingTime), initializationTime(initializationTime), firstPassTime(firstPassTime), localMoveTime(localMoveTime), aggregationTime(aggregationTime), affectedVertices(affectedVertices) {}
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
inline void louvainAllocateHashtablesW(vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, size_t S) {
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
inline void louvainFreeHashtablesW(vector<vector<K>*>& vcs, vector<vector<W>*>& vcout) {
  size_t N = vcs.size();
  for (size_t i=0; i<N; ++i) {
    delete vcs[i];
    delete vcout[i];
  }
}
#pragma endregion




#pragma region INITIALIZE
/**
 * Find the total edge weight of each vertex.
 * @param vtot total edge weight of each vertex (updated, must be initialized)
 * @param x original graph
 */
template <class G, class W>
inline void louvainVertexWeightsW(vector<W>& vtot, const G& x) {
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
inline void louvainVertexWeightsOmpW(vector<W>& vtot, const G& x) {
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
inline void louvainCommunityWeightsW(vector<W>& ctot, const G& x, const vector<K>& vcom, const vector<W>& vtot) {
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
inline void louvainCommunityWeightsOmpW(vector<W>& ctot, const G& x, const vector<K>& vcom, const vector<W>& vtot) {
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
 * Initialize communities such that each vertex is its own community.
 * @param vcom community each vertex belongs to (updated, must be initialized)
 * @param ctot total edge weight of each community (updated, must be initialized)
 * @param x original graph
 * @param vtot total edge weight of each vertex
 */
template <class G, class K, class W>
inline void louvainInitializeW(vector<K>& vcom, vector<W>& ctot, const G& x, const vector<W>& vtot) {
  x.forEachVertexKey([&](auto u) {
    vcom[u] = u;
    ctot[u] = vtot[u];
  });
}


#ifdef OPENMP
/**
 * Initialize communities such that each vertex is its own community.
 * @param vcom community each vertex belongs to (updated, must be initialized)
 * @param ctot total edge weight of each community (updated, must be initialized)
 * @param x original graph
 * @param vtot total edge weight of each vertex
 */
template <class G, class K, class W>
inline void louvainInitializeOmpW(vector<K>& vcom, vector<W>& ctot, const G& x, const vector<W>& vtot) {
  size_t S = x.span();
  #pragma omp parallel for schedule(static, 2048)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    vcom[u] = u;
    ctot[u] = vtot[u];
  }
}
#endif


/**
 * Initialize communities from given initial communities.
 * @param vcom community each vertex belongs to (updated, must be initialized)
 * @param ctot total edge weight of each community (updated, must be initialized)
 * @param y updated graph
 * @param vtot total edge weight of each vertex
 * @param q initial community each vertex belongs to
 */
template <class G, class K, class W>
inline void louvainInitializeFromW(vector<K>& vcom, vector<W>& ctot, const G& y, const vector<W>& vtot, const vector<K>& q) {
  y.forEachVertexKey([&](auto u) {
    K c = q[u];
    vcom[u]  = c;
    ctot[c] += vtot[u];
  });
}


#ifdef OPENMP
/**
 * Initialize communities from given initial communities.
 * @param vcom community each vertex belongs to (updated, must be initialized)
 * @param ctot total edge weight of each community (updated, must be initialized)
 * @param y updated graph
 * @param vtot total edge weight of each vertex
 * @param q initial community each vertex belongs to
 */
template <class G, class K, class W>
inline void louvainInitializeFromOmpW(vector<K>& vcom, vector<W>& ctot, const G& y, const vector<W>& vtot, const vector<K>& q) {
  size_t S = y.span();
  #pragma omp parallel for schedule(static, 2048)
  for (K u=0; u<S; ++u) {
    if (!y.hasVertex(u)) continue;
    K c = q[u];
    vcom[u]  = c;
    #pragma omp atomic
    ctot[c] += vtot[u];
  }
}
#endif


/**
 * Update weights using given edge deletions and insertions.
 * @param vtot total edge weight of each vertex (updated)
 * @param ctot total edge weight of each community (updated)
 * @param y updated graph
 * @param deletions edge deletions in batch update
 * @param insertions edge insertions in batch update
 * @param vcom community each vertex belongs to
 */
template <class G, class K, class V, class W>
inline void louvainUpdateWeightsFromU(vector<W>& vtot, vector<W>& ctot, const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom) {
  for (auto [u, v, w] : deletions) {
    K c = vcom[u];
    vtot[u] -= w;
    ctot[c] -= w;
  }
  for (auto [u, v, w] : insertions) {
    K c = vcom[u];
    vtot[u] += w;
    ctot[c] += w;
  }
}


#ifdef OPENMP
/**
 * Update weights using given edge deletions and insertions.
 * @param vtot total edge weight of each vertex (updated)
 * @param ctot total edge weight of each community (updated)
 * @param y updated graph
 * @param deletions edge deletions in batch update
 * @param insertions edge insertions in batch update
 * @param vcom community each vertex belongs to
 */
template <class G, class K, class V, class W>
inline void louvainUpdateWeightsFromOmpU(vector<W>& vtot, vector<W>& ctot, const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom) {
  #pragma omp parallel
  {
    for (auto [u, v, w] : deletions) {
      K c = vcom[u];
      if (belongsOmp(u)) vtot[u] -= w;
      if (belongsOmp(c)) ctot[c] -= w;
    }
    for (auto [u, v, w] : insertions) {
      K c = vcom[u];
      if (belongsOmp(u)) vtot[u] += w;
      if (belongsOmp(c)) ctot[c] += w;
    }
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
 */
template <bool SELF=false, class K, class V, class W>
inline void louvainScanCommunityW(vector<K>& vcs, vector<W>& vcout, K u, K v, V w, const vector<K>& vcom) {
  if (!SELF && u==v) return;
  K c = vcom[v];
  if (!vcout[c]) vcs.push_back(c);
  vcout[c] += w;
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
inline void louvainScanCommunitiesW(vector<K>& vcs, vector<W>& vcout, const G& x, K u, const vector<K>& vcom) {
  x.forEachEdge(u, [&](auto v, auto w) { louvainScanCommunityW<SELF>(vcs, vcout, u, v, w, vcom); });
}


/**
 * Clear communities scan data.
 * @param vcs total edge weight from vertex u to community C (updated)
 * @param vcout communities vertex u is linked to (updated)
 */
template <class K, class W>
inline void louvainClearScanW(vector<K>& vcs, vector<W>& vcout) {
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
inline auto louvainChooseCommunity(const G& x, K u, const vector<K>& vcom, const vector<W>& vtot, const vector<W>& ctot, const vector<K>& vcs, const vector<W>& vcout, double M, double R) {
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
 * Move vertex to another community C.
 * @param vcom community each vertex belongs to (updated)
 * @param ctot total edge weight of each community (updated)
 * @param x original graph
 * @param u given vertex
 * @param c community to move to
 * @param vtot total edge weight of each vertex
 */
template <class G, class K, class W>
inline void louvainChangeCommunityW(vector<K>& vcom, vector<W>& ctot, const G& x, K u, K c, const vector<W>& vtot) {
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
inline void louvainChangeCommunityOmpW(vector<K>& vcom, vector<W>& ctot, const G& x, K u, K c, const vector<W>& vtot) {
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
 * Louvain algorithm's local moving phase.
 * @param vcom community each vertex belongs to (initial, updated)
 * @param ctot total edge weight of each community (precalculated, updated)
 * @param vaff is vertex affected flag (updated)
 * @param vcs communities vertex u is linked to (temporary buffer, updated)
 * @param vcout total edge weight from vertex u to community C (temporary buffer, updated)
 * @param x original graph
 * @param vtot total edge weight of each vertex
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @param L max iterations
 * @param fc has local moving phase converged?
 * @param fa is vertex allowed to be updated?
 * @returns iterations performed (0 if converged already)
 */
template <class G, class K, class W, class B, class FC, class FA>
inline int louvainMoveW(vector<K>& vcom, vector<W>& ctot, vector<B>& vaff, vector<K>& vcs, vector<W>& vcout, const G& x, const vector<W>& vtot, double M, double R, int L, FC fc, FA fa) {
  int l = 0;
  W  el = W();
  for (; l<L;) {
    el = W();
    x.forEachVertexKey([&](auto u) {
      if (!fa(u) || !vaff[u]) return;
      louvainClearScanW(vcs, vcout);
      louvainScanCommunitiesW(vcs, vcout, x, u, vcom);
      auto [c, e] = louvainChooseCommunity(x, u, vcom, vtot, ctot, vcs, vcout, M, R);
      if (c)      { louvainChangeCommunityW(vcom, ctot, x, u, c, vtot); x.forEachEdgeKey(u, [&](auto v) { vaff[v] = B(1); }); }
      vaff[u] = B();
      el += e;  // l1-norm
    });
    if (fc(el, l++)) break;
  }
  return l>1 || el? l : 0;
}


/**
 * Louvain algorithm's local moving phase.
 * @param vcom community each vertex belongs to (initial, updated)
 * @param ctot total edge weight of each community (precalculated, updated)
 * @param vaff is vertex affected flag (updated)
 * @param vcs communities vertex u is linked to (temporary buffer, updated)
 * @param vcout total edge weight from vertex u to community C (temporary buffer, updated)
 * @param x original graph
 * @param vtot total edge weight of each vertex
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @param L max iterations
 * @param fc has local moving phase converged?
 * @returns iterations performed (0 if converged already)
 */
template <class G, class K, class W, class B, class FC>
inline int louvainMoveW(vector<K>& vcom, vector<W>& ctot, vector<B>& vaff, vector<K>& vcs, vector<W>& vcout, const G& x, const vector<W>& vtot, double M, double R, int L, FC fc) {
  auto fa = [](auto u) { return true; };
  return louvainMoveW(vcom, ctot, vaff, vcs, vcout, x, vtot, M, R, L, fc, fa);
}


#ifdef OPENMP
/**
 * Louvain algorithm's local moving phase.
 * @param vcom community each vertex belongs to (initial, updated)
 * @param ctot total edge weight of each community (precalculated, updated)
 * @param vaff is vertex affected flag (updated)
 * @param vcs communities vertex u is linked to (temporary buffer, updated)
 * @param vcout total edge weight from vertex u to community C (temporary buffer, updated)
 * @param x original graph
 * @param vtot total edge weight of each vertex
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @param L max iterations
 * @param fc has local moving phase converged?
 * @param fa is vertex allowed to be updated?
 * @returns iterations performed (0 if converged already)
 */
template <class G, class K, class W, class B, class FC, class FA>
inline int louvainMoveOmpW(vector<K>& vcom, vector<W>& ctot, vector<B>& vaff, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, const vector<W>& vtot, double M, double R, int L, FC fc, FA fa) {
  size_t S = x.span();
  int l = 0;
  W  el = W();
  for (; l<L;) {
    el = W();
    #pragma omp parallel for schedule(dynamic, 2048) reduction(+:el)
    for (K u=0; u<S; ++u) {
      int t = omp_get_thread_num();
      if (!x.hasVertex(u)) continue;
      if (!fa(u) || !vaff[u]) continue;
      louvainClearScanW(*vcs[t], *vcout[t]);
      louvainScanCommunitiesW(*vcs[t], *vcout[t], x, u, vcom);
      auto [c, e] = louvainChooseCommunity(x, u, vcom, vtot, ctot, *vcs[t], *vcout[t], M, R);
      if (c)      { louvainChangeCommunityOmpW(vcom, ctot, x, u, c, vtot); x.forEachEdgeKey(u, [&](auto v) { vaff[v] = B(1); }); }
      vaff[u] = B();
      el += e;  // l1-norm
    }
    if (fc(el, l++)) break;
  }
  return l>1 || el? l : 0;
}


/**
 * Louvain algorithm's local moving phase.
 * @param vcom community each vertex belongs to (initial, updated)
 * @param ctot total edge weight of each community (precalculated, updated)
 * @param vaff is vertex affected flag (updated)
 * @param vcs communities vertex u is linked to (temporary buffer, updated)
 * @param vcout total edge weight from vertex u to community C (temporary buffer, updated)
 * @param x original graph
 * @param vtot total edge weight of each vertex
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @param L max iterations
 * @param fc has local moving phase converged?
 * @returns iterations performed (0 if converged already)
 */
template <class G, class K, class W, class B, class FC>
inline int louvainMoveOmpW(vector<K>& vcom, vector<W>& ctot, vector<B>& vaff, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, const vector<W>& vtot, double M, double R, int L, FC fc) {
  auto fa = [](auto u) { return true; };
  return louvainMoveOmpW(vcom, ctot, vaff, vcs, vcout, x, vtot, M, R, L, fc, fa);
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
inline size_t louvainCommunityExistsW(vector<A>& a, const G& x, const vector<K>& vcom) {
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
inline size_t louvainCommunityExistsOmpW(vector<A>& a, const G& x, const vector<K>& vcom) {
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
inline void louvainCommunityTotalDegreeW(vector<A>& a, const G& x, const vector<K>& vcom) {
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
inline void louvainCommunityTotalDegreeOmpW(vector<A>& a, const G& x, const vector<K>& vcom) {
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
inline void louvainCountCommunityVerticesW(vector<A>& a, const G& x, const vector<K>& vcom) {
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
inline void louvainCountCommunityVerticesOmpW(vector<A>& a, const G& x, const vector<K>& vcom) {
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
inline void louvainCommunityVerticesW(vector<K>& coff, vector<K>& cdeg, vector<K>& cedg, const G& x, const vector<K>& vcom) {
  size_t C = coff.size() - 1;
  louvainCountCommunityVerticesW(coff, x, vcom);
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
inline void louvainCommunityVerticesOmpW(vector<K>& coff, vector<K>& cdeg, vector<K>& cedg, vector<K>& bufk, const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  size_t C = coff.size() - 1;
  louvainCountCommunityVerticesOmpW(coff, x, vcom);
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
inline void louvainLookupCommunitiesU(vector<K>& a, const vector<K>& vcom) {
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
inline void louvainLookupCommunitiesOmpU(vector<K>& a, const vector<K>& vcom) {
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
inline void louvainAggregateEdgesW(vector<K>& ydeg, vector<K>& yedg, vector<W>& ywei, vector<K>& vcs, vector<W>& vcout, const G& x, const vector<K>& vcom, const vector<K>& coff, const vector<K>& cedg, const vector<size_t>& yoff) {
  size_t C = coff.size() - 1;
  fillValueU(ydeg, K());
  for (K c=0; c<C; ++c) {
    K n = csrDegree(coff, c);
    if (n==0) continue;
    louvainClearScanW(vcs, vcout);
    csrForEachEdgeKey(coff, cedg, c, [&](auto u) {
      louvainScanCommunitiesW<true>(vcs, vcout, x, u, vcom);
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
inline void louvainAggregateEdgesOmpW(vector<K>& ydeg, vector<K>& yedg, vector<W>& ywei, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, const vector<K>& vcom, const vector<K>& coff, const vector<K>& cedg, const vector<size_t>& yoff) {
  size_t C = coff.size() - 1;
  fillValueOmpU(ydeg, K());
  #pragma omp parallel for schedule(dynamic, 2048)
  for (K c=0; c<C; ++c) {
    int t = omp_get_thread_num();
    K   n = csrDegree(coff, c);
    if (n==0) continue;
    louvainClearScanW(*vcs[t], *vcout[t]);
    csrForEachEdgeKey(coff, cedg, c, [&](auto u) {
      louvainScanCommunitiesW<true>(*vcs[t], *vcout[t], x, u, vcom);
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
inline size_t louvainRenumberCommunitiesW(vector<K>& vcom, vector<K>& cext, const G& x) {
  size_t C = exclusiveScanW(cext, cext);
  louvainLookupCommunitiesU(vcom, cext);
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
inline size_t louvainRenumberCommunitiesOmpW(vector<K>& vcom, vector<K>& cext, vector<K>& bufk, const G& x) {
  size_t C = exclusiveScanOmpW(cext, bufk, cext);
  louvainLookupCommunitiesOmpU(vcom, cext);
  return C;
}
#endif


/**
 * Louvain algorithm's community aggregation phase.
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
inline void louvainAggregateW(vector<size_t>& yoff, vector<K>& ydeg, vector<K>& yedg, vector<W>& ywei, vector<K>& vcs, vector<W>& vcout, const G& x, const vector<K>& vcom, vector<K>& coff, vector<K>& cedg) {
  size_t C = coff.size() - 1;
  louvainCommunityTotalDegreeW(yoff, x, vcom);
  yoff[C] = exclusiveScanW(yoff.data(), yoff.data(), C);
  louvainAggregateEdgesW(ydeg, yedg, ywei, vcs, vcout, x, vcom, coff, cedg, yoff);
}


#ifdef OPENMP
/**
 * Louvain algorithm's community aggregation phase.
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
inline void louvainAggregateOmpW(vector<size_t>& yoff, vector<K>& ydeg, vector<K>& yedg, vector<W>& ywei, vector<size_t>& bufs, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, const vector<K>& vcom, vector<K>& coff, vector<K>& cedg) {
  size_t C = coff.size() - 1;
  louvainCommunityTotalDegreeOmpW(yoff, x, vcom);
  yoff[C] = exclusiveScanOmpW(yoff.data(), bufs.data(), yoff.data(), C);
  louvainAggregateEdgesOmpW(ydeg, yedg, ywei, vcs, vcout, x, vcom, coff, cedg, yoff);
}
#endif
#pragma endregion




#pragma region ENVIRONMENT SETUP
/**
 * Setup and perform the Louvain algorithm.
 * @param x original graph
 * @param o louvain options
 * @param fi initializing community membership and total vertex/community weights (vcom, vtot, ctot)
 * @param fm marking affected vertices (vaff, vcs, vcout, vcom, vtot, ctot)
 * @param fa is vertex allowed to be updated? (u)
 * @returns louvain result
 */
template <bool DYNAMIC=false, class FLAG=char, class G, class FI, class FM, class FA>
inline auto louvainInvoke(const G& x, const LouvainOptions& o, FI fi, FM fm, FA fa) {
  using  K = typename G::key_type;
  using  W = LOUVAIN_WEIGHT_TYPE;
  using  B = FLAG;
  // Options.
  double R = o.resolution;
  int    L = o.maxIterations, l = 0;
  int    P = o.maxPasses, p = 0;
  // Get graph properties.
  size_t X = x.size();
  size_t S = x.span();
  double M = edgeWeight(x)/2;
  // Allocate buffers.
  vector<B> vaff(S);        // Affected vertex flag (any pass)
  vector<K> ucom, vcom(S);  // Community membership (first pass, current pass)
  vector<W> utot, vtot(S);  // Total vertex weights (first pass, current pass)
  vector<W> ctot;           // Total community weights (any pass)
  vector<K> vcs;       // Hashtable keys
  vector<W> vcout(S);  // Hashtable values
  if (!DYNAMIC) ucom.resize(S);
  if (!DYNAMIC) utot.resize(S);
  if (!DYNAMIC) ctot.resize(S);
  size_t Z = max(size_t(o.aggregationTolerance * X), X);
  size_t Y = max(size_t(o.aggregationTolerance * Z), Z);
  DiGraphCsr<K, None, None, K> cv(S, S);  // CSR for community vertices
  DiGraphCsr<K, None, W> y(S, Y);  // CSR for aggregated graph (input);  y(S, X)
  DiGraphCsr<K, None, W> z(S, Z);  // CSR for aggregated graph (output); z(S, X)
  // Perform Louvain algorithm.
  float tm = 0, ti = 0, tp = 0, tl = 0, ta = 0;  // Time spent in different phases
  float t  = measureDurationMarked([&](auto mark) {
    double E  = o.tolerance;
    auto   fc = [&](double el, int l) { return el<=E; };
    // Reset buffers, in case of multiple runs.
    fillValueU(vaff, B());
    fillValueU(ucom, K());
    fillValueU(vcom, K());
    fillValueU(utot, W());
    fillValueU(vtot, W());
    fillValueU(ctot, W());
    cv.respan(S);
    y .respan(S);
    z .respan(S);
    // Time the algorithm.
    mark([&]() {
      // Initialize community membership and total vertex/community weights.
      ti += measureDuration([&]() { fi(ucom, utot, ctot); });
      // Mark affected vertices.
      tm += measureDuration([&]() { fm(vaff, vcs, vcout, ucom, utot, ctot); });
      // Start timing first pass.
      auto t0 = timeNow(), t1 = t0;
      // Start local-moving, aggregation phases.
      // NOTE: In first pass, the input graph is a DiGraph.
      // NOTE: For subsequent passes, the input graph is a DiGraphCsr (optimization).
      for (l=0, p=0; M>0 && P>0;) {
        if (p==1) t1 = timeNow();
        bool isFirst = p==0;
        int m = 0;
        tl += measureDuration([&]() {
          if (isFirst) m = louvainMoveW(ucom, ctot, vaff, vcs, vcout, x, utot, M, R, L, fc, fa);
          else         m = louvainMoveW(vcom, ctot, vaff, vcs, vcout, y, vtot, M, R, L, fc);
        });
        l += max(m, 1); ++p;
        if (m<=1 || p>=P) break;
        size_t GN = isFirst? x.order() : y.order();
        size_t GS = isFirst? x.span()  : y.span();
        size_t CN = 0;
        if (isFirst) CN = louvainCommunityExistsW(cv.degrees, x, ucom);
        else         CN = louvainCommunityExistsW(cv.degrees, y, vcom);
        if (double(CN)/GN >= o.aggregationTolerance) break;
        if (isFirst) louvainRenumberCommunitiesW(ucom, cv.degrees, x);
        else         louvainRenumberCommunitiesW(vcom, cv.degrees, y);
        if (isFirst) {}
        else         louvainLookupCommunitiesU(ucom, vcom);
        cv.respan(CN); z.respan(CN);
        if (isFirst) louvainCommunityVerticesW(cv.offsets, cv.degrees, cv.edgeKeys, x, ucom);
        else         louvainCommunityVerticesW(cv.offsets, cv.degrees, cv.edgeKeys, y, vcom);
        ta += measureDuration([&]() {
          if (isFirst) louvainAggregateW(z.offsets, z.degrees, z.edgeKeys, z.edgeValues, vcs, vcout, x, ucom, cv.offsets, cv.edgeKeys);
          else         louvainAggregateW(z.offsets, z.degrees, z.edgeKeys, z.edgeValues, vcs, vcout, y, vcom, cv.offsets, cv.edgeKeys);
        });
        swap(y, z);
        // fillValueU(vcom.data(), CN, K());
        // fillValueU(ctot.data(), CN, W());
        fillValueU(vtot.data(), CN, W());
        fillValueU(vaff.data(), CN, B(1));
        louvainVertexWeightsW(vtot, y);
        louvainInitializeW(vcom, ctot, y, vtot);
        E /= o.toleranceDrop;
      }
      if (p<=1) {}
      else      louvainLookupCommunitiesU(ucom, vcom);
      if (p<=1) t1 = timeNow();
      tp += duration(t0, t1);
    });
  }, o.repeat);
  return LouvainResult<K, W>(ucom, utot, ctot, l, p, t, tm/o.repeat, ti/o.repeat, tp/o.repeat, tl/o.repeat, ta/o.repeat, countValue(vaff, B(1)));
}


#ifdef OPENMP
/**
 * Setup and perform the Louvain algorithm.
 * @param x original graph
 * @param o louvain options
 * @param fi initializing community membership and total vertex/community weights (vcom, vtot, ctot)
 * @param fm marking affected vertices (vaff, vcs, vcout, vcom, vtot, ctot)
 * @param fa is vertex allowed to be updated? (u)
 * @returns louvain result
 */
template <bool DYNAMIC=false, class FLAG=char, class G, class FI, class FM, class FA>
inline auto louvainInvokeOmp(const G& x, const LouvainOptions& o, FI fi, FM fm, FA fa) {
  using  K = typename G::key_type;
  using  W = LOUVAIN_WEIGHT_TYPE;
  using  B = FLAG;
  // Options.
  double R = o.resolution;
  int    L = o.maxIterations, l = 0;
  int    P = o.maxPasses, p = 0;
  // Get graph properties.
  size_t X = x.size();
  size_t S = x.span();
  double M = edgeWeightOmp(x)/2;
  // Allocate buffers.
  int    T = omp_get_max_threads();
  vector<B> vaff(S);        // Affected vertex flag (any pass)
  vector<K> ucom, vcom(S);  // Community membership (first pass, current pass)
  vector<W> utot, vtot(S);  // Total vertex weights (first pass, current pass)
  vector<W> ctot;           // Total community weights (any pass)
  vector<K> bufk(T);        // Buffer for exclusive scan
  vector<size_t> bufs(T);   // Buffer for exclusive scan
  vector<vector<K>*> vcs(T);    // Hashtable keys
  vector<vector<W>*> vcout(T);  // Hashtable values
  if (!DYNAMIC) ucom.resize(S);
  if (!DYNAMIC) utot.resize(S);
  if (!DYNAMIC) ctot.resize(S);
  louvainAllocateHashtablesW(vcs, vcout, S);
  size_t Z = max(size_t(o.aggregationTolerance * X), X);
  size_t Y = max(size_t(o.aggregationTolerance * Z), Z);
  DiGraphCsr<K, None, None, K> cv(S, S);  // CSR for community vertices
  DiGraphCsr<K, None, W> y(S, Y);         // CSR for aggregated graph (input);  y(S, X)
  DiGraphCsr<K, None, W> z(S, Z);         // CSR for aggregated graph (output); z(S, X)
  // Perform Louvain algorithm.
  float tm = 0, ti = 0, tp = 0, tl = 0, ta = 0;  // Time spent in different phases
  float t  = measureDurationMarked([&](auto mark) {
    double E  = o.tolerance;
    auto   fc = [&](double el, int l) { return el<=E; };
    // Reset buffers, in case of multiple runs.
    fillValueOmpU(vaff, B());
    fillValueOmpU(ucom, K());
    fillValueOmpU(vcom, K());
    fillValueOmpU(utot, W());
    fillValueOmpU(vtot, W());
    fillValueOmpU(ctot, W());
    cv.respan(S);
    y .respan(S);
    z .respan(S);
    // Time the algorithm.
    mark([&]() {
      // Initialize community membership and total vertex/community weights.
      ti += measureDuration([&]() { fi(ucom, utot, ctot); });
      // Mark affected vertices.
      tm += measureDuration([&]() { fm(vaff, vcs, vcout, ucom, utot, ctot); });
      // Start timing first pass.
      auto t0 = timeNow(), t1 = t0;
      // Start local-moving, aggregation phases.
      // NOTE: In first pass, the input graph is a DiGraph.
      // NOTE: For subsequent passes, the input graph is a DiGraphCsr (optimization).
      for (l=0, p=0; M>0 && P>0;) {
        if (p==1) t1 = timeNow();
        bool isFirst = p==0;
        int m = 0;
        tl += measureDuration([&]() {
          if (isFirst) m = louvainMoveOmpW(ucom, ctot, vaff, vcs, vcout, x, utot, M, R, L, fc, fa);
          else         m = louvainMoveOmpW(vcom, ctot, vaff, vcs, vcout, y, vtot, M, R, L, fc);
        });
        l += max(m, 1); ++p;
        if (m<=1 || p>=P) break;
        size_t GN = isFirst? x.order() : y.order();
        size_t GS = isFirst? x.span()  : y.span();
        size_t CN = 0;
        if (isFirst) CN = louvainCommunityExistsOmpW(cv.degrees, x, ucom);
        else         CN = louvainCommunityExistsOmpW(cv.degrees, y, vcom);
        if (double(CN)/GN >= o.aggregationTolerance) break;
        if (isFirst) louvainRenumberCommunitiesOmpW(ucom, cv.degrees, bufk, x);
        else         louvainRenumberCommunitiesOmpW(vcom, cv.degrees, bufk, y);
        if (isFirst) {}
        else         louvainLookupCommunitiesOmpU(ucom, vcom);
        cv.respan(CN); z.respan(CN);
        if (isFirst) louvainCommunityVerticesOmpW(cv.offsets, cv.degrees, cv.edgeKeys, bufk, x, ucom);
        else         louvainCommunityVerticesOmpW(cv.offsets, cv.degrees, cv.edgeKeys, bufk, y, vcom);
        ta += measureDuration([&]() {
          if (isFirst) louvainAggregateOmpW(z.offsets, z.degrees, z.edgeKeys, z.edgeValues, bufs, vcs, vcout, x, ucom, cv.offsets, cv.edgeKeys);
          else         louvainAggregateOmpW(z.offsets, z.degrees, z.edgeKeys, z.edgeValues, bufs, vcs, vcout, y, vcom, cv.offsets, cv.edgeKeys);
        });
        swap(y, z);
        // fillValueOmpU(vcom.data(), CN, K());
        // fillValueOmpU(ctot.data(), CN, W());
        fillValueOmpU(vtot.data(), CN, W());
        fillValueOmpU(vaff.data(), CN, B(1));
        louvainVertexWeightsOmpW(vtot, y);
        louvainInitializeOmpW(vcom, ctot, y, vtot);
        E /= o.toleranceDrop;
      }
      if (p<=1) {}
      else      louvainLookupCommunitiesOmpU(ucom, vcom);
      if (p<=1) t1 = timeNow();
      tp += duration(t0, t1);
    });
  }, o.repeat);
  louvainFreeHashtablesW(vcs, vcout);
  return LouvainResult<K, W>(ucom, utot, ctot, l, p, t, tm/o.repeat, ti/o.repeat, tp/o.repeat, tl/o.repeat, ta/o.repeat, countValueOmp(vaff, B(1)));
}
#endif
#pragma endregion




#pragma region REPEAT SETUP (DYNAMIC)
/**
 * Setup the Dynamic Louvain algorithm for multiple runs.
 * @param qs initial community membership for each run (updated)
 * @param qvtots initial total vertex weights for each run (updated)
 * @param qctots initial total community weights for each run (updated)
 * @param q initial community membership
 * @param qvtot initial total vertex weights
 * @param qctot initial total community weights
 * @param repeat number of runs
 */
template <class K, class W>
inline void louvainSetupInitialsW(vector2d<K>& qs, vector2d<W>& qvtots, vector2d<W>& qctots, const vector<K>& q, const vector<W>& qvtot, const vector<W>& qctot, int repeat) {
  qs    .resize(repeat);
  qvtots.resize(repeat);
  qctots.resize(repeat);
  for (int r=0; r<repeat; ++r) {
    qs[r]     = q;
    qvtots[r] = qvtot;
    qctots[r] = qctot;
  }
}
#pragma endregion




#pragma region STATIC APPROACH
/**
 * Obtain the community membership of each vertex with Static Louvain.
 * @param x original graph
 * @param o louvain options
 * @returns louvain result
 */
template <class FLAG=char, class G>
inline auto louvainStatic(const G& x, const LouvainOptions& o={}) {
  auto fi = [&](auto& vcom, auto& vtot, auto& ctot)  {
    louvainVertexWeightsW(vtot, x);
    louvainInitializeW(vcom, ctot, x, vtot);
  };
  auto fm = [ ](auto& vaff, const auto& vcom, const auto& vtot, const auto& ctot, auto& vcs,  auto& vcout) {
    fillValueU(vaff, FLAG(1));
  };
  auto fa = [ ](auto u) { return true; };
  return louvainInvoke<false, FLAG>(x, o, fi, fm, fa);
}


#ifdef OPENMP
/**
 * Obtain the community membership of each vertex with Static Louvain.
 * @param x original graph
 * @param o louvain options
 * @returns louvain result
 */
template <class FLAG=char, class G>
inline auto louvainStaticOmp(const G& x, const LouvainOptions& o={}) {
  auto fi = [&](auto& vcom, auto& vtot, auto& ctot)  {
    louvainVertexWeightsOmpW(vtot, x);
    louvainInitializeOmpW(vcom, ctot, x, vtot);
  };
  auto fm = [ ](auto& vaff, const auto& vcom, const auto& vtot, const auto& ctot, auto& vcs,  auto& vcout) {
    fillValueOmpU(vaff, FLAG(1));
  };
  auto fa = [ ](auto u) { return true; };
  return louvainInvokeOmp<false, FLAG>(x, o, fi, fm, fa);
}
#endif
#pragma endregion




#pragma region NAIVE-DYNAMIC APPROACH
/**
 * Obtain the community membership of each vertex with Naive-dynamic Louvain.
 * @param y updated graph
 * @param deletions edge deletions for this batch update (undirected)
 * @param insertions edge insertions for this batch update (undirected)
 * @param q initial community each vertex belongs to
 * @param qvtot initial total edge weight of each vertex
 * @param qctot initial total edge weight of each community
 * @param o louvain options
 * @returns louvain result
 */
template <class FLAG=char, class G, class K, class V, class W>
inline auto louvainNaiveDynamic(const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& q, const vector<W>& qvtot, const vector<W>& qctot, const LouvainOptions& o={}) {
  vector2d<K> qs;
  vector2d<W> qvtots, qctots;
  louvainSetupInitialsW(qs, qvtots, qctots, q, qvtot, qctot, o.repeat);
  int  r  = 0;
  auto fi = [&](auto& vcom, auto& vtot, auto& ctot)  {
    vcom = move(qs[r]);
    vtot = move(qvtots[r]);
    ctot = move(qctots[r]); ++r;
    louvainUpdateWeightsFromU(vtot, ctot, y, deletions, insertions, vcom);
  };
  auto fm = [ ](auto& vaff, const auto& vcom, const auto& vtot, const auto& ctot, auto& vcs,  auto& vcout) {
    fillValueU(vaff, FLAG(1));
  };
  auto fa = [ ](auto u) { return true; };
  return louvainInvoke<true, FLAG>(y, o, fi, fm, fa);
}


#ifdef OPENMP
/**
 * Obtain the community membership of each vertex with Naive-dynamic Louvain.
 * @param y updated graph
 * @param deletions edge deletions for this batch update (undirected)
 * @param insertions edge insertions for this batch update (undirected)
 * @param q initial community each vertex belongs to
 * @param qvtot initial total edge weight of each vertex
 * @param qctot initial total edge weight of each community
 * @param o louvain options
 * @returns louvain result
 */
template <class FLAG=char, class G, class K, class V, class W>
inline auto louvainNaiveDynamicOmp(const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& q, const vector<W>& qvtot, const vector<W>& qctot, const LouvainOptions& o={}) {
  vector2d<K> qs;
  vector2d<W> qvtots, qctots;
  louvainSetupInitialsW(qs, qvtots, qctots, q, qvtot, qctot, o.repeat);
  int  r  = 0;
  auto fi = [&](auto& vcom, auto& vtot, auto& ctot)  {
    vcom = move(qs[r]);
    vtot = move(qvtots[r]);
    ctot = move(qctots[r]); ++r;
    louvainUpdateWeightsFromOmpU(vtot, ctot, y, deletions, insertions, vcom);
  };
  auto fm = [ ](auto& vaff, const auto& vcom, const auto& vtot, const auto& ctot, auto& vcs,  auto& vcout) {
    fillValueOmpU(vaff, FLAG(1));
  };
  auto fa = [ ](auto u) { return true; };
  return louvainInvokeOmp<true, FLAG>(y, o, fi, fm, fa);
}
#endif
#pragma endregion




#pragma region DYNAMIC DELTA-SCREENING
/**
 * Find the vertices which should be processed upon a batch of edge insertions and deletions.
 * @param vertices vertex affected flags (output)
 * @param neighbors neighbor affected flags (output)
 * @param communities community affected flags (output)
 * @param vcs communities vertex u is linked to (temporary buffer, updated)
 * @param vcout total edge weight from vertex u to community C (temporary buffer, updated)
 * @param y updated graph
 * @param deletions edge deletions for this batch update (undirected, sorted by source vertex id)
 * @param insertions edge insertions for this batch update (undirected, sorted by source vertex id)
 * @param vcom community each vertex belongs to
 * @param vtot total edge weight of each vertex
 * @param ctot total edge weight of each community
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 */
template <class B, class G, class K, class V, class W>
inline auto louvainAffectedVerticesDeltaScreeningW(vector<B>& vertices, vector<B>& neighbors, vector<B>& communities, vector<K>& vcs, vector<W>& vcout, const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom, const vector<W>& vtot, const vector<W>& ctot, double M, double R=1) {
  fillValueU(vertices,    B());
  fillValueU(neighbors,   B());
  fillValueU(communities, B());
  for (const auto& [u, v] : deletions) {
    if (vcom[u] != vcom[v]) continue;
    vertices[u]  = 1;
    neighbors[u] = 1;
    communities[vcom[v]] = 1;
  }
  for (size_t i=0; i<insertions.size();) {
    K u = get<0>(insertions[i]);
    louvainClearScanW(vcs, vcout);
    for (; i<insertions.size() && get<0>(insertions[i])==u; ++i) {
      K v = get<1>(insertions[i]);
      V w = get<2>(insertions[i]);
      if (vcom[u] == vcom[v]) continue;
      louvainScanCommunityW(vcs, vcout, u, v, w, vcom);
    }
    auto [c, e] = louvainChooseCommunity(y, u, vcom, vtot, ctot, vcs, vcout, M, R);
    if (e<=0) continue;
    vertices[u]  = 1;
    neighbors[u] = 1;
    communities[c] = 1;
  }
  y.forEachVertexKey([&](auto u) {
    if (neighbors[u]) y.forEachEdgeKey(u, [&](auto v) { vertices[v] = 1; });
    if (communities[vcom[u]]) vertices[u] = 1;
  });
}


#ifdef OPENMP
/**
 * Find the vertices which should be processed upon a batch of edge insertions and deletions.
 * @param vertices vertex affected flags (output)
 * @param neighbors neighbor affected flags (output)
 * @param communities community affected flags (output)
 * @param vcs communities vertex u is linked to (temporary buffer, updated)
 * @param vcout total edge weight from vertex u to community C (temporary buffer, updated)
 * @param y updated graph
 * @param deletions edge deletions for this batch update (undirected, sorted by source vertex id)
 * @param insertions edge insertions for this batch update (undirected, sorted by source vertex id)
 * @param vcom community each vertex belongs to
 * @param vtot total edge weight of each vertex
 * @param ctot total edge weight of each community
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 */
template <class B, class G, class K, class V, class W>
inline auto louvainAffectedVerticesDeltaScreeningOmpW(vector<B>& vertices, vector<B>& neighbors, vector<B>& communities, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom, const vector<W>& vtot, const vector<W>& ctot, double M, double R=1) {
  size_t S = y.span();
  size_t D = deletions.size();
  size_t I = insertions.size();
  fillValueOmpU(vertices,    B());
  fillValueOmpU(neighbors,   B());
  fillValueOmpU(communities, B());
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<D; ++i) {
    K u = get<0>(deletions[i]);
    K v = get<1>(deletions[i]);
    if (vcom[u] != vcom[v]) continue;
    vertices[u]  = 1;
    neighbors[u] = 1;
    communities[vcom[v]] = 1;
  }
  #pragma omp parallel
  {
    int T = omp_get_num_threads();
    int t = omp_get_thread_num();
    K  u0 = I>0? get<0>(insertions[0]) : 0;
    for (size_t i=0, n=0; i<I;) {
      K u = get<0>(insertions[i]);
      if (u!=u0) { ++n; u0 = u; }
      if (n % T != t) { ++i; continue; }
      louvainClearScanW(*vcs[t], *vcout[t]);
      for (; i<I && get<0>(insertions[i])==u; ++i) {
        K v = get<1>(insertions[i]);
        V w = get<2>(insertions[i]);
        if (vcom[u] == vcom[v]) continue;
        louvainScanCommunityW(*vcs[t], *vcout[t], u, v, w, vcom);
      }
      auto [c, e] = louvainChooseCommunity(y, u, vcom, vtot, ctot, *vcs[t], *vcout[t], M, R);
      if (e<=0) continue;
      vertices[u]  = 1;
      neighbors[u] = 1;
      communities[c] = 1;
    }
  }
  #pragma omp parallel for schedule(auto)
  for (K u=0; u<S; ++u) {
    if (!y.hasVertex(u)) continue;
    if (neighbors[u]) y.forEachEdgeKey(u, [&](auto v) { vertices[v] = 1; });
    if (communities[vcom[u]]) vertices[u] = 1;
  }
}
#endif




/**
 * Obtain the community membership of each vertex with Dynamic Delta-screening Louvain.
 * @param y updated graph
 * @param deletions edge deletions in batch update (undirected, sorted by source vertex id)
 * @param insertions edge insertions in batch update (undirected, sorted by source vertex id)
 * @param q initial community each vertex belongs to
 * @param qvtot initial total edge weight of each vertex
 * @param qctot initial total edge weight of each community
 * @param o louvain options
 * @returns louvain result
 */
template <class FLAG=char, class G, class K, class V, class W>
inline auto louvainDynamicDeltaScreening(const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& q, const vector<W>& qvtot, const vector<W>& qctot, const LouvainOptions& o={}) {
  using  B = FLAG;
  size_t S = y.span();
  double R = o.resolution;
  double M = edgeWeight(y)/2;
  vector<B> vertices(S), neighbors(S), communities(S);
  vector2d<K> qs;
  vector2d<W> qvtots, qctots;
  louvainSetupInitialsW(qs, qvtots, qctots, q, qvtot, qctot, o.repeat);
  int  r  = 0;
  auto fi = [&](auto& vcom, auto& vtot, auto& ctot) {
    vcom = move(qs[r]);
    vtot = move(qvtots[r]);
    ctot = move(qctots[r]); ++r;
    louvainUpdateWeightsFromU(vtot, ctot, y, deletions, insertions, vcom);
  };
  auto fm = [&](auto& vaff, auto& vcs, auto& vcout, const auto& vcom, const auto& vtot, const auto& ctot) {
    louvainAffectedVerticesDeltaScreeningW(vertices, neighbors, communities, vcs, vcout, y, deletions, insertions, vcom, vtot, ctot, M, R);
    copyValuesW(vaff, vertices);
  };
  auto fa = [&](auto u) { return vertices[u] == B(1); };
  return louvainInvoke<true, FLAG>(y, o, fi, fm, fa);
}


#ifdef OPENMP
/**
 * Obtain the community membership of each vertex with Dynamic Delta-screening Louvain.
 * @param y updated graph
 * @param deletions edge deletions in batch update
 * @param insertions edge insertions in batch update
 * @param q initial community each vertex belongs to
 * @param qvtot initial total edge weight of each vertex
 * @param qctot initial total edge weight of each community
 * @param o louvain options
 * @returns louvain result
 */
template <class FLAG=char, class G, class K, class V, class W>
inline auto louvainDynamicDeltaScreeningOmp(const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& q, const vector<W>& qvtot, const vector<W>& qctot, const LouvainOptions& o={}) {
  using  B = FLAG;
  size_t S = y.span();
  double R = o.resolution;
  double M = edgeWeightOmp(y)/2;
  int    T = omp_get_max_threads();
  vector<B> vertices(S), neighbors(S), communities(S);
  vector2d<K> qs;
  vector2d<W> qvtots, qctots;
  louvainSetupInitialsW(qs, qvtots, qctots, q, qvtot, qctot, o.repeat);
  int  r  = 0;
  auto fi = [&](auto& vcom, auto& vtot, auto& ctot) {
    vcom = move(qs[r]);
    vtot = move(qvtots[r]);
    ctot = move(qctots[r]); ++r;
    louvainUpdateWeightsFromOmpU(vtot, ctot, y, deletions, insertions, vcom);
  };
  auto fm = [&](auto& vaff, auto& vcs, auto& vcout, const auto& vcom, const auto& vtot, const auto& ctot) {
    louvainAffectedVerticesDeltaScreeningOmpW(vertices, neighbors, communities, vcs, vcout, y, deletions, insertions, vcom, vtot, ctot, M, R);
    copyValuesOmpW(vaff, vertices);
  };
  auto fa = [&](auto u) { return vertices[u] == B(1); };
  return louvainInvokeOmp<true, FLAG>(y, o, fi, fm, fa);
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
inline void louvainAffectedVerticesFrontierW(vector<B>& vertices, const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom) {
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
 * @returns flags for each vertex marking whether it is affected
 */
template <class B, class G, class K, class V>
inline void louvainAffectedVerticesFrontierOmpW(vector<B>& vertices, const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom) {
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
 * Obtain the community membership of each vertex with Dynamic Frontier Louvain.
 * @param y updated graph
 * @param deletions edge deletions in batch update
 * @param insertions edge insertions in batch update
 * @param q initial community each vertex belongs to
 * @param qvtot initial total edge weight of each vertex
 * @param qctot initial total edge weight of each community
 * @param o louvain options
 * @returns louvain result
 */
template <class FLAG=char, class G, class K, class V, class W>
inline auto louvainDynamicFrontier(const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& q, const vector<W>& qvtot, const vector<W>& qctot, const LouvainOptions& o={}) {
  vector2d<K> qs;
  vector2d<W> qvtots, qctots;
  louvainSetupInitialsW(qs, qvtots, qctots, q, qvtot, qctot, o.repeat);
  int  r  = 0;
  auto fi = [&](auto& vcom, auto& vtot, auto& ctot) {
    vcom = move(qs[r]);
    vtot = move(qvtots[r]);
    ctot = move(qctots[r]); ++r;
    louvainUpdateWeightsFromU(vtot, ctot, y, deletions, insertions, vcom);
  };
  auto fm = [&](auto& vaff, auto& vcs, auto& vcout, const auto& vcom, const auto& vtot, const auto& ctot) {
    louvainAffectedVerticesFrontierW(vaff, y, deletions, insertions, vcom);
  };
  auto fa = [ ](auto u) { return true; };
  return louvainInvoke<true, FLAG>(y, o, fi, fm, fa);
}


#ifdef OPENMP
/**
 * Obtain the community membership of each vertex with Dynamic Frontier Louvain.
 * @param y updated graph
 * @param deletions edge deletions in batch update
 * @param insertions edge insertions in batch update
 * @param q initial community each vertex belongs to
 * @param qvtot initial total edge weight of each vertex
 * @param qctot initial total edge weight of each community
 * @param o louvain options
 * @returns louvain result
 */
template <class FLAG=char, class G, class K, class V, class W>
inline auto louvainDynamicFrontierOmp(const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& q, const vector<W>& qvtot, const vector<W>& qctot, const LouvainOptions& o={}) {
  vector2d<K> qs;
  vector2d<W> qvtots, qctots;
  louvainSetupInitialsW(qs, qvtots, qctots, q, qvtot, qctot, o.repeat);
  int  r  = 0;
  auto fi = [&](auto& vcom, auto& vtot, auto& ctot) {
    vcom = move(qs[r]);
    vtot = move(qvtots[r]);
    ctot = move(qctots[r]); ++r;
    louvainUpdateWeightsFromOmpU(vtot, ctot, y, deletions, insertions, vcom);
  };
  auto fm = [&](auto& vaff, auto& vcs, auto& vcout, const auto& vcom, const auto& vtot, const auto& ctot) {
    louvainAffectedVerticesFrontierOmpW(vaff, y, deletions, insertions, vcom);
  };
  auto fa = [ ](auto u) { return true; };
  return louvainInvokeOmp<true, FLAG>(y, o, fi, fm, fa);
}
#endif
#pragma endregion
#pragma endregion
