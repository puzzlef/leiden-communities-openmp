#pragma once
#include <utility>
#include <tuple>
#include <vector>
#include <cstdint>
#include "_main.hxx"
#ifdef OPENMP
#include <omp.h>
#endif

using std::pair;
using std::tuple;
using std::vector;
using std::make_pair;
using std::move;
using std::get;




#pragma region TYPES
/**
 * Options for RAK algorithm.
 */
struct RakOptions {
  #pragma region DATA
  /** Number of times to repeat the algorithm [1]. */
  int repeat;
  /** Tolerance for convergence [0.05]. */
  double tolerance;
  /** Maximum number of iterations [20]. */
  int maxIterations;
  #pragma endregion


  #pragma region CONSTRUCTORS
  /**
   * Define options for RAK algorithm.
   * @param repeat number of times to repeat the algorithm [1]
   * @param tolerance tolerance for convergence [0.05]
   * @param maxIterations maximum number of iterations [20]
   */
  RakOptions(int repeat=1, double tolerance=0.05, int maxIterations=20) :
  repeat(repeat), tolerance(tolerance), maxIterations(maxIterations) {}
  #pragma endregion
};


/** Weight to be used in hashtable. */
#define RAK_WEIGHT_TYPE double




/**
 * Result of RAK algorithm.
 * @tparam K key type (vertex-id)
 */
template <class K>
struct RakResult {
  #pragma region DATA
  /** Community membership each vertex belongs to. */
  vector<K> membership;
  /** Number of iterations performed. */
  int iterations;
  /** Time spent in milliseconds. */
  float time;
  /** Time spent in milliseconds for preprocessing. */
  float preprocessingTime;
  #pragma endregion


  #pragma region CONSTRUCTORS
  /**
   * Result of RAK algorithm.
   * @param membership community membership each vertex belongs to
   * @param iterations number of iterations performed
   * @param time time spent in milliseconds
   * @param preprocessingTime time spent in milliseconds for preprocessing
   */
  RakResult(vector<K>&& membership, int iterations=0, float time=0, float preprocessingTime=0) :
  membership(membership), iterations(iterations), time(time), preprocessingTime(preprocessingTime) {}


  /**
   * Result of RAK algorithm.
   * @param membership community membership each vertex belongs to (moved)
   * @param iterations number of iterations performed
   * @param time time spent in milliseconds
   * @param preprocessingTime time spent in milliseconds for preprocessing
   */
  RakResult(vector<K>& membership, int iterations=0, float time=0, float preprocessingTime=0) :
  membership(move(membership)), iterations(iterations), time(time), preprocessingTime(preprocessingTime) {}
  #pragma endregion
};
#pragma endregion




#pragma region METHODS
#pragma region HASHTABLES
/**
 * Allocate a number of hashtables.
 * @param vcs communities vertex u is linked to (temporary buffer, updated)
 * @param vcout total edge weight from vertex u to community C (temporary buffer, updated)
 * @param S size of each hashtable
 */
template <class K, class W>
inline void rakAllocateHashtablesW(vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, size_t S) {
  size_t N = vcs.size();
  for (size_t i=0; i<N; ++i) {
    vcs[i]   = new vector<K>();
    vcout[i] = new vector<W>(S);
  }
}


/**
 * Free a number of hashtables.
 * @param vcs communities vertex u is linked to (temporary buffer, updated)
 * @param vcout total edge weight from vertex u to community C (temporary buffer, updated)
 */
template <class K, class W>
inline void rakFreeHashtablesW(vector<vector<K>*>& vcs, vector<vector<W>*>& vcout) {
  size_t N = vcs.size();
  for (size_t i=0; i<N; ++i) {
    delete vcs[i];
    delete vcout[i];
  }
}
#pragma endregion




#pragma region INITIALIZE
/**
 * Initialize communities such that each vertex is its own community.
 * @param vcom community each vertex belongs to (updated, must be initialized)
 * @param x original graph
 */
template <class G, class K>
inline void rakInitializeW(vector<K>& vcom, const G& x) {
  x.forEachVertexKey([&](auto u) { vcom[u] = u; });
}


#ifdef OPENMP
/**
 * Initialize communities such that each vertex is its own community.
 * @param vcom community each vertex belongs to (updated, must be initialized)
 * @param x original graph
 */
template <class G, class K>
inline void rakInitializeOmpW(vector<K>& vcom, const G& x) {
  size_t S = x.span();
  #pragma omp parallel for schedule(auto)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    vcom[u] = u;
  }
}
#endif


/**
 * Initialize communities from given initial communities.
 * @param vcom community each vertex belongs to (updated, must be initialized)
 * @param x original graph
 * @param q initial community each vertex belongs to
 */
template <class G, class K>
inline void rakInitializeFromW(vector<K>& vcom, const G& x, const vector<K>& q) {
  x.forEachVertexKey([&](auto u) { vcom[u] = q[u]; });
}


#ifdef OPENMP
/**
 * Initialize communities from given initial communities.
 * @param vcom community each vertex belongs to (updated, must be initialized)
 * @param x original graph
 * @param q initial community each vertex belongs to
 */
template <class G, class K>
inline void rakInitializeFromOmpW(vector<K>& vcom, const G& x, const vector<K>& q) {
  size_t S = x.span();
  #pragma omp parallel for schedule(auto)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    vcom[u] = q[u];
  }
}
#endif
#pragma endregion




#pragma region CHOOSE COMMUNITY
/**
 * Scan an edge community connected to a vertex.
 * @param vcs communities vertex u is linked to (temporary buffer, updated)
 * @param vcout total edge weight from vertex u to community C (temporary buffer, updated)
 * @param u given vertex
 * @param v outgoing edge vertex
 * @param w outgoing edge weight
 * @param vcom community each vertex belongs to
 */
template <bool SELF=false, class K, class V, class W>
inline void rakScanCommunityW(vector<K>& vcs, vector<W>& vcout, K u, K v, V w, const vector<K>& vcom) {
  if (!SELF && u==v) return;
  K c = vcom[v];
  if (!vcout[c]) vcs.push_back(c);
  vcout[c] += w;
}


/**
 * Scan communities connected to a vertex.
 * @param vcs communities vertex u is linked to (temporary buffer, updated)
 * @param vcout total edge weight from vertex u to community C (temporary buffer, updated)
 * @param x original graph
 * @param u given vertex
 * @param vcom community each vertex belongs to
 */
template <bool SELF=false, class G, class K, class W>
inline void rakScanCommunitiesW(vector<K>& vcs, vector<W>& vcout, const G& x, K u, const vector<K>& vcom) {
  x.forEachEdge(u, [&](auto v, auto w) { rakScanCommunityW<SELF>(vcs, vcout, u, v, w, vcom); });
}


/**
 * Clear communities scan data.
 * @param vcs total edge weight from vertex u to community C (temporary buffer, updated)
 * @param vcout communities vertex u is linked to (updated)
 */
template <class K, class W>
inline void rakClearScanW(vector<K>& vcs, vector<W>& vcout) {
  for (K c : vcs)
    vcout[c] = W();
  vcs.clear();
}


/**
 * Choose connected community with most weight.
 * @param x original graph
 * @param u given vertex
 * @param vcom community each vertex belongs to
 * @param vcs communities vertex u is linked to
 * @param vcout total edge weight from vertex u to community C
 * @returns [best community, best edge weight to community]
 */
template <class G, class K, class W>
inline pair<K, W> rakChooseCommunity(const G& x, K u, const vector<K>& vcom, const vector<K>& vcs, const vector<W>& vcout) {
  K d = vcom[u];
  K cmax = K();
  W wmax = W();
  for (K c : vcs)
    if (vcout[c]>wmax) { cmax = c; wmax = vcout[c]; }
  return make_pair(cmax, wmax);
}
#pragma endregion




#pragma region MOVE ITERATION
/**
 * Move each vertex to its best community.
 * @param vcom community each vertex belongs to (updated)
 * @param vaff is vertex affected (updated)
 * @param vcs communities vertex u is linked to (temporary buffer, updated)
 * @param vcout total edge weight from vertex u to community C (temporary buffer, updated)
 * @param x original graph
 * @returns number of changed vertices
 */
template <class G, class K, class W, class F>
inline size_t rakMoveIterationW(vector<K>& vcom, vector<F>& vaff, vector<K>& vcs, vector<W>& vcout, const G& x) {
  size_t a = 0;
  x.forEachVertexKey([&](auto u) {
    if (!vaff[u]) return;
    K d = vcom[u];
    rakClearScanW(vcs, vcout);
    rakScanCommunitiesW(vcs, vcout, x, u, vcom);
    auto [c, w] = rakChooseCommunity(x, u, vcom, vcs, vcout);
    if (c && c!=d) { vcom[u] = c; ++a; x.forEachEdgeKey(u, [&](auto v) { vaff[v] = F(1); }); }
    vaff[u] = F(0);
  });
  return a;
}


#ifdef OPENMP
/**
 * Move each vertex to its best community.
 * @param vcom community each vertex belongs to (updated)
 * @param vaff is vertex affected (updated)
 * @param vcs communities vertex u is linked to (temporary buffer, updated)
 * @param vcout total edge weight from vertex u to community C (temporary buffer, updated)
 * @param x original graph
 * @returns number of changed vertices
 */
template <class G, class K, class W, class F>
inline size_t rakMoveIterationOmpW(vector<K>& vcom, vector<F>& vaff, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x) {
  size_t a = K();
  size_t S = x.span();
  #pragma omp parallel for schedule(dynamic, 2048) reduction(+:a)
  for (K u=0; u<S; ++u) {
    int t = omp_get_thread_num();
    if (!x.hasVertex(u)) continue;
    if (!vaff[u]) continue;
    K d = vcom[u];
    rakClearScanW(*vcs[t], *vcout[t]);
    rakScanCommunitiesW(*vcs[t], *vcout[t], x, u, vcom);
    auto [c, w] = rakChooseCommunity(x, u, vcom, *vcs[t], *vcout[t]);
    if (c && c!=d) { vcom[u] = c; ++a; x.forEachEdgeKey(u, [&](auto v) { vaff[v] = F(1); }); }
    vaff[u] = F(0);
  }
  return a;
}
#endif
#pragma endregion




#pragma region ENVIRONMENT SETUP
/**
 * Setup and perform the RAK algorithm.
 * @param x original graph
 * @param q initial community each vertex belongs to
 * @param o rak options
 * @param fm marking affected vertices / preprocessing to be performed (vaff)
 * @returns rak result
 */
template <class FLAG=char, class G, class K, class FM>
inline RakResult<K> rakInvoke(const G& x, const vector<K>* q, const RakOptions& o, FM fm) {
  using V = typename G::edge_value_type;
  using W = RAK_WEIGHT_TYPE;
  using F = FLAG;
  int l = 0;
  size_t S = x.span();
  size_t N = x.order();
  vector<K> vcom(S), vcs;
  vector<W> vcout(S);
  vector<F> vaff(S);
  float tm = 0;
  float t  = measureDuration([&]() {
    tm += measureDuration([&]() { fm(vaff); });
    if (q) rakInitializeFromW(vcom, x, *q);
    else   rakInitializeW(vcom, x);
    for (l=0; l<o.maxIterations;) {
      size_t n = rakMoveIterationW(vcom, vaff, vcs, vcout, x); ++l;
      if (double(n)/N <= o.tolerance) break;
    }
  }, o.repeat);
  return {vcom, l, t, tm/o.repeat};
}


#ifdef OPENMP
/**
 * Setup and perform the RAK algorithm.
 * @param x original graph
 * @param q initial community each vertex belongs to
 * @param o rak options
 * @param fm marking affected vertices / preprocessing to be performed (vaff)
 * @returns rak result
 */
template <class FLAG=char, class G, class K, class FM>
inline RakResult<K> rakInvokeOmp(const G& x, const vector<K>* q, const RakOptions& o, FM fm) {
  using V = typename G::edge_value_type;
  using W = RAK_WEIGHT_TYPE;
  using F = FLAG;
  int l = 0;
  int T = omp_get_max_threads();
  size_t S = x.span();
  size_t N = x.order();
  vector<K> vcom(S);
  vector<F> vaff(S);
  vector<vector<K>*> vcs(T);
  vector<vector<W>*> vcout(T);
  rakAllocateHashtablesW(vcs, vcout, S);
  float tm = 0;
  float t  = measureDuration([&]() {
    tm += measureDuration([&]() { fm(vaff); });
    if (q) rakInitializeFromOmpW(vcom, x, *q);
    else   rakInitializeOmpW(vcom, x);
    for (l=0; l<o.maxIterations;) {
      size_t n = rakMoveIterationOmpW(vcom, vaff, vcs, vcout, x); ++l;
      if (double(n)/N <= o.tolerance) break;
    }
  }, o.repeat);
  rakFreeHashtablesW(vcs, vcout);
  return {vcom, l, t, tm/o.repeat};
}
#endif
#pragma endregion




#pragma region STATIC/NAIVE-DYNAMIC
/**
 * Obtain the community membership of each vertex with Static/Naive-dynamic RAK.
 * @param x original graph
 * @param q initial community each vertex belongs to
 * @param o rak options
 * @returns rak result
 */
template <class FLAG=char, class G, class K>
inline RakResult<K> rakStatic(const G& x, const vector<K>* q=nullptr, const RakOptions& o={}) {
  auto fm = [](auto& vaff) { fillValueU(vaff, FLAG(1)); };
  return rakInvoke<FLAG>(x, q, o, fm);
}


#ifdef OPENMP
/**
 * Obtain the community membership of each vertex with Static/Naive-dynamic RAK.
 * @param x original graph
 * @param q initial community each vertex belongs to
 * @param o rak options
 * @returns rak result
 */
template <class FLAG=char, class G, class K>
inline RakResult<K> rakStaticOmp(const G& x, const vector<K>* q=nullptr, const RakOptions& o={}) {
  auto fm = [](auto& vaff) { fillValueOmpU(vaff, FLAG(1)); };
  return rakInvokeOmp<FLAG>(x, q, o, fm);
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
template <class G, class K, class V, class F>
inline void rakAffectedVerticesFrontierW(vector<F>& vertices, const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom) {
  fillValueU(vertices, F());
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
template <class G, class K, class V, class F>
inline void rakAffectedVerticesFrontierOmpW(vector<F>& vertices, const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom) {
  fillValueOmpU(vertices, F());
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
 * Obtain the community membership of each vertex with Dynamic Frontier RAK.
 * @param y updated graph
 * @param deletions edge deletions in batch update
 * @param insertions edge insertions in batch update
 * @param q initial community each vertex belongs to
 * @param o rak options
 * @returns rak result
 */
template <class FLAG=char, class G, class K, class V>
inline RakResult<K> rakDynamicFrontier(const G& y, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const RakOptions& o={}) {
  auto fm = [&](auto& vaff) { rakAffectedVerticesFrontierW(vaff, y, deletions, insertions, *q); };
  return rakInvoke<FLAG>(y, q, o, fm);
}


#ifdef OPENMP
/**
 * Obtain the community membership of each vertex with Dynamic Frontier RAK.
 * @param y updated graph
 * @param deletions edge deletions in batch update
 * @param insertions edge insertions in batch update
 * @param q initial community each vertex belongs to
 * @param o rak options
 * @returns rak result
 */
template <class FLAG=char, class G, class K, class V>
inline RakResult<K> rakDynamicFrontierOmp(const G& y, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const RakOptions& o={}) {
  auto fm = [&](auto& vaff) { rakAffectedVerticesFrontierOmpW(vaff, y, deletions, insertions, *q); };
  return rakInvokeOmp<FLAG>(y, q, o, fm);
}
#endif
#pragma endregion
#pragma endregion
