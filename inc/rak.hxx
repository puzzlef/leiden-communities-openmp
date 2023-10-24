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
  /** Time spent in milliseconds for initial marking of affected vertices. */
  float markingTime;
  /** Time spent in initializing community memberships. */
  float initializationTime;
  /** Number of vertices initially marked as affected. */
  size_t affectedVertices;
  #pragma endregion


  #pragma region CONSTRUCTORS
  /**
   * Result of RAK algorithm.
   * @param membership community membership each vertex belongs to
   * @param iterations number of iterations performed
   * @param time time spent in milliseconds
   * @param markingTime time spent in milliseconds for initial marking of affected vertices
   * @param initializationTime time spent in initializing community memberships
   * @param affectedVertices number of vertices initially marked as affected
   */
  RakResult(vector<K>&& membership, int iterations=0, float time=0, float markingTime=0, float initializationTime=0, size_t affectedVertices=0) :
  membership(membership), iterations(iterations), time(time), markingTime(markingTime), initializationTime(initializationTime), affectedVertices(affectedVertices) {}


  /**
   * Result of RAK algorithm.
   * @param membership community membership each vertex belongs to (moved)
   * @param iterations number of iterations performed
   * @param time time spent in milliseconds
   * @param markingTime time spent in milliseconds for initial marking of affected vertices
   * @param initializationTime time spent in initializing community memberships
   * @param affectedVertices number of vertices initially marked as affected
   */
  RakResult(vector<K>& membership, int iterations=0, float time=0, float markingTime=0, float initializationTime=0, size_t affectedVertices=0) :
  membership(move(membership)), iterations(iterations), time(time), markingTime(markingTime), initializationTime(initializationTime), affectedVertices(affectedVertices) {}
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
inline void rakScanCommunityU(vector<K>& vcs, vector<W>& vcout, K u, K v, V w, const vector<K>& vcom) {
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
  x.forEachEdge(u, [&](auto v, auto w) { rakScanCommunityU<SELF>(vcs, vcout, u, v, w, vcom); });
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
 * @param fa is vertex allowed to be updated? (u)
 * @returns number of changed vertices
 */
template <class G, class K, class W, class F, class FA>
inline size_t rakMoveIterationW(vector<K>& vcom, vector<F>& vaff, vector<K>& vcs, vector<W>& vcout, const G& x, FA fa) {
  size_t a = 0;
  x.forEachVertexKey([&](auto u) {
    if (!fa(u) || !vaff[u]) return;
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
 * @param fa is vertex allowed to be updated? (u)
 * @returns number of changed vertices
 */
template <class G, class K, class W, class F, class FA>
inline size_t rakMoveIterationOmpW(vector<K>& vcom, vector<F>& vaff, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, FA fa) {
  size_t a = K();
  size_t S = x.span();
  #pragma omp parallel for schedule(dynamic, 2048) reduction(+:a)
  for (K u=0; u<S; ++u) {
    int t = omp_get_thread_num();
    if (!x.hasVertex(u)) continue;
    if (!fa(u) || !vaff[u]) continue;
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
 * @param o rak options
 * @param fi initialzing community membership (vcom)
 * @param fm marking affected vertices (vaff, vcs, vcout, vcom)
 * @param fa is vertex allowed to be updated? (u)
 * @returns rak result
 */
template <bool DYNAMIC=false, class FLAG=char, class G, class FI, class FM, class FA>
inline auto rakInvoke(const G& x, const RakOptions& o, FI fi, FM fm, FA fa) {
  using K = typename G::key_type;
  using V = typename G::edge_value_type;
  using W = RAK_WEIGHT_TYPE;
  using F = FLAG;
  int l = 0;
  // Get graph properties.
  size_t S = x.span();
  size_t N = x.order();
  // Allocate buffers.
  vector<F> vaff(S);   // Affected vertex flag
  vector<K> vcom;      // Community membership
  vector<K> vcs;       // Hashtable keys
  vector<W> vcout(S);  // Hashtable values
  if (!DYNAMIC) vcom.resize(S);
  // Perform RAK algorithm.
  float tm = 0, ti = 0;  // Time spent in different phases
  float t  = measureDuration([&]() {
    // Initialize community membership.
    ti += measureDuration([&]() { fi(vcom); });
    // Mark affected vertices.
    tm += measureDuration([&]() { fm(vaff, vcs, vcout, vcom); });
    // Perform iterations.
    for (l=0; l<o.maxIterations;) {
      size_t n = rakMoveIterationW(vcom, vaff, vcs, vcout, x, fa); ++l;
      if (double(n)/N <= o.tolerance) break;
    }
  }, o.repeat);
  return RakResult<K>(vcom, l, t, tm/o.repeat, ti/o.repeat, countValue(vaff, F(1)));
}


#ifdef OPENMP
/**
 * Setup and perform the RAK algorithm.
 * @param x original graph
 * @param o rak options
 * @param fi initialzing community membership (vcom)
 * @param fm marking affected vertices (vaff, vcs, vcout, vcom)
 * @param fa is vertex allowed to be updated? (u)
 * @returns rak result
 */
template <bool DYNAMIC=false, class FLAG=char, class G, class FI, class FM, class FA>
inline auto rakInvokeOmp(const G& x, const RakOptions& o, FI fi, FM fm, FA fa) {
  using K = typename G::key_type;
  using V = typename G::edge_value_type;
  using W = RAK_WEIGHT_TYPE;
  using F = FLAG;
  int l = 0;
  int T = omp_get_max_threads();
  // Get graph properties.
  size_t S = x.span();
  size_t N = x.order();
  // Allocate buffers.
  vector<F> vaff(S);  // Affected vertex flag
  vector<K> vcom;     // Community membership
  vector<vector<K>*> vcs(T);    // Hashtable keys
  vector<vector<W>*> vcout(T);  // Hashtable values
  if (!DYNAMIC) vcom.resize(S);
  rakAllocateHashtablesW(vcs, vcout, S);
  // Perform RAK algorithm.
  float tm = 0, ti = 0;  // Time spent in different phases
  float t  = measureDuration([&]() {
    // Initialize community membership.
    ti += measureDuration([&]() { fi(vcom); });
    // Mark affected vertices.
    tm += measureDuration([&]() { fm(vaff, vcs, vcout, vcom); });
    // Perform iterations.
    for (l=0; l<o.maxIterations;) {
      size_t n = rakMoveIterationOmpW(vcom, vaff, vcs, vcout, x, fa); ++l;
      if (double(n)/N <= o.tolerance) break;
    }
  }, o.repeat);
  rakFreeHashtablesW(vcs, vcout);
  return RakResult<K>(vcom, l, t, tm/o.repeat, ti/o.repeat, countValueOmp(vaff, F(1)));
}
#endif
#pragma endregion




#pragma region REPEAT SETUP (DYNAMIC)
/**
 * Setup the Dynamic RAK algorithm for multiple runs.
 * @param qs initial community membership for each run (updated)
 * @param q initial community membership
 * @param repeat number of runs
 */
template <class K>
inline void rakSetupInitialsW(vector2d<K>& qs, const vector<K>& q, int repeat) {
  qs.resize(repeat);
  for (int r=0; r<repeat; ++r)
    qs[r] = q;
}
#pragma endregion




#pragma region STATIC
/**
 * Obtain the community membership of each vertex with Static RAK.
 * @param x original graph
 * @param o rak options
 * @returns rak result
 */
template <class FLAG=char, class G>
inline auto rakStatic(const G& x, const RakOptions& o={}) {
  auto fi = [&](auto& vcom) { rakInitializeW(vcom, x); };
  auto fm = [ ](auto& vaff, auto& vcs, auto& vcout, const auto& vcom) { fillValueU(vaff, FLAG(1)); };
  auto fa = [ ](auto u) { return true; };
  return rakInvoke<false, FLAG>(x, o, fi, fm, fa);
}


#ifdef OPENMP
/**
 * Obtain the community membership of each vertex with Static RAK.
 * @param x original graph
 * @param o rak options
 * @returns rak result
 */
template <class FLAG=char, class G>
inline auto rakStaticOmp(const G& x, const RakOptions& o={}) {
  auto fi = [&](auto& vcom) { rakInitializeOmpW(vcom, x); };
  auto fm = [ ](auto& vaff, auto& vcs, auto& vcout, const auto& vcom) { fillValueOmpU(vaff, FLAG(1)); };
  auto fa = [ ](auto u) { return true; };
  return rakInvokeOmp<false, FLAG>(x, o, fi, fm, fa);
}
#endif
#pragma endregion




#pragma region NAIVE-DYNAMIC
/**
 * Obtain the community membership of each vertex with Naive-dynamic RAK.
 * @param y updated graph
 * @param q initial community each vertex belongs to
 * @param o rak options
 * @returns rak result
 */
template <class FLAG=char, class G, class K>
inline RakResult<K> rakNaiveDynamic(const G& y, const vector<K>& q, const RakOptions& o={}) {
  vector2d<K> qs;
  rakSetupInitialsW(qs, q, o.repeat);
  int  r  = 0;
  auto fi = [&](auto& vcom) { vcom = move(qs[r]); ++r; };
  auto fm = [ ](auto& vaff, auto& vcs, auto& vcout, const auto& vcom) { fillValueU(vaff, FLAG(1)); };
  auto fa = [ ](auto u) { return true; };
  return rakInvoke<true, FLAG>(y, o, fi, fm, fa);
}


#ifdef OPENMP
/**
 * Obtain the community membership of each vertex with Naive-dynamic RAK.
 * @param y updated graph
 * @param q initial community each vertex belongs to
 * @param o rak options
 * @returns rak result
 */
template <class FLAG=char, class G, class K>
inline RakResult<K> rakNaiveDynamicOmp(const G& y, const vector<K>& q, const RakOptions& o={}) {
  vector2d<K> qs;
  rakSetupInitialsW(qs, q, o.repeat);
  int  r  = 0;
  auto fi = [&](auto& vcom) { vcom = move(qs[r]); ++r; };
  auto fm = [ ](auto& vaff, auto& vcs, auto& vcout, const auto& vcom) { fillValueOmpU(vaff, FLAG(1)); };
  auto fa = [ ](auto u) { return true; };
  return rakInvokeOmp<true, FLAG>(y, o, fi, fm, fa);
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
 * @returns flags for each vertex marking whether it is affected
 */
template <class B, class G, class K, class V, class W>
inline void rakAffectedVerticesDeltaScreeningW(vector<B>& vertices, vector<B>& neighbors, vector<B>& communities, vector<K>& vcs, vector<W>& vcout, const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom) {
  fillValueU(vertices,    B());
  fillValueU(neighbors,   B());
  fillValueU(communities, B());
  for (const auto& [u, v, w] : deletions) {
    if (vcom[u] != vcom[v]) continue;
    vertices[u]  = 1;
    neighbors[u] = 1;
    communities[vcom[v]] = 1;
  }
  for (size_t i=0; i<insertions.size();) {
    K u = get<0>(insertions[i]);
    rakClearScanW(vcs, vcout);
    for (; i<insertions.size() && get<0>(insertions[i])==u; ++i) {
      K v = get<1>(insertions[i]);
      V w = get<2>(insertions[i]);
      if (vcom[u] == vcom[v]) continue;
      rakScanCommunityU(vcs, vcout, u, v, w, vcom);
    }
    auto [c, w] = rakChooseCommunity(y, u, vcom, vcs, vcout);
    if (w<=0 || c==vcom[u]) continue;
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
 * @returns flags for each vertex marking whether it is affected
 */
template <class B, class G, class K, class V, class W>
inline void rakAffectedVerticesDeltaScreeningOmpW(vector<B>& vertices, vector<B>& neighbors, vector<B>& communities, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom) {
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
      rakClearScanW(*vcs[t], *vcout[t]);
      for (; i<I && get<0>(insertions[i])==u; ++i) {
        K v = get<1>(insertions[i]);
        V w = get<2>(insertions[i]);
        if (vcom[u] == vcom[v]) continue;
        rakScanCommunityU(*vcs[t], *vcout[t], u, v, w, vcom);
      }
      auto [c, w] = rakChooseCommunity(y, u, vcom, *vcs[t], *vcout[t]);
      if (w<=0 || c==vcom[u]) continue;
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
 * Obtain the community membership of each vertex with Dynamic Delta-screening RAK.
 * @param y updated graph
 * @param deletions edge deletions in batch update
 * @param insertions edge insertions in batch update
 * @param q initial community each vertex belongs to
 * @param o rak options
 * @returns rak result
 */
template <class FLAG=char, class G, class K, class V>
inline RakResult<K> rakDynamicDeltaScreening(const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& q, const RakOptions& o={}) {
  using  W = RAK_WEIGHT_TYPE;
  using  B = FLAG;
  size_t S = y.span();
  vector<B> vertices(S), neighbors(S), communities(S);
  vector2d<K> qs;
  rakSetupInitialsW(qs, q, o.repeat);
  int  r  = 0;
  auto fi = [&](auto& vcom) { vcom = move(qs[r]); ++r; };
  auto fm = [&](auto& vaff, auto& vcs, auto& vcout, const auto& vcom) {
    rakAffectedVerticesDeltaScreeningW(vertices, neighbors, communities, vcs, vcout, y, deletions, insertions, vcom);
    copyValuesW(vaff, vertices);
  };
  auto fa = [&](auto u) { return vertices[u]==B(1); };
  return rakInvoke<true, FLAG>(y, o, fi, fm, fa);
}


#ifdef OPENMP
/**
 * Obtain the community membership of each vertex with Dynamic Delta-screening RAK.
 * @param y updated graph
 * @param deletions edge deletions in batch update
 * @param insertions edge insertions in batch update
 * @param q initial community each vertex belongs to
 * @param o rak options
 * @returns rak result
 */
template <class FLAG=char, class G, class K, class V>
inline RakResult<K> rakDynamicDeltaScreeningOmp(const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& q, const RakOptions& o={}) {
  using  W = RAK_WEIGHT_TYPE;
  using  B = FLAG;
  size_t S = y.span();
  int    T = omp_get_max_threads();
  vector<B> vertices(S), neighbors(S), communities(S);
  vector2d<K> qs;
  rakSetupInitialsW(qs, q, o.repeat);
  int  r  = 0;
  auto fi = [&](auto& vcom) { vcom = move(qs[r]); ++r; };
  auto fm = [&](auto& vaff, auto& vcs, auto& vcout, const auto& vcom) {
    rakAffectedVerticesDeltaScreeningOmpW(vertices, neighbors, communities, vcs, vcout, y, deletions, insertions, vcom);
    copyValuesOmpW(vaff, vertices);
  };
  auto fa = [&](auto u) { return vertices[u]==B(1); };
  return rakInvokeOmp<true, FLAG>(y, o, fi, fm, fa);
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
inline void rakAffectedVerticesFrontierW(vector<B>& vertices, const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom) {
  fillValueU(vertices, B());
  for (const auto& [u, v, w] : deletions) {
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
inline void rakAffectedVerticesFrontierOmpW(vector<B>& vertices, const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom) {
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
 * Obtain the community membership of each vertex with Dynamic Frontier RAK.
 * @param y updated graph
 * @param deletions edge deletions in batch update
 * @param insertions edge insertions in batch update
 * @param q initial community each vertex belongs to
 * @param o rak options
 * @returns rak result
 */
template <class FLAG=char, class G, class K, class V>
inline auto rakDynamicFrontier(const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& q, const RakOptions& o={}) {
  vector2d<K> qs;
  rakSetupInitialsW(qs, q, o.repeat);
  int  r  = 0;
  auto fi = [&](auto& vcom) { vcom = move(qs[r]); ++r; };
  auto fm = [&](auto& vaff, auto& vcs, auto& vcout, const auto& vcom) {
    rakAffectedVerticesFrontierW(vaff, y, deletions, insertions, vcom);
  };
  auto fa = [ ](auto u) { return true; };
  return rakInvoke<true, FLAG>(y, o, fi, fm, fa);
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
inline auto rakDynamicFrontierOmp(const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& q, const RakOptions& o={}) {
  vector2d<K> qs;
  rakSetupInitialsW(qs, q, o.repeat);
  int  r  = 0;
  auto fi = [&](auto& vcom) { vcom = move(qs[r]); ++r; };
  auto fm = [&](auto& vaff, auto& vcs, auto& vcout, const auto& vcom) {
    rakAffectedVerticesFrontierOmpW(vaff, y, deletions, insertions, vcom);
  };
  auto fa = [ ](auto u) { return true; };
  return rakInvokeOmp<true, FLAG>(y, o, fi, fm, fa);
}
#endif
#pragma endregion
#pragma endregion
