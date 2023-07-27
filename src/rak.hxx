#pragma once
#include <utility>
#include <tuple>
#include <vector>
#include "_main.hxx"
#ifdef OPENMP
#include <omp.h>
#endif

using std::pair;
using std::tuple;
using std::vector;
using std::make_pair;
using std::swap;
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
  /** Community each vertex belongs to. */
  vector<K> membership;
  /** Number of iterations performed. */
  int iterations;
  /** Time spent in milliseconds. */
  float time;
  /** Time spent in milliseconds for preprocessing. */
  float preprocessingTime;
  /** Number of passes performed [1]. */
  int passes;
  #pragma endregion


  #pragma region CONSTRUCTORS
  /**
   * Result of RAK algorithm.
   * @param membership community each vertex belongs to
   * @param iterations number of iterations performed
   * @param time time spent in milliseconds
   * @param preprocessingTime time spent in milliseconds for preprocessing
   * @param passes number of passes performed [1]
   */
  RakResult(vector<K>&& membership, int iterations=0, float time=0, float preprocessingTime=0, int passes=1) :
  membership(membership), iterations(iterations), time(time), preprocessingTime(preprocessingTime), passes(passes) {}


  /**
   * Result of RAK algorithm.
   * @param membership community each vertex belongs to (moved)
   * @param iterations number of iterations performed
   * @param time time spent in milliseconds
   * @param preprocessingTime time spent in milliseconds for preprocessing
   * @param passes number of passes performed [1]
   */
  RakResult(vector<K>& membership, int iterations=0, float time=0, float preprocessingTime=0, int passes=1) :
  membership(move(membership)), iterations(iterations), time(time), preprocessingTime(preprocessingTime), passes(passes) {}
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
inline void rakAllocateHashtables(vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, size_t S) {
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
inline void rakFreeHashtables(vector<vector<K>*>& vcs, vector<vector<W>*>& vcout) {
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
 * @param vcom community each vertex belongs to (updated)
 * @param x original graph
 */
template <class G, class K>
inline void rakInitialize(vector<K>& vcom, const G& x) {
  x.forEachVertexKey([&](auto u) { vcom[u] = u; });
}


#ifdef OPENMP
/**
 * Initialize communities such that each vertex is its own community.
 * @param vcom community each vertex belongs to (updated)
 * @param x original graph
 */
template <class G, class K>
inline void rakInitializeOmp(vector<K>& vcom, const G& x) {
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
 * @param vcom community each vertex belongs to (updated)
 * @param x original graph
 * @param q initial community each vertex belongs to
 */
template <class G, class K>
inline void rakInitializeFrom(vector<K>& vcom, const G& x, const vector<K>& q) {
  x.forEachVertexKey([&](auto u) { vcom[u] = q[u]; });
}


#ifdef OPENMP
/**
 * Initialize communities from given initial communities.
 * @param vcom community each vertex belongs to (updated)
 * @param x original graph
 * @param q initial community each vertex belongs to
 */
template <class G, class K>
inline void rakInitializeFromOmp(vector<K>& vcom, const G& x, const vector<K>& q) {
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
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param u given vertex
 * @param v outgoing edge vertex
 * @param w outgoing edge weight
 * @param vcom community each vertex belongs to
 */
template <bool SELF=false, class K, class V, class W>
inline void rakScanCommunity(vector<K>& vcs, vector<W>& vcout, K u, K v, V w, const vector<K>& vcom) {
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
inline void rakScanCommunities(vector<K>& vcs, vector<W>& vcout, const G& x, K u, const vector<K>& vcom) {
  x.forEachEdge(u, [&](auto v, auto w) { rakScanCommunity<SELF>(vcs, vcout, u, v, w, vcom); });
}


/**
 * Clear communities scan data.
 * @param vcs total edge weight from vertex u to community C (updated)
 * @param vcout communities vertex u is linked to (updated)
 */
template <class K, class W>
inline void rakClearScan(vector<K>& vcs, vector<W>& vcout) {
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
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param x original graph
 * @param fa is vertex allowed to be updated?
 * @returns number of changed vertices
 */
template <class G, class K, class W, class B, class FA>
inline size_t rakMoveIteration(vector<K>& vcom, vector<B>& vaff, vector<K>& vcs, vector<W>& vcout, const G& x, FA fa) {
  size_t a = 0;
  x.forEachVertexKey([&](auto u) {
    if (!fa(u) || !vaff[u]) return;
    K d = vcom[u];
    rakClearScan(vcs, vcout);
    rakScanCommunities(vcs, vcout, x, u, vcom);
    auto [c, w] = rakChooseCommunity(x, u, vcom, vcs, vcout);
    if (c && c!=d) { vcom[u] = c; ++a; x.forEachEdgeKey(u, [&](auto v) { vaff[v] = B(1); }); }
    vaff[u] = B(0);
  });
  return a;
}

#ifdef OPENMP
/**
 * Move each vertex to its best community.
 * @param vcom community each vertex belongs to (updated)
 * @param vaff is vertex affected (updated)
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param x original graph
 * @param fa is vertex allowed to be updated?
 * @returns number of changed vertices
 */
template <class G, class K, class W, class B, class FA>
inline size_t rakMoveIterationOmp(vector<K>& vcom, vector<B>& vaff, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, FA fa) {
  size_t a = K();
  size_t S = x.span();
  #pragma omp parallel for schedule(dynamic, 2048) reduction(+:a)
  for (K u=0; u<S; ++u) {
    int t = omp_get_thread_num();
    if (!x.hasVertex(u)) continue;
    if (!fa(u) || !vaff[u]) continue;
    K d = vcom[u];
    rakClearScan(*vcs[t], *vcout[t]);
    rakScanCommunities(*vcs[t], *vcout[t], x, u, vcom);
    auto [c, w] = rakChooseCommunity(x, u, vcom, *vcs[t], *vcout[t]);
    if (c && c!=d) { vcom[u] = c; ++a; x.forEachEdgeKey(u, [&](auto v) { vaff[v] = B(1); }); }
    vaff[u] = B(0);
  }
  return a;
}
#endif
#pragma endregion




#pragma region MAIN
/**
 * Setup and perform the RAK algorithm.
 * @param x original graph
 * @param q initial communities
 * @param o RAK options
 * @param fm mark affected vertices ()
 * @param fa is vertex allowed to be updated?
 * @returns RAK result
 */
template <class FLAG=char, class G, class K, class FM, class FA>
inline RakResult<K> rakMain(const G& x, const vector<K>* q, const RakOptions& o, FM fm, FA fa) {
  using V = typename G::edge_value_type;
  using W = RAK_WEIGHT_TYPE;
  using B = FLAG;
  int l = 0;
  size_t S = x.span();
  size_t N = x.order();
  vector<K> vcom(S), vcs;
  vector<W> vcout(S);
  vector<B> vaff(S);
  float tm = 0;
  float t  = measureDuration([&]() {
    tm += measureDuration([&]() { fm(vaff); });
    if (q) rakInitializeFrom(vcom, x, *q);
    else   rakInitialize(vcom, x);
    for (l=0; l<o.maxIterations;) {
      size_t n = rakMoveIteration(vcom, vaff, vcs, vcout, x, fa); ++l;
      if (double(n)/N <= o.tolerance) break;
    }
  }, o.repeat);
  return {vcom, l, t, tm/o.repeat};
}


#ifdef OPENMP
/**
 * Setup and perform the RAK algorithm.
 * @param x original graph
 * @param q initial communities
 * @param o RAK options
 * @param fm mark affected vertices ()
 * @param fa is vertex allowed to be updated?
 * @returns RAK result
 */
template <class FLAG=char, class G, class K, class FM, class FA>
inline RakResult<K> rakMainOmp(const G& x, const vector<K>* q, const RakOptions& o, FM fm, FA fa) {
  using V = typename G::edge_value_type;
  using W = RAK_WEIGHT_TYPE;
  using B = FLAG;
  int l = 0;
  int T = omp_get_max_threads();
  size_t S = x.span();
  size_t N = x.order();
  vector<K> vcom(S);
  vector<B> vaff(S);
  vector<vector<K>*> vcs(T);
  vector<vector<W>*> vcout(T);
  rakAllocateHashtables(vcs, vcout, S);
  float tm = 0;
  float t  = measureDuration([&]() {
    tm += measureDuration([&]() { fm(vaff); });
    if (q) rakInitializeFromOmp(vcom, x, *q);
    else   rakInitializeOmp(vcom, x);
    for (l=0; l<o.maxIterations;) {
      size_t n = rakMoveIterationOmp(vcom, vaff, vcs, vcout, x, fa); ++l;
      if (double(n)/N <= o.tolerance) break;
    }
  }, o.repeat);
  rakFreeHashtables(vcs, vcout);
  return {vcom, l, t, tm/o.repeat};
}
#endif
#pragma endregion




#pragma region STATIC/NAIVE-DYNAMIC APPROACH
/**
 * Obtain the community membership of each vertex with Static/Naive-dynamic RAK.
 * @param x original graph
 * @param q initial communities
 * @param o RAK options
 * @returns RAK result
 */
template <class FLAG=char, class G, class K>
inline RakResult<K> rakStatic(const G& x, const vector<K>* q=nullptr, const RakOptions& o={}) {
  auto fm = [](auto& vaff) { fillValueU(vaff, FLAG(1)); };
  auto fa = [](auto u) { return true; };
  return rakMain<FLAG>(x, q, o, fm, fa);
}


#ifdef OPENMP
/**
 * Obtain the community membership of each vertex with Static/Naive-dynamic RAK.
 * @param x original graph
 * @param q initial communities
 * @param o RAK options
 * @returns RAK result
 */
template <class FLAG=char, class G, class K>
inline RakResult<K> rakStaticOmp(const G& x, const vector<K>* q=nullptr, const RakOptions& o={}) {
  auto fm = [](auto& vaff) { fillValueOmpU(vaff, FLAG(1)); };
  auto fa = [](auto u) { return true; };
  return rakMainOmp<FLAG>(x, q, o, fm, fa);
}
#endif
#pragma endregion




#pragma region DYNAMIC FRONTIER APPROACH
/**
 * Find the vertices which should be processed upon a batch of edge insertions and deletions.
 * @param vertices vertex affected flags (output)
 * @param x original graph
 * @param deletions edge deletions for this batch update (undirected, sorted by source vertex id)
 * @param insertions edge insertions for this batch update (undirected, sorted by source vertex id)
 * @param vcom community each vertex belongs to
 */
template <class B, class G, class K, class V>
inline void rakAffectedVerticesFrontier(vector<B>& vertices, const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom) {
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
 * @param x original graph
 * @param deletions edge deletions for this batch update (undirected, sorted by source vertex id)
 * @param insertions edge insertions for this batch update (undirected, sorted by source vertex id)
 * @param vcom community each vertex belongs to
 */
template <class B, class G, class K, class V>
inline void rakAffectedVerticesFrontierOmp(vector<B>& vertices, const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom) {
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
 * @param q initial communities
 * @param o RAK options
 * @returns RAK result
 */
template <class FLAG=char, class G, class K, class V>
inline RakResult<K> rakDynamicFrontier(const G& y, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const RakOptions& o={}) {
  using  B = FLAG;
  const vector<K>& vcom = *q;
  auto fm = [&](auto& vaff) { rakAffectedVerticesFrontier(vaff, y, deletions, insertions, vcom); };
  auto fa = [](auto u) { return true; };
  return rakMain<FLAG>(y, q, o, fm, fa);
}


#ifdef OPENMP
/**
 * Obtain the community membership of each vertex with Dynamic Frontier RAK.
 * @param y updated graph
 * @param deletions edge deletions in batch update
 * @param insertions edge insertions in batch update
 * @param q initial communities
 * @param o RAK options
 * @returns RAK result
 */
template <class FLAG=char, class G, class K, class V>
inline RakResult<K> rakDynamicFrontierOmp(const G& y, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const RakOptions& o={}) {
  using  B = FLAG;
  const vector<K>& vcom = *q;
  auto fm = [&](auto& vaff) { rakAffectedVerticesFrontierOmp(vaff, y, deletions, insertions, vcom); };
  auto fa = [](auto u) { return true; };
  return rakMainOmp<FLAG>(y, q, o, fm, fa);
}
#endif
#pragma endregion
#pragma endregion
