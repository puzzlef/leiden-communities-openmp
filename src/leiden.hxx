#pragma once
#include <utility>
#include <algorithm>
#include <vector>
#include "_main.hxx"
#include "Graph.hxx"
#include "duplicate.hxx"
#include "properties.hxx"
#include "modularity.hxx"

#ifdef OPENMP
#include <omp.h>
#endif

using std::pair;
using std::tuple;
using std::vector;
using std::make_pair;
using std::move;
using std::get;
using std::min;
using std::max;




// LEIDEN OPTIONS
// --------------

struct LeidenOptions {
  int    repeat;
  double resolution;
  double tolerance;
  double aggregationTolerance;
  double toleranceDecline;
  int    maxIterations;
  int    maxPasses;

  LeidenOptions(int repeat=1, double resolution=1, double tolerance=1e-2, double aggregationTolerance=0.8, double toleranceDecline=100, int maxIterations=20, int maxPasses=10) :
  repeat(repeat), resolution(resolution), tolerance(tolerance), aggregationTolerance(aggregationTolerance), toleranceDecline(toleranceDecline), maxIterations(maxIterations), maxPasses(maxPasses) {}
};

// Weight to be using in hashtable.
#define LEIDEN_WEIGHT_TYPE double




// LEIDEN RESULT
// -------------

template <class K>
struct LeidenResult {
  vector<K> membership;
  int   iterations;
  int   passes;
  float time;
  float preprocessingTime;

  LeidenResult(vector<K>&& membership, int iterations=0, int passes=0, float time=0, float preprocessingTime=0) :
  membership(membership), iterations(iterations), passes(passes), time(time), preprocessingTime(preprocessingTime) {}

  LeidenResult(vector<K>& membership, int iterations=0, int passes=0, float time=0, float preprocessingTime=0) :
  membership(move(membership)), iterations(iterations), passes(passes), time(time), preprocessingTime(preprocessingTime) {}
};




// LEIDEN HASHTABLES
// -----------------

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




// LEIDEN INITIALIZE
// -----------------

/**
 * Find the total edge weight of each vertex.
 * @param vtot total edge weight of each vertex (updated, should be initialized to 0)
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
template <class G, class W>
inline void leidenVertexWeightsOmpW(vector<W>& vtot, const G& x) {
  using  K = typename G::key_type;
  size_t S = x.span();
  #pragma omp parallel for schedule(auto)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    x.forEachEdge(u, [&](auto v, auto w) { vtot[u] += w; });
  }
}
#endif


/**
 * Find the total edge weight of each community.
 * @param ctot total edge weight of each community (updated, should be initialized to 0)
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
template <class G, class K, class W>
inline void leidenCommunityWeightsOmpW(vector<W>& ctot, const G& x, const vector<K>& vcom, const vector<W>& vtot) {
  size_t S = x.span();
  #pragma omp parallel for schedule(auto)
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
 * @param vcom community each vertex belongs to (updated, should be initialized to 0)
 * @param ctot total edge weight of each community (updated, should be initilized to 0)
 * @param x original graph
 * @param vtot total edge weight of each vertex
 */
template <class G, class K, class W>
inline void leidenInitializeW(vector<K>& vcom, vector<W>& ctot, const G& x, const vector<W>& vtot) {
  x.forEachVertexKey([&](auto u) {
    vcom[u] = u;
    ctot[u] = vtot[u];
  });
}

#ifdef OPENMP
template <class G, class K, class W>
inline void leidenInitializeOmpW(vector<K>& vcom, vector<W>& ctot, const G& x, const vector<W>& vtot) {
  size_t S = x.span();
  #pragma omp parallel for schedule(auto)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    vcom[u] = u;
    ctot[u] = vtot[u];
  }
}
#endif


/**
 * Initialize communities from given initial communities.
 * @param vcom community each vertex belongs to (updated, should be initialized to 0)
 * @param vcob community bound each vertex belongs to (updated, should be initialized to 0)
 * @param ctot total edge weight of each community (updated, should be initilized to 0)
 * @param x original graph
 * @param vtot total edge weight of each vertex
 * @param q initial community each vertex belongs to
 */
template <bool REFINE=false, class G, class K, class W>
inline void leidenInitializeFromW(vector<K>& vcom, vector<K>& vcob, vector<W>& ctot, const G& x, const vector<W>& vtot, const vector<K>& q) {
  x.forEachVertexKey([&](auto u) {
    K c = q[u];
    if (REFINE) {
      vcom[u] = u;
      vcob[u] = c;
      ctot[u] = vtot[u];
    }
    else {
      vcom[u]  = c;
      ctot[c] += vtot[u];
    }
  });
}

#ifdef OPENMP
template <bool REFINE=false, class G, class K, class W>
inline void leidenInitializeFromOmpW(vector<K>& vcom, vector<K>& vcob, vector<W>& ctot, const G& x, const vector<W>& vtot, const vector<K>& q) {
  size_t S = x.span();
  #pragma omp parallel for schedule(auto)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    K c = q[u];
    if (REFINE) {
      vcom[u] = u;
      vcob[u] = c;
      ctot[u] = vtot[u];
    }
    else {
      vcom[u]  = c;
      #pragma omp atomic
      ctot[c] += vtot[u];
    }
  }
}
#endif




// LEIDEN COMMUNITY VERTICES
// -------------------------

/**
 * Find the number of vertices in each community.
 * @param a number of vertices belonging to each community (updated)
 * @param x original graph
 * @param vcom community each vertex belongs to
 * @returns number of communities
 */
template <class G, class K>
inline size_t leidenCountCommunityVerticesW(K *a, const G& x, const K *vcom) {
  size_t S = x.span();
  size_t n = 0;
  fillValueU(a, S, K());
  x.forEachVertexKey([&](auto u) {
    K c = vcom[u];
    if (a[c]==0) ++n;
    ++a[c];
  });
  return n;
}
template <class G, class K>
inline size_t leidenCountCommunityVerticesW(vector<K>& a, const G& x, const vector<K>& vcom) {
  return leidenCountCommunityVerticesW(a.data(), x, vcom.data());
}


#ifdef OPENMP
template <class G, class K>
inline size_t leidenCountCommunityVerticesOmpW(K *a, const G& x, const K *vcom) {
  size_t S = x.span();
  size_t n = 0;
  fillValueOmpU(a, S, K());
  #pragma omp parallel for schedule(auto) reduction(+:n)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    K c = vcom[u], m = 0;
    #pragma omp atomic capture
    { m = a[c]; ++a[c]; }
    if (m==0) ++n;
  }
  return n;
}
template <class G, class K>
inline size_t leidenCountCommunityVerticesOmpW(vector<K>& a, const G& x, const vector<K>& vcom) {
  return leidenCountCommunityVerticesOmpW(a.data(), x, vcom.data());
}
#endif




/**
 * Find the vertices in each community.
 * @param co csr offsets for vertices belonging to each community (updated)
 * @param ce csr data vertices belonging to each community (updated)
 * @param cn number of vertices in each community (updated)
 * @param x original graph
 * @param vcom community each vertex belongs to
 */
template <class G, class K>
inline void leidenCommunityVerticesW(vector<K>& co, vector<K>& ce, vector<K>& cn, const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  co[S] = exclusiveScanW(co, cn);
  fillValueU(cn, K());
  x.forEachVertexKey([&](auto u) {
    K c = vcom[u];
    K i = cn[c]++;
    ce[co[c] + i] = u;
  });
}


#ifdef OPENMP
template <class G, class K>
inline void leidenCommunityVerticesOmpW(vector<K>& co, vector<K>& ce, vector<K>& cn, vector<K>& bufk, const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  co[S] = exclusiveScanOmpW(co, bufk, cn);
  fillValueOmpU(cn, K());
  #pragma omp parallel for schedule(auto)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    K c = vcom[u], i = 0;
    #pragma omp atomic capture
    { i = cn[c]; ++cn[c]; }
    ce[co[c] + i] = u;
  }
}
#endif




// LEIDEN LOOKUP COMMUNITIES
// -------------------------

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
template <class K>
inline void leidenLookupCommunitiesOmpU(vector<K>& a, const vector<K>& vcom) {
  size_t S = a.size();
  #pragma omp parallel for schedule(auto)
  for (size_t u=0; u<S; ++u)
    a[u] = vcom[a[u]];
}
#endif




// LEIDEN CHANGE COMMUNITY
// -----------------------

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
inline auto leidenChooseCommunity(const G& x, K u, const vector<K>& vcom, const vector<W>& vtot, const vector<W>& ctot, const vector<K>& vcs, const vector<W>& vcout, double M, double R) {
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
inline void leidenChangeCommunityW(vector<K>& vcom, vector<W>& ctot, const G& x, K u, K c, const vector<W>& vtot) {
  K d = vcom[u];
  ctot[d] -= vtot[u];
  ctot[c] += vtot[u];
  vcom[u] = c;
}

#ifdef OPENMP
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




// LEIDEN MOVE
// -----------

/**
 * Leiden algorithm's local moving phase.
 * @param vcom community each vertex belongs to (initial, updated)
 * @param ctot total edge weight of each community (precalculated, updated)
 * @param vaff is vertex affected flag (updated)
 * @param vpro has vertex been processed flag (updated)
 * @param vcs communities vertex u is linked to (temporary buffer, updated)
 * @param vcout total edge weight from vertex u to community C (temporary buffer, updated)
 * @param x original graph
 * @param vcob community bound each vertex belongs to
 * @param vtot total edge weight of each vertex
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @param L max iterations
 * @param fc has local moving phase converged?
 * @returns iterations performed (0 if converged already)
 */
template <bool REFINE=false, class G, class K, class W, class B, class FC>
inline int leidenMoveW(vector<K>& vcom, vector<W>& ctot, vector<B>& vaff, vector<B>& vpro, vector<K>& vcs, vector<W>& vcout, const G& x, const vector<K>& vcob, const vector<W>& vtot, double M, double R, int L, FC fc) {
  int l = 0;
  W  el = W();
  for (; l<L;) {
    el = W();
    x.forEachVertexKey([&](auto u) {
      if (!vaff[u]) return;
      leidenClearScanW(vcs, vcout);
      leidenScanCommunitiesW<false, REFINE>(vcs, vcout, x, u, vcom, vcob);
      auto [c, e] = leidenChooseCommunity(x, u, vcom, vtot, ctot, vcs, vcout, M, R);
      if (c)      { leidenChangeCommunityW(vcom, ctot, x, u, c, vtot); x.forEachEdgeKey(u, [&](auto v) { vaff[v] = B(1); }); }
      if (REFINE) vpro[u] = B(1);
      vaff[u] = B();
      el += e;  // l1-norm
    });
    if (fc(el, l++)) break;
  }
  // Move unprocessed vertices to their community bound
  if (REFINE) {
    x.forEachVertexKey([&](auto u) {
      if (!vpro[u]) vcom[u] = vcob[u];
    });
  }
  return l>1 || el? l : 0;
}

#ifdef OPENMP
template <bool REFINE=false, class G, class K, class W, class B, class FC>
inline int leidenMoveOmpW(vector<K>& vcom, vector<W>& ctot, vector<B>& vaff, vector<B>& vpro, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, const vector<K>& vcob, const vector<W>& vtot, double M, double R, int L, FC fc) {
  size_t S = x.span();
  int l = 0;
  W  el = W();
  for (; l<L;) {
    el = W();
    #pragma omp parallel for schedule(auto) reduction(+:el)
    for (K u=0; u<S; ++u) {
      int t = omp_get_thread_num();
      if (!x.hasVertex(u)) continue;
      if (!vaff[u]) continue;
      leidenClearScanW(*vcs[t], *vcout[t]);
      leidenScanCommunitiesW<false, REFINE>(*vcs[t], *vcout[t], x, u, vcom, vcob);
      auto [c, e] = leidenChooseCommunity(x, u, vcom, vtot, ctot, *vcs[t], *vcout[t], M, R);
      if (c)      { leidenChangeCommunityOmpW(vcom, ctot, x, u, c, vtot); x.forEachEdgeKey(u, [&](auto v) { vaff[v] = B(1); }); }
      if (REFINE) vpro[u] = B(1);
      vaff[u] = B();
      el += e;  // l1-norm
    }
    if (fc(el, l++)) break;
  }
  // Move unprocessed vertices to their community bound.
  if (REFINE) {
    #pragma omp parallel for schedule(auto)
    for (K u=0; u<S; ++u) {
      if (!x.hasVertex(u)) continue;
      if (!vpro[u]) vcom[u] = vcob[u];
    }
  }
  return l>1 || el? l : 0;
}
#endif




// LEIDEN AGGREGATE
// ----------------

/**
 * Leiden algorithm's community aggregation phase.
 * @param a output graph
 * @param vcs communities vertex u is linked to (temporary buffer, updated)
 * @param vcout total edge weight from vertex u to community C (temporary buffer, updated)
 * @param x original graph
 * @param vcom community each vertex belongs to
 * @param co csr offsets for vertices belonging to each community
 * @param ce csr data vertices belonging to each community
 */
template <class G, class K, class W>
inline void leidenAggregateW(G& a, vector<K>& vcs, vector<W>& vcout, const G& x, const vector<K>& vcom, const vector<K>& co, const vector<K>& ce) {
  size_t S = x.span();
  a.respan(S);
  for (K c=0; c<S; ++c) {
    K oc = co[c];
    K nc = co[c+1] - co[c];
    if (nc==0) continue;
    leidenClearScanW(vcs, vcout);
    for (K i=0; i<nc; ++i) {
      K u = ce[oc+i];
      leidenScanCommunitiesW<true>(vcs, vcout, x, u, vcom);
    }
    // a.reserveEdges(c, vcs.size());
    a.addVertex(c);
    for (auto d : vcs)
      a.addEdge(c, d, vcout[d]);
  }
  // Aggregated graph has unique edges, so an update may not be necessary.
  a.update();
}

#ifdef OPENMP
template <class G, class K, class W>
inline void leidenAggregateOmpW(G& a, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, const vector<K>& vcom, const vector<K>& co, const vector<K>& ce) {
  size_t S = x.span();
  a.respan(S);
  #pragma omp parallel for schedule(auto)
  for (K c=0; c<S; ++c) {
    int t = omp_get_thread_num();
    K oc = co[c];
    K nc = co[c+1] - co[c];
    if (nc==0) continue;
    leidenClearScanW(*vcs[t], *vcout[t]);
    for (K i=0; i<nc; ++i) {
      K u = ce[oc+i];
      leidenScanCommunitiesW<true>(*vcs[t], *vcout[t], x, u, vcom);
    }
    // a.reserveEdges(c, (*vcs[t]).size());
    for (auto d : *vcs[t])
      a.addEdge(c, d, (*vcout[t])[d]);
  }
  // Aggregated graph has unique edges, so an update may not be necessary.
  updateOmpU(a);
}
#endif


template <class G, class K, class W>
inline auto leidenAggregate(vector<K>& vcs, vector<W>& vcout, const G& x, const vector<K>& vcom, const vector<K>& co, const vector<K>& ce) {
  G a; leidenAggregateW(a, vcs, vcout, x, vcom, co, ce);
  return a;
}

#ifdef OPENMP
template <class G, class K, class W>
inline auto leidenAggregateOmp(vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, const vector<K>& vcom, const vector<K>& co, const vector<K>& ce) {
  G a; leidenAggregateOmpW(a, vcs, vcout, x, vcom, co, ce);
  return a;
}
#endif




// LEIDEN
// ------

/**
 * Find the community each vertex belongs to.
 * @param x original graph
 * @param q initial community each vertex belongs to
 * @param o leiden options
 * @param fm marking affected vertices / preprocessing to be performed (vaff)
 * @returns community each vertex belongs to
 */
template <bool JUMP=false, bool REFINE=false, class FLAG=char, class G, class K, class FM>
auto leidenSeq(const G& x, const vector<K> *q, const LeidenOptions& o, FM fm) {
  using  W = LOUVAIN_WEIGHT_TYPE;
  using  B = FLAG;
  double R = o.resolution;
  int    L = o.maxIterations, l = 0;
  int    P = o.maxPasses, p = 0, s = 0;
  size_t S = x.span();
  double M = edgeWeight(x)/2;
  vector<K> vcom(S), vcs, a(S);
  vector<W> vtot(S), ctot(S), vcout(S);
  vector<K> co(S+1), ce(S), cn(S);
  vector<B> vaff(S), vpro(S);
  vector<K>& vcob = a;
  float tm = 0;
  float t  = measureDurationMarked([&](auto mark) {
    double E  = o.tolerance;
    auto   fc = [&](double el, int l) {
      if (!JUMP) return el<=E;
      if (el>E)  return false;
      if (l>=1)  return true;
      for (; p<P-1 && el<=E; ++p)
        E /= o.toleranceDecline;
      return el<=E;
    };
    G y; y.respan(S);
    fillValueU(vpro, B());
    fillValueU(vcom, K());
    fillValueU(vtot, W());
    fillValueU(ctot, W());
    fillValueU(a, K());
    mark([&]() {
      tm = measureDuration([&]() { fm(vaff); });
      leidenVertexWeightsW(vtot, x);
      if (q) leidenInitializeFromW<REFINE>(vcom, vcob, ctot, x, vtot, *q);
      else   leidenInitializeW(vcom, ctot, x, vtot);
      for (l=0, p=0, s=0; M>0 && p<P;) {
        int m = 0;
        if (s==0) m = leidenMoveW<REFINE>(vcom, ctot, vaff, vpro, vcs, vcout, x, vcob, vtot, M, R, L, fc);
        else      m = leidenMoveW<false> (vcom, ctot, vaff, vpro, vcs, vcout, y, vcob, vtot, M, R, L, fc);
        if (s==0) copyValuesW(a, vcom);
        else      leidenLookupCommunitiesU(a, vcom);
        l += max(m, 1); ++p; ++s;
        if (m<=1 || p>=P) break;
        const G& g = s<=1? x : y;
        size_t gn = g.order();
        size_t yn = leidenCountCommunityVerticesW(cn, g, vcom);
        if (double(yn)/gn >= o.aggregationTolerance) break;
        leidenCommunityVerticesW(co, ce, cn, g, vcom);
        y = leidenAggregate(vcs, vcout, g, vcom, co, ce);
        fillValueU(vcom, K());
        fillValueU(vtot, W());
        fillValueU(ctot, W());
        fillValueU(vaff, B(1));
        leidenVertexWeightsW(vtot, y);
        leidenInitializeW(vcom, ctot, y, vtot);
        E /= o.toleranceDecline;
      }
    });
  }, o.repeat);
  return LeidenResult<K>(a, l, s, t, tm);
}

#ifdef OPENMP
template <bool JUMP=false, bool REFINE=false, class FLAG=char, class G, class K, class FM>
auto leidenOmp(const G& x, const vector<K> *q, const LeidenOptions& o, FM fm) {
  using  W = LOUVAIN_WEIGHT_TYPE;
  using  B = FLAG;
  double R = o.resolution;
  int    L = o.maxIterations, l = 0;
  int    P = o.maxPasses, p = 0, s = 0;
  size_t S = x.span();
  double M = edgeWeightOmp(x)/2;
  int    T = omp_get_max_threads();
  vector<K> bufk(T);
  vector<K> vcom(S), a(S);
  vector<W> vtot(S), ctot(S);
  vector<K> co(S+1), ce(S), cn(S);
  vector<B> vaff(S), vpro(S);
  vector<K>& vcob = a;
  vector<vector<K>*> vcs(T);
  vector<vector<W>*> vcout(T);
  leidenAllocateHashtablesW(vcs, vcout, S);
  float tm = 0;
  float t  = measureDurationMarked([&](auto mark) {
    double E  = o.tolerance;
    auto   fc = [&](double el, int l) {
      if (!JUMP) return el<=E;
      if (el>E)  return false;
      if (l>=1)  return true;
      for (; p<P-1 && el<=E; ++p)
        E /= o.toleranceDecline;
      return el<=E;
    };
    G y; y.respan(S);
    fillValueOmpU(vpro, B());
    fillValueOmpU(vcom, K());
    fillValueOmpU(vtot, W());
    fillValueOmpU(ctot, W());
    fillValueOmpU(a, K());
    mark([&]() {
      tm = measureDuration([&]() { fm(vaff); });
      leidenVertexWeightsOmpW(vtot, x);
      if (q) leidenInitializeFromOmpW<REFINE>(vcom, vcob, ctot, x, vtot, *q);
      else   leidenInitializeOmpW(vcom, ctot, x, vtot);
      for (l=0, p=0, s=0; M>0 && p<P;) {
        int m = 0;
        if (s==0) m = leidenMoveOmpW<REFINE>(vcom, ctot, vaff, vpro, vcs, vcout, x, vcob, vtot, M, R, L, fc);
        else      m = leidenMoveOmpW<false> (vcom, ctot, vaff, vpro, vcs, vcout, y, vcob, vtot, M, R, L, fc);
        if (s==0) copyValuesW(a, vcom);
        else      leidenLookupCommunitiesOmpU(a, vcom);
        l += max(m, 1); ++p; ++s;
        if (m<=1 || p>=P) break;
        const G& g = s<=1? x : y;
        size_t gn = g.order();
        size_t yn = leidenCountCommunityVerticesOmpW(cn, g, vcom);
        if (double(yn)/gn >= o.aggregationTolerance) break;
        leidenCommunityVerticesOmpW(co, ce, cn, bufk, g, vcom);
        y = leidenAggregateOmp(vcs, vcout, g, vcom, co, ce);
        fillValueOmpU(vcom, K());
        fillValueOmpU(vtot, W());
        fillValueOmpU(ctot, W());
        fillValueOmpU(vaff, B(1));
        leidenVertexWeightsOmpW(vtot, y);
        leidenInitializeOmpW(vcom, ctot, y, vtot);
        E /= o.toleranceDecline;
      }
    });
  }, o.repeat);
  leidenFreeHashtablesW(vcs, vcout);
  return LeidenResult<K>(a, l, s, t, tm);
}
#endif




// LEIDEN STATIC
// -------------

template <bool JUMP=false, bool REFINE=false, class FLAG=char, class G, class K>
inline auto leidenStaticSeq(const G& x, const vector<K>* q=nullptr, const LeidenOptions& o={}) {
  auto fm = [](auto& vertices) { fillValueU(vertices, FLAG(1)); };
  return leidenSeq<JUMP, REFINE, FLAG>(x, q, o, fm);
}

#ifdef OPENMP
template <bool JUMP=false, bool REFINE=false, class FLAG=char, class G, class K>
inline auto leidenStaticOmp(const G& x, const vector<K>* q=nullptr, const LeidenOptions& o={}) {
  auto fm = [](auto& vertices) { fillValueOmpU(vertices, FLAG(1)); };
  return leidenOmp<JUMP, REFINE, FLAG>(x, q, o, fm);
}
#endif
