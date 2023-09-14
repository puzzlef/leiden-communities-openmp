#pragma once
#include <vector>
#include <cstdint>
#include <cmath>
#include "_main.hxx"
#include "bfs.hxx"
#include "dfs.hxx"
#ifdef OPENMP
#include <omp.h>
#endif

using std::vector;
using std::pow;




#pragma region METHODS
#pragma region GRAPH DATA
/**
 * Obtain the vertex keys of a graph.
 * @param x given graph
 * @returns vertex keys
 */
template <class G>
inline auto vertexKeys(const G& x) {
  using  K = typename G::key_type;
  size_t N = x.order();
  vector<K> a;
  a.reserve(N);
  x.forEachVertexKey([&](auto u) { a.push_back(u); });
  return a;
}


/**
 * Obtain the vertex value of each vertex.
 * @param a vertex value of each vertex (output)
 * @param x given graph
 */
template <class G, class V>
inline void vertexValuesW(vector<V>& a, const G& x) {
  x.forEachVertex([&](auto u, auto d) { a[u] = d; });
}


/**
 * Obtain the outgoing degree of each vertex.
 * @param a degrees of each vertex (output)
 * @param x given graph
 */
template <class G, class K>
inline void degreesW(vector<K>& a, const G& x) {
  x.forEachVertexKey([&](auto u) { a[u] = x.degree(u); });
}
#pragma endregion




#pragma region EDGE WEIGHT
/**
 * Find the total outgoing edge weight of a vertex.
 * @param x given graph
 * @param u given vertex
 * @returns total outgoing weight of a vertex
 */
template <class G, class K>
inline double edgeWeight(const G& x, K u) {
  double a = 0;
  x.forEachEdge(u, [&](auto v, auto w) { a += w; });
  return a;
}


/**
 * Find the total edge weight of a graph.
 * @param x given graph
 * @returns total edge weight (undirected graph => each edge considered twice)
 */
template <class G>
inline double edgeWeight(const G& x) {
  double a = 0;
  x.forEachVertexKey([&](auto u) { a += edgeWeight(x, u); });
  return a;
}


#ifdef OPENMP
/**
 * Find the total edge weight of a graph.
 * @param x given graph
 * @returns total edge weight (undirected graph => each edge considered twice)
 */
template <class G>
inline double edgeWeightOmp(const G& x) {
  using K = typename G::key_type;
  double a = 0;
  size_t S = x.span();
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    a += edgeWeight(x, u);
  }
  return a;
}
#endif
#pragma endregion




#pragma region MODULARITY
/**
 * Find the modularity of a community C.
 * @param cin total weight of edges within community C
 * @param ctot total weight of edges of community C
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns modularity [-0.5, 1]
 * @see https://www.youtube.com/watch?v=0zuiLBOIcsw
 */
inline double modularityCommunity(double cin, double ctot, double M, double R=1) {
  ASSERT(cin>=0 && ctot>=0 && M>0 && R>0);
  return cin/(2*M) - R*pow(ctot/(2*M), 2);
}


/**
 * Find the modularity of a set of communities.
 * @param cin total weight of edges within each community
 * @param ctot total weight of edges of each community
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns modularity [-0.5, 1]
 */
template <class V>
inline double modularityCommunities(const vector<V>& cin, const vector<V>& ctot, double M, double R=1) {
  ASSERT(M>0 && R>0);
  double a = 0;
  for (size_t i=0, I=cin.size(); i<I; ++i)
    a += modularityCommunity(cin[i], ctot[i], M, R);
  return a;
}


#ifdef OPENMP
/**
 * Find the modularity of a set of communities.
 * @param cin total weight of edges within each community
 * @param ctot total weight of edges of each community
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns modularity [-0.5, 1]
 */
template <class V>
inline double modularityCommunitiesOmp(const vector<V>& cin, const vector<V>& ctot, double M, double R=1) {
  ASSERT(M>0 && R>0);
  double a = 0;
  size_t C = cin.size();
  #pragma omp parallel for schedule(static) reduction(+:a)
  for (size_t i=0; i<C; ++i)
    a += modularityCommunity(cin[i], ctot[i], M, R);
  return a;
}
#endif


/**
 * Find the modularity of a graph, based on community membership function.
 * @param x given graph
 * @param fc community membership function of each vertex (u)
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns modularity [-0.5, 1]
 */
template <class G, class FC>
inline double modularityBy(const G& x, FC fc, double M, double R=1) {
  using  K = typename G::key_type;
  ASSERT(M>0 && R>0);
  size_t S = x.span();
  vector<double> cin(S), ctot(S);
  x.forEachVertexKey([&](auto u) {
    K c = fc(u);
    x.forEachEdge(u, [&](auto v, auto w) {
      K d = fc(v);
      if (c==d) cin[c] += w;
      ctot[c] += w;
    });
  });
  return modularityCommunities(cin, ctot, M, R);
}


#ifdef OPENMP
/**
 * Find the modularity of a graph, based on community membership function.
 * @param x given graph
 * @param fc community membership function of each vertex (u)
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns modularity [-0.5, 1]
 */
template <class G, class FC>
inline double modularityByOmp(const G& x, FC fc, double M, double R=1) {
  using  K = typename G::key_type;
  ASSERT(M>0 && R>0);
  size_t S = x.span();
  vector<double> vin(S), vtot(S);
  vector<double> cin(S), ctot(S);
  // Compute the internal and total weight of each vertex.
  #pragma omp parallel for schedule(dynamic, 2048)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    K c = fc(u);
    x.forEachEdge(u, [&](auto v, auto w) {
      K d = fc(v);
      if (c==d) vin[u] += w;
      vtot[u] += w;
    });
  }
  // Compute the internal and total weight of each community.
  #pragma omp parallel for schedule(static, 2048)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    K c = fc(u);
    #pragma omp atomic
    cin[c]  += vin[u];
    #pragma omp atomic
    ctot[c] += vtot[u];
  }
  return modularityCommunitiesOmp(cin, ctot, M, R);
}
#endif
#pragma endregion




#pragma region DELTA MODULARITY
/**
 * Find the change in modularity when moving a vertex from community D to C.
 * @param vcout total weight of edges from vertex v to community C
 * @param vdout total weight of edges from vertex v to community D
 * @param vtot total weight of edges from vertex v
 * @param ctot total weight of edges from community C
 * @param dtot total weight of edges from community C
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns delta-modularity [-0.5, 1]
 * @see https://gist.github.com/wolfram77/a3c95cd94a38a100f9b075594a823928
 */
inline double deltaModularity(double vcout, double vdout, double vtot, double ctot, double dtot, double M, double R=1) {
  ASSERT(vcout>=0 && vdout>=0 && vtot>=0 && ctot>=0 && dtot>=0 && M>0 && R>0);
  return (vcout-vdout)/M - R*vtot*(vtot+ctot-dtot)/(2*M*M);
}
#pragma endregion




#pragma region COMMUNITIES
/**
 * Obtain the size of each community.
 * @param x given graph
 * @param vcom community each vertex belongs to
 * @returns size of each community
 */
template <class G, class K>
inline vector<K> communitySize(const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  vector<K> a(S);
  x.forEachVertexKey([&](auto u) {
    K c = vcom[u];
    ++a[c];
  });
  return a;
}


#ifdef OPENMP
/**
 * Obtain the size of each community.
 * @param x given graph
 * @param vcom community each vertex belongs to
 * @returns size of each community
 */
template <class G, class K>
inline vector<K> communitySizeOmp(const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  vector<K> a(S);
  #pragma omp parallel for schedule(static, 2048)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    K c = vcom[u];
    #pragma omp atomic
    ++a[c];
  }
  return a;
}
#endif


/**
 * Obtain the vertices belonging to each community.
 * @param x given graph
 * @param vcom community each vertex belongs to
 * @returns vertices belonging to each community
 */
template <class G, class K>
inline vector2d<K> communityVertices(const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  vector2d<K> a(S);
  x.forEachVertexKey([&](auto u) {
    K c = vcom[u];
    a[c].push_back(u);
  });
  return a;
}


#ifdef OPENMP
/**
 * Obtain the vertices belonging to each community.
 * @param x given graph
 * @param vcom community each vertex belongs to
 * @returns vertices belonging to each community
 */
template <class G, class K>
inline vector2d<K> communityVerticesOmp(const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  vector2d<K> a(S);
  #pragma omp parallel
  {
    for (K u=0; u<S; ++u) {
      if (!x.hasVertex(u)) continue;
      K c = vcom[u];
      if (belongsOmp(c)) a[c].push_back(u);
    }
  }
  return a;
}
#endif


/**
 * Obtain the community ids of vertices in a graph.
 * @param x given graph
 * @param vcom community each vertex belongs to
 * @returns community ids
 */
template <class G, class K>
inline vector<K> communities(const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  vector<K> a;
  vector<char> vis(S);
  x.forEachVertexKey([&](auto u) {
    K c = vcom[u];
    if (vis[c]) return;
    vis[c] = 1;
    a.push_back(c);
  });
  return a;
}
#pragma endregion




#pragma region DISCONNECTED COMMUNITIES
#ifdef OPENMP
/**
 * Examine if each community in a graph is disconnected (using single flag vector, BFS).
 * @param x given graph
 * @param vcom community each vertex belongs to
 * @returns whether each community is disconnected
 */
template <class G, class K>
inline vector<char> communitiesDisconnectedOmp(const G& x, const vector<K>& vcom) {
  size_t  S = x.span();
  int     T = omp_get_max_threads();
  auto coms = communitySizeOmp(x, vcom);
  vector<char> a(S), vis(S);
  vector2d<K>  us (T), vs(T);
  #pragma omp parallel
  {
    for (K u=0; u<S; ++u) {
      int t = omp_get_thread_num();
      K   c = vcom[u], reached = K();
      if (coms[c]==0 || !belongsOmp(c)) continue;
      auto ft = [&](auto v, auto d) { return vcom[v]==c; };
      auto fp = [&](auto v, auto d) { ++reached; };
      us[t].clear(); vs[t].clear(); us[t].push_back(u);
      bfsVisitedForEachU(vis, us[t], vs[t], x, ft, fp);
      if (reached < coms[c]) a[c] = 1;
      coms[c] = 0;
    }
  }
  return a;
}
#endif
#pragma endregion
#pragma endregion
