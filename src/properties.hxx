#pragma once
#include <vector>
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
#pragma region EDGE WEIGHT
/**
 * Find the total outgoing edge weight of a vertex.
 * @param x original graph
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
 * @param x original graph
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
 * @param x original graph
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




#pragma region DEGREES
/**
 * Find the outgoing degree of each vertex.
 * @param a degrees of each vertex (output)
 * @param x original graph
 * @returns outgoing degree of each vertex
 */
template <class G, class K>
inline void degreesW(vector<K>& a, const G& x) {
  x.forEachVertexKey([&](auto u) { a[u] = x.degree(u); });
}
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
 * @param cin total weight of edges within each community (updated, must be initialized to 0)
 * @param ctot total weight of edges of each community (updated, must be initialized to 0)
 * @param x original graph
 * @param fc community membership function of each vertex (u)
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns modularity [-0.5, 1]
 */
template <class G, class FC>
inline double modularityByW(vector<double>& cin, vector<double>& ctot, const G& x, FC fc, double M, double R=1) {
  ASSERT(M>0 && R>0);
  x.forEachVertexKey([&](auto u) {
    size_t c = fc(u);
    x.forEachEdge(u, [&](auto v, auto w) {
      size_t d = fc(v);
      if (c==d) cin[c] += w;
      ctot[c] += w;
    });
  });
  return modularityCommunities(cin, ctot, M, R);
}


#ifdef OPENMP
/**
 * Find the modularity of a graph, based on community membership function.
 * @param cin total weight of edges within each community (updated, must be initialized to 0)
 * @param ctot total weight of edges of each community (updated, must be initialized to 0)
 * @param x original graph
 * @param fc community membership function of each vertex (u)
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns modularity [-0.5, 1]
 */
template <class G, class FC>
inline double modularityByOmpW(vector2d<double>& cin, vector2d<double>& ctot, const G& x, FC fc, double M, double R=1) {
  using K = typename G::key_type;
  ASSERT(M>0 && R>0);
  size_t S = x.span();
  int    T = omp_get_max_threads();
  #pragma omp parallel for schedule(dynamic, 2048)
  for (K u=0; u<S; ++u) {
    int t = omp_get_thread_num();
    if (!x.hasVertex(u)) continue;
    size_t c = fc(u);
    x.forEachEdge(u, [&](auto v, auto w) {
      size_t d = fc(v);
      if (c==d) cin[t][c] += w;
      ctot[t][c] += w;
    });
  }
  #pragma omp parallel for schedule(auto)
  for (size_t c=0; c<S; ++c) {
    for (int t=1; t<T; ++t) {
      cin [0][c] += cin [t][c];
      ctot[0][c] += ctot[t][c];
    }
  }
  return modularityCommunitiesOmp(cin[0], ctot[0], M, R);
}
#endif


/**
 * Find the modularity of a graph, based on community membership function.
 * @param x original graph
 * @param fc community membership function of each vertex (u)
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns modularity [-0.5, 1]
 */
template <class G, class FC>
inline double modularityBy(const G& x, FC fc, double M, double R=1) {
  size_t S = x.span();
  vector<double> cin(S);
  vector<double> ctot(S);
  return modularityByW(cin, ctot, x, fc, M, R);
}


#ifdef OPENMP
/**
 * Find the modularity of a graph, based on community membership function.
 * @param x original graph
 * @param fc community membership function of each vertex (u)
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns modularity [-0.5, 1]
 */
template <class G, class FC>
inline double modularityByOmp(const G& x, FC fc, double M, double R=1) {
  size_t S = x.span();
  int    T = omp_get_max_threads();
  // Limit memory usage to 64GB.
  size_t VALUES = 64ULL*1024*1024*1024 / 8;
  int    TADJ   = int(max(VALUES / (2*S), size_t(1)));
  vector2d<double> cin (TADJ, vector<double>(S));
  vector2d<double> ctot(TADJ, vector<double>(S));
  // Run in parallel with limited threads
  omp_set_num_threads(TADJ);
  double Q = modularityByOmpW(cin, ctot, x, fc, M, R);
  omp_set_num_threads(T);
  return Q;
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
 * Obtain the vertices belonging to each community.
 * @param x original graph
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
 * @param x original graph
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
 * @param x original graph
 * @param vcom community each vertex belongs to
 * @returns community ids
 */
template <class G, class K>
inline vector<K> communities(const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  vector<K> a;
  vector<char> cflag(S);
  x.forEachVertexKey([&](auto u) {
    K c = vcom[u];
    if (cflag[c]) return;
    cflag[c] = 1;
    a.push_back(c);
  });
  return a;
}
#pragma endregion




#pragma region DISCONNECTED COMMUNITIES
/**
 * Obtain the community ids of vertices in a graph which are disconnected.
 * @param x original graph
 * @param vcom community each vertex belongs to
 * @returns community ids of disconnected communities
 */
template <class G, class K>
inline vector<K> disconnectedCommunities(const G& x, const vector<K>& vcom) {
  size_t  S = x.span();
  vector<char> vis(S);
  vector<K> a;
  auto comv = communityVertices(x, vcom);
  LOG("disconnectedCommunities(): communityVertices() done\n");
  for (K c=0; c<S; ++c) {
    if (comv[c].empty()) continue;
    K u = comv[c][0];
    size_t nvis = 0;
    fillValueU(vis, char());
    auto ft = [&](auto v, auto d) { return vcom[v]==c; };
    auto fp = [&](auto v, auto d) { ++nvis; };
    bfsVisitedForEachW(vis, x, u, ft, fp);
    if (nvis<comv[c].size()) a.push_back(c);
  }
  return a;
}


#ifdef OPENMP
/**
 * Obtain the community ids of vertices in a graph which are disconnected.
 * @param x original graph
 * @param vcom community each vertex belongs to
 * @returns community ids of disconnected communities
 */
template <class G, class K>
inline vector<K> disconnectedCommunitiesOmp(const G& x, const vector<K>& vcom) {
  size_t  S = x.span();
  int     T = omp_get_max_threads();
  vector2d<char> vis(T, vector<char>(S));
  vector2d<K> a(T);
  auto comv = communityVerticesOmp(x, vcom);
  #pragma omp parallel for schedule(auto)
  for (K c=0; c<S; ++c) {
    int t = omp_get_thread_num();
    if (comv[c].empty()) continue;
    K u = comv[c][0];
    size_t nvis = 0;
    fillValueU(vis[t], char());
    auto ft = [&](auto v, auto d) { return vcom[v]==c; };
    auto fp = [&](auto v, auto d) { ++nvis; };
    bfsVisitedForEachW(vis[t], x, u, ft, fp);
    if (nvis<comv[c].size()) a[t].push_back(c);
  }
  for (int t=1; t<T; ++t)
    a[0].insert(a[0].end(), a[t].begin(), a[t].end());
  return a[0];
}
#endif


/**
 * Obtain the community ids of vertices in a graph which are disconnected.
 * @param x original graph
 * @param vcom community each vertex belongs to
 * @returns community ids of disconnected communities
 */
template <class G, class K>
inline vector<K> disconnectedCommunitiesDfsOmp(const G& x, const vector<K>& vcom) {
  size_t  S = x.span();
  int     T = omp_get_max_threads();
  vector2d<char> vis(T, vector<char>(S));
  vector2d<K> a(T);
  auto comv = communityVerticesOmp(x, vcom);
  #pragma omp parallel for schedule(auto)
  for (K c=0; c<S; ++c) {
    int t = omp_get_thread_num();
    if (comv[c].empty()) continue;
    K u = comv[c][0];
    size_t nvis = 0;
    fillValueU(vis[t], char());
    auto ft = [&](auto v) { return vcom[v]==c; };
    auto fp = [&](auto v) { ++nvis; };
    dfsVisitedForEachW(vis[t], x, u, ft, fp);
    if (nvis<comv[c].size()) a[t].push_back(c);
  }
  for (int t=1; t<T; ++t)
    a[0].insert(a[0].end(), a[t].begin(), a[t].end());
  return a[0];
}
#pragma endregion
#pragma endregion
