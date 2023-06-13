#pragma once
#include <vector>
#include "_main.hxx"
#include "bfs.hxx"
#include "dfs.hxx"

using std::vector;




// WEIGHT
// ------

/**
 * Find the total outgoing edge weight of a vertex.
 * @param x original graph
 * @param u given vertex
 * @returns total outgoing weight of a vertex
 */
template <class G, class K>
inline double edgeWeight(const G& x, K u) {
  double a = 0;
  x.forEachEdgeValue(u, [&](auto w) { a += w; });
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
