#pragma once
#include <cstdint>
#include <vector>
#ifdef OPENMP
#include <omp.h>
#endif
#include "dfs.hxx"
#include "bfs.hxx"

using std::vector;




#pragma region SPLIT DISCONNECTED COMMUNITIES
#ifdef OPENMP
/**
 * Split disconnected communities using Label Propagation Algorithm (LPA).
 * @param vcom label/subcommunity each vertex belongs to (output)
 * @param vaff whether each vertex is affected, if pruning is enabled (scratch)
 * @param x given graph
 * @param vdom community each vertex belongs to
 */
template <bool PRUNE=false, class B, class G, class K>
inline void splitDisconnectedCommunitiesLpaOmpW(vector<K>& vcom, vector<B>& vaff, const G& x, const vector<K>& vdom) {
  size_t S = x.span();
  // Initialize each vertex to its own label/subcommunity.
  #pragma omp parallel for schedule(static, 2048)
  for (K u=0; u<S; ++u) {
    vcom[u] = u;
    if constexpr (PRUNE) vaff[u] = B(1);
  }
  // Perform label propagation within each community.
  while (true) {
    size_t ndel = 0;
    #pragma omp parallel for schedule(dynamic, 2048) reduction(+:ndel)
    for (K u=0; u<S; ++u) {
      if (!x.hasVertex(u)) continue;
      if (PRUNE && !vaff[u]) continue;
      K d = vdom[u];
      K c = vcom[u];
      // Find the minimum label of all neighbors in the same community.
      x.forEachEdgeKey(u, [&](auto v) {
        if (vdom[v]==d) c = min(c, vcom[v]);
      });
      if (c==vcom[u]) continue;
      // Update the label of this vertex.
      vcom[u] = c;
      ++ndel;
      if constexpr (!PRUNE) continue;
      vaff[u] = B();
      x.forEachEdgeKey(u, [&](auto v) { if (vdom[v]==d && !vaff[v]) vaff[v] = B(1); });
    }
    if (ndel==0) break;
  }
}


/**
 * Split disconnected communities using DFS.
 * @param vcom label/subcommunity each vertex belongs to (output)
 * @param vis vertex visited flags (scratch)
 * @param x given graph
 * @param vdom community each vertex belongs to
 */
template <class B, class G, class K>
inline void splitDisconnectedCommunitiesDfsOmpW(vector<K>& vcom, vector<B>& vis, const G& x, const vector<K>& vdom) {
  size_t S = x.span();
  // Initialize each vertex to its own label/subcommunity.
  #pragma omp parallel for schedule(static, 2048)
  for (K u=0; u<S; ++u) {
    vcom[u] = u;
    vis[u]  = B();
  }
  // Perform DFS from an untouched vertex, within each community (each thread picks a community atomically).
  #pragma omp parallel
  {
    int t = omp_get_thread_num();
    int T = omp_get_num_threads();
    for (K u=0; u<S; ++u) {
      if (!x.hasVertex(u)) continue;
      K d = vdom[u];
      K c = vcom[u];
      if (!belongsOmp(d, t, T) || vis[u]) continue;
      auto ft = [&](auto v) { return vdom[v]==d; };
      auto fp = [&](auto v) { vcom[v] = c; };
      dfsVisitedForEachU(vis, x, u, ft, fp);
    }
  }
}


/**
 * Split disconnected communities using BFS.
 * @param vcom label/subcommunity each vertex belongs to (output)
 * @param vis vertex visited flags (scratch)
 * @param us per-thread start vertices for BFS (scratch)
 * @param vs per-thread frontier vertices for BFS (scratch)
 * @param x given graph
 * @param vdom community each vertex belongs to
 */
template <class B, class G, class K>
inline void splitDisconnectedCommunitiesBfsOmpW(vector<K>& vcom, vector<B>& vis, vector<vector<K>*>& us, vector<vector<K>*>& vs, const G& x, const vector<K>& vdom) {
  size_t S = x.span();
  // Initialize each vertex to its own label/subcommunity.
  #pragma omp parallel for schedule(static, 2048)
  for (K u=0; u<S; ++u) {
    vcom[u] = u;
    vis[u]  = B();
  }
  // Perform DFS from an untouched vertex, within each community (each thread picks a community atomically).
  #pragma omp parallel
  {
    int t = omp_get_thread_num();
    int T = omp_get_num_threads();
    for (K u=0; u<S; ++u) {
      if (!x.hasVertex(u)) continue;
      K d = vdom[u];
      K c = vcom[u];
      if (!belongsOmp(d, t, T) || vis[u]) continue;
      auto ft = [&](auto v, auto _) { return vdom[v]==d; };
      auto fp = [&](auto v, auto _) { vcom[v] = c; };
      (*us[t]).clear(); (*vs[t]).clear(); (*us[t]).push_back(u);
      bfsVisitedForEachU(vis, *us[t], *vs[t], x, ft, fp);
    }
  }
}
#endif
#pragma endregion
