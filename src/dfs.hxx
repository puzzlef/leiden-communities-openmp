#pragma once
#include <vector>
#include <utility>

using std::vector;




// BFS
// ---
// Depth First Search (DFS) graph traversal.

/**
 * Find vertices visited with DFS.
 * @param vis vertex visited? (updated)
 * @param x original graph
 * @param u start vertex
 * @param ft should vertex be visited? (vertex)
 * @param fp action to perform on every visited vertex (vertex)
 */
template <class B, class G, class K, class FT, class FP>
inline void dfsVisitedForEachW(vector<B>& vis, const G& x, K u, FT ft, FP fp) {
  if (vis[u] || !ft(u)) return;
  vis[u] = B(1); fp(u);
  x.forEachEdgeKey(u, [&](K v) {
    dfsVisitedForEachW(vis, x, v, ft, fp);
  });
}
template <class B=bool, class G, class K, class FT, class FP>
inline vector<B> dfsVisitedForEach(const G& x, K u, FT ft, FP fp) {
  vector<B> vis(x.span());
  dfsVisitedForEachW(vis, x, u, ft, fp);
  return vis;
}
