#pragma once
#include <utility>
#include <vector>

using std::vector;
using std::swap;




#pragma region METHODS
/**
 * Find vertices visited with BFS.
 * @param vis vertex visited flags (updated)
 * @param us start vertices (updated)
 * @param vs frontier vertices (updated)
 * @param x original graph
 * @param ft should vertex be visited? (vertex, depth)
 * @param fp action to perform on every visited vertex (vertex, depth)
 */
template <class B, class G, class K, class FT, class FP>
inline void bfsVisitedForEachU(vector<B>& vis, vector<K>& us, vector<K>& vs, const G& x, FT ft, FP fp) {
  for (K u : us) {
    if (vis[u] || !ft(u, K())) continue;
    vis[u] = B(1);
    fp(u, K());
  }
  for (K d=1; !us.empty(); ++d) {
    vs.clear();
    for (K u : us) {
      x.forEachEdgeKey(u, [&](K v) {
        if (vis[v] || !ft(v, d)) return;
        vis[v] = B(1);
        vs.push_back(v);
        fp(v, d);
      });
    }
    swap(us, vs);
  }
}


/**
 * Find vertices visited with BFS.
 * @param vis vertex visited flags (updated)
 * @param x original graph
 * @param u start vertex
 * @param ft should vertex be visited? (vertex, depth)
 * @param fp action to perform on every visited vertex (vertex, depth)
 */
template <class B, class G, class K, class FT, class FP>
inline void bfsVisitedForEachU(vector<B>& vis, const G& x, K u, FT ft, FP fp) {
  vector<K> us {u}, vs;
  bfsVisitedForEachU(vis, us, vs, x, ft, fp);
}


/**
 * Find vertices visited with BFS.
 * @tparam FLAG visited flag type
 * @param x original graph
 * @param u start vertex
 * @param ft should vertex be visited? (vertex, depth)
 * @param fp action to perform on every visited vertex (vertex, depth)
 * @returns vertex visited flags
 */
template <class FLAG=bool, class G, class K, class FT, class FP>
inline vector<FLAG> bfsVisitedForEach(const G& x, K u, FT ft, FP fp) {
  vector<FLAG> vis(x.span());
  bfsVisitedForEachU(vis, x, u, ft, fp);
  return vis;
}
#pragma endregion
