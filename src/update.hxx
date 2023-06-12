#pragma once
#include <cstdint>
#include <utility>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

using std::pair;
using std::vector;




// ADD VERTICES IF
// ---------------
// Add a range of vertices.

template <class G, class K, class V, class FT>
inline void addVerticesIfU(G& a, K u, K U, V d, FT ft) {
  if (U<=1) return;
  a.respan(U);
  for (; u<U; ++u)
    if (ft(u, d)) a.addVertex(u, d);
}


template <class G, class K, class V=typename G::vertex_value_type>
inline void addVerticesU(G& a, K u, K U, V d=V()) {
  auto ft = [](auto u, auto d) { return true; };
  addVerticesIfU(a, u, U, d, ft);
}




// ADD EDGE
// --------
// Add an edge (in parallel).

template <class G, class K, class E, class FT>
inline void addEdgeU(G &a, K u, K v, E w, FT ft) {
  a.addEdge(u, v, w, ft);
}
template <class G, class K, class E=typename G::edge_value_type>
inline void addEdgeU(G& a, K u, K v, E w=E()) {
  a.addEdge(u, v, w);
}


#ifdef OPENMP
template <class G, class K, class E=typename G::edge_value_type>
inline void addEdgeOmpU(G& a, K u, K v, E w=E()) {
  auto ft = [](K u) { return belongsOmp(u); };
  a.addEdge(u, v, w, ft);
}
#endif




// REMOVE EDGE
// -----------
// Remove an edge (in parallel).

template <class G, class K, class FT>
inline void removeEdgeU(G &a, K u, K v, FT ft) {
  a.removeEdge(u, v, ft);
}
template <class G, class K>
inline void removeEdgeU(G& a, K u, K v) {
  a.removeEdge(u, v);
}


#ifdef OPENMP
template <class G, class K>
inline void removeEdgeOmpU(G& a, K u, K v) {
  auto ft = [](K u) { return belongsOmp(u); };
  a.removeEdge(u, v, ft);
}
#endif




// UPDATE
// ------
// Update changes made to a graph.

template <class G>
inline void updateU(G& a) {
  a.update();
}


#ifdef OPENMP
template <class G>
inline void updateOmpU(G& a) {
  using  K = typename G::key_type;
  using  E = typename G::edge_value_type;
  size_t S = a.span();
  // Create per-thread buffers for update operation.
  int THREADS = omp_get_max_threads();
  vector<vector<pair<K, E>>*> bufs(THREADS);
  for (int i=0; i<THREADS; ++i)
    bufs[i] = new vector<pair<K, E>>();
  // Update edges of each vertex individually.
  #pragma omp parallel for schedule(auto)
  for (K u=0; u<S; ++u) {
    int t = omp_get_thread_num();
    a.updateEdges(u, bufs[t]);
  }
  // Update the entire graph, find total vertices and edges.
  a.update();
  // Clean up.
  for (int i=0; i<THREADS; ++i)
    delete bufs[i];
}
#endif
