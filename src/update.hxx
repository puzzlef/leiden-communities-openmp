#pragma once
#include <utility>
#include <vector>
#include <cstdint>
#ifdef OPENMP
#include <omp.h>
#endif

using std::pair;
using std::vector;




#pragma region METHODS
#pragma region ADD VERTICES
/**
 * Add a range of vertices to a graph.
 * @param a graph to add vertices to
 * @param u begin vertex
 * @param U end vertex (excluding)
 * @param d vertex value
 * @param ft test function to determine if vertex should be added (vertex, vertex value)
 */
template <class G, class K, class V, class FT>
inline void addVerticesIfU(G& a, K u, K U, V d, FT ft) {
  if (U<=1) return;
  a.respan(U);
  for (; u<U; ++u)
    if (ft(u, d)) a.addVertex(u, d);
}

/**
 * Add a range of vertices to a graph.
 * @param a graph to add vertices to
 * @param u begin vertex
 * @param U end vertex (excluding)
 * @param d vertex value
 */
template <class G, class K, class V=typename G::vertex_value_type>
inline void addVerticesU(G& a, K u, K U, V d=V()) {
  auto ft = [](auto u, auto d) { return true; };
  addVerticesIfU(a, u, U, d, ft);
}
#pragma endregion




#pragma region ADD EDGE
/**
 * Add an edge to a graph.
 * @param a graph to add edge to
 * @param u source vertex
 * @param v destination vertex
 * @param w edge weight
 */
template <class G, class K, class E=typename G::edge_value_type>
inline void addEdgeU(G& a, K u, K v, E w=E()) {
  a.addEdge(u, v, w);
}


#ifdef OPENMP
/**
 * Add an edge to a graph in parallel.
 * @param a graph to add edge to
 * @param u source vertex
 * @param v destination vertex
 * @param w edge weight
 */
template <class G, class K, class E=typename G::edge_value_type>
inline void addEdgeOmpU(G& a, K u, K v, E w=E()) {
  auto ft = [](K u) { return belongsOmp(u); };
  a.addEdgeIf(u, v, w, ft);
}
#endif
#pragma endregion




#pragma region REMOVE EDGE
/**
 * Remove an edge from a graph.
 * @param a graph to remove edge from
 * @param u source vertex
 * @param v destination vertex
 */
template <class G, class K>
inline void removeEdgeU(G& a, K u, K v) {
  a.removeEdge(u, v);
}

#ifdef OPENMP
/**
 * Remove an edge from a graph in parallel.
 * @param a graph to remove edge from
 * @param u source vertex
 * @param v destination vertex
 */
template <class G, class K>
inline void removeEdgeOmpU(G& a, K u, K v) {
  auto ft = [](K u) { return belongsOmp(u); };
  a.removeEdgeIf(u, v, ft);
}
#endif
#pragma endregion




#pragma region UPDATE
/**
 * Update changes made to a graph.
 * @param a graph to update
 */
template <class G>
inline void updateU(G& a) {
  a.update();
}


#ifdef OPENMP
/**
 * Update changes made to a graph in parallel.
 * @param a graph to update
 */
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
  #pragma omp parallel for schedule(dynamic, 2048)
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
#pragma endregion
#pragma endregion
