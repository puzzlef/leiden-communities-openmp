#pragma once
#include "Graph.hxx"
#include "update.hxx"




#pragma region METHODS
#pragma region TRANSPOSE
/**
 * Transpose a graph.
 * @param a transposed graph (output)
 * @param x graph to transpose
 */
template <class H, class G>
inline void transposeW(H& a, const G& x) {
  a.reserve(x.span());
  x.forEachVertex([&](auto u, auto d) { a.addVertex(u, d); });
  x.forEachVertexKey([&](auto u) {
    x.forEachEdge(u, [&](auto v, auto w) { a.addEdge(v, u, w); });
  });
  a.update();
}

/**
 * Transpose a graph.
 * @param x graph to transpose
 * @returns transposed graph
 */
template <class G>
inline auto transpose(const G& x) {
  G a; transposeW(a, x);
  return a;
}


#ifdef OPENMP
/**
 * Transpose a graph in parallel.
 * @param a transposed graph (output)
 * @param x graph to transpose
 */
template <class H, class G>
inline void transposeOmpW(H& a, const G& x) {
  a.reserve(x.span());
  x.forEachVertex([&](auto u, auto d) { a.addVertex(u, d); });
  #pragma omp parallel
  {
    x.forEachVertexKey([&](auto u) {
      x.forEachEdge(u, [&](auto v, auto w) { addEdgeOmpU(a, v, u, w); });
    });
  }
  updateOmpU(a);
}

/**
 * Transpose a graph in parallel.
 * @param x graph to transpose
 * @returns transposed graph
 */
template <class G>
inline auto transposeOmp(const G& x) {
  G a; transposeOmpW(a, x);
  return a;
}
#endif
#pragma endregion




#pragma region TRANSPOSE WITH DEGREE
/**
 * Transpose a graph with degree.
 * @param a transposed graph with degree (output)
 * @param x graph to transpose
 */
template <class H, class G>
inline void transposeWithDegreeW(H& a, const G& x) {
  a.reserve(x.span());
  x.forEachVertexKey([&](auto u) { a.addVertex(u, x.degree(u)); });
  x.forEachVertexKey([&](auto u) {
    x.forEachEdge(u, [&](auto v, auto w) { a.addEdge(v, u, w); });
  });
  a.update();
}

/**
 * Transpose a graph with degree.
 * @param x graph to transpose
 * @returns transposed graph with degree
 */
template <class G>
inline auto transposeWithDegree(const G& x) {
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  using H = DiGraph<K, K, E>;
  H a; transposeWithDegreeW(a, x);
  return a;
}


#ifdef OPENMP
/**
 * Transpose a graph with degree in parallel.
 * @param a transposed graph with degree (output)
 * @param x graph to transpose
 */
template <class H, class G>
inline void transposeWithDegreeOmpW(H& a, const G& x) {
  a.reserve(x.span());
  x.forEachVertexKey([&](auto u) { a.addVertex(u, x.degree(u)); });
  #pragma omp parallel
  {
    x.forEachVertexKey([&](auto u) {
      x.forEachEdge(u, [&](auto v, auto w) { addEdgeOmpU(a, v, u, w); });
    });
  }
  updateOmpU(a);
}

/**
 * Transpose a graph with degree in parallel.
 * @param x graph to transpose
 * @returns transposed graph with degree
 */
template <class G>
inline auto transposeWithDegreeOmp(const G& x) {
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  using H = DiGraph<K, K, E>;
  H a; transposeWithDegreeOmpW(a, x);
  return a;
}
#endif
#pragma endregion
#pragma endregion
