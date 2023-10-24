#pragma once
#include "update.hxx"




#pragma region METHODS
/**
 * Obtain the symmetric version of a graph.
 * @param a output symmetric graph (empty, updated)
 * @param x input graph
 */
template <class H, class G>
inline void symmetricizeW(H& a, const G& x) {
  a.reserve(x.span());
  x.forEachVertex([&](auto u, auto d) { a.addVertex(u, d); });
  x.forEachVertexKey([&](auto u) {
    x.forEachEdge(u, [&](auto v, auto w) {
      a.addEdge(u, v, w);
      a.addEdge(v, u, w);
    });
  });
  a.update();
}

#ifdef OPENMP
/**
 * Obtain the symmetric version of a graph in parallel.
 * @param a output symmetric graph (empty, updated)
 * @param x input graph
 */
template <class H, class G>
inline void symmetricizeOmpW(H& a, const G& x) {
  a.reserve(x.span());
  x.forEachVertex([&](auto u, auto d) { a.addVertex(u, d); });
  #pragma omp parallel
  {
    x.forEachVertexKey([&](auto u) {
      x.forEachEdge(u, [&](auto v, auto w) {
        addEdgeOmpU(a, u, v, w);
        addEdgeOmpU(a, v, u, w);
      });
    });
  }
  updateOmpU(a);
}
#endif


/**
 * Obtain the symmetric version of a graph.
 * @param x input graph
 * @return symmetric graph
 */
template <class G>
inline auto symmetricize(const G& x) {
  G a = x;
  x.forEachVertexKey([&](auto u) {
    x.forEachEdge(u, [&](auto v, auto w) { a.addEdge(v, u, w); });
  });
  a.update();
  return a;
}

#ifdef OPENMP
/**
 * Obtain the symmetric version of a graph in parallel.
 * @param x input graph
 * @return symmetric graph
 */
template <class G>
inline auto symmetricizeOmp(const G& x) {
  G a = x;
  #pragma omp parallel
  {
    x.forEachVertexKey([&](auto u) {
      x.forEachEdge(u, [&](auto v, auto w) { addEdgeOmpU(a, v, u, w); });
    });
  }
  updateOmpU(a);
  return a;
}
#endif
#pragma endregion
