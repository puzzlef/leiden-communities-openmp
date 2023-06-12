#pragma once
#include "update.hxx"




// SYMMETRICIZE
// ------------

template <class G>
inline void symmetricizeU(G& a, const G& x) {
  x.forEachVertexKey([&](auto u) {
    x.forEachEdge(u, [&](auto v, auto w) { a.addEdge(v, u, w); });
  });
  a.update();
}

template <class G>
inline void symmetricizeOmpU(G& a, const G& x) {
  #pragma omp parallel
  {
    x.forEachVertexKey([&](auto u) {
      x.forEachEdge(u, [&](auto v, auto w) { addEdgeOmpU(a, v, u, w); });
    });
  }
  updateOmpU(a);
}


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


template <class G>
inline auto symmetricize(const G& x) {
  G a = x; symmetricizeU(a, x);
  return a;
}

template <class G>
inline auto symmetricizeOmp(const G& x) {
  G a = x; symmetricizeOmpU(a, x);
  return a;
}
