#pragma once
#include "update.hxx"




// DUPLICATE IF
// ------------
// Duplicate vertices/edges of graph if test passes.

/**
 * Duplicate vertices/edges of graph if test passes.
 * @param a output graph (updated)
 * @param x original graph
 * @param fv include vertex? (u, d)
 * @param fe include edge? (u, v, w)
 */
template <class H, class G, class FV, class FE>
inline void duplicateIfW(H& a, const G& x, FV fv, FE fe) {
  a.respan(x.span());
  x.forEachVertex   ([&](auto u, auto d) { if (fv(u, d)) a.addVertex(u, d); });
  x.forEachVertexKey([&](auto u) {
    x.forEachEdge(u, [&](auto v, auto w) { if (fe(u, v, w)) a.addEdge(u, v, w); });
  });
  a.update();
}
template <class G, class FV, class FE>
inline G duplicateIf(const G& x, FV fv, FE fe) {
  G a; duplicateIfW(a, x, fv, fe);
  return a;
}


#ifdef OPENMP
/**
 * Duplicate vertices/edges of graph if test passes.
 * @param a output graph (updated)
 * @param x original graph
 * @param fv include vertex? (u, d)
 * @param fe include edge? (u, v, w)
 */
template <class H, class G, class FV, class FE>
inline void duplicateIfOmpW(H& a, const G& x, FV fv, FE fe) {
  a.respan(x.span());
  x.forEachVertex([&](auto u, auto d) { if (fv(u, d)) a.addVertex(u, d); });
  #pragma omp parallel
  {
    x.forEachVertexKey([&](auto u) {
      x.forEachEdge(u, [&](auto v, auto w) { if (fe(u, v, w)) addEdgeOmpU(a, u, v, w); });
    });
  }
  updateOmpU(a);
}
template <class G, class FV, class FE>
inline G duplicateIfOmp(const G& x, FV fv, FE fe) {
  G a; duplicateIfOmpW(a, x, fv, fe);
  return a;
}
#endif




// DUPLICATE
// ---------
// Duplicate vertices/edges of graph (unconditonal love).

/**
 * Duplicate vertices/edges of graph.
 * @param a output graph (updated)
 * @param x original graph
 */
template <class H, class G>
inline void duplicateW(H& a, const G& x) {
  auto fv = [](auto u, auto d)         { return true; };
  auto fe = [](auto u, auto v, auto w) { return true; };
  duplicateIfW(a, x, fv, fe);
}
template <class G>
inline G duplicate(const G& x) {
  G a = x;  // Just use the copy constructor.
  return a;
}
