#pragma once
#include "update.hxx"




// HAS SELF-LOOP
// -------------

template <class G, class K>
inline bool hasSelfLoop(const G& x, K u) {
  return x.hasEdge(u, u);
}




// SELF-LOOPS
// ----------

template <class G, class F>
inline void selfLoopForEach(const G& x, F fn) {
  x.forEachVertexKey([&](auto u) { if (hasSelfLoop(x, u)) fn(u); });
}
template <class G>
inline auto selfLoops(const G& x) {
  using K = typename G::key_type; vector<K> a;
  selfLoopForEach(x, [&](auto u) { a.push_back(u); });
  return a;
}
template <class G>
inline auto selfLoopCount(const G& x) {
  using K = typename G::key_type; K a = 0;
  selfLoopForEach(x, [&](auto u) { ++a; });
  return a;
}




// SELF-LOOP
// ---------

template <class G, class E, class FT>
inline void selfLoopU(G& a, E w, FT ft) {
  a.forEachVertexKey([&](auto u) { if (ft(u)) a.addEdge(u, u, w); });
  a.update();
}
template <class G, class E, class FT>
inline auto selfLoop(const G& x, E w, FT ft) {
  G a = x; selfLoopU(a, w, ft);
  return a;
}


#ifdef OPENMP
template <class G, class E, class FT>
inline void selfLoopOmpU(G& a, E w, FT ft) {
  #pragma omp parallel
  {
    a.forEachVertexKey([&](auto u) { if (ft(u)) addEdgeOmpU(a, u, u, w); });
  }
  updateOmpU(a);
}
template <class G, class E, class FT>
inline auto selfLoopOmp(const G& x, E w, FT ft) {
  G a = x; selfLoopOmpU(a, w, ft);
  return a;
}
#endif
