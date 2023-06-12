#include <random>

using std::uniform_real_distribution;




// ADD RANDOM EDGE
// ---------------

template <class G, class R, class V, class FE>
inline bool addRandomEdge(const G& x, R& rnd, size_t i, size_t n, V w, FE fe) {
  using K = typename G::key_type;
  uniform_real_distribution<> dis(0.0, 1.0);
  K u = K(i + n*dis(rnd));
  K v = K(i + n*dis(rnd));
  return fe(u, v, w);
}

template <class G, class R, class V>
inline bool addRandomEdge(G& a, R& rnd, size_t i, size_t n, V w) {
  auto fe = [&](auto u, auto v, auto w) { a.addEdge(u, v, w); return true; };
  return addRandomEdge(a, rnd, i, n, w, fe);
}




// REMOVE RANDOM EDGE
// ------------------

template <class G, class R, class K, class FE>
inline bool removeRandomEdgeFrom(const G& x, R& rnd, K u, FE fe) {
  uniform_real_distribution<> dis(0.0, 1.0);
  if (x.degree(u) == 0) return false;
  K vi = K(dis(rnd) * x.degree(u)), i = 0;
  bool removed = false, skip = false;
  x.forEachEdgeKey(u, [&](auto v) {
    if (skip) return;
    if (i++ == vi) { removed = fe(u, v); skip = true; }
  });
  return removed;
}

template <class G, class R, class K>
inline bool removeRandomEdgeFrom(G& a, R& rnd, K u) {
  auto fe = [&](auto u, auto v) { a.removeEdge(u, v); return true; };
  return removeRandomEdgeFrom(a, rnd, u, fe);
}


template <class G, class R, class FE>
inline bool removeRandomEdge(const G& x, R& rnd, size_t i, size_t n, FE fe) {
  using K = typename G::key_type;
  uniform_real_distribution<> dis(0.0, 1.0);
  K u = K(i + n*dis(rnd));
  return removeRandomEdgeFrom(x, rnd, u, fe);
}

template <class G, class R>
inline bool removeRandomEdge(G& a, R& rnd, size_t i, size_t n) {
  auto fe = [&](auto u, auto v) { a.removeEdge(u, v); return true; };
  return removeRandomEdge(a, rnd, i, n, fe);
}
