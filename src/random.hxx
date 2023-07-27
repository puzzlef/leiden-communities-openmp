#include <random>

using std::uniform_real_distribution;




#pragma region METHODS
#pragma region ADD RANDOM EDGE
/**
 * Randomly decide upon a edge that may be added to a graph.
 * @param x input graph
 * @param rnd random number generator
 * @param i begin vertex
 * @param n number of vertices (range)
 * @param w edge weight
 * @param fe callback function to add edge
 * @returns true if edge was added
 */
template <class G, class R, class V, class FE>
inline bool addRandomEdge(const G& x, R& rnd, size_t i, size_t n, V w, FE fe) {
  using K = typename G::key_type;
  uniform_real_distribution<> dis(0.0, 1.0);
  K u = K(i + n*dis(rnd));
  K v = K(i + n*dis(rnd));
  return fe(u, v, w);
}

/**
 * Add a random edge to a graph.
 * @param a graph to add edge to (updated)
 * @param rnd random number generator
 * @param i begin vertex
 * @param n number of vertices (range)
 * @param w edge weight
 * @returns true if edge was added
 */
template <class G, class R, class V>
inline bool addRandomEdge(G& a, R& rnd, size_t i, size_t n, V w) {
  auto fe = [&](auto u, auto v, auto w) { a.addEdge(u, v, w); return true; };
  return addRandomEdge(a, rnd, i, n, w, fe);
}
#pragma endregion




#pragma region REMOVE RANDOM EDGE
/**
 * Randomly decide upon a edge that may be removed from a graph, from a given source vertex.
 * @param x input graph
 * @param rnd random number generator
 * @param u source vertex
 * @param fe callback function to remove edge
 * @returns true if edge was removed
 */
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

/**
 * Remove a random edge from a graph, from a given source vertex.
 * @param a graph to remove edge from (updated)
 * @param rnd random number generator
 * @param u source vertex
 * @returns true if edge was removed
 */
template <class G, class R, class K>
inline bool removeRandomEdgeFrom(G& a, R& rnd, K u) {
  auto fe = [&](auto u, auto v) { a.removeEdge(u, v); return true; };
  return removeRandomEdgeFrom(a, rnd, u, fe);
}

/**
 * Randomly decide upon a edge that may be removed from a graph.
 * @param x input graph
 * @param rnd random number generator
 * @param i begin vertex
 * @param n number of vertices (range)
 * @param fe callback function to remove edge
 * @returns true if edge was removed
 */
template <class G, class R, class FE>
inline bool removeRandomEdge(const G& x, R& rnd, size_t i, size_t n, FE fe) {
  using K = typename G::key_type;
  uniform_real_distribution<> dis(0.0, 1.0);
  K u = K(i + n*dis(rnd));
  return removeRandomEdgeFrom(x, rnd, u, fe);
}

/**
 * Remove a random edge from a graph.
 * @param a graph to remove edge from (updated)
 * @param rnd random number generator
 * @param i begin vertex
 * @param n number of vertices (range)
 * @returns true if edge was removed
 */
template <class G, class R>
inline bool removeRandomEdge(G& a, R& rnd, size_t i, size_t n) {
  auto fe = [&](auto u, auto v) { a.removeEdge(u, v); return true; };
  return removeRandomEdge(a, rnd, i, n, fe);
}
#pragma endregion
#pragma endregion
