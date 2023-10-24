#include <tuple>
#include <vector>
#include <random>
#include "_main.hxx"
#include "update.hxx"

using std::tuple;
using std::vector;
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




#pragma region GENERATE BATCH
/**
 * Add a batch of random edges to a graph.
 * @param a graph to add edges to (updated)
 * @param rnd random number generator
 * @param batchSize number of edges to add
 * @param i begin vertex
 * @param n number of vertices (range)
 * @param w edge weight
 * @returns inserted edges {u, v, w}
 */
template <class G, class R, class V>
inline auto addRandomEdges(G& a, R& rnd, size_t batchSize, size_t i, size_t n, V w) {
  using K = typename G::key_type;
  int retries = 5;
  vector<tuple<K, K, V>> insertions;
  auto fe = [&](auto u, auto v, auto w) {
    a.addEdge(u, v, w);
    a.addEdge(v, u, w);
    insertions.push_back(make_tuple(u, v, w));
    insertions.push_back(make_tuple(v, u, w));
    return true;
  };
  for (size_t l=0; l<batchSize; ++l)
    retry([&]() { return addRandomEdge(a, rnd, i, n, w, fe); }, retries);
  updateOmpU(a);
  return insertions;
}


/**
 * Remove a batch of random edges from a graph.
 * @param a graph to remove edges from (updated)
 * @param rnd random number generator
 * @param batchSize number of edges to remove
 * @param i begin vertex
 * @param n number of vertices (range)
 * @returns deleted edges {u, v}
 */
template <class G, class R>
inline auto removeRandomEdges(G& a, R& rnd, size_t batchSize, size_t i, size_t n) {
  using K = typename G::key_type;
  int retries = 5;
  vector<tuple<K, K>> deletions;
  auto fe = [&](auto u, auto v) {
    a.removeEdge(u, v);
    a.removeEdge(v, u);
    deletions.push_back(make_tuple(u, v));
    deletions.push_back(make_tuple(v, u));
    return true;
  };
  for (size_t l=0; l<batchSize; ++l)
    retry([&]() { return removeRandomEdge(a, rnd, i, n, fe); }, retries);
  updateOmpU(a);
  return deletions;
}
#pragma endregion
#pragma endregion
