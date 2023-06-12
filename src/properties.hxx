#pragma once




// WEIGHT
// ------

/**
 * Find the total outgoing edge weight of a vertex.
 * @param x original graph
 * @param u given vertex
 * @returns total outgoing weight of a vertex
 */
template <class G, class K>
inline double edgeWeight(const G& x, K u) {
  double a = 0;
  x.forEachEdgeValue(u, [&](auto w) { a += w; });
  return a;
}


/**
 * Find the total edge weight of a graph.
 * @param x original graph
 * @returns total edge weight (undirected graph => each edge considered twice)
 */
template <class G>
inline double edgeWeight(const G& x) {
  double a = 0;
  x.forEachVertexKey([&](auto u) { a += edgeWeight(x, u); });
  return a;
}

template <class G>
inline double edgeWeightOmp(const G& x) {
  using K = typename G::key_type;
  double a = 0;
  size_t S = x.span();
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    a += edgeWeight(x, u);
  }
  return a;
}


/**
 * Find the outgoing degree of each vertex.
 * @param a degrees of each vertex (output)
 * @param x original graph
 * @returns outgoing degree of each vertex
 */
template <class G, class K>
inline void degreesW(vector<K>& a, const G& x) {
  x.forEachVertexKey([&](auto u) { a[u] = x.degree(u); });
}
