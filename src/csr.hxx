#pragma once
#include <vector>
#include "_main.hxx"

using std::vector;




#pragma region PROPERTIES
/**
 * Get the degree of a vertex in the graph.
 * @param offsets offsets of the outgoing edges of vertices
 * @param u vertex id
 * @returns degree of the vertex
 */
template <class O, class K>
inline K csrDegree(const vector<O>& offsets, K u) {
  return offsets[u+1] - offsets[u];
}
#pragma endregion




#pragma region FOREACH
/**
 * Iterate over the target vertex ids of a source vertex in the graph.
 * @param offsets offsets of the outgoing edges of vertices
 * @param edgeKeys vertex ids of the outgoing edges of each vertex
 * @param u source vertex id
 * @param fp process function (target vertex id)
 */
template <class O, class K, class FP>
inline void csrForEachEdgeKey(const vector<O>& offsets, const vector<K>& edgeKeys, K u, FP fp) {
  O i = offsets[u];
  O I = offsets[u+1];
  for (; i<I; ++i)
    fp(edgeKeys[i]);
}


/**
 * Iterate over the target vertex ids of a source vertex in the graph.
 * @param offsets offsets of the outgoing edges of vertices
 * @param degrees degree of each vertex
 * @param edgeKeys vertex ids of the outgoing edges of each vertex
 * @param u source vertex id
 * @param fp process function (target vertex id)
 */
template <class O, class K, class FP>
inline void csrForEachEdgeKey(const vector<O>& offsets, const vector<K>& degrees, const vector<K>& edgeKeys, K u, FP fp) {
  O i = offsets[u];
  O I = offsets[u] + degrees[u];
  for (; i<I; ++i)
    fp(edgeKeys[i]);
}


/**
 * Iterate over the target vertex ids and edge weights of a source vertex in the graph.
 * @param offsets offsets of the outgoing edges of vertices
 * @param edgeKeys vertex ids of the outgoing edges of each vertex
 * @param edgeValues edge weights of the outgoing edges of each vertex
 * @param u source vertex id
 * @param fp process function (target vertex id, edge weight)
 */
template <class O, class K, class E, class FP>
inline void csrForEachEdge(const vector<O>& offsets, const vector<K>& edgeKeys, const vector<E>& edgeValues, K u, FP fp) {
  O i = offsets[u];
  O I = offsets[u+1];
  for (; i<I; ++i)
    fp(edgeKeys[i], edgeValues[i]);
}


/**
 * Iterate over the target vertex ids and edge weights of a source vertex in the graph.
 * @param offsets offsets of the outgoing edges of vertices
 * @param degrees degree of each vertex
 * @param edgeKeys vertex ids of the outgoing edges of each vertex
 * @param edgeValues edge weights of the outgoing edges of each vertex
 * @param u source vertex id
 * @param fp process function (target vertex id, edge weight)
 */
template <class O, class K, class E, class FP>
inline void csrForEachEdgeKey(const vector<O>& offsets, const vector<K>& degrees, const vector<K>& edgeKeys, const vector<E>& edgeValues, K u, FP fp) {
  O i = offsets[u];
  O I = offsets[u] + degrees[u];
  for (; i<I; ++i)
    fp(edgeKeys[i], edgeValues[i]);
}
#pragma endregion




#pragma region UPDATE
/**
 * Clear the graph.
 * @param offsets offsets of the outgoing edges of vertices
 */
template <class O>
inline void csrClearW(vector<O>& offsets) {
  fillValueU(offsets, O());
}

#ifdef OPENMP
/**
 * Clear the graph.
 * @param offsets offsets of the outgoing edges of vertices
 */
template <class O>
inline void csrClearOmpW(vector<O>& offsets) {
  fillValueOmpU(offsets, O());
}
#endif


/**
 * Clear the graph.
 * @param offsets offsets of the outgoing edges of vertices
 * @param degrees degree of each vertex
 */
template <class O, class K>
inline void csrClearW(vector<O>& offsets, vector<K>& degrees) {
  fillValueU(offsets, O());
  fillValueU(degrees, K());
}

#ifdef OPENMP
/**
 * Clear the graph.
 * @param offsets offsets of the outgoing edges of vertices
 * @param degrees degree of each vertex
 */
template <class O, class K>
inline void csrClearOmpW(vector<O>& offsets, vector<K>& degrees) {
  fillValueOmpU(offsets, O());
  fillValueOmpU(degrees, K());
}
#endif


/**
 * Add an edge to the graph.
 * @param degrees degree of each vertex
 * @param edgeKeys vertex ids of the outgoing edges of each vertex
 * @param offsets offsets of the outgoing edges of vertices
 * @param u source vertex id
 * @param v target vertex id
 * @note Does not check if the edge already exists, or is there is available space.
 */
template <class O, class K>
inline void csrAddEdgeU(vector<K>& degrees, vector<K>& edgeKeys, const vector<O>& offsets, K u, K v) {
  O n = degrees[u]++;
  O i = offsets[u] + n;
  edgeKeys[i] = v;
}

#ifdef OPENMP
/**
 * Add an edge to the graph.
 * @param degrees degree of each vertex
 * @param edgeKeys vertex ids of the outgoing edges of each vertex
 * @param offsets offsets of the outgoing edges of vertices
 * @param u source vertex id
 * @param v target vertex id
 * @note Does not check if the edge already exists, or is there is available space.
 */
template <class O, class K>
inline void csrAddEdgeOmpU(vector<K>& degrees, vector<K>& edgeKeys, const vector<O>& offsets, K u, K v) {
  O n = 0;
  #pragma omp atomic capture
  { n = degrees[u]; ++degrees[u]; }
  O i = offsets[u] + n;
  edgeKeys[i] = v;
}
#endif


/**
 * Add a weighted edge to the graph.
 * @param degrees degree of each vertex
 * @param edgeKeys vertex ids of the outgoing edges of each vertex
 * @param edgeValues edge weights of the outgoing edges of each vertex
 * @param offsets offsets of the outgoing edges of vertices
 * @param u source vertex id
 * @param v target vertex id
 * @param w associated weight of the edge
 * @note Does not check if the edge already exists, or is there is available space.
 */
template <class O, class K, class E>
inline void csrAddEdgeU(vector<K>& degrees, vector<K>& edgeKeys, vector<E>& edgeValues, const vector<O>& offsets, K u, K v, E w) {
  O n = degrees[u]++;
  O i = offsets[u] + n;
  edgeKeys[i]   = v;
  edgeValues[i] = w;
}

#ifdef OPENMP
/**
 * Add a weighted edge to the graph.
 * @param degrees degree of each vertex
 * @param edgeKeys vertex ids of the outgoing edges of each vertex
 * @param edgeValues edge weights of the outgoing edges of each vertex
 * @param offsets offsets of the outgoing edges of vertices
 * @param u source vertex id
 * @param v target vertex id
 * @param w associated weight of the edge
 * @note Does not check if the edge already exists, or is there is available space.
 */
template <class O, class K, class E>
inline void csrAddEdgeOmpU(vector<K>& degrees, vector<K>& edgeKeys, vector<E>& edgeValues, const vector<O>& offsets, K u, K v, E w) {
  O n = 0;
  #pragma omp atomic capture
  { n = degrees[u]; ++degrees[u]; }
  O i = offsets[u] + n;
  edgeKeys[i]   = v;
  edgeValues[i] = w;
}
#endif
#pragma endregion
