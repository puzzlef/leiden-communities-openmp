#pragma once
#include <vector>
#include <unordered_map>
#include "_main.hxx"

using std::vector;
using std::unordered_map;




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
 * @param edgeValues edge values of the outgoing edges of each vertex
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
 * @param edgeValues edge values of the outgoing edges of each vertex
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




#pragma region CREATE
/**
 * Obtain offsets of the outgoing edges of vertices.
 * @param offsets offsets of the outgoing edges of vertices (output)
 * @param x given graph
 */
template <class G, class O>
inline void csrCreateOffsetsW(vector<O>& offsets, const G& x) {
  size_t N = x.order();
  O o = O();
  offsets.clear();
  offsets.reserve(N+1);
  x.forEachVertexKey([&](auto u) {
    offsets.push_back(o);
    o += x.degree(u);
  });
  offsets.push_back(o);
}


/**
 * Obtain offsets of the outgoing edges of vertices, for a subset of vertices.
 * @param offsets offsets of the outgoing edges of vertices (output)
 * @param x given graph
 * @param ks vertex keys to include
 */
template <class G, class O, class K>
inline void csrCreateOffsetsW(vector<O>& offsets, const G& x, const vector<K>& ks) {
  size_t KS = ks.size();
  O o = O();
  offsets.clear();
  offsets.reserve(KS+1);
  for (auto u : ks) {
    offsets.push_back(o);
    o += x.degree(u);
  }
  offsets.push_back(o);
}


/**
 * Obtain degree of each vertex.
 * @param degrees degree of each vertex (output)
 * @param x given graph
 */
template <class G, class K>
inline void csrCreateDegreesW(vector<K>& degrees, const G& x) {
  size_t N = x.order();
  degrees.clear();
  degrees.reserve(N);
  x.forEachVertexKey([&](auto u) {
    degrees.push_back(x.degree(u));
  });
}


/**
 * Obtain degree of each vertex, for a subset of vertices.
 * @param degrees degree of each vertex (output)
 * @param x given graph
 * @param ks vertex keys to include
 */
template <class G, class K>
inline void csrCreateDegreesW(vector<K>& degrees, const G& x, const vector<K>& ks) {
  size_t KS = ks.size();
  degrees.clear();
  degrees.reserve(KS);
  for (auto u : ks)
    degrees.push_back(x.degree(u));
}


/**
 * Obtain vertex value of each vertex.
 * @param vertexValues vertex value of each vertex (output)
 * @param x given graph
 */
template <class G, class V>
inline void csrCreateVertexValuesW(vector<V>& vertexValues, const G& x) {
  size_t N = x.order();
  vertexValues.clear();
  vertexValues.reserve(N);
  x.forEachVertex([&](auto u, auto d) {
    vertexValues.push_back(d);
  });
}


/**
 * Obtain vertex value of each vertex, for a subset of vertices.
 * @param vertexValues vertex value of each vertex (output)
 * @param x given graph
 * @param ks vertex keys to include
 */
template <class G, class K, class V>
inline void csrCreateVertexValuesW(vector<V>& vertexValues, const G& x, const vector<K>& ks) {
  size_t KS = ks.size();
  vertexValues.clear();
  vertexValues.reserve(KS);
  for (auto u : ks)
    vertexValues.push_back(x.vertexValue(u));
}


/**
 * Obtain vertex ids of the outgoing edges of each vertex.
 * @param edgeKeys vertex ids of the outgoing edges of each vertex (output)
 * @param x given graph
 */
template <class G, class K>
inline void csrCreateEdgeKeysW(vector<K>& edgeKeys, const G& x) {
  size_t M = x.size();
  K i = K();
  unordered_map<K, K> ids;
  x.forEachVertexKey([&](auto u) { ids[u] = i++; });
  edgeKeys.clear();
  edgeKeys.reserve(M);
  x.forEachVertexKey([&](auto u) {
    x.forEachEdgeKey(u, [&](auto v) {
      edgeKeys.push_back(ids[v]);
    });
  });
}


/**
 * Obtain vertex ids of the outgoing edges of each vertex, for a subset of vertices.
 * @param edgeKeys vertex ids of the outgoing edges of each vertex (output)
 * @param x given graph
 * @param ks vertex keys to include
 */
template <class G, class K>
inline void csrCreateEdgeKeysW(vector<K>& edgeKeys, const G& x, const vector<K>& ks) {
  size_t M = 0;
  K i = K();
  unordered_map<K, K> ids;
  for (auto u : ks) {
    M += x.degree(u);
    ids[u] = i++;
  }
  edgeKeys.clear();
  edgeKeys.reserve(M);
  for (auto u : ks) {
    x.forEachEdgeKey(u, [&](auto v) {
      edgeKeys.push_back(ids[v]);
    });
  }
}


/**
 * Obtain edge values of the outgoing edges of each vertex.
 * @param edgeValues edge values of the outgoing edges of each vertex (output)
 * @param x given graph
 */
template <class G, class E>
inline void csrCreateEdgeValuesW(vector<E>& edgeValues, const G& x) {
  size_t M = x.size();
  edgeValues.clear();
  edgeValues.reserve(M);
  x.forEachVertexKey([&](auto u) {
    x.forEachEdge(u, [&](auto v, auto w) {
      edgeValues.push_back(w);
    });
  });
}


/**
 * Obtain edge values of the outgoing edges of each vertex, for a subset of vertices.
 * @param edgeValues edge values of the outgoing edges of each vertex (output)
 * @param x given graph
 * @param ks vertex keys to include
 */
template <class G, class K, class E>
inline void csrCreateEdgeValuesW(vector<E>& edgeValues, const G& x, const vector<K>& ks) {
  size_t M = 0;
  for (auto u : ks)
    M += x.degree(u);
  edgeValues.clear();
  edgeValues.reserve(M);
  for (auto u : ks) {
    x.forEachEdge(u, [&](auto v, auto w) {
      edgeValues.push_back(w);
    });
  }
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
 * @param edgeValues edge values of the outgoing edges of each vertex
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
 * @param edgeValues edge values of the outgoing edges of each vertex
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
