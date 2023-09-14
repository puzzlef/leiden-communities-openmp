#pragma once
#include <utility>
#include <vector>
#include <ostream>
#include <algorithm>
#include "_main.hxx"

using std::pair;
using std::vector;
using std::ostream;
using std::max;




#pragma region CLASSES
/**
 * Directed graph that memorizes only out-edges for each vertex.
 * @tparam K key type (vertex id)
 * @tparam V vertex value type (vertex data)
 * @tparam E edge value type (edge weight)
 */
template <class K=uint32_t, class V=None, class E=None>
class DiGraph {
  #pragma region TYPES
  public:
  /** Key type (vertex id). */
  using key_type = K;
  /** Vertex value type (vertex data). */
  using vertex_value_type = V;
  /** Edge value type (edge weight). */
  using edge_value_type   = E;
  #pragma endregion


  #pragma region DATA
  protected:
  /** Number of vertices. */
  size_t N = 0;
  /** Number of edges. */
  size_t M = 0;
  /** Vertex existence flags. */
  vector<bool> exists;
  /** Vertex values. */
  vector<V> values;
  /** Outgoing edges for each vertex (including edge weights). */
  vector<LazyBitset<K, E>> edges;
  #pragma endregion


  #pragma region METHODS
  #pragma region PROPERTIES
  public:
  /**
   * Get the size of buffer required to store data associated with each vertex
   * in the graph, indexed by its vertex-id.
   * @returns size of buffer required
   */
  inline size_t span() const noexcept {
    return exists.size();
  }

  /**
   * Get the number of vertices in the graph.
   * @returns |V|
   */
  inline size_t order() const noexcept {
    return N;
  }

  /**
   * Get the number of edges in the graph.
   * @returns |E|
   */
  inline size_t size() const noexcept {
    return M;
  }

  /**
   * Check if the graph is empty.
   * @returns is the graph empty?
   */
  inline bool empty() const noexcept {
    return N == 0;
  }

  /**
   * Check if the graph is directed.
   * @returns is the graph directed?
   */
  inline bool directed() const noexcept {
    return true;
  }
  #pragma endregion


  #pragma region FOREACH
  public:
  /**
   * Iterate over the vertices in the graph.
   * @param fp process function (vertex id, vertex data)
   */
  template <class FP>
  inline void forEachVertex(FP fp) const noexcept {
    for (K u=0; u<span(); ++u)
      if (exists[u]) fp(u, values[u]);
  }

  /**
   * Iterate over the vertex ids in the graph.
   * @param fp process function (vertex id)
   */
  template <class FP>
  inline void forEachVertexKey(FP fp) const noexcept {
    for (K u=0; u<span(); ++u)
      if (exists[u]) fp(u);
  }

  /**
   * Iterate over the outgoing edges of a source vertex in the graph.
   * @param u source vertex id
   * @param fp process function (target vertex id, edge weight)
   */
  template <class FP>
  inline void forEachEdge(K u, FP fp) const noexcept {
    edges[u].forEach(fp);
  }

  /**
   * Iterate over the target vertex ids of a source vertex in the graph.
   * @param u source vertex id
   * @param fp process function (target vertex id)
   */
  template <class FP>
  inline void forEachEdgeKey(K u, FP fp) const noexcept {
    edges[u].forEachKey(fp);
  }
  #pragma endregion


  #pragma region ACCESS
  public:
  /**
   * Check if a vertex exists in the graph.
   * @param u vertex id
   * @returns does the vertex exist?
   */
  inline bool hasVertex(K u) const noexcept {
    return u < span() && exists[u];
  }

  /**
   * Check if an edge exists in the graph.
   * @param u source vertex id
   * @param v target vertex id
   * @returns does the edge exist?
   */
  inline bool hasEdge(K u, K v) const noexcept {
    return u < span() && edges[u].has(v);
  }

  /**
   * Get the number of outgoing edges of a vertex in the graph.
   * @param u vertex id
   * @returns number of outgoing edges of the vertex
   */
  inline size_t degree(K u) const noexcept {
    return u < span()? edges[u].size() : 0;
  }

  /**
   * Get the vertex data of a vertex in the graph.
   * @param u vertex id
   * @returns associated data of the vertex
   */
  inline V vertexValue(K u) const noexcept {
    return u < span()? values[u] : V();
  }

  /**
   * Set the vertex data of a vertex in the graph.
   * @param u vertex id
   * @param d associated data of the vertex
   * @returns success?
   */
  inline bool setVertexValue(K u, V d) noexcept {
    if (!hasVertex(u)) return false;
    values[u] = d;
    return true;
  }

  /**
   * Get the edge weight of an edge in the graph.
   * @param u source vertex id
   * @param v target vertex id
   * @returns associated weight of the edge
   */
  inline E edgeValue(K u, K v) const noexcept {
    return u < span()? edges[u].get(v) : E();
  }

  /**
   * Set the edge weight of an edge in the graph.
   * @param u source vertex id
   * @param v target vertex id
   * @param w associated weight of the edge
   * @returns success?
   */
  inline bool setEdgeValue(K u, K v, E w) noexcept {
    if (!hasVertex(u) || !hasVertex(v)) return false;
    return edges[u].set(v, w);
  }
  #pragma endregion


  #pragma region UPDATE
  public:
  /**
   * Remove all vertices and edges from the graph.
   */
  inline void clear() noexcept {
    N = 0; M = 0;
    exists.clear();
    values.clear();
    edges.clear();
  }

  /**
   * Reserve space for outgoing edges of a vertex in the graph.
   * @param u source vertex id
   * @param deg expected degree of the vertex
   */
  inline void reserveEdges(K u, size_t deg) {
    if (u < span()) edges[u].reserve(deg);
  }


  /**
   * Reserve space for a number of vertices and edges in the graph.
   * @param n number of vertices to reserve space for
   * @param deg expected average degree of vertices
   */
  inline void reserve(size_t n, size_t deg=0) {
    size_t S = max(n, span());
    exists.resize(S);
    values.resize(S);
    edges.resize(S);
    if (deg==0) return;
    for (K u=0; u<S; ++u)
      edges[u].reserve(deg);
  }

  /**
   * Adjust the span of the graph.
   * @param n new span
   * @note This operation is lazy.
   */
  inline void respan(size_t n) {
    exists.resize(n);
    values.resize(n);
    edges.resize(n);
  }

  /**
   * Update the outgoing edges of a vertex in the graph to reflect the changes.
   * @param u source vertex id
   * @param buf scratch buffer for the update
   */
  inline void updateEdges(K u, vector<pair<K, E>> *buf=nullptr) {
    if (u < span()) edges[u].update(buf);
  }

  /**
   * Update the graph to reflect the changes.
   * @note This is an expensive operation.
   */
  inline void update() {
    vector<pair<K, E>> buf;
    N = 0; M = 0;
    forEachVertexKey([&](K u) {
      edges[u].update(&buf);
      M += degree(u); ++N;
    });
  }

  /**
   * Add a vertex to the graph.
   * @param u vertex id
   * @note This operation is lazy.
   */
  inline void addVertex(K u) {
    if (hasVertex(u)) return;
    if (u >= span()) respan(u+1);
    exists[u] = true;
  }

  /**
   * Add a vertex to the graph.
   * @param u vertex id
   * @param d associated data of the vertex
   * @note This operation is lazy.
   */
  inline void addVertex(K u, V d) {
    if (hasVertex(u)) { values[u] = d; return; }
    if (u >= span()) respan(u+1);
    exists[u] = true;
    values[u] = d;
  }

  /**
   * Add an outgoing edge to the graph if a condition is met.
   * @param u source vertex id
   * @param v target vertex id
   * @param w associated weight of the edge
   * @param ft test function (source vertex id)
   */
  template <class FT>
  inline void addEdgeIf(K u, K v, E w, FT ft) {
    addVertex(u);
    addVertex(v);
    if (ft(u)) edges[u].add(v, w);
  }

  /**
   * Add an outgoing edge to the graph.
   * @param u source vertex id
   * @param v target vertex id
   * @param w associated weight of the edge
   * @note This operation is lazy.
   */
  inline void addEdge(K u, K v, E w=E()) {
    auto ft = [](K u) { return true; };
    addEdgeIf(u, v, w, ft);
  }

  /**
   * Remove an outgoing edge from the graph if a condition is met.
   * @param u source vertex id
   * @param v target vertex id
   * @param ft test function (source vertex id)
   */
  template <class FT>
  inline void removeEdgeIf(K u, K v, FT ft) {
    if (!hasVertex(u) || !hasVertex(v)) return;
    if (ft(u)) edges[u].remove(v);
  }

  /**
   * Remove an outgoing edge from the graph.
   * @param u source vertex id
   * @param v target vertex id
   * @note This operation is lazy.
   */
  inline void removeEdge(K u, K v) {
    auto ft = [](K u) { return true; };
    removeEdgeIf(u, v, ft);
  }

  /**
   * Remove a vertex from the graph.
   * @param u vertex id
   * @note This operation is lazy.
   */
  inline void removeVertex(K u) {
    if (!hasVertex(u)) return;
    exists[u] = false;
    values[u] = V();
    edges[u].clear();
  }
  #pragma endregion
  #pragma endregion
};



/**
 * A directed graph with CSR representation.
 * @tparam K key type (vertex id)
 * @tparam V vertex value type (vertex data)
 * @tparam E edge value type (edge weight)
 * @tparam O offset type
 */
template <class K=uint32_t, class V=None, class E=None, class O=size_t>
class DiGraphCsr {
  #pragma region TYPES
  public:
  /** Key type (vertex id). */
  using key_type = K;
  /** Vertex value type (vertex data). */
  using vertex_value_type = V;
  /** Edge value type (edge weight). */
  using edge_value_type   = E;
  #pragma endregion


  #pragma region DATA
  public:
  /** Offsets of the outgoing edges of vertices. */
  vector<O> offsets;
  /** Degree of each vertex. */
  vector<K> degrees;
  /** Vertex values. */
  vector<V> values;
  /** Vertex ids of the outgoing edges of each vertex (lookup using offsets). */
  vector<K> edgeKeys;
  /** Edge weights of the outgoing edges of each vertex (lookup using offsets). */
  vector<E> edgeValues;
  #pragma endregion


  #pragma region METHODS
  #pragma region PROPERTIES
  public:
  /**
   * Get the size of buffer required to store data associated with each vertex
   * in the graph, indexed by its vertex-id.
   * @returns size of buffer required
   */
  inline size_t span() const noexcept {
    return degrees.size();
  }

  /**
   * Get the number of vertices in the graph.
   * @returns |V|
   */
  inline size_t order() const noexcept {
    return degrees.size();
  }

  /**
   * Obtain the number of edges in the graph.
   * @returns |E|
   */
  inline size_t size() const noexcept {
    size_t M = 0;
    for (auto d : degrees)
      M += d;
    return M;
  }

  /**
   * Check if the graph is empty.
   * @returns is the graph empty?
   */
  inline bool empty() const noexcept {
    return degrees.empty();
  }

  /**
   * Check if the graph is directed.
   * @returns is the graph directed?
   */
  inline bool directed() const noexcept {
    return true;
  }
  #pragma endregion


  #pragma region FOREACH
  public:
  /**
   * Iterate over the vertices in the graph.
   * @param fp process function (vertex id, vertex data)
   */
  template <class FP>
  inline void forEachVertex(FP fp) const noexcept {
    for (K u=0; u<span(); ++u)
      fp(u, values[u]);
  }

  /**
   * Iterate over the vertex ids in the graph.
   * @param fp process function (vertex id)
   */
  template <class FP>
  inline void forEachVertexKey(FP fp) const noexcept {
    for (K u=0; u<span(); ++u)
      fp(u);
  }

  /**
   * Iterate over the outgoing edges of a source vertex in the graph.
   * @param u source vertex id
   * @param fp process function (target vertex id, edge weight)
   */
  template <class FP>
  inline void forEachEdge(K u, FP fp) const noexcept {
    size_t i = offsets[u];
    size_t d = degrees[u];
    for (size_t I=i+d; i<I; ++i)
      fp(edgeKeys[i], edgeValues[i]);
  }

  /**
   * Iterate over the target vertex ids of a source vertex in the graph.
   * @param u source vertex id
   * @param fp process function (target vertex id)
   */
  template <class FP>
  inline void forEachEdgeKey(K u, FP fp) const noexcept {
    size_t i = offsets[u];
    size_t d = degrees[u];
    for (size_t I=i+d; i<I; ++i)
      fp(edgeKeys[i]);
  }
  #pragma endregion


  #pragma region OFFSET
  public:
  /**
   * Get the offset of an edge in the graph.
   * @param u source vertex id
   * @param v target vertex id
   * @returns offset of the edge, or -1 if it does not exist
   */
  inline size_t edgeOffset(K u, K v) const noexcept {
    if (!hasVertex(u) || !hasVertex(v)) return size_t(-1);
    size_t  i = offsets[u];
    size_t  d = degrees[u];
    auto   ib = edgeKeys.begin() + i;
    auto   ie = edgeKeys.begin() + i + d;
    auto   it = find(ib, ie, v);
    return it!=ie? it - edgeKeys.begin() : size_t(-1);
  }
  #pragma endregion


  #pragma region ACCESS
  public:
  /**
   * Check if a vertex exists in the graph.
   * @param u vertex id
   * @returns does the vertex exist?
   */
  inline bool hasVertex(K u) const noexcept {
    return u < span();
  }

  /**
   * Check if an edge exists in the graph.
   * @param u source vertex id
   * @param v target vertex id
   * @returns does the edge exist?
   */
  inline bool hasEdge(K u, K v) const noexcept {
    size_t o = edgeOffset(u, v);
    return o != size_t(-1);
  }

  /**
   * Get the number of outgoing edges of a vertex in the graph.
   * @param u vertex id
   * @returns number of outgoing edges of the vertex
   */
  inline size_t degree(K u) const noexcept {
    return u < span()? degrees[u] : 0;
  }

  /**
   * Get the vertex data of a vertex in the graph.
   * @param u vertex id
   * @returns associated data of the vertex
   */
  inline V vertexValue(K u) const noexcept {
    return u < span()? values[u] : V();
  }

  /**
   * Set the vertex data of a vertex in the graph.
   * @param u vertex id
   * @param d associated data of the vertex
   * @returns success?
   */
  inline bool setVertexValue(K u, V d) noexcept {
    if (!hasVertex(u)) return false;
    values[u] = d;
    return true;
  }

  /**
   * Get the edge weight of an edge in the graph.
   * @param u source vertex id
   * @param v target vertex id
   * @returns associated weight of the edge
   */
  inline E edgeValue(K u, K v) const noexcept {
    size_t o = edgeOffset(u, v);
    return o != size_t(-1)? edgeValues[o] : E();
  }

  /**
   * Set the edge weight of an edge in the graph.
   * @param u source vertex id
   * @param v target vertex id
   * @param w associated weight of the edge
   * @returns success?
   */
  inline bool setEdgeValue(K u, K v, E w) noexcept {
    size_t o = edgeOffset(u, v);
    if (o == size_t(-1)) return false;
    edgeValues[o] = w;
    return true;
  }
  #pragma endregion


  #pragma region UPDATE
  public:
  /**
   * Adjust the span of the graph (or the number of vertices).
   * @param n new span
   */
  inline void respan(size_t n) {
    offsets.resize(n+1);
    degrees.resize(n);
    values.resize(n);
  }
  #pragma endregion
  #pragma endregion


  #pragma region CONSTRUCTORS
  public:
  /**
   * Allocate space for CSR representation of a directed graph.
   * @param n number of vertices
   * @param m number of edges
   */
  DiGraphCsr(size_t n, size_t m) {
    offsets.resize(n+1);
    degrees.resize(n);
    values.resize(n);
    edgeKeys.resize(m);
    edgeValues.resize(m);
  }
  #pragma endregion
};
#pragma endregion




#pragma region METHODS
#pragma region WRITE
/**
 * Write the only the sizes of a graph to an output stream.
 * @tparam G graph type
 * @param a output stream
 * @param x graph
 */
template <class G>
inline void writeGraphSizes(ostream& a, const G& x) {
  a << "order: " << x.order() << " size: " << x.size();
  a << (x.directed()? " [directed]" : " [undirected]") << " {}";
}

/**
 * @brief Write the full details of a graph to an output stream.
 * @tparam G graph type
 * @param a output stream
 * @param x graph
 */
template <class G>
inline void writeGraphDetailed(ostream& a, const G& x) {
  a << "order: " << x.order() << " size: " << x.size();
  a << (x.directed()? " [directed]" : " [undirected]") << " {\n";
  x.forEachVertex([&](auto u, auto d) {
    a << u << ":" << d << " ->";
    x.forEachEdge(u, [&](auto v, auto w) {
      a << " " << v << ":" << w;
    });
    a << "\n";
  });
  a << "}";
}

/**
 * Write a graph to an output stream.
 * @tparam G graph type
 * @param a output stream
 * @param x graph
 * @param detailed write detailed information?
 */
template <class G>
inline void writeGraph(ostream& a, const G& x, bool detailed=false) {
  if (detailed) writeGraphDetailed(a, x);
  else writeGraphSizes(a, x);
}

/**
 * Write a graph to an output stream.
 * @tparam K vertex id type
 * @tparam V vertex data type
 * @tparam E edge weight type
 * @param a output stream
 * @param x graph
 * @param detailed write detailed information?
 */
template <class K, class V, class E>
inline void write(ostream& a, const DiGraph<K, V, E>& x, bool detailed=false) {
  writeGraph(a, x, detailed);
}

/**
 * Write only the sizes of a graph to an output stream.
 * @tparam K vertex id type
 * @tparam V vertex data type
 * @tparam E edge weight type
 * @param a output stream
 * @param x graph
 */
template <class K, class V, class E>
inline ostream& operator<<(ostream& a, const DiGraph<K, V, E>& x) {
  write(a, x);
  return a;
}
#pragma endregion
#pragma endregion
