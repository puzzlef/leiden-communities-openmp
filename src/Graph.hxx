#pragma once
#include <cstdint>
#include <utility>
#include <vector>
#include <ostream>
#include <iostream>
#include "_main.hxx"

using std::pair;
using std::vector;
using std::ostream;
using std::cout;




// HELPER MACROS
// -------------
// Helps create graphs.

#ifndef GRAPH_TYPES
#define GRAPH_TYPES(K, V, E) \
  using key_type = K; \
  using vertex_key_type   = K; \
  using vertex_value_type = V; \
  using vertex_pair_type  = pair<K, V>; \
  using edge_key_type     = K; \
  using edge_value_type   = E; \
  using edge_pair_type    = pair<K, E>;
#endif


#ifndef GRAPH_SIZE
#define GRAPH_SIZE(K, V, E, N, M, vexists)  \
  inline size_t span()  const noexcept { return vexists.size(); } \
  inline size_t order() const noexcept { return N; } \
  inline size_t size()  const noexcept { return M; } \
  inline bool   empty() const noexcept { return N==0; }
#endif


#ifndef GRAPH_DIRECTED
#define GRAPH_DIRECTED(K, V, E, de) \
  inline bool directed() const noexcept { return de; }
#endif


#ifndef GRAPH_VERTICES
#define GRAPH_VERTICES(K, V, E, vexists, vvalues) \
  inline auto vertexKeys() const noexcept { \
    auto vkeys = rangeIterable(span()); \
    return conditionalIterable(vkeys, vexists); \
  } \
  inline auto vertexValues() const noexcept { \
    return conditionalIterable(vvalues, vexists); \
  } \
  inline auto vertices() const noexcept { \
    auto vkeys = rangeIterable(span()); \
    auto pairs = pairIterable(vkeys, vvalues); \
    return conditionalIterable(pairs, vexists); \
  }

#define GRAPH_EDGES(K, V, E, eto, enone) \
  inline auto edgeKeys(K u) const noexcept { \
    return u<span()? eto[u].keys()   : enone.keys(); \
  } \
  inline auto edgeValues(K u) const noexcept { \
    return u<span()? eto[u].values() : enone.values(); \
  } \
  inline auto edges(K u) const noexcept { \
    return u<span()? eto[u].entries()  : enone.entries(); \
  }

#define GRAPH_EDGE_AT(K, V, E, eto, enone) \
  inline auto edgeKeyAt(K u, K i) const noexcept { \
    return u<span()? eto[u].keyAt(i)   : enone.keyAt(i); \
  } \
  inline auto edgeValueAt(K u, K i) const noexcept { \
    return u<span()? eto[u].valueAt(i) : enone.valueAt(i); \
  } \
  inline auto edgeAt(K u, K i) const noexcept { \
    return u<span()? eto[u].entryAt(i) : enone.entryAt(i); \
  }

#define GRAPH_INEDGES(K, V, E, efrom, enone) \
  inline auto inEdgeKeys(K v) const noexcept { \
    return v<span()? efrom[v].keys()   : enone.keys(); \
  } \
  inline auto inEdgeValues(K v) const noexcept { \
    return v<span()? efrom[v].values() : enone.values(); \
  } \
  inline auto inEdges(K v) const noexcept { \
    return v<span()? efrom[v].entries()  : enone.entries(); \
  }

#define GRAPH_INEDGES_SCAN(K, V, E, eto) \
  inline auto inEdgeKeys(K v) const noexcept { \
    auto vkeys = rangeIterable(span()); \
    auto fedge = [&](K u) { return eto[u].has(v); }; \
    return filterIterable(vkeys, fedge); \
  } \
  inline auto inEdgeValues(K v) const noexcept { \
    auto fvals = [&](K u) { return eto[u].get(v); }; \
    return transformIterable(inEdgeKeys(v), fvals); \
  } \
  inline auto inEdges(K v) const noexcept { \
    return pairIterable(inEdgeKeys(v), inEdgeValues(v)); \
  }

#define GRAPH_INEDGES_COPY(K, V, E) \
  inline auto inEdgeKeys(K v)   const noexcept { return edgeKeys(v); } \
  inline auto inEdgeValues(K v) const noexcept { return edgeValues(v); } \
  inline auto inEdges(K v)      const noexcept { return edges(v); }
#endif


#ifndef GRAPH_FOREACH_VERTEX
#define GRAPH_FOREACH_VERTEX(K, V, E, vexists, vvalues) \
  template <class F> \
  inline void forEachVertexKey(F fn) const noexcept { \
    for (K u=0; u<span(); ++u) \
      if (vexists[u]) fn(u); \
  } \
  template <class F> \
  inline void forEachVertexValue(F fn) const noexcept { \
    for (K u=0; u<span(); ++u) \
      if (vexists[u]) fn(vvalues[u]); \
  } \
  template <class F> \
  inline void forEachVertex(F fn) const noexcept { \
    for (K u=0; u<span(); ++u) \
      if (vexists[u]) fn(u, vvalues[u]); \
  }

#define GRAPH_FOREACH_XEDGE(K, V, E, name, eto, u, fn) \
  template <class F> \
  inline void forEach##name##Key(K u, F fn) const noexcept { \
    if (u<span()) eto[u].forEachKey(fn); \
  } \
  template <class F> \
  inline void forEach##name##Value(K u, F fn) const noexcept { \
    if (u<span()) eto[u].forEachValue(fn); \
  } \
  template <class F> \
  inline void forEach##name(K u, F fn) const noexcept { \
    if (u<span()) eto[u].forEach(fn); \
  }

#define GRAPH_FOREACH_EDGE(K, V, E, eto) \
  GRAPH_FOREACH_XEDGE(K, V, E, Edge,   eto,   u, fn)
#define GRAPH_FOREACH_INEDGE(K, V, E, efrom) \
  GRAPH_FOREACH_XEDGE(K, V, E, InEdge, efrom, v, fn)

#define GRAPH_FOREACH_INEDGE_SCAN(K, V, E, eto) \
  template <class F> \
  inline void forEachInEdgeKey(K v, F fn) const noexcept { \
    for (K u=0; u<span(); ++u) \
      if (eto[u].has(v)) fn(u); \
  } \
  template <class F> \
  inline void forEachInEdgeValue(K v, F fn) const noexcept { \
    for (K u=0; u<span(); ++u) \
      if (eto[u].has(v)) fn(eto[u].get(v)); \
  } \
  template <class F> \
  inline void forEachInEdge(K v, F fn) const noexcept { \
    for (K u=0; u<span(); ++u) \
      if (eto[u].has(v)) fn(u, eto[u].get(v)); \
  }

#define GRAPH_FOREACH_INEDGE_COPY(K, V, E) \
  template <class F> \
  inline void forEachInEdgeKey(K u, F fn)   const noexcept { forEachEdgeKey(u, fn); } \
  template <class F> \
  inline void forEachInEdgeValue(K u, F fn) const noexcept { forEachEdgeValue(u, fn); } \
  template <class F> \
  inline void forEachInEdge(K u, F fn)      const noexcept { forEachEdge(u, fn); }
#endif


#ifndef GRAPH_HAS
#define GRAPH_HAS(K, V, E, vexists, eto) \
  inline bool hasVertex(K u) const noexcept { \
    return u<span() && vexists[u]; \
  } \
  inline bool hasEdge(K u, K v) const noexcept { \
    return u<span() && eto[u].has(v); \
  }
#endif


#ifndef GRAPH_DEGREES
#define GRAPH_XDEGREE(K, V, E, name, eto) \
  inline K name(K u) const noexcept { \
    return u<span()? K(eto[u].size()) : 0; \
  }

#define GRAPH_INDEGREE_SCAN(K, V, E, eto) \
  inline K inDegree(K v) const noexcept { \
    auto fedge = [&](K u) { return eto[u].has(v); }; \
    return countIf(rangeIterable(span()), fedge); \
  }

#define GRAPH_DEGREES(K, V, E, eto, efrom) \
  GRAPH_XDEGREE(K, V, E, degree, eto) \
  GRAPH_XDEGREE(K, V, E, inDegree, efrom)
#define GRAPH_DEGREES_SCAN(K, V, E, eto) \
  GRAPH_XDEGREE(K, V, E, degree, eto) \
  GRAPH_INDEGREE_SCAN(K, V, E, eto)
#endif


#ifndef GRAPH_VALUES
#define GRAPH_VERTEX_VALUE(K, V, E, vvalues) \
  inline V vertexValue(K u) const noexcept { \
    return u<span()? vvalues[u] : V(); \
  }

#define GRAPH_EDGE_VALUE(K, V, E, eto) \
  inline E edgeValue(K u, K v) const noexcept { \
    return u<span()? eto[u].get(v) : E(); \
  }

#define GRAPH_VALUES(K, V, E, vvalues, eto) \
  GRAPH_VERTEX_VALUE(K, V, E, vvalues) \
  GRAPH_EDGE_VALUE(K, V, E, eto)

#define GRAPH_SET_VERTEX_VALUE(K, V, E, vvalues) \
  inline void setVertexValue(K u, V d) noexcept { \
    if (!hasVertex(u)) return; \
    vvalues[u] = d; \
  }

#define GRAPH_SET_EDGE_VALUE_X(K, V, E, u, v, w, ft, e0, e1) \
  template <class FT> \
  inline void setEdgeValue(K u, K v, E w, FT ft) noexcept { \
    if (!hasVertex(u) || !hasVertex(v)) return; \
    e0; \
    e1; \
  } \
  inline void setEdgeValue(K u, K v, E w) noexcept { \
    auto ft = [](K u) { return true; }; \
    setEdgeValue(u, v, w, ft); \
  }

#define GRAPH_SET_VALUES(K, V, E, vvalues, eto, efrom) \
  GRAPH_SET_VERTEX_VALUE(K, V, E, vvalues) \
  GRAPH_SET_EDGE_VALUE_X(K, V, E, u, v, w, ft, if (ft(u)) eto[u].set(v, w), if (ft(v)) efrom[v].set(u, w))

#define GRAPH_SET_VALUES_SCAN(K, V, E, vvalues, eto) \
  GRAPH_SET_VERTEX_VALUE(K, V, E, vvalues) \
  GRAPH_SET_EDGE_VALUE_X(K, V, E, u, v, w, ft, if (ft(u)) eto[u].set(v, w), false)
#endif


#ifndef GRAPH_RESERVE
#define GRAPH_RESERVE_EDGES_X(K, V, E, u, deg, e0, e1) \
  inline void reserveEdges(K u, size_t deg) { \
    e0; \
    e1; \
  }

#define GRAPH_RESERVE_EDGES(K, V, E, eto, efrom) \
  GRAPH_RESERVE_EDGES_X(K, V, E, u, deg, eto[u].reserve(deg), efrom[u].reserve(deg))
#define GRAPH_RESERVE_EDGES_SCAN(K, V, E, eto) \
  GRAPH_RESERVE_EDGES_X(K, V, E, u, deg, eto[u].reserve(deg), false)

#define GRAPH_RESERVE_X(K, V, E, vexists, vvalues, n, deg, e0, e1) \
  inline void reserve(size_t n, size_t deg=0) { \
    vexists.resize(max(n, span())); \
    vvalues.resize(max(n, span())); \
    e0; \
    e1; \
    if (deg==0) return; \
    for (K u=0; u<span(); ++u) \
      reserveEdges(u, deg); \
  }

#define GRAPH_RESERVE(K, V, E, vexists, vvalues, eto, efrom) \
  GRAPH_RESERVE_X(K, V, E, vexists, vvalues, n, deg, eto.resize(max(n, span())), efrom.resize(max(n, span())))
#define GRAPH_RESERVE_SCAN(K, V, E, vexists, vvalues, eto) \
  GRAPH_RESERVE_X(K, V, E, vexists, vvalues, n, deg, eto.resize(max(n, span())), false)
#endif


#ifndef GRAPH_UPDATE
#define GRAPH_UPDATE_EDGES_X(K, V, E, u, buf, e0, e1) \
  inline void updateEdges(K u, vector<pair<K, E>> *buf=nullptr) { \
    e0; \
    e1; \
  }

#define GRAPH_UPDATE_EDGES(K, V, E, eto, efrom) \
  GRAPH_UPDATE_EDGES_X(K, V, E, u, buf, eto[u].update(buf), efrom[u].update(buf))
#define GRAPH_UPDATE_EDGES_SCAN(K, V, E, eto) \
  GRAPH_UPDATE_EDGES_X(K, V, E, u, buf, eto[u].update(buf), false)

#define GRAPH_UPDATE(K, V, E, N, M) \
  inline void update() { \
    N = 0; M = 0; \
    vector<pair<K, E>> buf; \
    forEachVertexKey([&](K u) { \
      updateEdges(u, &buf); \
      M += degree(u); ++N; \
    }); \
  }
#endif


#ifndef GRAPH_RESPAN
#define GRAPH_RESPAN_X(K, V, E, vexists, vvalues, n, e0, e1) \
  inline void respan(size_t n) { \
    vexists.resize(n); \
    vvalues.resize(n); \
    e0; \
    e1; \
  }

#define GRAPH_RESPAN(K, V, E, vexists, vvalues, eto, efrom) \
  GRAPH_RESPAN_X(K, V, E, vexists, vvalues, n, eto.resize(n), efrom.resize(n))
#define GRAPH_RESPAN_SCAN(K, V, E, vexists, vvalues, eto) \
  GRAPH_RESPAN_X(K, V, E, vexists, vvalues, n, eto.resize(n), false)
#endif


#ifndef GRAPH_CLEAR
#define GRAPH_CLEAR_X(K, V, E, N, M, vexists, vvalues, e0, e1) \
  inline void clear() noexcept { \
    N = 0; M = 0; \
    vexists.clear(); \
    vvalues.clear(); \
    e0; \
    e1; \
  }

#define GRAPH_CLEAR(K, V, E, N, M, vexists, vvalues, eto, efrom) \
  GRAPH_CLEAR_X(K, V, E, N, M, vexists, vvalues, eto.clear(), efrom.clear())
#define GRAPH_CLEAR_SCAN(K, V, E, N, M, vexists, vvalues, eto) \
  GRAPH_CLEAR_X(K, V, E, N, M, vexists, vvalues, eto.clear(), false)
#endif


#ifndef GRAPH_ADD_VERTEX
#define GRAPH_ADD_VERTEX(K, V, E, vexists, vvalues) \
  inline void addVertex(K u, V d=V()) { \
    if (hasVertex(u)) return; \
    if (u>=span()) respan(u+1); \
    vexists[u] = true; \
    vvalues[u] = d; \
  }
#endif


#ifndef GRAPH_ADD_EDGE
#define GRAPH_ADD_EDGE_X(K, V, E, u, v, w, ft, e0, e1) \
  template <class FT> \
  inline void addEdge(K u, K v, E w, FT ft) { \
    addVertex(u); addVertex(v); \
    e0; \
    e1; \
  } \
  inline void addEdge(K u, K v, E w=E()) { \
    auto ft = [](K u) { return true; }; \
    addEdge(u, v, w, ft); \
  }

#define GRAPH_ADD_EDGE(K, V, E, eto, efrom) \
  GRAPH_ADD_EDGE_X(K, V, E, u, v, w, ft, if (ft(u)) eto[u].add(v, w), if (ft(v)) efrom[v].add(u, w))
#define GRAPH_ADD_EDGE_SCAN(K, V, E, eto) \
  GRAPH_ADD_EDGE_X(K, V, E, u, v, w, ft, if (ft(u)) eto[u].add(v, w), false)
#endif


#ifndef GRAPH_REMOVE_EDGE
#define GRAPH_REMOVE_EDGE_X(K, V, E, u, v, ft, e0, e1) \
  template <class FT> \
  inline void removeEdge(K u, K v, FT ft) { \
    if (!hasVertex(u) || !hasVertex(v)) return; \
    e0; \
    e1; \
  } \
  inline void removeEdge(K u, K v) { \
    auto ft = [](K u) { return true; }; \
    removeEdge(u, v, ft); \
  }

#define GRAPH_REMOVE_EDGE(K, V, E, eto, efrom) \
  GRAPH_REMOVE_EDGE_X(K, V, E, u, v, ft, if (ft(u)) eto[u].remove(v), if (ft(v)) efrom[v].remove(u))
#define GRAPH_REMOVE_EDGE_SCAN(K, V, E, eto) \
  GRAPH_REMOVE_EDGE_X(K, V, E, u, v, ft, if (ft(u)) eto[u].remove(v), false)
#endif


#ifndef GRAPH_REMOVE_EDGES
#define GRAPH_REMOVE_EDGES(K, V, E, eto, efrom) \
  template <class FT> \
  inline void removeEdges(K u, FT ft) { \
    if (!hasVertex(u)) return; \
    eto[u].forEachKey([&](K v) { if (ft(v)) efrom[v].remove(u); }); \
    if (ft(u)) eto[u].clear(); \
  } \
  inline void removeEdges(K u) { \
    auto ft = [](K u) { return true; }; \
    removeEdges(u, ft); \
  }

#define GRAPH_REMOVE_EDGES_SCAN(K, V, E, eto) \
  template <class FT> \
  inline void removeEdges(K u, FT ft) { \
    if (!hasVertex(u)) return; \
    if (ft(u)) eto[u].clear(); \
  } \
  inline void removeEdges(K u) { \
    auto ft = [](K u) { return true; }; \
    removeEdges(u, ft); \
  }
#endif


#ifndef GRAPH_REMOVE_INEDGES
#define GRAPH_REMOVE_INEDGES(K, V, E, eto, efrom) \
  template <class FT> \
  inline void removeInEdges(K v, FT ft) { \
    if (!hasVertex(v)) return; \
    efrom[v].forEachKey([&](K u) { if (ft(u)) eto[u].remove(v); }); \
    if (ft(v)) efrom[v].clear(); \
  } \
  inline void removeInEdges(K v) { \
    auto ft = [](K u) { return true; }; \
    removeInEdges(v, ft); \
  }

#define GRAPH_REMOVE_INEDGES_SCAN(K, V, E, eto) \
  template <class FT> \
  inline void removeInEdges(K v, FT ft) { \
    if (!hasVertex(v)) return; \
    forEachVertexKey([&](K u) { if (ft(u)) eto[u].remove(v); }); \
  } \
  inline void removeInEdges(K v) { \
    auto ft = [](K u) { return true; }; \
    removeInEdges(v, ft); \
  }
#endif


#ifndef GRAPH_REMOVE_VERTEX
#define GRAPH_REMOVE_VERTEX(K, V, E, vexists, vvalues) \
  template <class FT> \
  inline void removeVertex(K u, FT ft) { \
    if (!hasVertex(u)) return; \
    removeEdges(u, ft); \
    removeInEdges(u, ft); \
    vexists[u] = false; \
    vvalues[u] = V(); \
  } \
  inline void removeVertex(K u) { \
    auto ft = [](K u) { return true; }; \
    removeVertex(u, ft); \
  }
#endif


#ifndef GRAPH_WRITE
#define GRAPH_WRITE(K, V, E, Bitset, Graph) \
  template <class K, class V, class E, tclass2 Bitset> \
  inline void write(ostream& a, const Graph<K, V, E, Bitset>& x, bool detailed=false) { writeGraph(a, x, detailed); } \
  template <class K, class V, class E, tclass2 Bitset> \
  inline ostream& operator<<(ostream& a, const Graph<K, V, E, Bitset>& x) { write(a, x); return a; }
#endif




// DI-GRAPH
// --------
// Directed graph that memorizes in- and out-edges for each vertex.

template <class K=uint32_t, class V=NONE, class E=NONE, tclass2 Bitset=LazyBitset>
class DiGraph {
  // Data.
  protected:
  size_t N = 0, M = 0;
  vector<bool> vexists;
  vector<V>    vvalues;
  vector<Bitset<K, E>> eto;
  vector<Bitset<K, E>> efrom;
  Bitset<K, E> enone;

  // Types.
  public:
  GRAPH_TYPES(K, V, E)

  // Property operations.
  public:
  GRAPH_SIZE(K, V, E, N, M, vexists)
  GRAPH_DIRECTED(K, V, E, true)

  // Scan operations.
  public:
  GRAPH_VERTICES(K, V, E, vexists, vvalues)
  GRAPH_EDGES(K, V, E, eto, enone)
  GRAPH_EDGE_AT(K, V, E, eto, enone)
  GRAPH_INEDGES(K, V, E, efrom, enone)
  GRAPH_FOREACH_VERTEX(K, V, E, vexists, vvalues)
  GRAPH_FOREACH_EDGE(K, V, E, eto)
  GRAPH_FOREACH_INEDGE(K, V, E, efrom)

  // Access operations.
  public:
  GRAPH_HAS(K, V, E, vexists, eto)
  GRAPH_DEGREES(K, V, E, eto, efrom)
  GRAPH_VALUES(K, V, E, vvalues, eto)
  GRAPH_SET_VALUES(K, V, E, vvalues, eto, efrom)

  // Update operations.
  public:
  GRAPH_RESERVE_EDGES(K, V, E, eto, efrom)
  GRAPH_RESERVE(K, V, E, vexists, vvalues, eto, efrom)
  GRAPH_UPDATE_EDGES(K, V, E, eto, efrom)
  GRAPH_UPDATE(K, V, E, N, M)
  GRAPH_RESPAN(K, V, E, vexists, vvalues, eto, efrom)
  GRAPH_CLEAR(K, V, E, N, M, vexists, vvalues, eto, efrom)
  GRAPH_ADD_VERTEX(K, V, E, vexists, vvalues)
  GRAPH_ADD_EDGE(K, V, E, eto, efrom)
  GRAPH_REMOVE_EDGE(K, V, E, eto, efrom)
  GRAPH_REMOVE_EDGES(K, V, E, eto, efrom)
  GRAPH_REMOVE_INEDGES(K, V, E, eto, efrom)
  GRAPH_REMOVE_VERTEX(K, V, E, vexists, vvalues)
};

template <class K=uint32_t, class V=NONE, class E=NONE>
using UnorderedDiGraph = DiGraph<K, V, E, LazyBitset>;




// OUT DI-GRAPH
// ------------
// Directed graph that memorizes only out-edges for each vertex.

template <class K=uint32_t, class V=NONE, class E=NONE, tclass2 Bitset=LazyBitset>
class OutDiGraph {
  // Data.
  protected:
  size_t N = 0, M = 0;
  vector<bool> vexists;
  vector<V>    vvalues;
  vector<Bitset<K, E>> eto;
  Bitset<K, E> enone;

  // Types.
  public:
  GRAPH_TYPES(K, V, E)

  // Property operations.
  public:
  GRAPH_SIZE(K, V, E, N, M, vexists)
  GRAPH_DIRECTED(K, V, E, true)

  // Scan operations.
  public:
  GRAPH_VERTICES(K, V, E, vexists, vvalues)
  GRAPH_EDGES(K, V, E, eto, enone)
  GRAPH_EDGE_AT(K, V, E, eto, enone)
  GRAPH_INEDGES_SCAN(K, V, E, eto)
  GRAPH_FOREACH_VERTEX(K, V, E, vexists, vvalues)
  GRAPH_FOREACH_EDGE(K, V, E, eto)
  GRAPH_FOREACH_INEDGE_SCAN(K, V, E, eto)

  // Access operations.
  public:
  GRAPH_HAS(K, V, E, vexists, eto)
  GRAPH_DEGREES_SCAN(K, V, E, eto)
  GRAPH_VALUES(K, V, E, vvalues, eto)
  GRAPH_SET_VALUES_SCAN(K, V, E, vvalues, eto)

  // Update operations.
  public:
  GRAPH_RESERVE_EDGES_SCAN(K, V, E, eto)
  GRAPH_RESERVE_SCAN(K, V, E, vexists, vvalues, eto)
  GRAPH_UPDATE_EDGES_SCAN(K, V, E, eto)
  GRAPH_UPDATE(K, V, E, N, M)
  GRAPH_RESPAN_SCAN(K, V, E, vexists, vvalues, eto)
  GRAPH_CLEAR_SCAN(K, V, E, N, M, vexists, vvalues, eto)
  GRAPH_ADD_VERTEX(K, V, E, vexists, vvalues)
  GRAPH_ADD_EDGE_SCAN(K, V, E, eto)
  GRAPH_REMOVE_EDGE_SCAN(K, V, E, eto)
  GRAPH_REMOVE_EDGES_SCAN(K, V, E, eto)
  GRAPH_REMOVE_INEDGES_SCAN(K, V, E, eto)
  GRAPH_REMOVE_VERTEX(K, V, E, vexists, vvalues)
};

template <class K=uint32_t, class V=NONE, class E=NONE>
using LazyOutDiGraph = OutDiGraph<K, V, E, LazyBitset>;




// GRAPH
// -----
// Undirected graph.

template <class K=uint32_t, class V=NONE, class E=NONE, tclass2 Bitset=LazyBitset>
class Graph : public OutDiGraph<K, V, E, Bitset> {
  using G = OutDiGraph<K, V, E, Bitset>;

  // Property operations.
  public:
  inline size_t size() const noexcept { return G::size()/2; }
  GRAPH_DIRECTED(K, V, E, false)

  // Scan operations.
  public:
  GRAPH_INEDGES_COPY(K, V, E)
  GRAPH_FOREACH_INEDGE_COPY(K, V, E)

  // Access operations.
  public:
  inline K inDegree(K v) const noexcept {
    return degree(v);
  }

  template <class FT>
  inline void setEdgeValue(K u, K v, E w, FT ft) {
    G::setEdgeValue(u, v, w, ft);
    G::setEdgeValue(v, u, w, ft);
  }
  inline void setEdgeValue(K u, K v, E w) {
    auto ft = [](K u) { return true; };
    setEdgeValue(u, v, w, ft);
  }

  // Update operations.
  public:
  template <class FT>
  inline void addEdge(K u, K v, E w, FT ft) {
    G::addEdge(u, v, w, ft);
    G::addEdge(v, u, w, ft);
  }
  inline void addEdge(K u, K v, E w=E()) {
    auto ft = [](K u) { return true; };
    addEdge(u, v, w, ft);
  }

  template <class FT>
  inline void removeEdge(K u, K v, FT ft) {
    G::removeEdge(u, v, ft);
    G::removeEdge(v, u, ft);
  }
  inline void removeEdge(K u, K v) {
    auto ft = [](K u) { return true; };
    removeEdge(u, v, ft);
  }

  template <class FT>
  inline void removeEdges(K u, FT ft) {
    forEachEdgeKey(u, [&](K v) { G::removeEdge(v, u, ft); });
    G::removeEdges(u, ft);
  }
  inline void removeEdges(K u) {
    auto ft = [](K u) { return true; };
    removeEdges(u, ft);
  }

  template <class FT>
  inline void removeInEdges(K v, FT ft) {
    removeEdges(v, ft);
  }
  inline void removeInEdges(K v) {
    auto ft = [](K u) { return true; };
    removeInEdges(v, ft);
  }
};

template <class K=uint32_t, class V=NONE, class E=NONE>
using LazyGraph = Graph<K, V, E, LazyBitset>;




// RETYPE
// ------

template <class K, class V, class E, tclass2 B, class KA=K, class VA=V, class EA=E>
constexpr auto retype(const DiGraph<K, V, E, B>& x, KA _k=KA(), VA _v=VA(), EA _e=E()) {
  return DiGraph<KA, VA, EA, B>();
}
template <class K, class V, class E, tclass2 B, class KA=K, class VA=V, class EA=E>
constexpr auto retype(const OutDiGraph<K, V, E, B>& x, KA _k=KA(), VA _v=VA(), EA _e=E()) {
  return OutDiGraph<KA, VA, EA, B>();
}
template <class K, class V, class E, tclass2 B, class KA=K, class VA=V, class EA=E>
constexpr auto retype(const Graph<K, V, E, B>& x, KA _k=KA(), VA _v=VA(), EA _e=E()) {
  return Graph<KA, VA, EA, B>();
}




// WTITE
// -----

template <class G>
void writeGraphSizes(ostream& a, const G& x) {
  a << "order: " << x.order() << " size: " << x.size();
  a << (x.directed()? " [directed]" : " [undirected]") << " {}";
}

template <class G>
void writeGraphDetailed(ostream& a, const G& x) {
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

template <class G>
inline void writeGraph(ostream& a, const G& x, bool detailed=false) {
  if (detailed) writeGraphDetailed(a, x);
  else writeGraphSizes(a, x);
}

GRAPH_WRITE(K, V, E, Bitset, DiGraph)
GRAPH_WRITE(K, V, E, Bitset, OutDiGraph)
GRAPH_WRITE(K, V, E, Bitset, Graph)
