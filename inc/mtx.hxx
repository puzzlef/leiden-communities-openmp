#pragma once
#include <utility>
#include <string>
#include <istream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include "_main.hxx"
#include "Graph.hxx"
#include "update.hxx"
#ifdef OPENMP
#include <omp.h>
#endif

using std::tuple;
using std::string;
using std::istream;
using std::istringstream;
using std::ifstream;
using std::ofstream;
using std::move;
using std::max;
using std::getline;




#pragma region METHODS
#pragma region READ MTX HEADER
/**
 * Read header of MTX file.
 * @param s input stream
 * @param symmetric is graph symmetric (updated)
 * @param rows number of rows (updated)
 * @param cols number of columns (updated)
 * @param size number of lines/edges (updated)
 */
inline void readMtxHeader(istream& s, bool& symmetric, size_t& rows, size_t& cols, size_t& size) {
  string line, h0, h1, h2, h3, h4;
  // Skip past the comments and read the graph type.
  while (true) {
    getline(s, line);
    if (line.find('%')!=0) break;
    if (line.find("%%")!=0) continue;
    istringstream sline(line);
    sline >> h0 >> h1 >> h2 >> h3 >> h4;
  }
  if (h1!="matrix" || h2!="coordinate") { symmetric = false; rows = 0; cols = 0; size = 0; return; }
  symmetric = h4=="symmetric" || h4=="skew-symmetric";
  // Read rows, cols, size.
  istringstream sline(line);
  sline >> rows >> cols >> size;
}
inline void readMtxHeader(const char* pth, bool& symmetric, size_t& rows, size_t& cols, size_t& size) {
  ifstream s(pth);
  return readMtxHeader(s, symmetric, rows, cols, size);
}


/**
 * Read order of graph in MTX file.
 * @param s input stream
 * @returns number of vertices (1-N)
 */
inline size_t readMtxOrder(istream& s) {
  bool symmetric; size_t rows, cols, size;
  readMtxHeader(s, symmetric, rows, cols, size);
  return max(rows, cols);
}
inline size_t readMtxOrder(const char* pth) {
  ifstream s(pth);
  return readMtxOrder(s);
}


/**
 * Read size of graph in MTX file.
 * @param s input stream
 * @returns number of edges
 */
inline size_t readMtxSize(istream& s) {
  bool symmetric; size_t rows, cols, size;
  readMtxHeader(s, symmetric, rows, cols, size);
  return size;
}
inline size_t readMtxSize(const char* pth) {
  ifstream s(pth);
  return readMtxSize(s);
}


/**
 * Read span of graph in MTX file.
 * @param s input stream
 * @returns number of vertices + 1
 */
inline size_t readMtxSpan(istream& s) {
  return readMtxOrder(s) + 1;
}
inline size_t readMtxSpan(const char* pth) {
  ifstream s(pth);
  return readMtxSpan(s);
}
#pragma endregion




#pragma region READ MTX DO
/**
 * Read contents of MTX file.
 * @param s input stream
 * @param weighted is it weighted?
 * @param fh on header (symmetric, rows, cols, size)
 * @param fb on body line (u, v, w)
 */
template <class FH, class FB>
inline void readMtxDo(istream& s, bool weighted, FH fh, FB fb) {
  bool symmetric; size_t rows, cols, size;
  readMtxHeader(s, symmetric, rows, cols, size);
  fh(symmetric, rows, cols, size);
  size_t n = max(rows, cols);
  if (n==0) return;
  // Process body lines sequentially.
  string line;
  while (getline(s, line)) {
    size_t u, v; double w = 1;
    istringstream sline(line);
    if (!(sline >> u >> v)) break;
    if (weighted) sline >> w;
    fb(u, v, w);
    if (symmetric) fb(v, u, w);
  }
}
template <class FH, class FB>
inline void readMtxDo(const char *pth, bool weighted, FH fh, FB fb) {
  ifstream s(pth);
  readMtxDo(s, weighted, fh, fb);
}


#ifdef OPENMP
/**
 * Read contents of MTX file.
 * @param s input stream
 * @param weighted is it weighted?
 * @param fh on header (symmetric, rows, cols, size)
 * @param fb on body line (u, v, w)
 */
template <class FH, class FB>
inline void readMtxDoOmp(istream& s, bool weighted, FH fh, FB fb) {
  bool symmetric; size_t rows, cols, size;
  readMtxHeader(s, symmetric, rows, cols, size);
  fh(symmetric, rows, cols, size);
  size_t n = max(rows, cols);
  if (n==0) return;
  // Process body lines in parallel.
  const int LINES   = 131072;
  vector<string> lines(LINES);
  vector<tuple<size_t, size_t, double>> edges(LINES);
  while (true) {
    // Read several lines from the stream.
    int READ = 0;
    for (int i=0; i<LINES; ++i, ++READ)
      if (!getline(s, lines[i])) break;
    if (READ==0) break;
    // Parse lines using multiple threads.
    #pragma omp parallel for schedule(dynamic, 1024)
    for (int i=0; i<READ; ++i) {
      char *line = (char*) lines[i].c_str();
      size_t u = strtoull(line, &line, 10);
      size_t v = strtoull(line, &line, 10);
      double w = weighted? strtod(line, &line) : 0;
      edges[i] = {u, v, w? w : 1};
    }
    // Notify parsed lines.
    #pragma omp parallel
    {
      for (int i=0; i<READ; ++i) {
        const auto& [u, v, w] = edges[i];
        fb(u, v, w);
        if (symmetric) fb(v, u, w);
      }
    }
  }
}
template <class FH, class FB>
inline void readMtxDoOmp(const char *pth, bool weighted, FH fh, FB fb) {
  ifstream s(pth);
  readMtxDoOmp(s, weighted, fh, fb);
}
#endif
#pragma endregion




#pragma region READ MTX IF
/**
 * Read MTX file as graph if test passes.
 * @param a output graph (updated)
 * @param s input stream
 * @param weighted is it weighted?
 * @param fv include vertex? (u, d)
 * @param fe include edge? (u, v, w)
 */
template <class G, class FV, class FE>
inline void readMtxIfW(G &a, istream& s, bool weighted, FV fv, FE fe) {
  using K = typename G::key_type;
  using V = typename G::vertex_value_type;
  using E = typename G::edge_value_type;
  auto fh = [&](auto symmetric, auto rows, auto cols, auto size) { addVerticesIfU(a, K(1), K(max(rows, cols)+1), V(), fv); };
  auto fb = [&](auto u, auto v, auto w) { if (fe(K(u), K(v), K(w))) a.addEdge(K(u), K(v), E(w)); };
  readMtxDo(s, weighted, fh, fb);
  a.update();
}
template <class G, class FV, class FE>
inline void readMtxIfW(G &a, const char *pth, bool weighted, FV fv, FE fe) {
  ifstream s(pth);
  readMtxIfW(a, s, weighted, fv, fe);
}


#ifdef OPENMP
/**
 * Read MTX file as graph if test passes.
 * @param a output graph (updated)
 * @param s input stream
 * @param weighted is it weighted?
 * @param fv include vertex? (u, d)
 * @param fe include edge? (u, v, w)
 */
template <class G, class FV, class FE>
inline void readMtxIfOmpW(G &a, istream& s, bool weighted, FV fv, FE fe) {
  using K = typename G::key_type;
  using V = typename G::vertex_value_type;
  using E = typename G::edge_value_type;
  auto fh = [&](auto symmetric, auto rows, auto cols, auto size) { addVerticesIfU(a, K(1), K(max(rows, cols)+1), V(), fv); };
  auto fb = [&](auto u, auto v, auto w) { if (fe(K(u), K(v), K(w))) addEdgeOmpU(a, K(u), K(v), E(w)); };
  readMtxDoOmp(s, weighted, fh, fb);
  updateOmpU(a);
}
template <class G, class FV, class FE>
inline void readMtxIfOmpW(G &a, const char *pth, bool weighted, FV fv, FE fe) {
  ifstream s(pth);
  readMtxIfOmpW(a, s, weighted, fv, fe);
}
#endif
#pragma endregion




#pragma region READ MTX
/**
 * Read MTX file as graph.
 * @param a output graph (updated)
 * @param s input stream
 * @param weighted is it weighted?
 */
template <class G>
inline void readMtxW(G& a, istream& s, bool weighted=false) {
  auto fv = [](auto u, auto d)         { return true; };
  auto fe = [](auto u, auto v, auto w) { return true; };
  readMtxIfW(a, s, weighted, fv, fe);
}
template <class G>
inline void readMtxW(G& a, const char *pth, bool weighted=false) {
  ifstream s(pth);
  readMtxW(a, s, weighted);
}


#ifdef OPENMP
/**
 * Read MTX file as graph.
 * @param a output graph (updated)
 * @param s input stream
 * @param weighted is it weighted?
 */
template <class G>
inline void readMtxOmpW(G& a, istream& s, bool weighted=false) {
  auto fv = [](auto u, auto d)         { return true; };
  auto fe = [](auto u, auto v, auto w) { return true; };
  readMtxIfOmpW(a, s, weighted, fv, fe);
}
template <class G>
inline void readMtxOmpW(G& a, const char *pth, bool weighted=false) {
  ifstream s(pth);
  readMtxOmpW(a, s, weighted);
}
#endif
#pragma endregion
#pragma endregion
