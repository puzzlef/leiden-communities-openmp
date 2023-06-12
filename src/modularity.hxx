#pragma once
#include <cmath>
#include <vector>
#include "_main.hxx"

using std::pow;
using std::vector;




// MODULARITY
// ----------

/**
 * Find the modularity of a community C.
 * @param cin total weight of edges within community C
 * @param ctot total weight of edges of community C
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns modularity [-0.5, 1]
 * @see https://www.youtube.com/watch?v=0zuiLBOIcsw
 */
inline double modularityCommunity(double cin, double ctot, double M, double R=1) {
  ASSERT(cin>=0 && ctot>=0 && M>0 && R>0);
  return cin/(2*M) - R*pow(ctot/(2*M), 2);
}


/**
 * Find the modularity of a set of communities.
 * @param cin total weight of edges within each community
 * @param ctot total weight of edges of each community
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns modularity [-0.5, 1]
 */
template <class T>
inline double modularityCommunities(const vector<T>& cin, const vector<T>& ctot, double M, double R=1) {
  ASSERT(M>0 && R>0);
  double a = 0;
  for (size_t i=0, I=cin.size(); i<I; ++i)
    a += modularityCommunity(cin[i], ctot[i], M, R);
  return a;
}


template <class T>
inline double modularityCommunitiesOmp(const vector<T>& cin, const vector<T>& ctot, double M, double R=1) {
  ASSERT(M>0 && R>0);
  double a = 0;
  size_t C = cin.size();
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (size_t i=0; i<C; ++i)
    a += modularityCommunity(cin[i], ctot[i], M, R);
  return a;
}




/**
 * Find the modularity of a graph, based on community membership function.
 * @param x original graph
 * @param fc community membership function of each vertex (u)
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns modularity [-0.5, 1]
 */
template <class G, class FC>
inline double modularityBy(const G& x, FC fc, double M, double R=1) {
  ASSERT(M>0 && R>0);
  size_t S = x.span();
  vector<double> cin(S), ctot(S);
  x.forEachVertexKey([&](auto u) {
    size_t c = fc(u);
    x.forEachEdge(u, [&](auto v, auto w) {
      size_t d = fc(v);
      if (c==d) cin[c] += w;
      ctot[c] += w;
    });
  });
  return modularityCommunities(cin, ctot, M, R);
}

template <class G, class FC>
inline double modularityByOmp(const G& x, FC fc, double M, double R=1) {
  using K = typename G::key_type;
  ASSERT(M>0 && R>0);
  size_t S = x.span();
  vector<double> cin(S), ctot(S);
  #pragma omp parallel for schedule(auto)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    size_t c = fc(u);
    x.forEachEdge(u, [&](auto v, auto w) {
      size_t d = fc(v);
      if (c==d) {
        #pragma omp atomic
        cin[c] += w;
      }
      #pragma omp atomic
      ctot[c] += w;
    });
  }
  return modularityCommunities(cin, ctot, M, R);
}


/**
 * Find the modularity of a graph, where each vertex is its own community.
 * @param x original graph
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns modularity [-0.5, 1]
 */
template <class G>
inline double modularity(const G& x, double M, double R=1) {
  ASSERT(M>0 && R>0 && R<=1);
  auto fc = [](auto u) { return u; };
  return modularityBy(x, fc, M, R);
}

template <class G>
inline double modularityOmp(const G& x, double M, double R=1) {
  ASSERT(M>0 && R>0 && R<=1);
  auto fc = [](auto u) { return u; };
  return modularityByOmp(x, fc, M, R);
}




// DELTA MODULARITY
// ----------------

/**
 * Find the change in modularity when moving a vertex from community D to C.
 * @param vcout total weight of edges from vertex v to community C
 * @param vdout total weight of edges from vertex v to community D
 * @param vtot total weight of edges from vertex v
 * @param ctot total weight of edges from community C
 * @param dtot total weight of edges from community C
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns delta-modularity [-0.5, 1]
 * @see https://gist.github.com/wolfram77/a3c95cd94a38a100f9b075594a823928
 */
inline double deltaModularity(double vcout, double vdout, double vtot, double ctot, double dtot, double M, double R=1) {
  ASSERT(vcout>=0 && vdout>=0 && vtot>=0 && ctot>=0 && dtot>=0 && M>0 && R>0);
  return (vcout-vdout)/M - R*vtot*(vtot+ctot-dtot)/(2*M*M);
}
