Design of OpenMP-based Dynamic [Louvain algorithm] for [community detection].

[Louvain algorithm] is an **agglomerative-hierarchical** community detection
method that **greedily optimizes** for [modularity]. Given an *undirected*
*weighted graph*, all vertices are first considered to be *their own*
*communities*. In the **first phase**, each vertex greedily decides to move to
the community of one of its neighbors which gives greatest increase in
modularity. If moving to no neighbor's community leads to an increase in
modularity, the vertex chooses to stay with its own community. This is done
sequentially for all the vertices. If the total change in modularity is more
than a certain threshold, this phase is repeated. Once this **local-moving**
**phase** is complete, all vertices have formed their first hierarchy of
communities. The **next phase** is called the **aggregation phase**, where all
the *vertices belonging to a community* are *collapsed* into a single
**super-vertex**, such that edges between communities are represented as edges
between respective super-vertices (edge weights are combined), and edges within
each community are represented as self-loops in respective super-vertices
(again, edge weights are combined). Together, the local-moving and the
aggregation phases constitute a **pass**. This super-vertex graph is then used
as input for the next pass. This process continues until the increase in
modularity is below a certain threshold. As a result from each pass, we have a
*hierarchy of community memberships* for each vertex as a **dendrogram**. We
generally consider the *top-level hierarchy* as the *final result* of community
detection process.

*Louvain* algorithm is a hierarchical algorithm, and thus has two different
tolerance parameters: `tolerance` and `passTolerance`. **tolerance** defines the
minimum amount of increase in modularity expected, until the local-moving phase
of the algorithm is considered to have converged. We compare the increase in
modularity in each iteration of the local-moving phase to see if it is below
`tolerance`. **passTolerance** defines the minimum amount of increase in
modularity expected, until the entire algorithm is considered to have converged.
We compare the increase in modularity across all iterations of the local-moving
phase in the current pass to see if it is below `passTolerance`. `passTolerance`
is normally set to `0` (we want to maximize our modularity gain), but the same
thing does not apply for `tolerance`. Adjusting values of `tolerance` between
each pass have been observed to impact the runtime of the algorithm, without
significantly affecting the modularity of obtained communities. In this
experiment, we compare the performance of *three different types* of OpenMP-based
**dynamic Louvain** with respect to the *static* version.

**Naive-dynamic**:
- We start with previous community membership of each vertex (instead of each vertex its own community).

**Dynamic Delta-screening**:
- All edge batches are undirected, and sorted by source vertex-id.
- For edge deletions within the same community `i` and `j`,
  `i`'s neighbors and `j`'s community is marked as affected.
- For edge insertions across communities with source vertex `i` and highest modularity changing edge vertex `j*`,
  `i`'s neighbors and `j*`'s community is marked as affected.

**Dynamic Frontier**:
- All edge batches are undirected.
- For edge deletions within the same community `i` and `j`,
  `i` is marked as affected.
- For edge insertions across communities with source vertex `i` and destination vertex `j`,
  `i` is marked as affected.
- Vertices whose communities change in local-moving phase have their neighbors marked as affected.

The input data used for below experiments is available from the [SuiteSparse Matrix Collection].
The experiments were done with guidance from [Prof. Kishore Kothapalli] and
[Prof. Dip Sankar Banerjee].

[Louvain algorithm]: https://en.wikipedia.org/wiki/Louvain_method
[community detection]: https://en.wikipedia.org/wiki/Community_search
[modularity]: https://en.wikipedia.org/wiki/Modularity_(networks)
[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu

<br>


### Comparision on large graphs

In this experiment ([input-large]), we first compute the community membership of
each vertex using the static Louvain algorithm. We then generate random batch
updates consisting of an equal mix of *deletions (-)* and  *insertions (+)* of
edges of size `10^-7 |E|` to `0.1 |E|` in multiples of `10` (where `|E|` is the
number of edges in the original graph after making it undirected). For each
batch size, we generate *five* different batches for the purpose of *averaging*.
Each batch of edges (insertion / deletion) is generated randomly such that the
selection of each vertex (as endpoint) is *equally probable*. We choose the
Louvain *parameters* as `resolution = 1.0`, `tolerance = 1e-2` (for local-moving
phase) with *tolerance* decreasing after every pass by a factor of
`toleranceDeclineFactor = 10`, and a `passTolerance = 0.0` (when passes stop).
In addition we limit the maximum number of iterations in a single local-moving
phase with `maxIterations = 20`, and limit the maximum number of passes with
`maxPasses = 20`. We run the Louvain algorithm until convergence (or until the
maximum limits are exceeded), and measure the **time taken** for the
*computation* and *pre-processing* (for dynamic approaches), the **modularity**
**score**, the **total number of iterations** (in the *local-moving phase*), and
the number of **passes**. This is repeated for each input graph.

From the results, we make make the following observations. **Dynamic Frontier**
based **Louvain** converges the fastest, which obtaining communities with
equivalent modularity. We also observe that **Dynamic Delta-screening** based
**Louvain** has the same performance as that of the Naive-dynamic approach, but
has poorer performance for larger batch sizes (due to its high pre-processing
cost/overhead for large batch sizes). Therefore, **Dynamic Frontier based**
**Louvain** would be the **best choice**. We also not that **Louvain** algorithm
does not scaled too well with an increase in the number of threads. This is
likely due to higher pressure on cache coherence system as well as the algorithm
becoming closer to a synchronous approach, which is inherently slower than an
asynchronous approach. Trying to avoid community swaps with parallel approach
does not seem to improve performance by any significant amount. However, it is
possible that if synchronous approach is used with OpenMP, then its performance
may be a bit better. All outputs are saved in a [gist] and a small part of the
output is listed here. Some [charts] are also included below, generated from
[sheets].

[![](https://i.imgur.com/FDFFa4F.png)][sheetp]
[![](https://i.imgur.com/cb8M5dO.png)][sheetp]
[![](https://i.imgur.com/uls0R3W.png)][sheetp]
[![](https://i.imgur.com/SlWkSBc.png)][sheetp]
[![](https://i.imgur.com/goTmd1W.png)][sheetp]
[![](https://i.imgur.com/eW0rCWo.png)][sheetp]

[input-large]: https://github.com/puzzlef/louvain-communities-openmp-dynamic/tree/input-large
[gist]: https://gist.github.com/wolfram77/2d64f933f6524ba15ee7593f7e3b10f5
[charts]: https://imgur.com/a/Gbc8WgO
[sheets]: https://docs.google.com/spreadsheets/d/1F6Z-lWNDYynm6m2PTsIN_nxMu8Y9CrkIQagCU0Nr2LU/edit?usp=sharing
[sheetp]: https://docs.google.com/spreadsheets/d/e/2PACX-1vS5LH03ALzgcv6QNV9I9Wl1_000Vl9BNZKnMdF04d4qeG5dqQ60fFHL4xynG_8LnVFbsyaJAucWuen6/pubhtml

<br>


### Measure communities

In this experiment ([measure-communities]), we **measure** the **properties of**
**communities obtained** with *Static Louvain* algorithm. These include the
*number of communities*, the *size distribution of communities* (*gini*
*coefficient*), and the *overall modularity score*.

[measure-communities]: https://github.com/puzzlef/louvain-communities-openmp-dynamic/tree/measure-communities

<br>


### Measure affected vertices

In this experiment ([measure-affected]), we **measure** the number of **affected**
**vertices** with *Dynamic* *Delta-screening* and *Dynamic Frontier* based
*Louvain* for random batch updates consisting of edge insertions, with the size
of batch update varying from `10^-6 |E|` to `0.1 |E|`.

Results show that *Dynamic Delta-screening* marks `15000x`, `2000x`, `440x`,
`44x`, `6.4x`, and `1.7x` the number of affected vertices as *Dynamic Frontier*
based approach on batch updates of size `10^-6 |E|` to `0.1 |E|`.

[measure-affected]: https://github.com/puzzlef/louvain-communities-openmp-dynamic/tree/measure-affected

<br>


### Multi-batch updates

In this experiment ([multi-batch]), we generate `5000` random **multi-batch updates** consisting
of *edge insertions* of size `10^-3 |E|` one after the other on graphs
`web-Stanford` and `web-BerkStan` and observe the performance and modularity of
communities obtained with *Static*, *Naive-dynamic*, *Dynamic Delta-screening*,
and *Dynamic Frontier* based *Louvain*. We do this to measure after how many
batch updates do we need to re-run the static algorithm.

Our results indicate that we need to rerun the static algorithm after `~1300`
batch updates with *Dynamic Delta-screening* based *Louvain*, and after `~2800`
batch updates with *Dynamic Frontier* based *Louvain*.

[multi-batch]: https://github.com/puzzlef/louvain-communities-openmp-dynamic/tree/multi-batch

<br>
<br>


## Build instructions

To run the [input-large] experiment, download this repository and run the
following. Note that input graphs must be placed in `~/Data` directory, and
output logs will be written to `~/Logs` directory.

```bash
# Perform comparision on large graphs
$ DOWNLOAD=0 ./mains.sh

# Perform comparision on large graphs with custom number of threads
$ DOWNLOAD=0 MAX_THREADS=4 ./mains.sh
```


<br>
<br>


## References

- [Fast unfolding of communities in large networks; Vincent D. Blondel et al. (2008)](https://arxiv.org/abs/0803.0476)
- [Community Detection on the GPU; Md. Naim et al. (2017)](https://arxiv.org/abs/1305.2006)
- [Scalable Static and Dynamic Community Detection Using Grappolo; Mahantesh Halappanavar et al. (2017)](https://ieeexplore.ieee.org/document/8091047)
- [From Louvain to Leiden: guaranteeing well-connected communities; V.A. Traag et al. (2019)](https://www.nature.com/articles/s41598-019-41695-z)
- [CS224W: Machine Learning with Graphs | Louvain Algorithm; Jure Leskovec (2021)](https://www.youtube.com/watch?v=0zuiLBOIcsw)
- [The University of Florida Sparse Matrix Collection; Timothy A. Davis et al. (2011)](https://doi.org/10.1145/2049662.2049663)
- [Fetch-and-add using OpenMP atomic operations](https://stackoverflow.com/a/7918281/1413259)

<br>
<br>


[![](https://i.imgur.com/UGB0g2L.jpg)](https://www.youtube.com/watch?v=pIF3wOet-zw)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)
