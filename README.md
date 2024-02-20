Design of OpenMP-based Parallel [Leiden algorithm][Leiden] for [community detection], \
that prevents internally disconnected communities.

> [!NOTE]
> For the code of [GVE-Leiden][report1], refer to the [arXiv-2312.13936] branch.

<br>

Community detection entails the identification of clusters of vertices that exhibit stronger connections within themselves compared to the wider network. The Louvain method, a commonly utilized heuristic for this task, employs a two-step process comprising a local-moving phase and an aggregation phase. This process iteratively optimizes the modularity metric, a measure of community quality. Despite its popularity, the Louvain method has been noted for producing internally fragmented and weakly connected communities. In response to these limitations, Traag et al. propose the Leiden algorithm, which incorporates a refinement phase between the local-moving and aggregation phases. This refinement step enables vertices to explore and potentially establish sub-communities within the identified communities from the local-moving phase.

However, the Leiden algorithm is not guaranteed to avoid internally disconnected communities, a flaw that has largely escaped attention. We illustrate this through both a [counterexample][report1] and [empirical findings][report1]. In our experimental evaluation, we note that approximately `1.3×10^−4` fraction of the communities identified using the original Leiden implementation exhibit this issue. Although this fraction is small, addressing the presence of disconnected communities is crucial for ensuring the accuracy and dependability of community detection algorithms. Several studies have addressed internally disconnected communities as a post-processing step. However, this may exacerbate the problem of poorly connected communities. Furthermore, the surge in data volume and their graph representations in recent years has been unprecedented. Nonetheless, applying the original Leiden algorithm to massive graphs has posed computational hurdles, primarily due to its inherently sequential nature, akin to the Louvain method. To tackle these challenged, we propose two new *parallel algorithms*: **[GSP-Leiden]** and **[GSP-Louvain]**, based on the [Leiden] and [Louvain] algorithms, respectively.

Below we plot the time taken by the [original Leiden], [igraph] Leiden, [NetworKit] Leiden, GSP-Leiden, and GSP-Louvain on 13 different graphs. GSP-Leiden surpasses the original Leiden, igraph Leiden, and NetworKit Leiden by `190×`, `46×`, and `3.4×` respectively, achieving a processing rate of `195M` edges/s on a `3.8𝐵` edge graph.

[![](https://i.imgur.com/bgTuZsm.png)][sheets-o1]

Below we plot the speedup of GSP-Leiden and GSP-Louvain wrt original Leiden, igraph Leiden, and NetworKit Leiden.

[![](https://i.imgur.com/8jtfe7p.png)][sheets-o1]

Next, we compare the modularity of communities identified by the original Leiden algorithm, igraph Leiden, NetworKit Leiden, GSP-Leiden, and GSP-Leiden. On average, GSP-Leiden achieves `0.07%` and `0.02%` lower modularity than the original Leiden and igraph Leiden, respectively, and `26%` higher modularity than NetworKit Leiden, particularly evident on road networks and protein k-mer graphs.

[![](https://i.imgur.com/gKKH1dg.png)][sheets-o1]

Finally, we plot the fraction of disconnected communities identified by each implementation. Absence of bars indicates the absence of disconnected communities. As anticipated, both GSP-Leiden and GSP-Louvain detect no disconnected communities. However, on average, the original Leiden, igraph Leiden, and NetworKit Leiden exhibit fractions of disconnected communities amounting to `1.3×10^−4`, `7.9×10^−5`, and `1.5×10^−2`, respectively, particularly on web graphs (and especially on social networks with NetworKit Leiden).

[![](https://i.imgur.com/FgI5GT9.png)][sheets-o1]

Refer to our technical reports for more details: \
[GVE-Leiden: Fast Leiden Algorithm for Community Detection in Shared Memory Setting][report1]. \
[Addressing Internally-Disconnected Communities in Leiden and Louvain Community Detection Algorithms][report2].

<br>

> [!NOTE]
> You can just copy `main.sh` to your system and run it. \
> For the code, refer to `main.cxx`.

[Leiden]: https://www.nature.com/articles/s41598-019-41695-z
[Louvain]: https://arxiv.org/abs/0803.0476
[original Leiden]: https://github.com/vtraag/libleidenalg
[igraph]: https://github.com/igraph/igraph
[NetworKit]: https://github.com/networkit/networkit
[community detection]: https://en.wikipedia.org/wiki/Community_search
[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
[sheets-o1]: https://docs.google.com/spreadsheets/d/1oyx44kRewQmk9y23V-lJWwx51qYlGRxNR2kiT_tiMdo/edit?usp=sharing
[sheets-o2]: https://docs.google.com/spreadsheets/d/12CzNfXe3yO4NsOvs7sbmcwC6qBHTi6DmBF7yK2eXf7I/edit?usp=sharing
[report1]: https://arxiv.org/abs/2312.13936
[report2]: https://arxiv.org/abs/2402.11454
[GSP-Leiden]: https://github.com/puzzlef/leiden-communities-openmp
[GSP-Louvain]: https://github.com/puzzlef/louvain-communities-openmp
[arXiv-2312.13936]: https://github.com/puzzlef/leiden-communities-openmp/tree/arXiv-2312.13936

<br>
<br>


### Code structure

The code structure of GVE-Leiden is as follows:

```bash
- inc/_algorithm.hxx: Algorithm utility functions
- inc/_bitset.hxx: Bitset manipulation functions
- inc/_cmath.hxx: Math functions
- inc/_ctypes.hxx: Data type utility functions
- inc/_cuda.hxx: CUDA utility functions
- inc/_debug.hxx: Debugging macros (LOG, ASSERT, ...)
- inc/_iostream.hxx: Input/output stream functions
- inc/_iterator.hxx: Iterator utility functions
- inc/_main.hxx: Main program header
- inc/_mpi.hxx: MPI (Message Passing Interface) utility functions
- inc/_openmp.hxx: OpenMP utility functions
- inc/_queue.hxx: Queue utility functions
- inc/_random.hxx: Random number generation functions
- inc/_string.hxx: String utility functions
- inc/_utility.hxx: Runtime measurement functions
- inc/_vector.hxx: Vector utility functions
- inc/batch.hxx: Batch update generation functions
- inc/bfs.hxx: Breadth-first search algorithms
- inc/csr.hxx: Compressed Sparse Row (CSR) data structure functions
- inc/dfs.hxx: Depth-first search algorithms
- inc/duplicate.hxx: Graph duplicating functions
- inc/Graph.hxx: Graph data structure functions
- inc/leiden.hxx: Leiden algorithm functions
- inc/leidenSplit.hxx: Leiden with no disconnected communities
- inc/louvian.hxx: Louvian algorithm functions
- inc/louvainSplit.hxx: Louvain with no disconnected communities
- inc/main.hxx: Main header
- inc/mtx.hxx: Graph file reading functions
- inc/properties.hxx: Graph Property functions
- inc/selfLoop.hxx: Graph Self-looping functions
- inc/symmetricize.hxx: Graph Symmetricization functions
- inc/transpose.hxx: Graph transpose functions
- inc/update.hxx: Update functions
- main.cxx: Experimentation code
- process.js: Node.js script for processing output logs
```

Note that each branch in this repository contains code for a specific experiment. The `main` branch contains code for the final experiment. If the intention of a branch in unclear, or if you have comments on our technical report, feel free to open an issue.

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


[![](https://i.imgur.com/atJbkL1.png)](https://www.youtube.com/watch?v=yqO7wVBTuLw&pp)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)
[![DOI](https://zenodo.org/badge/652482935.svg)](https://zenodo.org/doi/10.5281/zenodo.10428321)


[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
