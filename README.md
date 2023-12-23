Design of OpenMP-based [Leiden algorithm] for [community detection].

Community detection involves identifying subsets of vertices that display higher connectivity within themselves than with the rest of the network. The widely used Louvain method, a heuristic-based approach for community detection, employs a two-phase process consisting of a local-moving phase and an aggregation phase. This iterative optimization targets the modularity metric, a measure of community quality. Despite its popularity, the Louvain method has been observed to generate internally-disconnected and poorly connected communities. In response to these limitations, Traag et al. propose the Leiden algorithm, which introduces a refinement phase between the local-moving and aggregation phases. This refinement phase allows vertices to explore and potentially form sub-communities within the identified communities from the local-moving phase, enabling the Leiden algorithm to identify well-connected communities.

Nevertheless, the original Leiden algorithm encounters computational bottlenecks when applied to massive graphs, primarily due to its inherently sequential nature, akin to the Louvain method. In scenarios where scalability is crucial, the development of an optimized parallel Leiden algorithm becomes essential, especially in the multicore/shared memory setting, given its energy efficiency and the prevalence of hardware with large memory sizes. Despite existing studies proposing various parallelization techniques for the Leiden algorithm, they do not address optimization for the aggregation phase, which emerges as a bottleneck after optimizing the local-moving phase. Additionally, several optimization techniques applicable to the Louvain method are also relevant to the Leiden algorithm. To tackle these challenges, we present **GVE-Leiden**, an optimized *parallel implementation of the Leiden algorithm* designed for shared memory multicores.

Below we plot the time taken by the [original Leiden], [igraph] Leiden, [NetworKit] Leiden, and GVE-Leiden on 13 different graphs. GVE-Leiden surpasses the original Leiden, igraph Leiden, and NetworKit Leiden by `373×`, `86×`, and `7.2×` respectively, achieving a processing rate of `1.4𝐵` edges/s on a `3.8𝐵` edge graph.

[![](https://i.imgur.com/fiZPpQy.png)][sheets-o1]

Below we plot the speedup of GVE-Leiden wrt original Leiden, igraph Leiden, and NetworKit Leiden.

[![](https://i.imgur.com/le5WOJk.png)][sheets-o1]

Next, we plot the modularity of communities identified by original Leiden, igraph Leiden, NetworKit Leiden, and GVE-Leiden. GVE-Leiden on average obtains `0.1%` lower modularity than original Leiden and igraph Leiden, and `26%` higher modularity than NetworKit Leiden (especially on road networks and protein k-mer graphs).

[![](https://i.imgur.com/h6vOSE9.png)][sheets-o1]

Then, we plot the fraction of disconnected communities obtained with each implementation. Here, the absence of bars indicates the absence of disconnected communities. Communities identified by GVE-Leiden on average have `88×`, `145×`, and `0.76×` disconnected communities than the original Leiden, igraph Leiden, and NetworKit Leiden respectively. While this compares unfavorably with the original Leiden and igraph Leiden (especially on social networks, road networks, and protein k-mer graphs), it may be simpler to split the disconnected communities obtained from GVE-Leiden as a post-processing step. We would like to address this issue some time in the future.

[![](https://i.imgur.com/BLsirqM.png)][sheets-o1]

Finally, we plot the strong scaling behaviour of GVE-Leiden. With doubling of threads, GVE-Leiden exhibits an average performance scaling of `1.6×`.

[![](https://i.imgur.com/WxkvAUm.png)][sheets-o2]

Refer to our technical report for more details:
[GVE-Leiden: Fast Leiden Algorithm for Community Detection in Shared Memory Setting][report].

<br>

> [!NOTE]
> You can just copy `main.sh` to your system and run it. \
> For the code, refer to `main.cxx`.

[Leiden algorithm]: https://www.nature.com/articles/s41598-019-41695-z
[original Leiden]: https://github.com/vtraag/libleidenalg
[igraph]: https://github.com/igraph/igraph
[NetworKit]: https://github.com/networkit/networkit
[community detection]: https://en.wikipedia.org/wiki/Community_search
[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
[sheets-o1]: https://docs.google.com/spreadsheets/d/1oyx44kRewQmk9y23V-lJWwx51qYlGRxNR2kiT_tiMdo/edit?usp=sharing
[sheets-o2]: https://docs.google.com/spreadsheets/d/12CzNfXe3yO4NsOvs7sbmcwC6qBHTi6DmBF7yK2eXf7I/edit?usp=sharing
[report]: https://arxiv.org/abs/2312.13936

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
- inc/louvian.hxx: Louvian algorithm functions
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


[![](https://i.imgur.com/Z0g3W0u.jpg)](https://www.youtube.com/watch?v=yqO7wVBTuLw&pp)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)
[![DOI](https://zenodo.org/badge/652482935.svg)](https://zenodo.org/doi/10.5281/zenodo.10428321)


[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
