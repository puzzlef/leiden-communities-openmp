I am trying out some experiments for a dynamic Leiden algoritm. I had optimized modularity earlier, and now i am trying to optimize the method for examining disconnected communities.

The base algorithm is as follows:
- Find the size of each community
- Pick one vertex from each community
- Traverse from that vertex within the community
- If all vertices within the community could not be reached, we know its disconnected

I try four different variations of this algorithm which differ in whether we use:
- Parallel DFS or BFS
- Per-thread or shared "visited" flags (shared flags is called light DFS/BFS)

If it is shared, each thread should scan all vertices but only process that community which belongs to it -> this is determined by community-id. It appears BFS traversal with a shared flag vector is the fastest. This is not a heuristic algorithm, so all algorithms return the same result.

> See
> [code](https://github.com/puzzlef/leiden-communities-openmp/tree/measure-disconnected-performance),
> [output](https://gist.github.com/wolfram77/8fa39d57a26833693cff51b4c2bd8f02), or
> [sheets].

<br>

[![](https://i.imgur.com/JqNtEaZ.png)][sheets]
[![](https://i.imgur.com/BWAd7Nr.png)][sheets]

[sheets]: https://docs.google.com/spreadsheets/d/13UZm9FJguOo1sMzsD_pb1NyuWJTR7iQ68_J7EXll7Es/edit?usp=sharing

<br>
<br>


## References

- [Fast unfolding of communities in large networks; Vincent D. Blondel et al. (2008)](https://arxiv.org/abs/0803.0476)
- [Community Detection on the GPU; Md. Naim et al. (2017)](https://arxiv.org/abs/1305.2006)
- [From Louvain to Leiden: guaranteeing well-connected communities; V.A. Traag et al. (2019)](https://www.nature.com/articles/s41598-019-41695-z)
- [CS224W: Machine Learning with Graphs | Louvain Algorithm; Jure Leskovec (2021)](https://www.youtube.com/watch?v=0zuiLBOIcsw)
- [The University of Florida Sparse Matrix Collection; Timothy A. Davis et al. (2011)](https://doi.org/10.1145/2049662.2049663)

<br>
<br>


[![](https://img.youtube.com/vi/dk8pwE3IByg/maxresdefault.jpg)](https://www.youtube.com/watch?v=dk8pwE3IByg)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)
