Design of OpenMP-based Leiden algorithm for community detection.

I tried two variations of Leiden algorithm (like before):
- **Greedy Leiden**: Where the refinement phase is greedy, just like local-moving phase
- **Random Leiden**: Where refinement phase is proportionally random to delta-modularity, as proposed by Traag et al.

For random leiden, i use xorshift32 random number generator, C++ default RNG too slow. Further i try 3 different variations of Greedy/Random Leiden:
- **Normal**: Uses same paramter config as our optimized Louvain (tolerance=10^-2, 20 max iterations, 5 max passes)
- **Medium**: Uses lower start tolerance of 10^-6, 100 max iterations, 100 max passes
- **Heavy**: Uses even lower start tolerance of 10^-10, 100 max iterations, 100 max passes

The idea is to see if my normal Leiden configs are good or not. Please ignore road networks and k-mer protien graphs, as Leiden does not work there, will try to explain why so later. Here are the observations:
- Greedy Leiden is faster than Random Leiden
- Medium/Heady variations generally dont do well, so my normal configs are good
- Louvain is genrally faster than Leiden

[![](https://i.imgur.com/Sx3P0JP.png)][sheets]

Louvain generally obtain communities of higher modularity than Leiden.

[![](https://i.imgur.com/UZGMsNB.png)][sheets]

However, a large fraction of communities obtained by Louvain are disconnected.

[![](https://i.imgur.com/b3shYWt.png)][sheets]
[![](https://i.imgur.com/0qpSdhF.png)][sheets]

- Q> Why does Leiden have some disconnected communities still?
- Q> Why does Leiden not work with road networks?

I think the answer for both questions has to do with parallelism. There is a race between threads to pick a suitable community. If the size of each community is large, this generally does not badly affect modularity. But with small community bounds (in refinement phase of Leiden) the race can lead to bad community memberships. I will try to come up a few solutions to this. This is not an issue with sequential, so Traag et al. dont observe this issue in their [original paper][com-traag19].

> See
> [code](https://github.com/puzzlef/leiden-communities-openmp/tree/measure-performance),
> [output](https://gist.github.com/wolfram77/6e9a7772d1c71563572293b9044aeee4), or
> [sheets].

[com-traag19]: https://www.nature.com/articles/s41598-019-41695-z
[sheets]: https://docs.google.com/spreadsheets/d/1khv6NIhY1dWuXiuQE2qn2u8Om28saqH983MM13S6LCQ/edit?usp=sharing

<br>
<br>


## References

- [Fast unfolding of communities in large networks; Vincent D. Blondel et al. (2008)](https://arxiv.org/abs/0803.0476)
- [Community Detection on the GPU; Md. Naim et al. (2017)](https://arxiv.org/abs/1305.2006)
- [Scalable Static and Dynamic Community Detection Using Grappolo; Mahantesh Halappanavar et al. (2017)](https://ieeexplore.ieee.org/document/8091047)
- [From Louvain to Leiden: guaranteeing well-connected communities; V.A. Traag et al. (2019)](https://www.nature.com/articles/s41598-019-41695-z)
- [CS224W: Machine Learning with Graphs | Louvain Algorithm; Jure Leskovec (2021)](https://www.youtube.com/watch?v=0zuiLBOIcsw)
- [The University of Florida Sparse Matrix Collection; Timothy A. Davis et al. (2011)](https://doi.org/10.1145/2049662.2049663)

<br>
<br>


[![](https://img.youtube.com/vi/dk8pwE3IByg/maxresdefault.jpg)](https://www.youtube.com/watch?v=dk8pwE3IByg)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)
