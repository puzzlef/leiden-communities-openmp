import os
import igraph as ig
import leidenalg


# Load graph
file = os.path.expanduser('~/Data/asia_osm.edgelist')
x = ig.read(file, format='edge', directed=False)

# Compute the best partition
ans = leidenalg.find_partition(x, leidenalg.ModularityVertexPartition)

# Print the modularity score
print(ans.modularity)
