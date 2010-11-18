#!/usr/bin/env python

from sys import argv
import conan

g = conan.undirected_graph()
g.read_adj_matrix_file(argv[1])

for v in g.vertices:
  v.name = str( v.id() )

# Select the main component from network
c = g.components()
max_num_vertices = 0
index = 0
main_component_index = 0
for subg in c:
  if subg.num_vertices() > max_num_vertices:
    main_component_index = index
    max_num_vertices = subg.num_vertices()
  index += 1

main_component = c[main_component_index]

m = main_component.newman_communities()
