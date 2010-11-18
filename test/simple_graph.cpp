#include <conan/graphs.hpp>
#include <conan/properties.hpp>

int main()
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;
  Graph g(2);
  Graph::edge_descriptor e = boost::add_edge(0, 1, g).first;
  g[e].weight = 1.0;
  // Dijkstra have a problem with this graph
  conan::decimal asp = conan::graph_avg_shortest_path(g);
  std::cout << "num_vertices = " << conan::num_vertices(g) << std::endl
            << "num_edges = " << conan::num_edges(g) << std::endl
            << "asp = " << asp << std::endl
            << "entropy = " << conan::graph_entropy(g) << std::endl
            << "clustering = " << conan::graph_avg_clustering(g) << std::endl;
  if (asp != 1.0)
    std::cerr << "asp != 1.0, but it must be equal" << std::endl;
  return(0);
}
