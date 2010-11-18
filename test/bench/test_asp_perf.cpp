#include <conan/config.hpp>
#include <conan/graphs.hpp>
#include <conan/graph_models.hpp>
#include <conan/properties/asp.hpp>

int main()
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;
  const int V = 1500, edges_by_vertex = 4;
  Graph g = conan::generate_scale_free_network<Graph>(V, edges_by_vertex);
  std::cout << conan::graph_avg_shortest_path(g) << std::endl;
  return(0);
}
