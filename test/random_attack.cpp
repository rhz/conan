#include <conan/graphs.hpp>
#include <conan/graph_models.hpp>
#include <conan/transformations.hpp>
#include <conan/properties.hpp>

int main()
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;
  Graph original_graph(conan::generate_scale_free_network<Graph>(40, 3));
  for (size_t i = 0; i <= 10; ++i)
  {
    double p = i * .1;
    Graph g(conan::random_attack(original_graph, p));
    std::cout << "num vertices = " << conan::num_vertices(g) << std::endl
              << "num edges = " << conan::num_edges(g) << std::endl
              << "asp = " << conan::graph_avg_shortest_path(g) << std::endl;
  }
  return(0);
}
