#include <conan/graphs.hpp>
#include <conan/graph_models/random_graphs.hpp>
#include <conan/properties/color.hpp>

int main()
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;

  Graph g( conan::generate_erdos_renyi_graph<Graph>(30, .4) );

  std::cout << "Erdos-Renyi graph:" << std::endl
            << "\tnum vertices = " << conan::num_vertices(g) << std::endl
            << "\tnum edges = " << conan::num_edges(g) << std::endl
            << "\tchromatic number = " << conan::chromatic_number(g) << std::endl;

  return(0);
}
