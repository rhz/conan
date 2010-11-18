#include <conan/graphs.hpp>
#include <conan/properties.hpp>
#include <conan/graph_models.hpp>
#include <conan/utils.hpp>

int main()
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;
  int num_vertices = 100;
  Graph random_graph = conan::generate_random_graph<Graph>(num_vertices, 2 * num_vertices); // (V, E)
  std::cout << "random_graph:" << std::endl
            << "\tnum_vertices: " << boost::num_vertices(random_graph) << std::endl
            << "\tnum_edges: " << boost::num_edges(random_graph) << std::endl
            << "\tentropy: " << conan::graph_entropy(random_graph) << std::endl;
  conan::write_dotfile (random_graph,
                       "random_V" + conan::to_string<int>(boost::num_vertices(random_graph)) + \
                       "_E" + conan::to_string<int>(boost::num_edges(random_graph)) + ".dot");

  return EXIT_SUCCESS;
}
