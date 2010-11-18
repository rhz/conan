#include <conan/graphs.hpp>
#include <conan/graph_models.hpp>
#include <conan/io.hpp>
#include <conan/utils.hpp>

template <class Graph>
Graph generate_random_graph(std::size_t num_vertices)
{
  return conan::generate_erdos_renyi_graph<Graph>(num_vertices, .65);
}

int main()
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;
  typedef conan::graph_traits<Graph>::edge_iterator edge_iter;

  Graph g = conan::generate_scale_free_network<Graph>(140, 5, &generate_random_graph);
  conan::write_dotfile(g, "scale_free.dot");

  // iterate over all edges and print them
  edge_iter ei, ei_end;
  for (tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
  {
    std::cout << *ei << " ";
  }
  std::cout << std::endl;

  return(0);
}
