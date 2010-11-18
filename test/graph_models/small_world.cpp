#include <conan/graphs.hpp>
#include <conan/graph_models/watts_strogatz.hpp>
#include <conan/io.hpp>
#include <conan/utils.hpp>

int main()
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;

  Graph g = conan::generate_small_world_network<Graph>(140, 4, .3);
  conan::write_dotfile(g, "small_world.dot");

  return(0);
}
