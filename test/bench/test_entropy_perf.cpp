#include <conan/config.hpp>
#include <conan/graphs.hpp>
#include <conan/graph_models.hpp>
#include <conan/properties/entropy.hpp>

int main()
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;
  const int V = 1500, E = V * 6;
  Graph g = conan::generate_random_graph<Graph>(V, E);
  std::cout << conan::graph_entropy(g) << std::endl;
  return(0);
}
