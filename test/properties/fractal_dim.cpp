#include <conan/graphs.hpp>
#include <conan/graph_models.hpp>
#include <conan/properties/fractal_dim.hpp>

int main()
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;
  Graph g = conan::generate_scale_free_network<Graph>(100, 4);
  std::cout << conan::graph_fractal_dimension(g) << std::endl;
  return(0);
}
