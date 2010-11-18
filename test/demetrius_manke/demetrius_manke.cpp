#include <conan/graphs.hpp>
#include <conan/graph_models/demetrius_manke.hpp>
#include <conan/io.hpp>
#include <conan/properties.hpp>

int main()
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;

  double T = -100;
  Graph g = conan::generate_Demetrius_Manke_network<Graph>(30, T, 2, false);
  conan::write_dotfile(g, "demetrius_manke.dot");
  std::cout << "Entropy = " << conan::graph_entropy(g) << std::endl;

  return(0);
}
