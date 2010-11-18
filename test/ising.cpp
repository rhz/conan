#include <conan/graphs.hpp>
#include <conan/graph_models.hpp>
#include <conan/dynamics/ising.hpp>

int main()
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;
  typedef conan::ising<Graph> Ising;

  double T = 1000;
  Graph g(conan::generate_small_world_network<Graph>(30, 4, .3));
  std::cout << "Small-World network created" << std::endl;
  Ising i(g);
  std::cout << "Ising object created" << std::endl;
  i.equilibrate_system(T, 100);
  Graph new_g(i.run_montecarlo(T, 100));
  std::cout << "Simulation ran" << std::endl;

  return(0);
}
