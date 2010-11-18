#include <conan/graphs.hpp>
#include <conan/graph_models/fibonacci.hpp>
#include <conan/properties.hpp>

int main(
    int argc,
    char* argv[]
    )
{
  //const int num_iterations = 10;
  if (argc != 2)
  {
    std::cout << "Usage: " << argv[0] << " <num_iterations>" << std::endl;
    exit(1);
  }

  typedef conan::undirected_graph<conan::adj_listS> Graph;
  int num_iterations = conan::from_string<int>(argv[1]);

  Graph g = conan::generate_fibonacci_graph<Graph>(num_iterations);
  conan::write_dotfile(g, "fibonacci_" + conan::to_string(argv[1]) + ".dot");

  std::cout << "Number of vertices = " << conan::num_vertices(g) << std::endl;
  std::cout << "Number of edges = " << conan::num_edges(g) << std::endl;
  std::cout << "ASP = " << conan::graph_avg_shortest_path(g) << std::endl;
  std::cout << "Clustering coef. = " << conan::graph_avg_clustering(g) << std::endl;
  std::cout << "Entropy = " << conan::graph_entropy(g) << std::endl;

  return(0);
}
