#include <conan/config.hpp>
#include <conan/graphs.hpp>
#include <conan/graph_models.hpp>
#include <conan/properties.hpp>

int main(
    int argc,
    char* argv[]
    )
{
  typedef conan::decimal decimal;
  typedef conan::undirected_graph<conan::adj_listS> Graph;
  
  long num_vertices = 5000,
       avg_degree = 6,
       num_edges = num_vertices * avg_degree;
  std::cout << "Generating graph" << std::endl;
  system("date");

  std::cout << "V = " << num_vertices << std::endl
            << "E = " << num_edges    << std::endl;
  Graph big_network = conan::generate_scale_free_network<Graph>(num_vertices, avg_degree / 2);

  std::cout << "Graph generated" << std::endl;
  system("date");

  std::cout << "Calculating network properties" << std::endl;
  std::cout << "ASP = " << conan::graph_avg_shortest_path(big_network) << std::endl;
  system("date");
  std::cout << "Clustering coef. = " << conan::graph_avg_clustering(big_network) << std::endl;
  system("date");
#if 0
  std::vector<int> dd = conan::degree_distribution(big_network);
  std::cout << "Degree distribution not shown" << std::endl;
  system("date");
#endif
  std::cout << "Network properties calculated" << std::endl;

  return(0);
}
