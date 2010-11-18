#include <conan/config.hpp>
#include <conan/properties.hpp>
#include <conan/topodb.hpp>
#include <conan/graphs.hpp>

int main(
    int argc,
    char* argv[]
    )
{
  typedef conan::decimal decimal;
  typedef conan::topodb<conan::undirectedS> TopoDB;
  typedef conan::undirected_graph<conan::adj_listS> Graph;
  
  if (argc != 2)
  {
    std::cout << "Usage: " << argv[0] << " <input_file>" << std::endl;
    return(1);
  }

  std::cout << "Reading input file " << argv[1] << std::endl;
  system("date");

  TopoDB big_network_db;
  big_network_db.read_tab_data_file(argv[1]);
  Graph big_network = big_network_db.make_graph<Graph>();

  std::cout << "File readed" << std::endl;
  system("date");

  std::cout << "Calculating network properties" << std::endl;
  std::cout << "ASP = " << conan::graph_avg_shortest_path(big_network) << std::endl;
  system("date");
  std::cout << "Clustering coef. = " << conan::graph_avg_clustering(big_network) << std::endl;
  system("date");
#if 0
  /*std::vector<int> dd = */conan::degree_distribution(big_network);
  std::cout << "Degree distribution not shown" << std::endl;
  system("date");
#endif
  std::cout << "Network properties calculated" << std::endl;

  return(0);
}
