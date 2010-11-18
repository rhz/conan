#include <conan/graphs.hpp>
#include <conan/graph_models.hpp>
#include <conan/properties.hpp>

int main()
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;
  const int num_vertices = 800,
            num_edges_by_vertex = 5;
  int i = 0;
  
  // Preferential Attachment
  for (i = 0; i < 5; ++i)
  {
    Graph g = conan::generate_scale_free_network<Graph>(num_vertices, num_edges_by_vertex);
    std::cout << "V = " << boost::num_vertices(g) << '\t'
              << "E = " << boost::num_edges(g) << '\t'
              << "Clustering = " << conan::graph_avg_clustering(g) << '\t'
              << "ASP = " << conan::graph_avg_shortest_path(g) << '\t'
              << "Entropy = " << conan::graph_entropy(g) << std::endl;
  }
 
  std::cout << std::endl;
  // Small-World
  for (i = 0; i < 5; ++i)
  {
    Graph g = conan::generate_small_world_network<Graph>(num_vertices, 3, .2);
    std::cout << "V = " << boost::num_vertices(g) << '\t'
              << "E = " << boost::num_edges(g) << '\t'
              << "Clustering = " << conan::graph_avg_clustering(g) << '\t'
              << "ASP = " << conan::graph_avg_shortest_path(g) << '\t'
              << "Entropy = " << conan::graph_entropy(g) << std::endl;
  }
  return(0);
}
