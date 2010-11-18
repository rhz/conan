#include <conan/graphs.hpp>
#include <conan/properties.hpp>
#include <conan/utils.hpp>

int main(
    int argc,
    char* argv[]
    )
{
  if (argc != 2)
  {
    std::cout << "Usage: " << argv[0] << " <num_vertices>" << std::endl;
    exit(1);
  }

  const int V = conan::from_string<int>(argv[1]);
  double** m;
  m = (double**) malloc(sizeof(double*) * V);
  for (int i = 0; i < V; ++i)
    m[i] = (double*) malloc(sizeof(double) * V);

  for (int i = 0; i < V; ++i)
  {
    m[i][i] = 0;
    for (int j = i + 1; j < V; ++j)
    {
      m[i][j] = (i + j) % 2;
      m[j][i] = (i + j) % 2;
    }
  }

#if 0
  for (int i = 0; i < V; ++i)
  {
    for (int j = 0; j < V; ++j)
      std::cout << m[i][j] << '\t';
    std::cout << std::endl;
  }
  std::cout << std::endl;
#endif

  typedef conan::undirected_graph<conan::adj_matrixS> Graph;
  Graph g = conan::make_graph_from_adj_matrix<Graph>(m, V);
  conan::write_dotfile(g, "test_adj_matrix2.dot");
  double asp = conan::graph_avg_shortest_path(g);
  std::cout << "ASP = " << asp << std::endl;
  std::cout << "Clustering coef. = " << conan::graph_avg_clustering(g) << std::endl;
  std::cout << "Entropy = " << conan::graph_entropy(g) << std::endl;

  for (int i = 0; i < V; ++i)
    free(m[i]);
  free(m);

  return(0);
}
