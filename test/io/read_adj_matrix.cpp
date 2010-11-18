#include <conan/graphs.hpp>
#include <conan/io.hpp>

int main(
    int argc,
    char * argv[]
    )
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;

  if (argc != 2)
    return(1);

  Graph g = conan::read_adj_matrix_file<Graph>(argv[1], ',');

  std::cout << "V = " << conan::num_vertices(g) << std::endl
            << "E = " << conan::num_edges(g) << std::endl;

  conan::write_adj_matrix_file(g, "adj_matrix.txt", ',');

  return(0);
}
