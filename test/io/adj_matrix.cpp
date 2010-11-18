#include <conan/graphs.hpp>
#include <conan/graph_models.hpp>
#include <conan/io.hpp>

int main()
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;
  Graph g1 = conan::generate_random_graph<Graph>(20, 30);
  conan::write_dotfile(g1, "adj_matrix_1.dot");
  conan::write_adj_matrix_file(g1, "adj_matrix_1.txt", ',');
  Graph g2 = conan::read_adj_matrix_file<Graph>("adj_matrix_1.txt", ',');
  conan::write_dotfile(g2, "adj_matrix_2.dot");
  conan::write_adj_matrix_file(g2, "adj_matrix_2.txt", ',');
  return(0);
}
