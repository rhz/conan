#include <conan/graphs.hpp>
#include <conan/graph_models.hpp>
#include <conan/transformations.hpp>
#include <conan/io.hpp>

int main()
{
  typedef conan::directed_graph<conan::adj_listS> Graph;
  Graph g = conan::generate_erdos_renyi_graph<Graph>(14, .5);
  conan::write_adj_matrix_file(g, "g.txt");
  Graph g_transpose = conan::transpose(g);
  conan::write_adj_matrix_file(g_transpose, "g_transpose.txt");
  return(0);
}
