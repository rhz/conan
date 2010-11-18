#include <conan/graphs.hpp>
#include <conan/graph_models.hpp>
#include <conan/transformations.hpp>
#include <conan/io.hpp>

int main()
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;
  Graph g = conan::generate_erdos_renyi_graph<Graph>(14, .5);
  conan::write_adj_matrix_file(g, "g.txt");
  Graph g_complement = conan::complement(g);
  conan::write_adj_matrix_file(g_complement, "g_complement.txt");
  return(0);
}
