#include <conan/graphs.hpp>
#include <conan/io.hpp>
#include <conan/graph_models/canonical_graphs.hpp>

int main(
    int argc,
    char *argv[]
    )
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;

  Graph graph_to_store = conan::generate_erdos_renyi_graph<Graph>(10, .2);
  conan::write_gml(graph_to_store, "g1.gml");
  conan::write_dotfile(graph_to_store, "g1.dot");
  Graph restored_graph = conan::read_gml<Graph>("g1.gml");

  return(0);
}
