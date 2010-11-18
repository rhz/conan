#include <conan/graphs.hpp>
#include <conan/io.hpp>
#include <conan/graph_models/random_graphs.hpp>

int main(
    int argc,
    char *argv[]
    )
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;

  Graph graph_to_store = conan::generate_erdos_renyi_graph<Graph>(10, .2);
  conan::write_graphml(graph_to_store, "g1.graphml");
  conan::write_dotfile(graph_to_store, "g1.dot");
  Graph restored_graph = conan::read_graphml<Graph>("g1.graphml");

  return(0);
}
