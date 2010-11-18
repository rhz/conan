#include <conan/graphs.hpp>
#include <conan/io.hpp>
#include <conan/graph_models/random_graphs.hpp>

int main(
    int argc,
    char *argv[]
    )
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;

  if (argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " <input_filename>" << std::endl;
    return(1);
  }

  Graph g = conan::read_graphml<Graph>(argv[1]);
  conan::write_dotfile(g, "g1.dot");

  return(0);
}
