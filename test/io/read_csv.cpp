#include <conan/graphs.hpp>
#include <conan/io.hpp>

int main(
    int argc,
    char *argv[]
    )
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;

  if (argc != 2)
    return(1);

  Graph g = conan::read_csv<Graph>(conan::to_string(argv[1]), "stw");
  conan::write_dotfile(g, "test_read_csv.dot");

  return(0);
}
