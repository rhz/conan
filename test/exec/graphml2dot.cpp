#include <conan/graphs.hpp>
#include <conan/io.hpp>

int main(int argc, char* argv[])
{
  if (argc != 2)
  {
    std::cout << "Usage: graphml2dot example.graphml" << std::endl;
    exit(1);
  }

  typedef conan::undirected_graph<conan::adj_listS> Graph;
  std::string infile = argv[argc - 1];
  Graph g = conan::read_graphml<Graph>(infile);
  size_t last_dot = infile.find_last_of(".");
  conan::write_dot(g, infile.substr(0, last_dot) + ".dot");
  return(0);
}
