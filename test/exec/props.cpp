#include <conan/graphs.hpp>
#include <conan/properties.hpp>
#include <conan/io.hpp>
#include <string.h>

template <typename T>
void print_property(const std::string & name, const T & value)
{
  std::cout << "  <property name=\"" << name << "\" value=\"" << value << "\" />" << std::endl;
}

int main(int argc, char* argv[])
{
  if (argc < 3 || argc > 9)
  {
    std::cout << "Usage: props [-V] [-E] [-asp] [-cc] [-S] example.graphml" << std::endl;
    exit(1);
  }

  typedef conan::undirected_graph<conan::adj_listS> Graph;
  Graph g = conan::read_graphml<Graph>(argv[argc - 1]);

  std::cout << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << std::endl;
  std::cout << "<graph name=\"" << argv[argc - 1] << "\">" << std::endl;
  for (int i = 1; i < argc - 1; ++i)
  {
    if ( !strcmp(argv[i], "-V") )
      print_property("number of vertices", boost::num_vertices(g));
    else if ( !strcmp(argv[i], "-E") )
      print_property("number of edges", boost::num_edges(g));
    else if ( !strcmp(argv[i], "-asp") )
      print_property("average shortest path", conan::graph_avg_shortest_path(g));
    else if ( !strcmp(argv[i], "-cc") )
      print_property("clustering coefficient", conan::graph_avg_clustering(g));
    else if ( !strcmp(argv[i], "-S") )
      print_property("topological entropy", conan::graph_entropy(g));
  }
  std::cout << "</graph>" << std::endl;
  return(0);
}
