#include <conan/graphs.hpp>
#include <conan/subgraph.hpp>
#include <conan/io.hpp>
#include <conan/properties/modularity.hpp>

int main(
    int argc,
    char * argv[]
    )
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;
  typedef conan::graph_traits<Graph> GraphTraits;
  typedef GraphTraits::vertex_descriptor vertex;

  if (argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " <input_gml_file>" << std::endl;
    exit(1);
  }

  Graph g(conan::read_gml<Graph>(argv[1]));

  conan::write_dotfile(g, "g.dot");

  size_t V = conan::num_vertices(g);

  std::cout << V << std::endl;

  std::vector<uint> module(V);
  conan::split_in_two(g, module);

  std::cout << "Graph splitted in two modules" << std::endl;

  std::vector<vertex> group[2]; // an array of vectors with the vertex_descriptors in each group
  for (size_t i = 0; i < V; ++i)
  {
    group[module[i]].push_back(i);
    std::cout << "module[" << i << "] = " << module[i] << std::endl;
  }

  for (size_t i = 0; i < 2; ++i)
  {
    Graph module_i( conan::subgraph(&g, group[i]) );
    conan::write_dotfile(module_i, "module" + conan::to_string(i) + ".dot");
    std::cout << conan::num_vertices(module_i) << std::endl;
  }

  std::cout << "Finished" << std::endl;

  return(0);
}
