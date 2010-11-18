#include <conan/graphs.hpp>
#include <conan/subgraph.hpp>
#include <conan/io.hpp>
#include <conan/graph_models/random_graphs.hpp>
#include <conan/properties/modularity.hpp>

int main()
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;
  typedef conan::graph_traits<Graph>::vertex_descriptor vertex;
  typedef conan::graph_traits<Graph>::vertex_iterator vertex_iter;

  uint const V = 10;
  double const p = .3;
  Graph g(conan::generate_erdos_renyi_graph<Graph>(V, p));

  //conan::write_dot(g, "g.dot");

  vertex_iter vi, viend;
  for (tie(vi, viend) = conan::vertices(g); vi != viend; ++vi)
    g[*vi].name = conan::to_string(*vi);

  std::cout << "Erdos-Renyi graph generated" << std::endl;

  std::vector< std::vector<size_t> > modules;
  conan::newman_communities(g, modules);

  std::cout << "Graph splitted in " << modules.size() << " modules" << std::endl;

  /*
  std::vector<vertex> group[num_modules];
  for (size_t i = 0; i < V; ++i)
  {
    group[module[i]].push_back(i);
    std::cout << "module[" << i << "] = " << module[i] << std::endl;
  }

  for (size_t i = 0; i < num_modules; ++i)
  {
    Graph module_i( conan::subgraph(&g, group[i]) );
    conan::write_dot(module_i, "module" + conan::to_string(i) + ".dot");
  }
  */

  std::cout << "Finished" << std::endl;

  return(0);
}
