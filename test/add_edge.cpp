#include <conan/graphs.hpp>

int main()
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;
  typedef boost::graph_traits<Graph>::edge_iterator edge_iter;

  Graph g(2);
  boost::add_edge(0, 1, 1.0, g);

  // print all edges... mmmmm... there's only one =P
  edge_iter ei, eiend;
  for (tie(ei, eiend) = boost::edges(g); ei != eiend; ++ei)
    std::cout << *ei << '\t' << g[*ei].weight << std::endl;

  return(0);
}
