#include <conan/graphs.hpp>
#include <conan/io.hpp>
#include <conan/graph_models/random_graphs.hpp>

using namespace conan;

int main()
{
  typedef undirected_graph<> Graph;
  Graph g1( generate_erdos_renyi_graph<Graph>(20, .2) );
  write_dot(g1, "g1.dot");
  Graph g2( generate_erdos_renyi_graph<Graph>(20, .2) );
  write_dot(g2, "g2.dot");
  std::list<Graph> graph_list;
  graph_list.push_back(g1);
  graph_list.push_back(g2);
  Graph g3( merge_graphs(graph_list) );
  write_dot(g3, "g3.dot");
  return(0);
}
