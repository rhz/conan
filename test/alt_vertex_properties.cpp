#include <boost/graph/adjacency_list.hpp>

struct Node {
  int a;
  std::string b;
};

int main()
{
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Node> Graph;
  Graph g(1);
  Graph::vertex_iterator vi = vertices(g).first;
  g[*vi].a = 2;
  g[*vi].b = "hola";
  return(0);
}
