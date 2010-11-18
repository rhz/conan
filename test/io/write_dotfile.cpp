#include <conan/graphs.hpp>
#include <conan/io.hpp>
#include <conan/graph_models.hpp>

struct CustomVertexProperties { };

int main()
{
  typedef conan::undirected_graph<conan::adj_listS> GraphWithDefaultVertexProperties;
  typedef conan::undirected_graph<conan::adj_listS, conan::no_property> GraphWithoutVertexProperties;
  typedef conan::graph_traits<GraphWithDefaultVertexProperties>::vertex_iterator vertex_iter;

  GraphWithDefaultVertexProperties g1(
      conan::generate_erdos_renyi_graph<GraphWithDefaultVertexProperties>(5, .3));
  vertex_iter vi, viend;
  char chr = 'A';
  for (tie(vi, viend) = conan::vertices(g1); vi != viend; ++vi)
  {
    g1[*vi].name = chr++;
    std::cout << g1[*vi].name << std::endl;
  }
  conan::write_dotfile(g1, "g1.dot");

  GraphWithoutVertexProperties g2(
      conan::generate_erdos_renyi_graph<GraphWithoutVertexProperties>(5, .3));
  conan::write_dotfile(g2, "g2.dot");
  return(0);
}
