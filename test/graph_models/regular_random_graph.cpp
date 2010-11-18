#include <conan/graphs.hpp>
#include <conan/graph_models/regular_random.hpp>
#include <conan/properties/distributions.hpp>
#include <conan/io.hpp>

int main()
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;
  typedef conan::degree_distribution<Graph>::DegreeList DegreeList;
  typedef conan::degree_distribution<Graph>::VertexList VertexList;

  Graph g = conan::generate_random_regular_graph<Graph>(300, 8);
  conan::write_adj_matrix_file(g, "regular_random.txt");
  conan::write_dotfile(g, "regular_random.dot");

  conan::degree_distribution<Graph> d(g);
  DegreeList dl = d.non_zero_degrees();
  for (DegreeList::iterator li = dl.begin(); li != dl.end(); ++li)
  {
    std::cout << "P(" << *li << ") = " << d.P(*li) << std::endl;
    if (d.P(*li) < .1) // dirty workaround
    {
      std::cout << "vertices(" << *li << ") = ";
      VertexList vl = d.vertices(*li);
      for (VertexList::iterator vli = vl.begin(); vli != vl.end(); ++vli)
        std::cout << *vli << "  ";
      std::cout << std::endl;
    }
  }

  return(0);
}
