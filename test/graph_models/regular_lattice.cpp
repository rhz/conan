#include <conan/graphs.hpp>
#include <conan/graph_models/canonical_graphs.hpp>
#include <conan/properties/distributions.hpp>
#include <conan/io.hpp>

int main(
    int argc,
    char* argv[]
    )
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;
  typedef conan::degree_distribution<Graph>::DegreeList DegreeList;

  if (argc != 3)
  {
    std::cerr << "Usage: " << argv[0] << " <num_vertices> <degree>" << std::endl;
    exit(1);
  }
  std::size_t num_vertices = conan::from_string<std::size_t>(argv[1]),
              degree = conan::from_string<std::size_t>(argv[2]);

  Graph g = conan::generate_regular_lattice<Graph>(num_vertices, degree);
  conan::write_dotfile(g, "regular_lattice.dot");

  conan::degree_distribution<Graph> d(g);
  DegreeList dl = d.non_zero_degrees();
  for (DegreeList::iterator li = dl.begin(); li != dl.end(); ++li)
  {
    std::cout << "P(" << *li << ") = " << d.P(*li) << std::endl;
  }

  return(0);
}
