#include <conan/graphs.hpp>
#include <conan/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

int main(
    int argc,
    char * argv[]
    )
{
  if (argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " <input_dotfile>" << std::endl;
    return(1);
  }

  typedef conan::undirected_graph<> Graph;
  typedef conan::graph_traits<Graph>::vertex_iterator vertex_iter;
  typedef boost::numeric::ublas::matrix<double> matrix;

  Graph g(conan::read_dotfile<Graph>(argv[1]));

  vertex_iter vi, viend;
  for (tie(vi, viend) = conan::vertices(g); vi != viend; ++vi)
  {
    std::cout << *vi << " " << g[*vi].name << std::endl;
  }

  matrix m = conan::get_adj_matrix<Graph, matrix>(g);

  for (uint i = 0; i < m.size1(); ++i)
  {
    std::cout << m(i, 0);
    for (uint j = 1; j < m.size2(); ++j)
      std::cout << ' ' << m(i, j);

    std::cout << std::endl;
  }

  return(0);
}
