#include <conan/graphs.hpp>
#include <conan/graph_models/barabasi_albert.hpp>
#include <conan/io.hpp>
#include <conan/utils.hpp>
#include <boost/numeric/ublas/matrix.hpp>

int main(
    int argc,
    char* argv[]
    )
{
  typedef conan::undirected_graph<> Graph;
  typedef boost::numeric::ublas::matrix<double> matrix;

  Graph g(conan::generate_scale_free_network<Graph>(10, 2));
  conan::write_adj_matrix_file(g, "g.adj_matrix");

  Graph g_copy(conan::read_adj_matrix_file<Graph>("g.adj_matrix"));

  matrix m = conan::get_adj_matrix<Graph, matrix>(g_copy);

  for (uint i = 0; i < m.size1(); ++i)
  {
    std::cout << m(i, 0);
    for (uint j = 1; j < m.size2(); ++j)
      std::cout << ' ' << m(i, j);

    std::cout << std::endl;
  }

  return(0);
}
