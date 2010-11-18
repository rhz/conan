#include <conan/graphs.hpp>
#include <conan/utils.hpp>
#include <conan/network_inference/maximum_entropy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <cstdlib>

int main(
    int argc,
    char * argv[]
    )
{
  typedef conan::undirected_graph<> Graph;
  typedef boost::numeric::ublas::matrix<double> matrix;

  double a[] = { 1.0, 2.0,
                 3.0, 4.0,
                 5.0, 6.0 };

  gsl_matrix_view A = gsl_matrix_view_array(a, 3, 2);

  // Graph g(conan::inference::maximum_entropy_with_fixed_num_edges<Graph>(&A.matrix, 2));
  Graph g(conan::inference::maximum_entropy_with_bootstrapping<Graph>(&A.matrix, 10, 1,
                                                                      conan::inference::covarianceS(), true, false));

  matrix M = conan::get_adj_matrix<Graph, matrix>(g);

  for (size_t r = 0; r < M.size1(); ++r)
  {
    std::cout << M(r, 0);
    for (size_t c = 1; c < M.size2(); ++c)
      std::cout << ' ' << M(r, c);

    std::cout << std::endl;
  }

  return(0);
}
