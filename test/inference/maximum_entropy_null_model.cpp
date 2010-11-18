#include <conan/graphs.hpp>
#include <conan/utils.hpp>
#include <conan/io.hpp>
#include <conan/network_inference/maximum_entropy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <cstdlib>

#define NUM_GENES 10000
#define NUM_OBS   10

int main(
    int argc,
    char * argv[]
    )
{
  typedef conan::undirected_graph<> Graph;
  typedef boost::numeric::ublas::matrix<double> matrix;

  double a[NUM_GENES * NUM_OBS];

  boost::mt19937 rng;
  boost::uniform_real<> dist(0, 1);
  boost::variate_generator< boost::mt19937&, boost::uniform_real<> >
    random_number(rng, dist);

  for (size_t i = 0; i < NUM_GENES * NUM_OBS; ++i)
    a[i] = random_number();

  gsl_matrix_view A = gsl_matrix_view_array(a, NUM_GENES, NUM_OBS);

  std::cout << "Data matrix:" << std::endl;
  for (size_t r = 0; r < A.matrix.size1; ++r)
  {
    std::cout << gsl_matrix_get(&A.matrix, r, 0);

    for (size_t c = 1; c < A.matrix.size2; ++c)
      std::cout << ' ' << gsl_matrix_get(&A.matrix, r, c);

    std::cout << std::endl;
  }
  std::cout << std::endl;

  Graph g(conan::inference::maximum_entropy_with_fixed_num_edges<Graph>(&A.matrix, 100));

  conan::write_dot(g, "me_null_model.dot");

  matrix M = conan::get_adj_matrix<Graph, matrix>(g);

  return(0);
}
