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
  using namespace conan::inference;
  typedef conan::undirected_graph<> Graph;
  typedef boost::numeric::ublas::matrix<double> matrix;

  // Read data matrix file
  std::ifstream data_matrix_file("data_matrix_4.txt");
  size_t data_matrix_size1, data_matrix_size2;
  data_matrix_file >> data_matrix_size1 >> data_matrix_size2;

  gsl_matrix * data_matrix = gsl_matrix_calloc(data_matrix_size1, data_matrix_size2);

  size_t r = 0, c = 0;
  std::string line;
  while(!std::getline(data_matrix_file, line).eof())
  {
    if (line.empty())
      continue;

    double tmp = conan::from_string<double>(line);
    gsl_matrix_set(data_matrix, r, c++, tmp);

    if (c == data_matrix->size2)
    {
      c = 0;
      ++r;
    }
  }

  data_matrix_file.close();

  // Read matrix mask file
  std::ifstream matrix_mask_file("matrix_mask_4.txt");
  size_t matrix_mask_size1, matrix_mask_size2;
  matrix_mask_file >> matrix_mask_size1 >> matrix_mask_size2;
  matrix_mask_file.close();

  matrix_mask Amask(matrix_mask_size1, matrix_mask_size2);
  Amask.read("matrix_mask_4.txt", false);

  // Create the graph
  size_t min_num_complete_obs = 5;
  covariance_pairwise_complete_obs_S filter(Amask, min_num_complete_obs, true);
  Graph g(maximum_entropy_with_bootstrapping<Graph>(data_matrix, 100, 3, filter, false, false));

  gsl_matrix_free(data_matrix);

#if 0
  // Show the result
  matrix M = conan::get_adj_matrix<Graph, matrix>(g);

  for (r = 0; r < M.size1(); ++r)
  {
    std::cout << M(r, 0);
    for (c = 1; c < M.size2(); ++c)
      std::cout << ' ' << M(r, c);

    std::cout << std::endl;
  }
#endif

  return(0);
}
