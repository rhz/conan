#include <conan/graphs.hpp>
#include <conan/network_inference/base.hpp>

int main()
{
  double a[] = { 1.0, 2.0,
                 3.0, 4.0,
                 5.0, 6.0 };

  gsl_matrix_view A = gsl_matrix_view_array(a, 3, 2);

  conan::inference::matrix_mask Amask(3, 2);
  Amask(0, 0) = true;
  Amask(2, 1) = true;

  gsl_matrix * covar = gsl_matrix_alloc(3, 3);

  calc_covariance_matrix_from_pairwise_complete_obs(&A.matrix, Amask, 1, covar);

  gsl_matrix_free(covar);

  return(0);
}
