#include <conan/network_inference/pinv.hpp>

int main()
{
  gsl_matrix * Ainv;
  {
    double a[] = { 1.0, 2.0,
                   3.0, 4.0 };

    gsl_matrix_view A = gsl_matrix_view_array(a, 2, 2);
    Ainv = conan::linalg::pinv(&A.matrix);
  }

  for (size_t r = 0; r < 2; ++r)
  {
    std::cout << gsl_matrix_get(Ainv, r, 0);
    for (size_t c = 1; c < 2; ++c)
    {
      std::cout << '\t' << gsl_matrix_get(Ainv, r, c);
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  gsl_matrix * I;
  {
    double a[] = { 1.0, 2.0,
                   3.0, 4.0 };

    gsl_matrix_view A = gsl_matrix_view_array(a, 2, 2);
    I = gsl_matrix_calloc(2, 2);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Ainv, &A.matrix, 0.0, I);
  }

  for (size_t r = 0; r < 2; ++r)
  {
    std::cout << gsl_matrix_get(I, r, 0);
    for (size_t c = 1; c < 2; ++c)
    {
      std::cout << '\t' << gsl_matrix_get(I, r, c);
    }
    std::cout << std::endl;
  }

  gsl_matrix_free(Ainv);
  gsl_matrix_free(I);

  return(0);
}
