#include <conan/utils.hpp>
#include <conan/network_inference/pinv.hpp>

int main(
    int argc,
    char * argv[]
    )
{
  if (argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " <input_matrix_file>" << std::endl;
    return 1;
  }

  size_t num_rows, num_cols;
  num_rows = num_cols = 389;

  gsl_matrix * A = gsl_matrix_calloc(num_rows, num_cols);
  FILE * infile = fopen(argv[1], "rb");
  if (!infile)
  {
    std::cout << "Error opening input file " << argv[1] << std::endl;
    return(1);
  }

  int status = gsl_matrix_fread(infile, A);
  if (status != 0)
  {
    std::cout << "Error reading input matrix" << std::endl;
    return(1);
  }
  fclose(infile);

  gsl_matrix * Ainv = gsl_matrix_calloc(A->size2, A->size1);
  conan::linalg::pinv(A, Ainv);

  std::cout << "Printing inverse matrix" << std::endl;
  for (size_t r = 0; r < Ainv->size1; ++r)
  {
    std::cout << gsl_matrix_get(Ainv, r, 0);

    for (size_t c = 1; c < Ainv->size2; ++c)
      std::cout << '\t' << gsl_matrix_get(Ainv, r, c);

    std::cout << std::endl;
  }

  gsl_matrix_free(A);
  gsl_matrix_free(Ainv);

  return(0);
}
