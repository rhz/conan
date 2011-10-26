#include <conan/graphs.hpp>
#include <conan/io.hpp>
#include <conan/inference/mutual_information.hpp>

int main(
    int argc,
    char * argv[]
    )
{
  if (argc < 2 || argc > 8)
  {
    std::cout << "Usage: aracne [-eps e] [-threshold t] [-sigma s] <aracne_input_filename>" << std::endl;
    exit(1);
  }

  // Parameter default values
  double eps = 0.1, // DPI tolerance
         threshold = 0.04, // MI threshold
         sigma = 0.15; // kernel width
  for (int i = 1; i < argc; ++i)
  {
    if ( !strcmp(argv[i], "-eps") )
      eps = conan::from_string<double>(argv[++i]);
    else if ( !strcmp(argv[i], "-threshold") )
      threshold = conan::from_string<double>(argv[++i]);
    else if ( !strcmp(argv[i], "-sigma") )
      sigma = conan::from_string<double>(argv[++i]);
  }

  aracne::Parameter p;
  p.infile = argv[argc - 1];

  p.eps = eps; // DPI tolerance
  p.threshold = threshold; // MI threshold
  p.sigma = sigma; // kernel width

  typedef conan::undirected_graph<conan::adj_listS> Graph;
  Graph g( conan::inference::mi::mutual_information<Graph>(p) );
  conan::write_graphml(g, p.outfile.substr(0, p.outfile.size() - 4) + ".graphml");
  return(0);
}
