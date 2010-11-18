#include <conan/graphs.hpp>
#include <conan/io.hpp>
#include <conan/network_inference/mutual_information.hpp>

int main(
    int argc,
    char * argv[]
    )
{
  using conan::decimal;
  using aracne::Parameter;
  typedef conan::undirected_graph<conan::adj_listS> Graph;

  if (argc != 2)
  {
    std::cout << "Usage: " << argv[0] << " <aracne_input_filename>" << std::endl;
    exit(1);
  }

  Parameter p;
  p.infile = argv[1];

  p.eps = 0.1; // DPI tolerance
  p.threshold = 0.04; // MI threshold
  p.sigma = 0.15; // kernel width

  Graph g( conan::inference::mi::mutual_information<Graph>(p) );

  conan::write_csv(g, p.outfile.substr(0, p.outfile.size() - 4) + ".csv");

  return(0);
}
