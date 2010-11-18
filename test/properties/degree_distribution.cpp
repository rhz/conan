#include <conan/graphs.hpp>
#include <conan/graph_models.hpp>
#include <conan/properties/distributions.hpp>

int main()
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;
  typedef conan::degree_distribution<Graph>::DegreeList DegreeList;
  typedef conan::degree_distribution<Graph>::VertexList VertexList;

  //Graph g = conan::generate_fibonacci_graph<Graph>(6);
  Graph g = conan::generate_scale_free_network<Graph>(100, 2);
  std::cout << "Number of vertices = " << conan::num_vertices(g) << std::endl;
  std::cout << "Number of edges = " << conan::num_edges(g) << std::endl;
 
  system("date");
  conan::degree_distribution<Graph> d(g);
  system("date");

  /*
  DegreeList dl = d.non_zero_degrees();
  for (DegreeList::iterator li = dl.begin(); li != dl.end(); ++li)
  {
    std::cout << "P(" << *li << ") = " << d.P(*li) << std::endl;
  }

  VertexList::iterator vi, viend;
  for (tie(vi, viend) = d.vertices(3); vi != viend; ++vi)
    std::cout << *vi << std::endl;

  for (tie(vi, viend) = d.vertices(1); vi != viend; ++vi)
    std::cout << *vi << std::endl;
  */

  std::cout << "The degree distribution was best fitted by a ";
  if (d.best_fit == conan::linear)
    std::cout << "linear ";
  else if (d.best_fit == conan::power)
    std::cout << "power ";
  else if (d.best_fit == conan::exponential)
    std::cout << "exponential ";
  else if (d.best_fit == conan::logaritmic)
    std::cout << "logaritmic ";
  else
    std::cout << "unknown ";
  std::cout << "regression" << std::endl;

  double* best_fit_params = d._fit_parameters(d.best_fit);
  std::cout << "r^2 = " << best_fit_params[0] * best_fit_params[0] << std::endl;
  std::cout << "c_1 = " << best_fit_params[1] << std::endl;
  std::cout << "c_0 = " << best_fit_params[2] << std::endl << std::endl;

  std::cout << "entropy = " << d.entropy() << std::endl;

  return(0);
}
