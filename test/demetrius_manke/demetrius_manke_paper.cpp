#include <conan/graphs.hpp>
#include <conan/graph_models.hpp>
#include <conan/properties.hpp>

int main()
{
  typedef conan::undirected_graph<> Graph;
  double T[] = {-100, 0, .25, 1};
  std::cout << "T\tV\tE\tL\tS\tdegree distribution entropy" << std::endl;
  for (size_t i = 0; i < 4; ++i)
  {
    Graph g(conan::generate_Demetrius_Manke_network<Graph>(7, /*2,*/ T[i]));
    for (size_t V = 10; V <= 300; V += 3)
    {
      conan::add_vertex_as_Demetrius_Manke(g, /*2,*/ T[i]);
      conan::add_vertex_as_Demetrius_Manke(g, /*2,*/ T[i]);
      conan::add_vertex_as_Demetrius_Manke(g, /*2,*/ T[i]);

      std::cout << T[i] << '\t' << conan::num_vertices(g) << '\t' << conan::num_edges(g) << '\t'
                << conan::graph_avg_shortest_path(g) << '\t' << conan::graph_entropy(g) << '\t';
      conan::degree_distribution<Graph> d(g);
      std::cout << d.entropy() << std::endl;
    }
  }
  return(0);
}
