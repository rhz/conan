#include <conan/graphs.hpp>
#include <conan/io.hpp>
#include <conan/subgraph/newman_communities.hpp>

using namespace conan;
using namespace std;

int main(int argc,
         char * argv[])
{
  if (argc != 2)
  {
    std::cout << "Usage: communities example.graphml" << std::endl;
    exit(1);
  }

  typedef undirected_graph<> Graph;
  typedef vector<size_t> Module;
  typedef vector<Module> ModuleVector;

  Graph g = conan::read_graphml<Graph>(argv[1]);
  ModuleVector modules;
  newman_communities( g, modules );

  BOOST_FOREACH( Module m, modules )
  {
    bool first = true;
    BOOST_FOREACH( size_t v, m )
    {
      if (!first)
        cout << "\t";

      if (!g[v].name.empty())
        cout << g[v].name;
      else
        cout << v;
      
      first = false;
    }
    cout << endl;
  }

  return(0);
}
