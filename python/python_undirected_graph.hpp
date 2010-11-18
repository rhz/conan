#ifndef PYTHON_UNDIRECTED_GRAPH_HPP
#define PYTHON_UNDIRECTED_GRAPH_HPP

#include <conan/graphs.hpp>
#include <conan/subgraph.hpp>
#include <conan/subgraph/roles.hpp>
#include <conan/properties.hpp>
#include <conan/io.hpp>
#include <conan/utils.hpp>
#include <conan/transformations.hpp>
#include "vertex_and_edge_wrap.hpp"

#include <boost/python.hpp>

using namespace boost::python;


// undirected_graph class
struct python_undirected_graph
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;
  typedef boost::graph_traits<Graph> GraphTraits;
  typedef GraphTraits::vertex_descriptor vertex;
  typedef GraphTraits::vertex_iterator vertex_iterator;
  typedef GraphTraits::edge_descriptor edge;
  typedef GraphTraits::edge_iterator edge_iterator;

  typedef VertexWrap<Graph> Vertex;
  typedef EdgeWrap<Graph> Edge;

  typedef std::vector<vertex> Module;
  typedef std::vector<Module> ModuleVector;

  typedef conan::decimal decimal;


  inline python_undirected_graph()
    : g() { }

  inline python_undirected_graph(const Graph& x)
    : g(x) { }

  inline python_undirected_graph(const python_undirected_graph& x)
    : g(x.g) { }

  inline python_undirected_graph(size_t num_vertices)
    : g(num_vertices) { }

  inline python_undirected_graph(list adj_matrix)
  {
    size_t num_vertices = extract<size_t>(adj_matrix.attr("__len__")());

    double **c_adj_matrix;
    c_adj_matrix = (double**) malloc(sizeof(double*) * num_vertices);
    for (size_t i = 0; i < num_vertices; ++i)
      c_adj_matrix[i] = (double*) malloc(sizeof(double) * num_vertices);

    for (size_t i = 0; i < num_vertices; ++i)
    {
      list row = extract<list>(adj_matrix[i]);
      size_t row_length = extract<size_t>(row.attr("__len__")());
      if (row_length != num_vertices)
        throw std::runtime_error("In undirected_graph constructor: input adjacency matrix format error");

      for (size_t j = 0; j < num_vertices; ++j)
      {
        c_adj_matrix[i][j] = extract<double>(row[j]);
      }
    }

    g = conan::make_graph_from_adj_matrix<Graph>(c_adj_matrix, num_vertices);

    for (size_t i = 0; i < num_vertices; ++i)
      free(c_adj_matrix[i]);
    free(c_adj_matrix);
    return;
  }

  list vertices()
  {
    list l;

    vertex_iterator vi, viend;
    for (tie(vi, viend) = boost::vertices(g); vi != viend; ++vi)
      l.append( Vertex(*vi, g) );

    return l;
  }

  list edges()
  {
    list l;

    edge_iterator ei, eiend;
    for (tie(ei, eiend) = boost::edges(g); ei != eiend; ++ei)
      l.append( Edge(*ei, g) );

    return l;
  }

  size_t num_vertices()
  { return boost::num_vertices(g); }

  size_t num_edges()
  { return boost::num_edges(g); }

  Vertex add_vertex()
  { return Vertex(boost::add_vertex(g), g); }

  Vertex add_vertex(std::string name)
  { return Vertex(boost::add_vertex(name, g), g); }

  Edge add_edge(Vertex a, Vertex b)
  {
    if (boost::edge(a.v, b.v, g).second)
      throw std::runtime_error("Edge between vertices " + conan::to_string(a.v) + " and " + conan::to_string(b.v) + " already exists");

    bool edge_created;
    edge e;
    tie(e, edge_created) = boost::add_edge(a.v, b.v, g);
    if (!edge_created)
      throw std::runtime_error("Edge between vertices " + conan::to_string(a.v) + " and " + conan::to_string(b.v) + " could not be added to the graph");

    g[e].weight = 1.0;
    return Edge(e, g);
  }

  Edge add_edge(vertex a, vertex b)
  {
    if (boost::edge(a, b, g).second)
      throw std::runtime_error("Edge between vertices " + conan::to_string(a) + " and " + conan::to_string(b) + " already exists");

    bool edge_created;
    edge e;
    tie(e, edge_created) = boost::add_edge(a, b, g);
    if (!edge_created)
      throw std::runtime_error("Edge between vertices " + conan::to_string(a) + " and " + conan::to_string(b) + " could not be added to the graph");

    g[e].weight = 1.0;
    return Edge(e, g);
  }

  void remove_vertex(Vertex a)
  {
    boost::clear_vertex(a.v, g);
    return boost::remove_vertex(a.v, g);
  }

  void remove_vertex(vertex a)
  {
    boost::clear_vertex(a, g);
    return boost::remove_vertex(a, g);
  }

  void remove_edge(Edge e)
  {
    if (& e.parent_graph != & g)
      throw std::runtime_error("Given edge does not belong to graph");

    vertex source_vertex = boost::source(e.e, g),
           target_vertex = boost::target(e.e, g);

    if (boost::edge(source_vertex, target_vertex, g).second)
      return boost::remove_edge(source_vertex, target_vertex, g);
    else
      throw std::runtime_error("Edge between vertices " + conan::to_string(source_vertex) + " and " + conan::to_string(target_vertex) + " does not exist");
  }

  void remove_edge(Vertex a, Vertex b)
  {
    if (boost::edge(a.v, b.v, g).second)
      return boost::remove_edge(a.v, b.v, g);
    else
      throw std::runtime_error("Edge between vertices " + conan::to_string(a.v) + " and " + conan::to_string(b.v) + " does not exists");
  }

  void remove_edge(vertex a, vertex b)
  {
    if (boost::edge(a, b, g).second)
      return boost::remove_edge(a, b, g);
    else
      throw std::runtime_error("Edge between vertices " + conan::to_string(a) + " and " + conan::to_string(b) + " does not exists");
  }

  Vertex get_vertex(vertex a)
  { return Vertex(a, g); }

  Vertex get_vertex(std::string name)
  {
    vertex a = find_vertex_by_name(name, g);
    if (a == boost::graph_traits<Graph>::null_vertex())
      throw std::runtime_error("Vertex \"" + name + "\" not found");
    return Vertex(a, g);
  }

  Edge get_edge(vertex a, vertex b)
  { return Edge(a, b, g); }

  Edge get_edge(Vertex a, Vertex b)
  { return Edge(a.v, b.v, g); }

  bool has_edge(Edge e)
  { return conan::edge(e.get_source_vertex().v, e.get_target_vertex().v, g).second; }

  bool has_edge(Vertex a, Vertex b)
  { return conan::edge(a.v, b.v, g).second; }

  bool has_edge(vertex a, vertex b)
  { return conan::edge(a, b, g).second; }

  list get_adj_matrix()
  {
    typedef ublas::matrix<decimal> matrix;
    matrix m(conan::get_adj_matrix<Graph, matrix>(g));
    list adj_matrix_list;
    for (size_t i = 0; i < m.size1(); ++i)
    {
      std::vector<decimal> row;
      for (size_t j = 0; j < m.size2(); ++j)
      {
        row.push_back(m(i, j));
      }
      adj_matrix_list.append(list(row));
    }
    return adj_matrix_list;
  }

  python_undirected_graph subgraph(list l)
  {
    std::set<vertex> vertex_set;
    std::set<edge> edge_set;

    size_t vertex_list_length = extract<size_t>( l.attr( "__len__" )() );

    for (size_t i = 0; i < vertex_list_length; ++i)
    {
      extract<Vertex> get_vertex( l[i] );

      if ( get_vertex.check() )
      {
        Vertex v = get_vertex();
        vertex_set.insert( v.v );
      }
      else
      {
        extract<Edge> get_edge( l[i] );

        if ( get_edge.check() )
        {
          Edge e = get_edge();
          edge_set.insert( e.e );
        }
        else
        {
          throw std::runtime_error("Input list should contain only vertices or edges.");
        }
      }
    }

    if (vertex_set.size() != 0 && edge_set.size() != 0)
    {
      throw std::runtime_error("Input list should contain only vertices or edges, not both.");
    }
    else if (vertex_set.size() != 0)
    {
      return python_undirected_graph( conan::induced_subgraph< Graph, std::set<vertex> >( g, vertex_set ) );
    }
    else if (edge_set.size() != 0)
    {
      return python_undirected_graph( conan::make_subgraph_from_edge_set< Graph, std::set<edge> >( g, edge_set ) );
    }
    else
    {
      throw std::runtime_error("There were neither vertices nor edges in input list.");
    }
  }

  list components()
  {
    typedef std::vector<vertex> Module;
    typedef std::vector<Module> ModuleVector;

    ModuleVector all_components;

    conan::components(g, all_components);

    list components_list;
    BOOST_FOREACH( Module & m, all_components )
    {
      list component;
      BOOST_FOREACH( size_t & v, m )
      {
        component.append( Vertex(v, g) );
      }
      components_list.append( component );
    }

    return components_list;
  }

  list newman_communities(bool check_quality_function = true)
  {
    typedef std::vector<vertex> Module;
    typedef std::vector<Module> ModuleVector;

    ModuleVector modules;

    conan::newman_communities(g, modules, check_quality_function);

    list communities;
    BOOST_FOREACH( Module & m, modules )
    {
      list community;
      BOOST_FOREACH( size_t & v, m )
      {
        community.append( Vertex(v, g) );
      }
      communities.append( community );
    }

    return communities;
  }

  list newman_two_communities()
  {
    typedef std::vector<vertex> Module;
    typedef std::vector<Module> ModuleVector;

    ModuleVector modules;

    conan::newman_two_communities(g, modules);

    list communities;
    BOOST_FOREACH( Module & m, modules )
    {
      list community;
      BOOST_FOREACH( size_t & v, m )
      {
        community.append( Vertex(v, g) );
      }
      communities.append( community );
    }

    return communities;
  }

  list guimera_amaral_communities(double initial_temp = .0001, double final_temp = .00006, double f = 1.0)
  {
    typedef std::vector<vertex> Module;
    typedef std::vector<Module> ModuleVector;

    ModuleVector modules;

    conan::guimera_amaral_communities(g, modules, initial_temp, final_temp, f);

    list communities;
    BOOST_FOREACH( Module & m, modules )
    {
      list community;
      BOOST_FOREACH( size_t & v, m )
      {
        community.append( Vertex(v, g) );
      }
      communities.append( community );
    }

    return communities;
  }

  ModuleVector get_module_vector(list modules) // FIXME: missing in python_directed_graph.hpp
  {
    ModuleVector module_vector;

    size_t num_modules = extract<size_t>( modules.attr("__len__")() );

    for (size_t i = 0; i < num_modules; ++i)
    {
      size_t num_vertices = extract<size_t>( modules[ i ].attr("__len__")() );

      Module module;

      for (size_t j = 0; j < num_vertices; ++j)
      {
        Vertex v = extract<Vertex>( modules[ i ][ j ] );
        module.push_back( v.v );
      }

      module_vector.push_back( module );
    }

    return module_vector;
  }

  double modularity(list modules) // FIXME: missing in python_directed_graph.hpp
  {
    ModuleVector module_vector( get_module_vector( modules ) );

    return conan::modularity( g, module_vector );
  }

  double z_score(Vertex a, list modules) // FIXME: missing in python_directed_graph.hpp
  {
    ModuleVector module_vector( get_module_vector( modules ) );

    return conan::z_score(g, module_vector, a.v);
  }

  double participation_coefficient(Vertex a, list modules) // FIXME: missing in python_directed_graph.hpp
  {
    ModuleVector module_vector( get_module_vector( modules ) );

    return conan::participation_coefficient(g, module_vector, a.v);
  }

  // Transformations
  void booleanize()
  {
    Graph g_copy(g);
    g = conan::booleanize_graph(g_copy);
    return;
  }

  void randomize_edges(conan::decimal p, bool maintain_connected = false)
  {
    Graph g_copy(g);
    g = conan::randomize_graph_edges(g_copy, p, maintain_connected);
    return;
  }

  python_undirected_graph complement()
  {
    return python_undirected_graph( conan::complement(g) );
  }

  python_undirected_graph clusterize(list clusters, bool include_all_vertices = false)
  {
    ModuleVector cluster_vector( get_module_vector( clusters ) );

    return python_undirected_graph( conan::clusterize(g, cluster_vector, include_all_vertices) );
  }

  // Edge set extractor functions
  list nearest_neighbour_ring(vertex start_vertex = 0)
  {
    typedef std::vector<edge> EdgeVector;

    EdgeVector edge_set;
    conan::nearest_neighbour_ring( g, start_vertex, edge_set );

    list edge_list;
    for (EdgeVector::iterator evi = edge_set.begin(); evi != edge_set.end(); ++evi)
    {
      edge_list.append( Edge( *evi, g ) );
    }

    return edge_list;
  }

  list minimum_spanning_tree()
  {
    typedef std::vector<edge> EdgeVector;

    EdgeVector edge_set;
    conan::minimum_spanning_tree(g, edge_set);

    list edge_list;
    for (EdgeVector::iterator evi = edge_set.begin(); evi != edge_set.end(); ++evi)
    {
      edge_list.append( Edge( *evi, g ) );
    }

    return edge_list;
  }
  
  // IO write functions
  void write_dot(std::string filename)
  { return conan::write_dot(g, filename); }

  void write_graphml(std::string filename)
  { return conan::write_graphml(g, filename); }

  void write_gml(std::string filename)
  { return conan::write_gml(g, filename); }

  void write_pajek(std::string filename)
  { return conan::write_pajek(g, filename); }

  void write_csv(std::string filename, std::string fields_order = "stw", char sep_char = '\t')
  { return conan::write_csv(g, filename, fields_order, sep_char); }

  void write_adj_matrix_file(std::string filename, char sep_char = ' ')
  { return conan::write_adj_matrix_file(g, filename); }

  // IO read functions
  void read_dot(std::string filename)
  { g = conan::read_dot<Graph>(filename); }

  void read_graphml(std::string filename)
  { g = conan::read_graphml<Graph>(filename); }

  void read_gml(std::string filename)
  { g = conan::read_gml<Graph>(filename); }

  void read_pajek(std::string filename)
  { g = conan::read_pajek<Graph>(filename); }

  void read_csv(std::string filename, std::string fields_order = "stw", char sep_char = '\t')
  { g = conan::read_csv<Graph>(filename, fields_order, sep_char); }

  void read_adj_matrix_file(std::string filename, char sep_char = ' ')
  { g = conan::read_adj_matrix_file<Graph>(filename, sep_char); }

  // Network properties
  conan::decimal asp()
  { return conan::graph_avg_shortest_path(g); }

  conan::decimal asp_normalized()
  { return conan::graph_avg_shortest_path_normalized(g); }

  conan::decimal degree_centrality()
  { return conan::graph_degree_centrality(g); }

  conan::decimal avg_clustering()
  { return conan::graph_avg_clustering(g); }

  conan::decimal avg_clustering_normalized()
  { return conan::graph_avg_clustering_normalized(g); }

  conan::decimal entropy()
  { return conan::graph_entropy(g); }

  conan::decimal connectance()
  { return conan::graph_connectance(g); }

  conan::decimal fractal_dimension(decimal start = 2.0, decimal end = 9.0, decimal step = 1.0, size_t num_iterations = 10)
  { return conan::graph_fractal_dimension(g, start, end, step, num_iterations); }

  python_undirected_graph boxcounting(decimal boxsize)
  { return python_undirected_graph(conan::detail::boxcounting(g, boxsize)); }

  conan::decimal distance(Vertex source, Vertex target)
  { return conan::distance(source.v, target.v, g); }

  Graph g;
};

#endif // PYTHON_UNDIRECTED_GRAPH_HPP
