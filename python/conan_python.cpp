#include "python_undirected_graph.hpp"
#include "python_directed_graph.hpp"
#include "misc_classes.hpp"
#include "debug.hpp"
#include <ctime>

#include <boost/python/tuple.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/raw_function.hpp>


/**
 * This file contains the '_conan' module interface definition.
 * It also defines all function overloads needed by python_undirected_graph and python_directed_graph and
 * the deep copy functions for them.
 */

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(read_csv_overloads, read_csv, 1, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(write_csv_overloads, write_csv, 1, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(read_adj_matrix_file_overloads, read_adj_matrix_file, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(write_adj_matrix_file_overloads, write_adj_matrix_file, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(fractal_dimension_overloads, fractal_dimension, 0, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(randomize_edges_overloads, randomize_edges, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(nearest_neighbour_ring_overloads, nearest_neighbour_ring, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(newman_communities_overloads, newman_communities, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(guimera_amaral_communities_overloads, guimera_amaral_communities, 0, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(clusterize_overloads, clusterize, 1, 2)


// add_vertex overloads
python_undirected_graph::Vertex (python_undirected_graph:: *undirected_graph_add_vertex1)()
  = &python_undirected_graph::add_vertex;
python_undirected_graph::Vertex (python_undirected_graph:: *undirected_graph_add_vertex2)(std::string)
  = &python_undirected_graph::add_vertex;

python_directed_graph::Vertex (python_directed_graph:: *directed_graph_add_vertex1)()
  = &python_directed_graph::add_vertex;
python_directed_graph::Vertex (python_directed_graph:: *directed_graph_add_vertex2)(std::string)
  = &python_directed_graph::add_vertex;


// add_edge overloads
python_undirected_graph::Edge (python_undirected_graph:: *undirected_graph_add_edge1)(python_undirected_graph::Vertex, python_undirected_graph::Vertex)
  = &python_undirected_graph::add_edge;
python_undirected_graph::Edge (python_undirected_graph:: *undirected_graph_add_edge2)(python_undirected_graph::vertex, python_undirected_graph::vertex)
  = &python_undirected_graph::add_edge;
//python_undirected_graph::Edge (python_undirected_graph:: *undirected_graph_add_edge2)(python_undirected_graph::Edge)
//  = &python_undirected_graph::add_edge;

python_directed_graph::Edge (python_directed_graph:: *directed_graph_add_edge1)(python_directed_graph::Vertex, python_directed_graph::Vertex)
  = &python_directed_graph::add_edge;
python_directed_graph::Edge (python_directed_graph:: *directed_graph_add_edge2)(python_directed_graph::vertex, python_directed_graph::vertex)
  = &python_directed_graph::add_edge;
//python_directed_graph::Edge (python_directed_graph:: *directed_graph_add_edge2)(python_directed_graph::Edge)
//  = &python_directed_graph::add_edge;


// remove_vertex overloads
void (python_undirected_graph:: *undirected_graph_remove_vertex1)(python_undirected_graph::Vertex)
  = &python_undirected_graph::remove_vertex;
void (python_undirected_graph:: *undirected_graph_remove_vertex2)(python_undirected_graph::vertex)
  = &python_undirected_graph::remove_vertex;

void (python_directed_graph:: *directed_graph_remove_vertex1)(python_directed_graph::Vertex)
  = &python_directed_graph::remove_vertex;
void (python_directed_graph:: *directed_graph_remove_vertex2)(python_directed_graph::vertex)
  = &python_directed_graph::remove_vertex;


// remove_edge overloads
void (python_undirected_graph:: *undirected_graph_remove_edge1)(python_undirected_graph::Edge)
  = &python_undirected_graph::remove_edge;
void (python_undirected_graph:: *undirected_graph_remove_edge2)(python_undirected_graph::Vertex, python_undirected_graph::Vertex)
  = &python_undirected_graph::remove_edge;
void (python_undirected_graph:: *undirected_graph_remove_edge3)(python_undirected_graph::vertex, python_undirected_graph::vertex)
  = &python_undirected_graph::remove_edge;

void (python_directed_graph:: *directed_graph_remove_edge1)(python_directed_graph::Edge)
  = &python_directed_graph::remove_edge;
void (python_directed_graph:: *directed_graph_remove_edge2)(python_directed_graph::Vertex, python_directed_graph::Vertex)
  = &python_directed_graph::remove_edge;
void (python_directed_graph:: *directed_graph_remove_edge3)(python_directed_graph::vertex, python_directed_graph::vertex)
  = &python_directed_graph::remove_edge;


// has_edge overloads
bool (python_undirected_graph:: *undirected_graph_has_edge1) (python_undirected_graph::Edge)
  = &python_undirected_graph::has_edge;
bool (python_undirected_graph:: *undirected_graph_has_edge2) (python_undirected_graph::Vertex, python_undirected_graph::Vertex)
  = &python_undirected_graph::has_edge;
bool (python_undirected_graph:: *undirected_graph_has_edge3) (python_undirected_graph::vertex, python_undirected_graph::vertex)
  = &python_undirected_graph::has_edge;


// get_vertex overloads
python_undirected_graph::Vertex (python_undirected_graph:: *undirected_graph_get_vertex1)(size_t)
  = &python_undirected_graph::get_vertex;
python_undirected_graph::Vertex (python_undirected_graph:: *undirected_graph_get_vertex2)(std::string)
  = &python_undirected_graph::get_vertex;

python_directed_graph::Vertex (python_directed_graph:: *directed_graph_get_vertex1)(size_t)
  = &python_directed_graph::get_vertex;
python_directed_graph::Vertex (python_directed_graph:: *directed_graph_get_vertex2)(std::string)
  = &python_directed_graph::get_vertex;


// get_edge overloads
python_undirected_graph::Edge (python_undirected_graph:: *undirected_graph_get_edge1)(size_t, size_t)
  = &python_undirected_graph::get_edge;
python_undirected_graph::Edge (python_undirected_graph:: *undirected_graph_get_edge2)(python_undirected_graph::Vertex, python_undirected_graph::Vertex)
  = &python_undirected_graph::get_edge;

python_directed_graph::Edge (python_directed_graph:: *directed_graph_get_edge1)(size_t, size_t)
  = &python_directed_graph::get_edge;
python_directed_graph::Edge (python_directed_graph:: *directed_graph_get_edge2)(python_directed_graph::Vertex, python_directed_graph::Vertex)
  = &python_directed_graph::get_edge;


// merge_graphs
python_undirected_graph merge_undirected_graphs(
    tuple args,
    dict kw
    )
{
  std::list<python_undirected_graph::Graph> graph_list;
  size_t l = len(args);
  for (size_t i = 0; i < l; ++i)
  {
    python_undirected_graph g = extract<python_undirected_graph>( args[i] );
    graph_list.push_back( g.g );
  }

  bool overlap = false;
  try
  { 
    extract<bool> get_overlap(kw["overlap"]);
    if (get_overlap.check())
      overlap = get_overlap();
  }
  catch (...) { }

  return python_undirected_graph(
      conan::merge_graphs(graph_list, overlap) );
}

//BOOST_PYTHON_FUNCTION_OVERLOADS(merge_undirected_graphs_overloads, merge_undirected_graphs, 1, 2)


python_directed_graph merge_directed_graphs(
    tuple args,
    dict kw
    )
{
  std::list<python_directed_graph::Graph> graph_list;
  size_t l = len(args);
  for (size_t i = 0; i < l; ++i)
  {
    python_directed_graph g = extract<python_directed_graph>( args[i] );
    graph_list.push_back( g.g );
  }

  bool overlap = false;
  try
  { 
    extract<bool> get_overlap(kw["overlap"]);
    if (get_overlap.check())
      overlap = get_overlap();
  }
  catch (...) { }

  return python_directed_graph(
      conan::merge_graphs(graph_list, overlap) );
}

//BOOST_PYTHON_FUNCTION_OVERLOADS(merge_directed_graphs_overloads, merge_directed_graphs, 1, 2)


// deep copy
python_undirected_graph undirected_graph_deep_copy(
    const python_undirected_graph& x
    )
{
  python_undirected_graph x_copy(x);
  return x_copy;
}

python_directed_graph directed_graph_deep_copy(
    const python_directed_graph& x
    )
{
  python_directed_graph x_copy(x);
  return x_copy;
}


BOOST_PYTHON_MODULE(_conan)
{
  docstring_options doc_options;
  //doc_options.disable_py_signatures();
  doc_options.disable_cpp_signatures();

  // Free functions
  def("version", version);
  def("compilation_info", compilation_info);
  def("random_graph_average_degree", conan::random_graph_average_degree);
  def("random_graph_avg_shortest_path", conan::random_graph_avg_shortest_path);
  def("random_graph_avg_clustering", conan::random_graph_avg_clustering);

  def("undirected_graph_deep_copy", undirected_graph_deep_copy);
  def("directed_graph_deep_copy", directed_graph_deep_copy);

  def("merge_undirected_graphs", raw_function(merge_undirected_graphs, 2));
  def("merge_directed_graphs", raw_function(merge_directed_graphs, 2));

  // Classes
  // Cpp
  class_< python_undirected_graph::Graph >("cpp_graph", no_init)
    ;
  class_< python_undirected_graph::edge >("cpp_edge_descriptor")
    ;
  class_< std::list<python_undirected_graph::edge> >("cpp_edge_list")
    .def("__iter__", iterator< std::list<python_undirected_graph::edge> >())
    ;
  class_< std::list<python_undirected_graph::vertex> >("cpp_vertex_list")
    .def("__iter__", iterator< std::list<python_undirected_graph::vertex> >())
    ;
  class_< std::vector<python_undirected_graph::vertex> >("cpp_vertex_vector")
    .def("__iter__", iterator< std::vector<python_undirected_graph::vertex> >())
    ;
  class_< std::vector<conan::decimal> >("cpp_decimal_vector")
    .def("__iter__", iterator< std::vector<conan::decimal> >())
    ;
  class_< std::list<python_undirected_graph> >("cpp_graph_list")
    .def("__iter__", iterator< std::list<python_undirected_graph> >())
    ;

  // Python

  // Undirected graph vertex and edge wraps
  class_< python_undirected_graph::Vertex >("undirected_vertex", no_init)
    .def("_id", &python_undirected_graph::Vertex::get_id)
    .add_property("name", &python_undirected_graph::Vertex::get_name, &python_undirected_graph::Vertex::set_name)
    .def("degree", &python_undirected_graph::Vertex::get_degree)
    .def("clustering", &python_undirected_graph::Vertex::get_clustering)
    .def("asp", &python_undirected_graph::Vertex::get_asp)
    .def("degree_centrality", &python_undirected_graph::Vertex::get_degree_centrality)
    .def("closeness", &python_undirected_graph::Vertex::get_closeness)
    .def("incident_edges", &python_undirected_graph::Vertex::get_incident_edges)
    .def("adjacent_vertices", &python_undirected_graph::Vertex::get_adjacent_vertices)
    .def("__str__", &python_undirected_graph::Vertex::get_str)
    .def("__repr__", &python_undirected_graph::Vertex::get_str)
    ;
  class_< python_undirected_graph::Edge >("undirected_edge", no_init)
    .def("_id", &python_undirected_graph::Edge::get_id)
    .add_property("weight", &python_undirected_graph::Edge::get_weight, &python_undirected_graph::Edge::set_weight)
    .add_property("type", &python_undirected_graph::Edge::get_type, &python_undirected_graph::Edge::set_type)
    .def("source_vertex", &python_undirected_graph::Edge::get_source_vertex)
    .def("target_vertex", &python_undirected_graph::Edge::get_target_vertex)
    .def("__str__", &python_undirected_graph::Edge::get_str)
    .def("__repr__", &python_undirected_graph::Edge::get_str)
    ;

  // Undirected graph
  class_< python_undirected_graph >("undirected_graph")
    // Ctors.
    .def(init<size_t>())
    .def(init<python_undirected_graph::Graph>())
    .def(init<list>())
    // Base methods
    .def("vertices", &python_undirected_graph::vertices, "Returns the set of vertices of the graph.\n")
    .def("edges", &python_undirected_graph::edges, "Returns the set of edges of the graph.\n")
    .def("num_vertices", &python_undirected_graph::num_vertices, "Returns the total number of vertices in the graph.\n")
    .def("num_edges", &python_undirected_graph::num_edges, "Returns the total number of edges in the graph.\n")
    .def("add_vertex", undirected_graph_add_vertex1, "Add an unnamed vertex to the graph.\n")
    .def("add_vertex", undirected_graph_add_vertex2, "Add a named vertex to the graph.\n")
    .def("add_edge", undirected_graph_add_edge1, "Add an edge to the graph.\n")
    .def("add_edge", undirected_graph_add_edge2, "Add an edge to the graph.\n")
    .def("remove_vertex", undirected_graph_remove_vertex1, "Remove a vertex from the graph.\n")
    .def("remove_vertex", undirected_graph_remove_vertex2, "Remove a vertex from the graph.\n")
    .def("remove_edge", undirected_graph_remove_edge1, "Remove an edge from the graph.\n")
    .def("remove_edge", undirected_graph_remove_edge2, "Remove an edge from the graph.\n")
    .def("remove_edge", undirected_graph_remove_edge3, "Remove an edge from the graph.\n")
    .def("get_vertex", undirected_graph_get_vertex1, "Get a vertex object by giving the vertex id.\n")
    .def("get_vertex", undirected_graph_get_vertex2, "Get a vertex object by giving the vertex id.\n")
    .def("get_edge", undirected_graph_get_edge1, "Get an edge object by giving the id of the vertices it connects.\n")
    .def("get_edge", undirected_graph_get_edge2, "Get an edge object by giving the id of the vertices it connects.\n")
    .def("get_adj_matrix", &python_undirected_graph::get_adj_matrix, "Get the adjancency matrix in the form of a list of lists.\n")
    .def("has_edge", undirected_graph_has_edge1, "Check if an edge exists.\n")
    .def("has_edge", undirected_graph_has_edge2, "Check if an edge exists.\n")
    .def("has_edge", undirected_graph_has_edge3, "Check if an edge exists.\n")
    // IO
    .def("write_dot", &python_undirected_graph::write_dot, "Write an output DOT file to visualize the graph with GraphViz.\n")
    .def("write_graphml", &python_undirected_graph::write_graphml, "Write the graph in GraphML format.\n")
    .def("write_gml", &python_undirected_graph::write_gml, "Write the graph in GML format.\n")
    .def("write_pajek", &python_undirected_graph::write_pajek, "Write the graph in Pajek format.\n")
    .def("write_csv", &python_undirected_graph::write_csv,
        write_csv_overloads(args("filename", "print_edge_weight", "print_edge_type"),
          "Write the graph in a CSV-like format.\n"))
    .def("write_adj_matrix_file", &python_undirected_graph::write_adj_matrix_file,
        write_adj_matrix_file_overloads(args("filename", "sep_char"),
          "Write the graph in a raw adjacency matrix file.\n"))
    .def("read_dot", &python_undirected_graph::read_dot)
    .def("read_graphml", &python_undirected_graph::read_graphml)
    .def("read_gml", &python_undirected_graph::read_gml)
    .def("read_pajek", &python_undirected_graph::read_pajek)
    .def("read_csv", &python_undirected_graph::read_csv,
        read_csv_overloads(args("filename", "fields_order", "sep_char"),
          "Read a graph from a CSV-like formatted file.\n"))
    .def("read_adj_matrix_file", &python_undirected_graph::read_adj_matrix_file, read_adj_matrix_file_overloads())
    // Properties
    .def("asp", &python_undirected_graph::asp)
    .def("asp_normalized", &python_undirected_graph::asp_normalized)
    .def("degree_centrality", &python_undirected_graph::degree_centrality)
    .def("avg_clustering", &python_undirected_graph::avg_clustering)
    .def("avg_clustering_normalized", &python_undirected_graph::avg_clustering_normalized)
    .def("connectance", &python_undirected_graph::connectance)
    .def("entropy", &python_undirected_graph::entropy)
    .def("fractal_dimension", &python_undirected_graph::fractal_dimension, fractal_dimension_overloads())
    .def("distance", &python_undirected_graph::distance)
    // Division
    .def("subgraph", &python_undirected_graph::subgraph)
    .def("components", &python_undirected_graph::components)
    .def("newman_communities", &python_undirected_graph::newman_communities,
        newman_communities_overloads(args("check_quality_function"),
          "Find communities within the graph using the algorithm designed by Newman (ref.)\n"))
    .def("newman_two_communities", &python_undirected_graph::newman_two_communities)
    .def("guimera_amaral_communities", &python_undirected_graph::guimera_amaral_communities,
        guimera_amaral_communities_overloads(args("initial_temp", "final_temp", "f"),
          "Find communities within the graph using the algorithm designed by Guimera and Amaral (ref.)\n"))
    .def("modularity", &python_undirected_graph::modularity)
    .def("z_score", &python_undirected_graph::z_score)
    .def("participation_coefficient", &python_undirected_graph::participation_coefficient)
    // Transformations and generator functions
    .def("booleanize", &python_undirected_graph::booleanize)
    .def("randomize_edges", &python_undirected_graph::randomize_edges, randomize_edges_overloads())
    .def("complement", &python_undirected_graph::complement)
    .def("clusterize", &python_undirected_graph::clusterize,
        clusterize_overloads(args("clusters", "include_all_vertices"),
          "Return a graph where the clusters are nodes.\n"))
    .def("minimum_spanning_tree", &python_undirected_graph::minimum_spanning_tree)
    .def("nearest_neighbour_ring", &python_undirected_graph::nearest_neighbour_ring, nearest_neighbour_ring_overloads())
    .def("boxcounting", &python_undirected_graph::boxcounting)
    ;
  typedef degree_distribution< python_undirected_graph > undirected_degree_distribution;
  class_< undirected_degree_distribution >("undirected_degree_distribution", init<python_undirected_graph &>())
    .def(init<python_undirected_graph &, bool>())
    .def("best_fit", &undirected_degree_distribution::best_fit)
    .def("P", &undirected_degree_distribution::P)
    .def("vertices", &undirected_degree_distribution::vertices)
    .def("non_zero_degrees", &undirected_degree_distribution::non_zero_degrees)
    .def("gamma", &undirected_degree_distribution::gamma)
    .def("exp_decay", &undirected_degree_distribution::lambda)
    .def("slope", &undirected_degree_distribution::slope)
    .def("r", &undirected_degree_distribution::r)
    .def("rr", &undirected_degree_distribution::rr)
    .def("entropy", &undirected_degree_distribution::entropy)
    ;
  typedef shortest_paths< python_undirected_graph > undirected_shortest_paths;
  class_< undirected_shortest_paths >("undirected_shortest_paths", init<python_undirected_graph &>())
    .def("vertex_asp", &undirected_shortest_paths::vertex_asp)
    .def("graph_asp", &undirected_shortest_paths::graph_asp)
    .def("vertex_eccentricity", &undirected_shortest_paths::vertex_eccentricity)
    .def("graph_diameter", &undirected_shortest_paths::graph_diameter)
    .def("graph_radius", &undirected_shortest_paths::graph_radius)
    .def("graph_center", &undirected_shortest_paths::graph_center)
    ;
  typedef betweenness_centrality< python_undirected_graph > undirected_betweenness_centrality;
  class_< undirected_betweenness_centrality >("undirected_betweenness_centrality", init<python_undirected_graph &>())
    .def("vertex_centrality", &undirected_betweenness_centrality::vertex_centrality)
    .def("vertex_relative_centrality", &undirected_betweenness_centrality::vertex_relative_centrality)
    .def("edge_centrality", &undirected_betweenness_centrality::edge_centrality)
    .def("central_point_dominance", &undirected_betweenness_centrality::central_point_dominance)
    ;
  typedef eigenvector_centrality< python_undirected_graph > undirected_eigenvector_centrality;
  class_< undirected_eigenvector_centrality >("undirected_eigenvector_centrality", init<python_undirected_graph &>())
    .def("vertex_eigenvector_centrality", &undirected_eigenvector_centrality::vertex_eigenvector_centrality)
    ;

  // Directed graph vertex and edge wraps
  class_< python_directed_graph::Vertex >("directed_vertex", no_init)
    .def("_id", &python_directed_graph::Vertex::get_id)
    .add_property("name", &python_directed_graph::Vertex::get_name, &python_directed_graph::Vertex::set_name)
    .def("degree", &python_directed_graph::Vertex::get_degree)
    .def("clustering", &python_directed_graph::Vertex::get_clustering)
    .def("asp", &python_directed_graph::Vertex::get_asp)
    .def("degree_centrality", &python_directed_graph::Vertex::get_degree_centrality)
    .def("closeness", &python_directed_graph::Vertex::get_closeness)
    .def("incident_edges", &python_directed_graph::Vertex::get_incident_edges)
    .def("adjacent_vertices", &python_directed_graph::Vertex::get_adjacent_vertices)
    .def("__str__", &python_directed_graph::Vertex::get_str)
    .def("__repr__", &python_directed_graph::Vertex::get_str)
    ;
  class_< python_directed_graph::Edge >("directed_edge", no_init)
    .def("_id", &python_directed_graph::Edge::get_id)
    .add_property("weight", &python_directed_graph::Edge::get_weight, &python_directed_graph::Edge::set_weight)
    .add_property("type", &python_directed_graph::Edge::get_type, &python_directed_graph::Edge::set_type)
    .def("source_vertex", &python_directed_graph::Edge::get_source_vertex)
    .def("target_vertex", &python_directed_graph::Edge::get_target_vertex)
    .def("__str__", &python_directed_graph::Edge::get_str)
    .def("__repr__", &python_directed_graph::Edge::get_str)
    ;

  // Directed graph
  class_< python_directed_graph >("directed_graph")
    .def(init<size_t>())
    .def(init<python_directed_graph::Graph>())
    .def(init<list>())
    // Base methods
    .def("vertices", &python_directed_graph::vertices, "Returns the set of vertices of the graph.\n")
    .def("edges", &python_directed_graph::edges, "Returns the set of edges of the graph.\n")
    .def("num_vertices", &python_directed_graph::num_vertices, "Returns the total number of vertices in the graph.\n")
    .def("num_edges", &python_directed_graph::num_edges, "Returns the total number of edges in the graph.\n")
    .def("add_vertex", directed_graph_add_vertex1, "Add an unnamed vertex to the graph.\n")
    .def("add_vertex", directed_graph_add_vertex2, "Add a named vertex to the graph.\n")
    .def("add_edge", directed_graph_add_edge1, "Add an edge to the graph.\n")
    .def("add_edge", directed_graph_add_edge2, "Add an edge to the graph.\n")
    .def("remove_vertex", directed_graph_remove_vertex1, "Remove a vertex from the graph.\n")
    .def("remove_vertex", directed_graph_remove_vertex2, "Remove a vertex from the graph.\n")
    .def("remove_edge", directed_graph_remove_edge1, "Remove an edge from the graph.\n")
    .def("remove_edge", directed_graph_remove_edge2, "Remove an edge from the graph.\n")
    .def("remove_edge", directed_graph_remove_edge3, "Remove an edge from the graph.\n")
    .def("get_vertex", directed_graph_get_vertex1, "Get a vertex object by giving the vertex id.\n")
    .def("get_vertex", directed_graph_get_vertex2, "Get a vertex object by giving the vertex id.\n")
    .def("get_edge", directed_graph_get_edge1, "Get an edge object by giving the id of the vertices it connects.\n")
    .def("get_edge", directed_graph_get_edge2, "Get an edge object by giving the id of the vertices it connects.\n")
    .def("get_adj_matrix", &python_directed_graph::get_adj_matrix, "Get the adjancency matrix in the form of a list of lists.\n")
    // IO
    .def("write_dot", &python_directed_graph::write_dot, "Write an output DOT file to visualize the graph with GraphViz.\n")
    .def("write_graphml", &python_directed_graph::write_graphml, "Write the graph in GraphML format.\n")
    .def("write_gml", &python_directed_graph::write_gml, "Write the graph in GML format.\n")
    .def("write_pajek", &python_directed_graph::write_pajek, "Write the graph in Pajek format.\n")
    .def("write_csv", &python_directed_graph::write_csv,
        write_csv_overloads(args("filename", "print_edge_weight", "print_edge_type"),
          "Write the graph in a CSV-like format.\n"))
    .def("write_adj_matrix_file", &python_directed_graph::write_adj_matrix_file, "Write the graph in a raw adjacency matrix file.\n")
    .def("read_dot", &python_directed_graph::read_dot)
    .def("read_graphml", &python_directed_graph::read_graphml)
    .def("read_gml", &python_directed_graph::read_gml)
    .def("read_pajek", &python_directed_graph::read_pajek)
    .def("read_csv", &python_directed_graph::read_csv,
        read_csv_overloads(args("filename", "fields_order", "sep_char"),
          "Read a graph from a CSV-like formatted file.\n"))
    .def("read_adj_matrix_file", &python_directed_graph::read_adj_matrix_file, read_adj_matrix_file_overloads())
    // Properties
    .def("asp", &python_directed_graph::asp)
    .def("asp_normalized", &python_directed_graph::asp_normalized)
    .def("degree_centrality", &python_directed_graph::degree_centrality)
    .def("avg_clustering", &python_directed_graph::avg_clustering)
    .def("avg_clustering_normalized", &python_directed_graph::avg_clustering_normalized)
    .def("connectance", &python_directed_graph::connectance)
    .def("fractal_dimension", &python_directed_graph::fractal_dimension, fractal_dimension_overloads())
    .def("distance", &python_directed_graph::distance)
    // Division
    .def("subgraph", &python_directed_graph::subgraph)
    .def("newman_communities", &python_directed_graph::newman_communities,
        newman_communities_overloads(args("check_quality_function"),
          "Find communities within the graph using the algorithm disegned by Newman (ref.)\n"))
    .def("newman_two_communities", &python_directed_graph::newman_two_communities)
    .def("modularity", &python_directed_graph::modularity)
    //.def("components", &python_directed_graph::components)
    // Transformations and generator functions
    .def("booleanize", &python_directed_graph::booleanize)
    .def("transpose", &python_directed_graph::transpose)
    .def("boxcounting", &python_directed_graph::boxcounting)
    //.def("randomize_edges", &python_directed_graph::randomize_edges, randomize_edges_overloads()) // this function uses connected_components
    ;
  typedef degree_distribution< python_directed_graph > directed_degree_distribution;
  class_< directed_degree_distribution >("directed_degree_distribution", init<python_directed_graph &>())
    .def(init<python_directed_graph &, bool>())
    .def("best_fit", &directed_degree_distribution::best_fit)
    .def("P", &directed_degree_distribution::P)
    .def("vertices", &directed_degree_distribution::vertices)
    .def("non_zero_degrees", &directed_degree_distribution::non_zero_degrees)
    .def("gamma", &directed_degree_distribution::gamma)
    .def("lambda", &directed_degree_distribution::lambda)
    .def("slope", &directed_degree_distribution::slope)
    .def("r", &directed_degree_distribution::r)
    .def("rr", &directed_degree_distribution::rr)
    .def("entropy", &directed_degree_distribution::entropy)
    ;
  typedef shortest_paths< python_directed_graph > directed_shortest_paths;
  class_< directed_shortest_paths >("directed_shortest_paths", init<python_directed_graph &>())
    .def("vertex_asp", &directed_shortest_paths::vertex_asp)
    .def("graph_asp", &directed_shortest_paths::graph_asp)
    .def("vertex_eccentricity", &directed_shortest_paths::vertex_eccentricity)
    .def("graph_diameter", &directed_shortest_paths::graph_diameter)
    .def("graph_radius", &directed_shortest_paths::graph_radius)
    .def("graph_center", &directed_shortest_paths::graph_center)
    ;
  typedef betweenness_centrality< python_directed_graph > directed_betweenness_centrality;
  class_< directed_betweenness_centrality >("directed_betweenness_centrality", init<python_directed_graph &>())
    .def("vertex_centrality", &directed_betweenness_centrality::vertex_centrality)
    .def("vertex_relative_centrality", &directed_betweenness_centrality::vertex_relative_centrality)
    .def("edge_centrality", &directed_betweenness_centrality::edge_centrality)
    .def("central_point_dominance", &directed_betweenness_centrality::central_point_dominance)
    ;
  typedef eigenvector_centrality< python_directed_graph > directed_eigenvector_centrality;
  class_< directed_eigenvector_centrality >("directed_eigenvector_centrality", init<python_directed_graph &>())
    .def("vertex_eigenvector_centrality", &directed_eigenvector_centrality::vertex_eigenvector_centrality)
    ;
}
