#ifndef PYTHON_GRAPH_MODELS_HPP
#define PYTHON_GRAPH_MODELS_HPP

#include "python_undirected_graph.hpp"
#include <conan/graph_models.hpp>

// Generator functions
python_undirected_graph generate_random_graph(
    size_t num_vertices,
    size_t num_edges)
{
  return python_undirected_graph(
      conan::generate_random_graph<python_undirected_graph::Graph>(num_vertices, num_edges));
}

python_undirected_graph generate_erdos_renyi_graph(
    size_t num_vertices,
    conan::decimal p)
{
  return python_undirected_graph(
      conan::generate_erdos_renyi_graph<python_undirected_graph::Graph>(num_vertices, p));
}

python_undirected_graph generate_star_graph(
    size_t num_vertices,
    size_t num_centers)
{
  return python_undirected_graph(
      conan::generate_star_graph<python_undirected_graph::Graph>(num_vertices, num_centers));
}

python_undirected_graph generate_regular_lattice(
    size_t num_vertices,
    size_t degree)
{
  return python_undirected_graph(
      conan::generate_regular_lattice<python_undirected_graph::Graph>(num_vertices, degree));
}

python_undirected_graph generate_complete_graph(
    size_t num_vertices)
{
  return python_undirected_graph(
      conan::generate_complete_graph<python_undirected_graph::Graph>(num_vertices));
}

python_undirected_graph generate_scale_free_network(
    size_t num_vertices,
    size_t edges_by_vertex
    )
{
  // how to pass a custom function to generate the base graph?
  // or alternatively, how to pass arguments to conan::generate_scale_free_network in order
  // to call a base graph generator function that recieves more arguments than just num_vertices?
  return python_undirected_graph(
      conan::generate_scale_free_network<python_undirected_graph::Graph>(num_vertices, edges_by_vertex));
}

void expand_scale_free_network(
    python_undirected_graph& g,
    size_t num_vertices,
    size_t edges_by_vertex)
{
  return conan::expand_scale_free_network(g.g, num_vertices, edges_by_vertex);
}

python_undirected_graph generate_Demetrius_Manke_network(
    size_t num_vertices,
    conan::decimal T,
    size_t edges_by_vertex)
{
  return python_undirected_graph(
      conan::generate_Demetrius_Manke_network<python_undirected_graph::Graph>(num_vertices, T, edges_by_vertex));
}

void add_vertex_as_Demetrius_Manke(
    python_undirected_graph& g,
    conan::decimal T,
    size_t edges_by_vertex)
{
  return conan::add_vertex_as_Demetrius_Manke(g.g, T, edges_by_vertex);
}

python_undirected_graph generate_small_world_network(
    size_t num_vertices,
    size_t degree,
    conan::decimal percent_replaced_edges)
{
  return python_undirected_graph(
      conan::generate_small_world_network<python_undirected_graph::Graph>(num_vertices, degree, percent_replaced_edges));
}

python_undirected_graph generate_path(
    size_t num_vertices)
{
  return python_undirected_graph(
      conan::generate_path<python_undirected_graph::Graph>(num_vertices));
}

python_undirected_graph generate_cycle(
    size_t num_vertices)
{
  return python_undirected_graph(
      conan::generate_cycle<python_undirected_graph::Graph>(num_vertices));
}


BOOST_PYTHON_MODULE(_graph_models)
{
  def("generate_random_graph", generate_random_graph);
  def("generate_erdos_renyi_graph", generate_erdos_renyi_graph);
  def("generate_star_graph", generate_star_graph);
  def("generate_regular_lattice", generate_regular_lattice);
  def("generate_complete_graph", generate_complete_graph);
  def("generate_small_world_network", generate_small_world_network);
  def("generate_Demetrius_Manke_network", generate_Demetrius_Manke_network);
  def("add_vertex_as_Demetrius_Manke", add_vertex_as_Demetrius_Manke);
  def("generate_scale_free_network", generate_scale_free_network);
  def("expand_scale_free_network", expand_scale_free_network);
  def("generate_path", generate_path);
  def("generate_cycle", generate_cycle);
}

#endif //PYTHON_GRAPH_MODELS_HPP
