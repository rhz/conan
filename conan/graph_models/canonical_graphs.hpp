/*
 * Conan - COmplex Network ANalisys
 * Copyright (C) 2008-2009  Ricardo Honorato Zimmer [rikardo.horo@gmail.com]
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef CANONICAL_GRAPHS_HPP
#define CANONICAL_GRAPHS_HPP
#include <conan/config.hpp>

namespace conan {

  // ***** Declarations *****
  
  /**
   * @brief Return a graph in which all vertices are linked to the centers of the star.
   *
   * @param V number of vertices of the graph.
   * @param num_centers number of center in the star. Default value = 1.
   * @return a star graph.
   */
  template <class Graph>
  Graph generate_star_graph(
      size_t V,
      size_t num_centers = 1
      );

  /**
   * @brief Return a graph in which all vertices have the same degree.
   *
   * @param V number of vertices.
   * @param degree the number of neighbours to which each vertex will be linked to.
   * @return a regular lattice graph.
   */
  template <class Graph>
  Graph generate_regular_lattice(
      size_t V,
      size_t degree
      );

  /**
   * @brief Return a graph where all vertices are connected to each other.
   *
   * @param V number of vertices.
   * @return a complete graph.
   */
  template <class Graph>
  Graph generate_complete_graph(
      size_t V
      );

  /**
   * @brief Return a path (graph in which each vertex is connected to the next one).
   * @param V number of vertices.
   */
  template <class Graph>
  Graph generate_path(
      size_t V
      );

  /**
   * @brief Return a cycle (graph in which each vertex is connected to the next one and the last vertex
   *   is connected to the first one).
   * @param V number of vertices
   */
  template <class Graph>
  Graph generate_cycle(
      size_t V
      );


  // ***** Definitions ******

  template <class Graph>
  Graph generate_star_graph(
      size_t V,
      size_t num_centers
      )
  { // O (V)
    typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor edge;
    typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iter;

    vertex first_vertices[num_centers];
    Graph g(V);
    bool on_first_vertices = true;
    size_t count = 0;
    vertex_iter vi, vi_end;
    for (tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi)
    {
      if (on_first_vertices)
      {
        first_vertices[count] = *vi;
        if (++count >= num_centers)
          on_first_vertices = false;
      }
      else
      {
        for (count = 0; count < num_centers; ++count)
        {
          edge e = boost::add_edge(*vi, first_vertices[count], g).first;
          g[e].weight = 1.0;
        }
      } // if on_first_vertices
    } // for vertex in vertices(g)
    return g;
  }


  template <class Graph>
  Graph generate_regular_lattice(
      size_t V,
      size_t degree
      )
  { // O (V * degree)
    // if degree < V this function returns an empty graph (ie. num_vertices(returned_graph) = 0)
    typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor edge;
    typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iter;
    
    if (degree > V)
      return Graph(0);

    bool is_degree_odd = false;
    if (degree % 2 != 0)
    { // maybe regular lattices with odd degree exists and we must take them into account
      is_degree_odd = true;
      throw std::runtime_error("a regular lattice can not have an odd degree");
    }
    size_t half_degree = degree / 2;

    Graph g(V);

    size_t vertex_index = 0;
    vertex first_ones[half_degree + (is_degree_odd? 1 : 0)];
    vertex_iter vi, vi_end;
    for (tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi, ++vertex_index)
    { // iterate over all vertices
      if (vertex_index < half_degree + (is_degree_odd? 1 : 0))
        first_ones[vertex_index] = *vi;

      for (size_t i = 1; i <= half_degree; ++i)
      { // connect the current vertex with the next 'half_degree' ones
        edge e;
        if (vertex_index < V - i)
          e = boost::add_edge(*vi, *(vi + i), g).first;
        else
          e = boost::add_edge(*vi, first_ones[vertex_index + i - V], g).first;
        g[e].weight = 1.0;
      }
    }
    return g;
  }


  template <class Graph>
  Graph generate_complete_graph(
      size_t V
      )
  {
    typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iterator;
    typedef typename boost::graph_traits<Graph>::edge_descriptor edge;

    Graph g(V);
    vertex_iterator vi, vi_end, dest_vi;
    for (tie(vi, vi_end) = boost::vertices(g); vi != vi_end - 1; ++vi)
    {
      for (dest_vi = vi + 1; dest_vi != vi_end; ++dest_vi)
      {
        edge e = boost::add_edge(*vi, *dest_vi, g).first;
        g[e].weight = 1.0;
      }
    }
    return g;
  }

  template <class Graph>
  Graph generate_path(
      size_t V
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_iterator vertex_iter;
    typedef typename GraphTraits::vertex_descriptor vertex;

    Graph out_g(V);
    vertex_iter vi, viend;
    tie(vi, viend) = boost::vertices(out_g);
    vertex last_vertex = *vi;
    for (++vi; vi != viend; ++vi)
    {
      boost::add_edge(last_vertex, *vi, 1.0, out_g);
      last_vertex = *vi;
    }

    return out_g;
  }

  template <class Graph>
  Graph generate_cycle(
      size_t V
      )
  { // This function is the same as generate_path, but connects the first and the last vertices, making a cycle
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_iterator vertex_iter;
    typedef typename GraphTraits::vertex_descriptor vertex;

    Graph out_g(V);
    vertex_iter vi, viend;
    tie(vi, viend) = boost::vertices(out_g);
    vertex last_vertex = *vi;
    for (++vi; vi != viend; ++vi)
    {
      boost::add_edge(last_vertex, *vi, 1.0, out_g);
      last_vertex = *vi;
    }
    boost::add_edge(last_vertex, 0, 1.0, out_g);

    return out_g;
  }

} //conan

#endif //CANONICAL_GRAPHS_HPP
