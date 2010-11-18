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

#ifndef BARABASI_ALBERT_HPP
#define BARABASI_ALBERT_HPP
#include <conan/config.hpp>
#include <conan/graph_models/canonical_graphs.hpp>
#include <conan/utils.hpp>


namespace conan {

  // ***** Declarations *****
  
  /**
   * @brief Return a graph created based on the Preferential attachment model, proposed by
   * Barabasi et al.
   *
   * @param V number of vertices
   * @param edges_by_vertex the number of edges created by each new vertex.
   * @return a scale-free graph
   */
  template <class Graph>
  Graph generate_scale_free_network(
      size_t V,
      int edges_by_vertex,
      Graph (*base_graph_generator) (size_t) = &generate_complete_graph
      );

  /**
   * @brief Expand a graph based on Preferential attachment model's rules.
   * (see conan::generate_scale_free_network).
   *
   * Parameters:
   * @param g base graph.
   * @param number_edges the number of edges to add to the base graph.
   * @param edges_by_vertex the number of edges created by each new vertex.
   */
  template <class Graph>
  void expand_scale_free_network(
      Graph & g,
      size_t number_vertices,
      int edges_by_vertex
      );


  // ***** Definitions ******

  template <class Graph>
  Graph generate_scale_free_network(
      size_t V,
      int edges_by_vertex,
      Graph (*base_graph_generator) (size_t) // = &generate_complete_graph
      )
  {
    if (V < 2)
      return Graph(1); // shortcut

    Graph g(base_graph_generator(edges_by_vertex));
    expand_scale_free_network(g, V - edges_by_vertex, edges_by_vertex);
    return g;
  }


  template <class Graph>
  void expand_scale_free_network (
      Graph & g,
      size_t number_vertices,
      int edges_by_vertex
      )
  {
    typedef typename Graph::edge_property_type Weight;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor edge;
    typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iter;

    size_t V = boost::num_vertices(g);

    if ( edges_by_vertex <= 0 )
      throw std::runtime_error("conan::expand_scale_free_network: edges_by_vertex must be a positive integer.");
    else if ( V < size_t(edges_by_vertex) )
      throw std::runtime_error("conan::expand_scale_free_network: According to the Barabasi-Albert model, m must be equal or greater than m_0.");

    for (size_t i = 0; i < number_vertices; ++i)
    {
      vertex v = boost::add_vertex(g); ++V;

      std::vector<vertex> vertex_vector(V);
      std::vector<size_t> degree_vector(V);
      
      // Copy graph vertices to vertex_vector
      vertex_iter vi, vi_end;
      tie(vi, vi_end) = boost::vertices(g);
      std::copy(vi, vi_end, vertex_vector.begin());  //FIXME: can a prior vertex_vector elem appear in the next iteration?
      
      if (V == 1)
      {
        continue;
      }
      if (V == 2)
      {
        boost::add_edge(vertex_vector[ 0 ], vertex_vector[ 1 ], 1.0, g);
        continue;
      }
      
      for (int j = 0; j < edges_by_vertex; ++j)
      {
        size_t sum = 0, vertex_index = 0;

        vertex dest_vertex;
        
        // (re)make degree_vector
        for (tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi, ++vertex_index)
          sum += degree_vector[vertex_index] = boost::out_degree(*vi, g);
        
        // select any dest_vertex different from the current one
        while ( (dest_vertex = vertex_vector[ detail::weighted_choice(degree_vector, sum) ]) == v )
          { }

        if ( ! boost::edge(v, dest_vertex, g).second ) // check if an edge between v and dest_vertex already exists
        {
          edge e = boost::add_edge(v, dest_vertex, g).first;
          g[e].weight = 1.0;
        }
        else
        {
          --j;
        }
      } // for j
    } // for i
  }

} //conan

#endif //BARABASI_ALBERT_HPP
