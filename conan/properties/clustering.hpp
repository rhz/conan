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

#ifndef CLUSTERING_HPP
#define CLUSTERING_HPP
#include <conan/config.hpp>
#include <conan/graphs.hpp>
#include <conan/utils/combinatoria.hpp>
#include <conan/properties/shortest_paths.hpp>

namespace conan {

  /**
   * @brief Calculate the clustering coefficient of a vertex.
   *
   * @param g a graph.
   * @param v vertex of g.
   * @return the vertex's clustering coefficient value
   */
  template <class Graph>
  decimal vertex_clustering(
      typename boost::graph_traits<Graph>::vertex_descriptor v,
      const Graph & g
      );
  
  /**
   * @brief Return the clustering coefficient value for the graph.
   * @param g a graph.
   */
  template <class Graph>
  decimal graph_avg_clustering(
      const Graph & g
      );
  
  /**
   * @brief Return the clustering coefficient value for the graph normalized by the value for
   *   a random graph with the same number of vertices and edges.
   * @param g a graph.
   */
  template <class Graph>
  decimal graph_avg_clustering_normalized(
      const Graph & g
      );


  template <class Graph>
  inline decimal graph_avg_clustering_normalized(
      const Graph & g
      )
  {
    return graph_avg_clustering(g) / random_graph_avg_clustering(boost::num_vertices(g), boost::num_edges(g));
  }
  
  
  template <class Graph>
  decimal vertex_clustering(
      typename boost::graph_traits<Graph>::vertex_descriptor v,
      const Graph & g
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::adjacency_iterator adjacency_iter;
    typedef typename GraphTraits::vertex_descriptor vertex;
    
    // calculate total_triagles and num_adj_vertices
    size_t num_adj_vertices = boost::out_degree(v, g);
    if (num_adj_vertices < 2)
      return 0.0;  // by definition

    size_t total_triangles = (num_adj_vertices * (num_adj_vertices - 1)) / 2;

    // make an array with the adjacent vertices to vertex
    vertex adj_vertices[num_adj_vertices];
    adjacency_iter vi, vi_end; size_t i = 0;
    for (tie(vi, vi_end) = boost::adjacent_vertices(v, g); vi != vi_end; ++vi)
      adj_vertices[i++] = *vi;

    size_t existent_triangles = 0;

    for (i = 0; i < num_adj_vertices; ++i)
      for (size_t j = i + 1; j < num_adj_vertices; ++j)
      {
        if (boost::edge(adj_vertices[i], adj_vertices[j], g).second)
          ++existent_triangles;

        if (boost::is_directed(g) && boost::edge(adj_vertices[j], adj_vertices[i], g).second)
          ++existent_triangles; // CHECK THIS !!
      }

    return decimal(existent_triangles) / decimal(total_triangles);
  }


  template <class Graph>
  decimal extended_vertex_clustering(
      typename boost::graph_traits<Graph>::vertex_descriptor v,
      size_t order,
      const Graph & g
      )
  {
    typedef typename Graph::implementation_type implementation_type;
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::adjacency_iterator adjacency_iter;
    typedef typename GraphTraits::vertex_descriptor vertex;

    size_t V = boost::num_vertices(g);

    decimal **D;
    D = (double**) malloc(sizeof(double*) * V);
    for (int i = 0; i < V; ++i)
      D[i] = (double*) malloc(sizeof(double) * V);
    
    // FIXME: conan::adj_matrixS should be replaced by implementation_type
    all_pairs_shortest_path_dispatch<Graph>(g, D, conan::adj_matrixS());

    size_t num_adj_vertices = boost::out_degree(v, g);
    if (num_adj_vertices < 2)
      return 0;  //FIXME: Do the clustering coef of a vertex with a number of adjacent vertices less than 2 equal to 0? Where is this defined?

    size_t total_poligons = detail::binomial_coefficient(num_adj_vertices, 2);

    size_t i = 0;
    
    // make an array with the adjacent vertices to vertex
    vertex adj_vertices[num_adj_vertices];
    adjacency_iter vi, vi_end;
    for (tie(vi, vi_end) = boost::adjacent_vertices(v, g); vi != vi_end; ++vi)
      adj_vertices[i++] = *vi;

    size_t existent_poligons = 0;

    bool is_directed = boost::is_directed(g);

    for (i = 0; i < num_adj_vertices; ++i)
    {
      for (size_t j = i + 1; j < num_adj_vertices; ++j)
      {
        if ( size_t(D[ adj_vertices[i] ][ adj_vertices[j] ]) == order )
          ++existent_poligons;

        if ( is_directed and size_t(D[ adj_vertices[j] ][ adj_vertices[i] ]) == order )
          ++existent_poligons; // CHECK THIS !!
      }
    }

    for (int i = 0; i < V; ++i)
      free(D[i]);
    free(D);

    return decimal(existent_poligons) / decimal(total_poligons);
  }


  template <class Graph>
  decimal graph_avg_clustering(
      const Graph & g
      )
  {
    typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iter;
    
    int V = boost::num_vertices(g);
    decimal sum = 0;
    vertex_iter vi, vi_end;
    for (tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi)
      sum += vertex_clustering(*vi, g);

    return sum / (decimal)V;
  }

} // conan

#endif //CLUSTERING_HPP
