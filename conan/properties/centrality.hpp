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

#ifndef CENTRALITY_HPP
#define CENTRALITY_HPP

#include <conan/config.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#define HAVE_INLINE
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>


namespace conan {

  template <class Graph>
  struct betweenness_centrality
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex;
    typedef typename GraphTraits::edge_descriptor edge;
    typedef typename std::vector<decimal> VertexCentralityMap;
    typedef typename std::map<edge, decimal> EdgeCentralityMap;

    betweenness_centrality(Graph & g)
      :
      parent_graph(g),
      absolute_vertex_centrality( boost::num_vertices(g) ),
      relative_vertex_centrality(),
      edge_centrality_(),
      central_point_dominance_(0.0)
    {
      boost::associative_property_map< EdgeCentralityMap > edge_centrality_m(edge_centrality_);

      boost::brandes_betweenness_centrality(parent_graph,
          boost::centrality_map( boost::make_iterator_property_map(absolute_vertex_centrality.begin(), boost::get(boost::vertex_index, g), decimal()) ).
          edge_centrality_map( edge_centrality_m ).
          weight_map( get(&DefaultEdgeProperties::weight, g) )
          );

      relative_vertex_centrality = absolute_vertex_centrality;

      boost::relative_betweenness_centrality(parent_graph,
          boost::make_iterator_property_map(relative_vertex_centrality.begin(), boost::get(boost::vertex_index, g), decimal())
          );

      central_point_dominance_ = boost::central_point_dominance(parent_graph,
          boost::make_iterator_property_map(relative_vertex_centrality.begin(), boost::get(boost::vertex_index, g), decimal())
          );

      return;
    }

    decimal vertex_centrality(vertex v)
    { return absolute_vertex_centrality[v]; }

    decimal vertex_relative_centrality(vertex v)
    { return relative_vertex_centrality[v]; }

    decimal edge_centrality(edge e)
    { return edge_centrality_[e]; }

    decimal central_point_dominance()
    { return central_point_dominance_; }

    Graph & parent_graph;
    VertexCentralityMap absolute_vertex_centrality;
    VertexCentralityMap relative_vertex_centrality;
    EdgeCentralityMap edge_centrality_;
    decimal central_point_dominance_;
  };

  template <class Graph>
  decimal vertex_degree_centrality(
      typename boost::graph_traits<Graph>::vertex_descriptor v,
      const Graph & g
      )
  {
    return decimal( boost::out_degree(v, g) ) / decimal( boost::num_vertices(g) - 1 );
  }

  template <class Graph>
  decimal graph_degree_centrality(
      const Graph & g
      )
  {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex;
    decimal sum = 0.0, max_vdc = 0.0;
    size_t V = boost::num_vertices(g);
    std::vector<decimal> degree_centrality;
    BOOST_FOREACH( vertex v, boost::vertices(g) )
    {
      decimal vdc = vertex_degree_centrality(v, g);
      degree_centrality.push_back( vdc );
      if ( vdc > max_vdc )
        max_vdc = vdc;
    }

    BOOST_FOREACH( decimal vdc, degree_centrality )
      sum += (max_vdc - vdc);

    if (sum == 0.0)
      return 0.0;
    if (V <= 2)
      throw std::runtime_error("Graph is too small (|V| <= 2) for computing its degree centrality.");

    return sum / decimal( V - 2 );
  }

  template <class Graph>
  decimal vertex_closeness(
      typename boost::graph_traits<Graph>::vertex_descriptor v,
      Graph & g
      )
  {
    int V = boost::num_vertices(g);
    decimal sum = 0;
    std::vector<decimal> d(V);

    boost::dijkstra_shortest_paths(g, v,
        boost::weight_map(get(&DefaultEdgeProperties::weight, g)).distance_map(&d[0]));

    BOOST_FOREACH( decimal distance, d )
    {
      if ( distance >= std::numeric_limits<decimal>::max() )
      {
        --V;
        continue;
      }
      sum += distance;
    }

    if (V < 2)
      throw std::runtime_error("Vertex closeness could be computed because the size of the connected components reacheable from v is less than 2.");
    
    return decimal(V - 1) / sum;
  }

  template <class Graph>
  struct eigenvector_centrality
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::edge_descriptor edge;
    typedef typename GraphTraits::vertex_descriptor vertex;
    typedef typename GraphTraits::directed_category directed_category;

    eigenvector_centrality(const Graph & g)
      : parent_graph(g), eigencentrality()
    {
      size_t V = boost::num_vertices(g);
      bool is_undirected = boost::is_undirected(g);

      // Create adjacency matrix
      gsl_matrix * A = gsl_matrix_calloc(V, V);

      BOOST_FOREACH( const edge & e, boost::edges(g) )
      {
        gsl_matrix_set(A, boost::source(e, g), boost::target(e, g), 1.0);
        if ( is_undirected )
          gsl_matrix_set(A, boost::target(e, g), boost::source(e, g), 1.0);
      }

      // Compute eigenvalues
      gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(V);
      gsl_matrix * evec = gsl_matrix_alloc(V, V);
      gsl_vector * eval = gsl_vector_alloc(V);

      gsl_eigen_symmv(A, eval, evec, w);

      gsl_eigen_symmv_free(w);
      gsl_matrix_free(A);

      size_t max_eigenvalue_index = 0;
      double max_eigenvalue = 0;
      for (size_t i = 0; i < eval->size; ++i)
      {
        if (gsl_vector_get(eval, i) > max_eigenvalue)
        {
          max_eigenvalue = gsl_vector_get(eval, i);
          max_eigenvalue_index = i;
        }
      }

      gsl_vector_view max_eval_eigenvector = gsl_matrix_column(evec, max_eigenvalue_index);

      gsl_vector_free(eval);

      // search for the maximum eigenvector centrality to normalize values
      double max_eigencentrality = 0.0, max_eigencentrality_abs = 0.0;
      for (size_t i = 0; i < max_eval_eigenvector.vector.size; ++i)
        if ( max_eigencentrality_abs < fabs(gsl_vector_get(&max_eval_eigenvector.vector, i)) )
          max_eigencentrality_abs = fabs( max_eigencentrality = gsl_vector_get(&max_eval_eigenvector.vector, i) );

      // store the normalized eigenvector centrality
      for (size_t i = 0; i < max_eval_eigenvector.vector.size; ++i)
        eigencentrality.push_back( gsl_vector_get(&max_eval_eigenvector.vector, i) / max_eigencentrality );

      gsl_matrix_free(evec);
      
      return;
    }

    decimal operator()(vertex v)
    { return eigencentrality[v]; }

    const Graph & parent_graph;
    std::vector<decimal> eigencentrality;
  };

};

#endif // CENTRALITY_HPP
