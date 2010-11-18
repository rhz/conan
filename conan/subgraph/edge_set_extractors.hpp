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

#ifndef EDGE_SET_EXTRACTORS_HPP
#define EDGE_SET_EXTRACTORS_HPP
#include <conan/graphs.hpp>
#include <limits>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

namespace conan {
  /**
   * nearest_neighbour_ring (NN = Nearest Neighbor)
   * @brief Create a ring-shaped graph in which each vertex has only two edges, no more, no less.
   *
   * It returns a ring-graph constructed by the iterative process of taking a vertex (in the first
   * step we take start_vertex) and link it to its nearest vertex. Then we take the this one and
   * apply the same process on it and so on.
   * 
   * @param g base graph
   * @param start_vertex the vertex from which start the chain (ring actually).
   * @param nn_ring the ring graph
   */
  template <class Graph>
  void nearest_neighbour_ring(
      const Graph &,
      typename boost::graph_traits<Graph>::vertex_descriptor,
      std::vector< typename boost::graph_traits<Graph>::edge_descriptor > &
      );
 
  /// Simpler variant of conan::nearest_neighbour_ring that begins at vertex 0.
  template <class Graph>
  void nearest_neighbour_ring(
      const Graph &,
      std::vector< typename boost::graph_traits<Graph>::edge_descriptor > &
      );

  /**
   * minimum_spanning_tree
   * @brief Return the minimum spanning tree graph for the given undirected graph.
   * 
   * For reference on what the minimum spanning tree is, I recommend you to read the
   * Boost Graph Library documentation website or a good graph-theoretic book :) \n
   * @warning This algorithm works only for %undirected graphs.
   * Used to get an initial graph over which add the real graph (phantom graph).
   * So we'll always have only one connected component.
   * 
   * @param g Base graph
   * @param mst Vector which at the end of the function will contain the edge set of the minimum spanning tree
   */
  template <class Graph>
  void minimum_spanning_tree(
      Graph &,
      std::vector< typename boost::graph_traits<Graph>::edge_descriptor > &
      );


  // ***** Declarations *****


  namespace detail {

    template <class Edge>
    int is_vertex_in_edge_list(
        std::vector<Edge> & edge_list,
        size_t vertex
        )
    {
      typedef typename std::vector<Edge> EdgeVector;
      typedef typename EdgeVector::iterator EdgeVectorIter;

      for (EdgeVectorIter evi = edge_list.begin(); evi != edge_list.end(); ++evi)
      {
        if (evi->m_source == vertex || evi->m_target == vertex)
          return true;
      }

      return false;
    }

  } // detail


  template <class Graph>
  void nearest_neighbour_ring(
      const Graph & g,
      std::vector< typename boost::graph_traits<Graph>::edge_descriptor > & nn_ring
      )
  {
    return nearest_neighbour_ring(g, boost::vertex(0, g), nn_ring);
  }
  

  template <class Graph>
  void nearest_neighbour_ring(
      const Graph & g, // connected weighted graph
      typename boost::graph_traits<Graph>::vertex_descriptor start_vertex,
      std::vector< typename boost::graph_traits<Graph>::edge_descriptor > & nn_ring
      )
  {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor edge;
    typedef typename boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;

    typedef typename std::vector<edge> EdgeVector;
    typedef typename EdgeVector::iterator EdgeVectorIter;

    nn_ring.clear();

    out_edge_iterator ei, ei_end;
    vertex s = start_vertex;

    // least_weighted_edge and min_weight store the edge_descriptor and the weight of the least weighted
    edge least_weighted_edge;
    decimal min_weight;

    // Slow algorithm... many edges are viewed more than once
    while (true)
    {
      min_weight = std::numeric_limits<decimal>::max();

      // select next vertex
      for (tie(ei, ei_end) = boost::out_edges(s, g); ei != ei_end; ++ei)
      {
        // check if the target vertex is already in nn_ring
        if ( detail::is_vertex_in_edge_list( nn_ring, boost::target(*ei, g) ) == true )
          continue;

        // store the edge with the smallest weight coming from start_vertex
        if (min_weight > g[ *ei ].weight)
        {
          min_weight = g[ *ei ].weight;
          least_weighted_edge = *ei;
        }
      }

      if ( min_weight == std::numeric_limits<decimal>::max() ) // there are no more vertices to connect to
        break;

      nn_ring.push_back(least_weighted_edge);

      s = boost::target(least_weighted_edge, g);
    }
    
    // bind the last vertex with the first one to close the circle
    edge e;
    bool edge_found;
    tie(e, edge_found) = boost::edge(s, start_vertex, g);
    if (edge_found)
    {
      nn_ring.push_back(e);
    }
    else
    {
      throw std::runtime_error("Input graph is not complete and therefore Nearest Neighbour algorithm cannot be carried out on it.");
    }

    return;
  }


  template <class Graph>
  void minimum_spanning_tree(
      Graph & g, // connected weighted graph
      std::vector< typename boost::graph_traits<Graph>::edge_descriptor > & mst // minimum spanning tree
      )
  {
    typedef typename boost::graph_traits<Graph>::directed_category directed_category;
    
    BOOST_STATIC_ASSERT( ( boost::is_same< directed_category, boost::undirected_tag >::value ) ); // check if g is directed

    mst.clear();
    
    boost::kruskal_minimum_spanning_tree ( g ,
            std::back_inserter ( mst ) ,
            boost::weight_map ( get ( &DefaultEdgeProperties::weight, g ) )
            );

    return;
  }
}

#endif // EDGE_SET_EXTRACTORS_HPP
