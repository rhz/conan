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

#ifndef BOXCOUNTING_HPP
#define BOXCOUNTING_HPP
#include <conan/graphs.hpp>
#include <conan/utils/clusterize.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>


namespace conan { namespace detail {

  /**
   * @brief Boxcounting is a method to calculate the fractal dimension of a graph.
   *
   * It pick randomly a vertex and join it to all vertices near than boxsize, forming a cluster-vertex.
   * This process is repeated for each vertex that isn't already part of a cluster.
   * The result is a graph containing a set of clusters randomly selected from the input graph.
   * 
   * In order to assess the graph's fractal dimension, you must graph the number of vertices versus
   * the boxsize. The fractal's dimension is the slope of this graph.
   * 
   * @param g the graph over which apply the boxcounting process.
   * @param boxsize the box size.
   * @return the renormalized graph (it's the same type of graph as the input graph).
   */
  template <class Graph> Graph boxcounting (Graph& g, decimal boxsize);


  template <typename T1>
  void erase_in_vector(
      std::vector<T1> & l,
      const T1 & elem
      )
  {
    for (typename std::vector<T1>::iterator li = l.begin(); li != l.end(); ++li)
    {
      if (*li == elem)
      {
        l.erase(li);
        return;
      }
    }

    return;
  }


  template <class Graph, class RNG>
  void create_clusters(
      Graph & g,
      std::vector< std::vector<size_t> > & clusters,
      decimal boxsize,
      RNG & rng
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex;
    typedef typename GraphTraits::vertex_iterator vertex_iter;

    typedef typename std::vector<vertex> VertexList;
    typedef typename std::vector<vertex> Cluster;
    typedef typename std::vector<Cluster> ClusterVector;
    typedef typename std::vector<decimal> DistanceVector;

    typedef typename VertexList::iterator VertexListIter;
    
    typename boost::uniform_real<decimal> percent_dist(0, 1);
    typename boost::variate_generator< RNG &, boost::uniform_real<decimal> >
      percent(rng, percent_dist);

    VertexList not_counted_vertices( boost::num_vertices(g) );

    { // Copy all vertices to non_counted_vertices
      vertex_iter vi, viend;
      tie(vi, viend) = boost::vertices(g);
      std::copy(vi, viend, not_counted_vertices.begin());
    }
    
    size_t length = size_t(boxsize);
    decimal fraction = boxsize - decimal(length);

    size_t current_cluster = 0;

    while (not_counted_vertices.size())
    {
      typename boost::uniform_int<size_t> dist(0, not_counted_vertices.size() - 1);
      typename boost::variate_generator< RNG &, boost::uniform_int<size_t> >
        dice(rng, dist);

      clusters.push_back( Cluster() );

      vertex current_vertex = dice();
      clusters[current_cluster].push_back(current_vertex);

      // Compute distances
      DistanceVector d( boost::num_vertices(g) );

      VertexListIter li = not_counted_vertices.begin();
      li += current_vertex;

      boost::dijkstra_shortest_paths(g, *li,
          boost::weight_map(get(&DefaultEdgeProperties::weight, g)).distance_map(&d[0]));

      // Add to current cluster all vertices within a distance of lenght from the current_vertex
      size_t target_vertex_index = 0;

      BOOST_FOREACH( decimal & distance, d )
      { // iterate over all asp values for current_vertex
        if (distance < length || ((distance < length + 1) && (percent() < fraction))) // check if the target is in our box
          clusters[current_cluster].push_back(target_vertex_index);

        ++target_vertex_index;
      }

      // Delete added vertices from g and not_counted_vertices
      BOOST_FOREACH( size_t & vertex_to_remove, clusters[current_cluster] )
      {
        boost::clear_vertex(vertex_to_remove, g);
        erase_in_vector(not_counted_vertices, vertex_to_remove);
      }

      ++current_cluster;
    }

    return;
  }


  template <class Graph>
  Graph boxcounting(
      Graph & g,
      decimal boxsize
  )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex;
    typedef typename GraphTraits::edge_iterator edge_iter;

    typedef typename std::vector<vertex> VertexCluster;
    typedef typename std::vector<VertexCluster> ClusterList;
    
    size_t V = boost::num_vertices(g);

    if (V <= size_t(boxsize))
      return Graph(1); // shortcut

    // Copy the original graph to g_copy
    Graph g_copy(V);

    edge_iter ei, ei_end;
    for (tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
      boost::add_edge(boost::source(*ei, g), boost::target(*ei, g), g[*ei].weight, g_copy);

    // Initialize random number generators
#ifdef DETERMINISTIC_RNG
    typename boost::mt19937 rng(time(NULL) - V - int(boxsize * 10));
#else
    typename boost::random_device rng;
#endif
    
    ClusterList clusters;
    create_clusters(g_copy, clusters, boxsize, rng);

    Graph renormalized_graph( clusterize( g, clusters ) );

    return renormalized_graph;
  } // boxcounting

}} // conan::detail

#endif //BOXCOUNTING_HPP
