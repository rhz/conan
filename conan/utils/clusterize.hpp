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
#ifndef CLUSTERIZE_HPP
#define CLUSTERIZE_HPP
#include <conan/graphs.hpp>


namespace conan {
  
  namespace detail {

    template <class Graph>
    bool are_clusters_linked(
        Graph & g,
        const std::vector<size_t> & source_cluster,
        const std::vector<size_t> & target_cluster
        )
    {
      BOOST_FOREACH( const size_t & source_vertex, source_cluster )
      {
        BOOST_FOREACH( const size_t & target_vertex, target_cluster )
        {
          if ( boost::edge(source_vertex, target_vertex, g).second )
            return true;
        }
      }

      return false;
    }


    bool is_in_cluster_vector(
        size_t v,
        const std::vector< std::vector<size_t> > & clusters
        )
    {
      typedef std::vector<size_t> Cluster;

      BOOST_FOREACH( const Cluster & c, clusters )
      {
        BOOST_FOREACH( const size_t & v_c, c )
        {
          if ( v == v_c )
            return true;
        }
      }

      return false;
    }

  } // detail


  template <class Graph>
  Graph clusterize(
      const Graph & original_graph,
      std::vector< std::vector<size_t> > & clusters,
      bool include_all_vertices = false
      )
  {
    using namespace detail;
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex;

    typedef typename std::vector<vertex> Cluster;
    typedef typename std::vector<Cluster> ClusterVector;
    typedef typename ClusterVector::iterator ClusterVectorIter;

    Graph renormalized_graph( clusters.size() );

    ClusterVector * clusters_ptr = &clusters;
    if ( include_all_vertices )
    {
      clusters_ptr = new ClusterVector( clusters );
      BOOST_FOREACH( const vertex & v, boost::vertices( original_graph ) )
      {
        if ( ! is_in_cluster_vector( v, clusters ) )
        {
          Cluster c;
          c.push_back( v );
          clusters_ptr->push_back( c );
        }
      }
    }
    
    size_t source_cluster_index = 0;
    bool is_directed = boost::is_directed( renormalized_graph );

    ClusterVectorIter source_cluster_end, target_cluster_end;
    if (is_directed)
      source_cluster_end = target_cluster_end = clusters_ptr->end();
    else
      source_cluster_end = ( target_cluster_end = clusters_ptr->end() ) - 1;

    // If it is found any edge between source_cluster and target_cluster, then this two cluster will be linked in renormalized_graph
    for (ClusterVectorIter source_cluster = clusters_ptr->begin(); source_cluster != source_cluster_end; ++source_cluster, ++source_cluster_index)
    {
      size_t target_cluster_index = (is_directed? 0 : source_cluster_index + 1);

      for (ClusterVectorIter target_cluster = clusters_ptr->begin() + ( is_directed? 0 : target_cluster_index );
           target_cluster != target_cluster_end; ++target_cluster, ++target_cluster_index)
      {
        bool edge_found = are_clusters_linked(original_graph, *source_cluster, *target_cluster);

        if (edge_found)
          boost::add_edge(source_cluster_index, target_cluster_index, 1.0, renormalized_graph);
      }
    }

    if ( clusters_ptr != &clusters )
      delete clusters_ptr;

    return renormalized_graph;
  }

} // conan

#endif // CLUSTERIZE_HPP
