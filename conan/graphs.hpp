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

#ifndef GRAPHS_HPP
#define GRAPHS_HPP
#include <conan/config.hpp>

#include <boost/version.hpp>
#if BOOST_VERSION > 104100
  #include <boost/property_map/property_map.hpp>
#else
  #include <boost/property_map.hpp> // deprecated for Boost >= 1.41
#endif

// Boost Graph Library
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>

#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>

namespace conan {

  using boost::add_vertex;
  using boost::add_edge;
  using boost::graph_traits;
  using boost::edge;

  // Structs

  /**
   * @brief Adjacency List implementation tag.
   *
   * It's used to define the implementation a graph will use.
   */
  struct adj_listS { };

  /**
   * @brief Adjacency Matrix implementation tag.
   *
   * It's used to define the implementation a graph will use.
   */
  struct adj_matrixS { };

  /**
   * @defgroup graph Graph classes
   * @{
   */

  /// @cond
  // Null variant
  template <class Impl = conan::adj_listS,
            class VertexProperties = conan::DefaultVertexProperties,
            class EdgeProperties = conan::DefaultEdgeProperties>
  class undirected_graph
  {
    BOOST_STATIC_ASSERT(( boost::is_same<Impl, conan::adj_listS>::value ||
                          boost::is_same<Impl, conan::adj_matrixS>::value ));
  };
  /// @endcond
  
  /**
   * undirected_graph class
   * @brief Stores an %undirected graph in an adjacency list.
   */
  // Adjacency List variant
  template <class VertexProperties,
            class EdgeProperties>
  class undirected_graph<conan::adj_listS, VertexProperties, EdgeProperties>
    : public boost::adjacency_list<boost::vecS,
                                   boost::vecS,
                                   boost::undirectedS,
                                   VertexProperties,
                                   EdgeProperties>
  {
  public:
    // Typedefs
    typedef boost::vecS OutEdgeList;
    typedef boost::vecS VertexList;
    typedef typename boost::adjacency_list<OutEdgeList,
                                           VertexList,
                                           boost::undirectedS,
                                           VertexProperties,
                                           EdgeProperties>
      BaseGraph;
    typedef typename conan::adj_listS                implementation_type;
    typedef typename boost::graph_traits<BaseGraph>  Traits;
    typedef typename BaseGraph::vertex_descriptor    vertex_descriptor;
    typedef typename BaseGraph::vertex_iterator      vertex_iterator;
    typedef typename BaseGraph::adjacency_iterator   adjacency_iterator;
    typedef typename BaseGraph::edge_descriptor      edge_descriptor;
    typedef typename BaseGraph::edge_iterator        edge_iterator;
    typedef typename BaseGraph::out_edge_iterator    out_edge_iterator;
    typedef typename BaseGraph::in_edge_iterator     in_edge_iterator;
    typedef typename BaseGraph::vertices_size_type   vertices_size_type;
    typedef typename BaseGraph::edges_size_type      edges_size_type;
    
    // These typedefs are for compatibility with boost::graph_traits
    typedef typename Traits::directed_category       directed_category;
    typedef typename Traits::edge_parallel_category  edge_parallel_category;
    typedef typename Traits::traversal_category      traversal_category;
    typedef typename Traits::degree_size_type        degree_size_type;

    // Properties
    typedef typename BaseGraph::vertex_property_type vertex_property_type;
    typedef typename BaseGraph::edge_property_type   edge_property_type;

    // The types that are actually bundled
    typedef typename BaseGraph::vertex_bundled       vertex_bundled;
    typedef typename BaseGraph::edge_bundled         edge_bundled;

    // Constructors
    inline undirected_graph()
      : BaseGraph() { }

    inline undirected_graph(const BaseGraph& x)
      : BaseGraph(x) { }

    inline undirected_graph(vertices_size_type num_vertices)
      : BaseGraph(num_vertices) { }

    template <class EdgeIterator>
    inline undirected_graph(EdgeIterator first, EdgeIterator last,
                            vertices_size_type n)
      : BaseGraph(first, last, n) { }

    template <class EdgeIterator>
    inline undirected_graph(EdgeIterator first, EdgeIterator last,
                            decimal weight,
                            vertices_size_type n)
      : BaseGraph(first, last, n)
    {
      BOOST_FOREACH( const edge_descriptor & e, boost::edges(*this) )
      { (*this)[ e ].weight = weight; }
    }

    template <class EdgeIterator, class EdgePropertyIterator>
    inline undirected_graph(EdgeIterator first, EdgeIterator last,
                            EdgePropertyIterator ep_iter,
                            vertices_size_type n)
      : BaseGraph(first, last, ep_iter, n) { }
  };
  
  /**
   * undirected_graph class
   * @brief Stores an %undirected graph in an adjacency matrix.
   */
  // Adjacency Matrix variant
  template <class VertexProperties,
            class EdgeProperties>
  class undirected_graph<conan::adj_matrixS, VertexProperties, EdgeProperties>
    : public boost::adjacency_matrix<boost::undirectedS,
                                     VertexProperties,
                                     EdgeProperties>
  {
  public:
    // Typedefs
    typedef typename boost::adjacency_matrix<boost::undirectedS,
                                             VertexProperties,
                                             EdgeProperties>
      BaseGraph;
    typedef typename conan::adj_matrixS              implementation_type;
    typedef typename boost::graph_traits<BaseGraph>  Traits;
    typedef typename BaseGraph::vertex_descriptor    vertex_descriptor;
    typedef typename BaseGraph::vertex_iterator      vertex_iterator;
    typedef typename BaseGraph::adjacency_iterator   adjacency_iterator;
    typedef typename BaseGraph::edge_descriptor      edge_descriptor;
    typedef typename BaseGraph::edge_iterator        edge_iterator;
    typedef typename BaseGraph::out_edge_iterator    out_edge_iterator;
    typedef typename BaseGraph::in_edge_iterator     in_edge_iterator;
    typedef typename BaseGraph::vertices_size_type   vertices_size_type;
    typedef typename BaseGraph::edges_size_type      edges_size_type;
    
    // These typedefs are for compatibility with boost::graph_traits
    typedef typename Traits::directed_category       directed_category;
    typedef typename Traits::edge_parallel_category  edge_parallel_category;
    typedef typename Traits::traversal_category      traversal_category;
    typedef typename Traits::degree_size_type        degree_size_type;

    // Properties
    typedef typename BaseGraph::vertex_property_type vertex_property_type;
    typedef typename BaseGraph::edge_property_type   edge_property_type;

    // The types that are actually bundled
    typedef typename BaseGraph::vertex_bundled       vertex_bundled;
    typedef typename BaseGraph::edge_bundled         edge_bundled;

    // Constructors
    inline undirected_graph()
      : BaseGraph(0) { }

    inline undirected_graph(vertices_size_type n)
      : BaseGraph(n) { }

    template <class EdgeIterator>
    inline undirected_graph(EdgeIterator first, EdgeIterator last,
                            vertices_size_type n)
      : BaseGraph(first, last, n) { }

    template <class EdgeIterator>
    inline undirected_graph(EdgeIterator first, EdgeIterator last,
                            decimal weight,
                            vertices_size_type n)
      : BaseGraph(first, last, n)
    {
      BOOST_FOREACH( edge_descriptor & e, boost::edges(*this) )
      { (*this)[ e ].weight = weight; }
    }

    template <class EdgeIterator, class EdgePropertyIterator>
    inline undirected_graph(EdgeIterator first, EdgeIterator last,
                            EdgePropertyIterator ep_iter,
                            vertices_size_type n)
      : BaseGraph(first, last, ep_iter, n) { }
  };


  /// @cond
  // Null variant
  template <class Impl = conan::adj_listS,
            class VertexProperties = conan::DefaultVertexProperties,
            class EdgeProperties = conan::DefaultEdgeProperties>
  class directed_graph
  {
    BOOST_STATIC_ASSERT(( boost::is_same<Impl, conan::adj_listS>::value ||
                          boost::is_same<Impl, conan::adj_matrixS>::value ));
  };
  /// @endcond

  /**
   * directed_graph class
   * @brief Stores an %directed graph in an adjacency list.
   */
  // Adjacency List variant
  template <class VertexProperties,
            class EdgeProperties>
  class directed_graph<conan::adj_listS, VertexProperties, EdgeProperties>
    : public boost::adjacency_list<boost::vecS,
                                   boost::vecS,
                                   boost::directedS,
                                   VertexProperties,
                                   EdgeProperties>
  {
  public:
    // Typedefs
    typedef boost::vecS OutEdgeList;
    typedef boost::vecS VertexList;
    typedef typename boost::adjacency_list<OutEdgeList,
                                           VertexList,
                                           boost::directedS,
                                           VertexProperties,
                                           EdgeProperties>
      BaseGraph;
    typedef typename conan::adj_listS                implementation_type;
    typedef typename boost::graph_traits<BaseGraph>  Traits;
    typedef typename BaseGraph::vertex_descriptor    vertex_descriptor;
    typedef typename BaseGraph::vertex_iterator      vertex_iterator;
    typedef typename BaseGraph::adjacency_iterator   adjacency_iterator;
    typedef typename BaseGraph::edge_descriptor      edge_descriptor;
    typedef typename BaseGraph::edge_iterator        edge_iterator;
    typedef typename BaseGraph::out_edge_iterator    out_edge_iterator;
    typedef typename BaseGraph::in_edge_iterator     in_edge_iterator;
    typedef typename BaseGraph::vertices_size_type   vertices_size_type;
    typedef typename BaseGraph::edges_size_type      edges_size_type;
    
    // These typedefs are for compatibility with boost::graph_traits
    typedef typename Traits::directed_category       directed_category;
    typedef typename Traits::edge_parallel_category  edge_parallel_category;
    typedef typename Traits::traversal_category      traversal_category;
    typedef typename Traits::degree_size_type        degree_size_type;

    // Properties
    typedef typename BaseGraph::vertex_property_type vertex_property_type;
    typedef typename BaseGraph::edge_property_type   edge_property_type;

    // The types that are actually bundled
    typedef typename BaseGraph::vertex_bundled       vertex_bundled;
    typedef typename BaseGraph::edge_bundled         edge_bundled;

    // Constructors
    inline directed_graph()
      : BaseGraph() { }

    inline directed_graph(const BaseGraph& x)
      : BaseGraph(x) { }

    inline directed_graph(vertices_size_type num_vertices)
      : BaseGraph(num_vertices) { }

    template <class EdgeIterator>
    inline directed_graph(EdgeIterator first, EdgeIterator last,
                          vertices_size_type n)
      : BaseGraph(first, last, n) { }

    template <class EdgeIterator>
    inline directed_graph(EdgeIterator first, EdgeIterator last,
                          decimal weight,
                          vertices_size_type n)
      : BaseGraph(first, last, n)
    {
      BOOST_FOREACH( edge_descriptor & e, boost::edges(*this) )
      { (*this)[ e ].weight = weight; }
    }

    template <class EdgeIterator, class EdgePropertyIterator>
    inline directed_graph(EdgeIterator first, EdgeIterator last,
                          EdgePropertyIterator ep_iter,
                          vertices_size_type n)
      : BaseGraph(first, last, ep_iter, n) { }
  };
  
  /**
   * directed_graph class
   * @brief Stores an %directed graph in an adjacency matrix.
   */
  // Adjacency Matrix variant
  template <class VertexProperties,
            class EdgeProperties>
  class directed_graph<conan::adj_matrixS, VertexProperties, EdgeProperties>
    : public boost::adjacency_matrix<boost::directedS,
                                     VertexProperties,
                                     EdgeProperties>
  {
  public:
    // Typedefs
    typedef typename boost::adjacency_matrix<boost::directedS,
                                             VertexProperties,
                                             EdgeProperties>
      BaseGraph;
    typedef typename conan::adj_matrixS              implementation_type;
    typedef typename boost::graph_traits<BaseGraph>  Traits;
    typedef typename BaseGraph::vertex_descriptor    vertex_descriptor;
    typedef typename BaseGraph::vertex_iterator      vertex_iterator;
    typedef typename BaseGraph::adjacency_iterator   adjacency_iterator;
    typedef typename BaseGraph::edge_descriptor      edge_descriptor;
    typedef typename BaseGraph::edge_iterator        edge_iterator;
    typedef typename BaseGraph::out_edge_iterator    out_edge_iterator;
    typedef typename BaseGraph::in_edge_iterator     in_edge_iterator;
    typedef typename BaseGraph::vertices_size_type   vertices_size_type;
    typedef typename BaseGraph::edges_size_type      edges_size_type;
    
    // These typedefs are for compatibility with boost::graph_traits
    typedef typename Traits::directed_category       directed_category;
    typedef typename Traits::edge_parallel_category  edge_parallel_category;
    typedef typename Traits::traversal_category      traversal_category;
    typedef typename Traits::degree_size_type        degree_size_type;

    // Properties
    typedef typename BaseGraph::vertex_property_type vertex_property_type;
    typedef typename BaseGraph::edge_property_type   edge_property_type;

    // The types that are actually bundled
    typedef typename BaseGraph::vertex_bundled       vertex_bundled;
    typedef typename BaseGraph::edge_bundled         edge_bundled;

    // Constructors
    inline directed_graph()
      : BaseGraph(0) { }

    inline directed_graph(vertices_size_type n)
      : BaseGraph(n) { }

    template <class EdgeIterator>
    inline directed_graph(EdgeIterator first, EdgeIterator last,
                          vertices_size_type n)
      : BaseGraph(first, last, n) { }

    template <class EdgeIterator>
    inline directed_graph(EdgeIterator first, EdgeIterator last,
                          decimal weight,
                          vertices_size_type n)
      : BaseGraph(first, last, n)
    {
      BOOST_FOREACH( edge_descriptor & e, boost::edges(*this) )
      { (*this)[ e ].weight = weight; }
    }

    template <class EdgeIterator, class EdgePropertyIterator>
    inline directed_graph(EdgeIterator first, EdgeIterator last,
                          EdgePropertyIterator ep_iter,
                          vertices_size_type n)
      : BaseGraph(first, last, ep_iter, n) { }
  };

  /** @} */


  /**
   * @brief Find a vertex in a graph by giving its name.
   *
   * @param g The graph in which the search will occur.
   * @param name The vertex name you are searching for.
   * @return The descriptor of the vertex found, if any. Otherwise, a null descriptor
   *   (boost::graph_traits<Graph>::null_vertex).
   */
  template <class Graph>
  typename boost::graph_traits<Graph>::vertex_descriptor find_vertex_by_name(
      std::string name,
      const Graph & g
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_iterator vertex_iter;

    vertex_iter vi, viend;
    for (tie(vi, viend) = boost::vertices(g); vi != viend; ++vi)
    {
      if (g[*vi].name == name)
        return *vi;
    }
    
    return boost::graph_traits<Graph>::null_vertex();
  }


  /**
   * @brief Calculate the average degree of a graph.
   *
   * @param g A graph.
   * @return The value of the average degree.
   */
  template <class Graph>
  decimal avg_degree(
      const Graph & g
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_iterator vertex_iter;

    decimal sum_degree = 0.0;

    vertex_iter vi, viend;
    for (tie(vi, viend) = boost::vertices(g); vi != viend; ++vi)
      sum_degree += boost::out_degree(*vi, g);

    return sum_degree / decimal(boost::num_vertices(g));
  }

  namespace detail
  {
    struct found_target { }; // exception for termination

    template <class DistanceMap, class Graph, class Tag>
    struct my_distance_recorder
      : public boost::base_visitor<my_distance_recorder<DistanceMap, Graph, Tag> >
    {
      typedef Tag event_filter;

      my_distance_recorder(
          DistanceMap pa,
          typename boost::graph_traits<Graph>::vertex_descriptor target
          )
        : m_distance(pa), target_vertex(target) { }

      template <class Edge>
      void operator()(Edge e, const Graph & g)
      {
        typename boost::graph_traits<Graph>::vertex_descriptor
          u = source(e, g), v = target(e, g);

        put(m_distance, v, get(m_distance, u) + g[ e ].weight);

        if (v == target_vertex)
          throw found_target();
      }

      DistanceMap m_distance;
      typename boost::graph_traits<Graph>::vertex_descriptor target_vertex;
    };
  } // detail

  /// @brief Find the distance between two vertices.
  template <class Graph>
  decimal distance(
      typename boost::graph_traits<Graph>::vertex_descriptor source_vertex,
      typename boost::graph_traits<Graph>::vertex_descriptor target_vertex,
      Graph & g
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertices_size_type vertices_size_type;

    size_t V = boost::num_vertices(g);
    if ( V == 0 )
      return 0.0;

    decimal d[ V ];
    std::fill_n(d, V, 0.0);

    try {
      boost::breadth_first_search( g, source_vertex,
          boost::visitor(
            boost::make_bfs_visitor(
              detail::my_distance_recorder<decimal *, Graph, boost::on_tree_edge>(&d[0], target_vertex)
              ) ) );
    }
    catch ( detail::found_target & e ) { }

    return d[ target_vertex ];
  }

} // conan

namespace boost {

  /// @brief Overload function for adding an edge with a given weight.
  template <class G>
  inline typename std::pair< typename graph_traits<G>::edge_descriptor, bool >
  add_edge(
      typename graph_traits<G>::vertex_descriptor u,
      typename graph_traits<G>::vertex_descriptor v,
      conan::decimal weight,
      G & g
      )
  {
    typedef typename graph_traits<G>::edge_descriptor edge;
    bool inserted;
    edge e;
    tie(e, inserted) = add_edge(u, v, g);
    g[e].weight = weight;
    return std::make_pair(e, inserted);
  }

  /// @brief Overload function for adding a vertex with a given name.
  template <class G>
  inline typename graph_traits<G>::vertex_descriptor
  add_vertex(
      std::string name,
      G & g
      )
  {
    typedef typename graph_traits<G>::vertex_descriptor vertex;
    vertex v = add_vertex(g);
    g[v].name = name;
    return v;
  }

} // boost

namespace conan {

  template <class GraphList>
  typename GraphList::value_type
  merge_graphs(
      GraphList & graphs,
      bool overlap = false
      )
  {
    typedef typename GraphList::value_type
        Graph;

    if (graphs.size() == 0)
      return Graph();

    typedef typename boost::graph_traits<Graph>::vertex_descriptor
      vertex;
    typedef typename boost::graph_traits<Graph>::vertex_iterator
      vertex_iter;
    typedef typename boost::graph_traits<Graph>::edge_descriptor
      edge;
    typedef typename std::unordered_map<Graph *, std::unordered_map<vertex, vertex> >
      VertexMap;

    size_t N = 0;
    if (overlap)
    {
      BOOST_FOREACH( Graph & g, graphs )
        if (N < boost::num_vertices(g))
          N = boost::num_vertices(g);
    }
    else
    {
      BOOST_FOREACH( Graph & g, graphs )
        N += boost::num_vertices(g);
    }
    
    Graph out(N);
    VertexMap vertex_map;

    if (overlap)
    {
      BOOST_FOREACH( Graph & g, graphs )
      {
        vertex_iter out_v = boost::vertices(out).first;
        BOOST_FOREACH( vertex v, boost::vertices(g) )
        {
          vertex_map[&g][v] = *out_v;
          out[*out_v] = g[v];
          ++out_v;
        }
      }
    }
    else
    {
      vertex_iter out_v = boost::vertices(out).first;
      BOOST_FOREACH( Graph & g, graphs )
      {
        BOOST_FOREACH( vertex v, boost::vertices(g) )
        {
          vertex_map[&g][v] = *out_v;
          out[*out_v] = g[v];
          ++out_v;
        }
      }
    }

    BOOST_FOREACH( Graph & g, graphs )
      BOOST_FOREACH( edge e, boost::edges(g) )
      {
        vertex s = vertex_map[&g][ boost::source(e, g) ],
               t = vertex_map[&g][ boost::target(e, g) ];

        if (not boost::edge(s, t, out).second)
          boost::add_edge(s, t, g[e].weight, out);
      }

    return out;
  }

} // conan

#include <conan/subgraph.hpp>

#endif //GRAPHS_HPP
