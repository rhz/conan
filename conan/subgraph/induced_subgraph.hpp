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

#ifndef INDUCED_SUBGRAPH_HPP
#define INDUCED_SUBGRAPH_HPP
#include <conan/graphs.hpp>
#include <limits>
#include <utility> // for std::pair and std::make_pair

namespace conan {

  /**
   * @brief Class designed for making successives/incremental subgraphs from a graph.
   * The graph and the selected vertices can change between one call to
   * incremental_induced_subgraph::make_graph() and the next.
   */
  template <class Graph>
  class incremental_induced_subgraph
  {
    //const int null_edge = 2 * std::numeric_limits<int>::max() - 1;

  public:
    typedef typename boost::graph_traits<Graph> Traits;
    
    typedef typename Traits::vertex_descriptor   vertex_descriptor; // 0-based unsigned int
    typedef typename Traits::edge_descriptor     edge_descriptor;
    typedef typename Traits::vertices_size_type  vertices_size_type;
    typedef typename Traits::edges_size_type     edges_size_type;

    typedef typename std::vector<vertex_descriptor> VertexList;
    typedef typename std::vector<edge_descriptor> EdgeList;

    typedef typename VertexList::iterator vertex_iterator;
    typedef typename VertexList::iterator adjacency_iterator;
    typedef typename EdgeList::iterator edge_iterator;
    typedef typename EdgeList::iterator out_edge_iterator;
    typedef typename EdgeList::iterator in_edge_iterator;

    // These are for compatibility with boost::graph_traits only
    typedef typename Traits::directed_category directed_category;
    typedef typename Traits::edge_parallel_category edge_parallel_category;
    typedef typename Traits::traversal_category traversal_category;
    typedef typename Traits::degree_size_type degree_size_type;

    // Properties
    typedef typename Graph::vertex_property_type vertex_property_type;
    typedef typename Graph::edge_property_type edge_property_type;

    // The types that are actually bundled
    typedef typename Graph::vertex_bundled vertex_bundled;
    typedef typename Graph::edge_bundled edge_bundled;

  protected:
    Graph & _parent_graph;
    VertexList _vertex_list;
    mutable EdgeList _tmp_edge_list;

  public:
    // Ctor.
    incremental_induced_subgraph(Graph & g)
      : _parent_graph(g) { }

  protected:
    void construct_edge_list();

  public:
    inline Graph parent_graph()
    { return _parent_graph; }

    inline Graph parent_graph() const
    { return _parent_graph; }

    vertex_descriptor add_vertex(vertex_descriptor v)
    { _vertex_list.push_back(v); return v; }

    void remove_vertex(vertex_descriptor v)
    { }

    template <class OutGraph>
    OutGraph make_graph();

    vertices_size_type num_vertices() const
    { return _vertex_list.size(); }
  };

} // conan

namespace boost {

  template <typename G>
  typename conan::incremental_induced_subgraph<G>::vertex_descriptor
  add_vertex(typename conan::incremental_induced_subgraph<G>::vertex_descriptor v,
             conan::incremental_induced_subgraph<G> & g)
  { return g.add_vertex( v ); }

  template <typename G>
  typename conan::incremental_induced_subgraph<G>::vertex_descriptor
  remove_vertex(typename conan::incremental_induced_subgraph<G>::vertex_descriptor v,
                conan::incremental_induced_subgraph<G> & g)
  { return g.remove_vertex( v ); }

} // boost

namespace conan {

  template <class Graph>
  void incremental_induced_subgraph<Graph>::construct_edge_list()
  {
    if (_vertex_list.empty())
      return;

    _tmp_edge_list.clear(); // clear edge list

    bool undirected = boost::is_undirected(_parent_graph);

    vertex_iterator source = _vertex_list.begin(),
                    list_end = _vertex_list.end();

    size_t cnt = 1;
    for (; source != list_end - 1; ++source, ++cnt)
    {
      vertex_iterator target = _vertex_list.begin() + (undirected? cnt : 0);

      for (; target != list_end; ++target)
      {
        if (source == target)
          continue;

        if (boost::edge(*source, *target, _parent_graph).second)
          _tmp_edge_list.push_back( boost::edge(*source, *target, _parent_graph).first );
      }
    }

    return;
  }

  template <class Graph>
  template <class OutGraph>
  OutGraph incremental_induced_subgraph<Graph>::make_graph()
  {
    OutGraph g(_vertex_list.size());

    std::map<vertex_descriptor, typename OutGraph::vertex_descriptor> vertex_index_map;

    // Copy vertex properties
    vertex_iterator vi = _vertex_list.begin(),
                    vi_end = _vertex_list.end();

    size_t vertex_index = 0;

    for (; vi != vi_end; ++vi)
    {
      g[vertex_index].name = _parent_graph[*vi].name;
      vertex_index_map[*vi] = vertex_index++;
    }

    // Create edges and copy their properties
    construct_edge_list();

    edge_iterator ei = _tmp_edge_list.begin(),
                  ei_end = _tmp_edge_list.end();

    for (; ei != ei_end; ++ei)
    {
      typename OutGraph::edge_descriptor e =
          boost::add_edge(vertex_index_map[boost::source(*ei, _parent_graph)],
                          vertex_index_map[boost::target(*ei, _parent_graph)],
                          g).first;
      g[e].weight = _parent_graph[*ei].weight;
    }

    return g;
  }

  // End of incremental_induced_subgraph member functions


  /**
   * @brief Handy wrapper for making just one subgraph from a graph.
   *
   * @param g Pointer to a graph.
   * @param vertices STL list (or vector) of vertex descriptors. Usually you must specify
   *   the template type as std::vector< conan::graph_traits< conan::undirected_graph<> >::vertex_descriptor >.
   * @return The subgraph.
   */
  template <class Graph, class List>
  Graph induced_subgraph(
      Graph & g,
      List vertices // elements of this list must be of type boost::graph_traits<Graph>::vertex_descriptor
      )
  {
    typedef typename List::iterator ListIter;

    incremental_induced_subgraph<Graph> subg(g);

    for (ListIter li = vertices.begin(), liend = vertices.end(); li != liend; ++li)
      subg.add_vertex(*li);

    return subg.template make_graph<Graph>();
  }

  /**
   * @brief Handy wrapper for making just one subgraph from a graph.
   *
   * @param g Pointer to a graph.
   * @param vertex_descriptors C array of vertex descriptors.
   * @param n Size of the C array.
   * @return The subgraph.
   */
  template <class Graph>
  Graph induced_subgraph(
      Graph & g,
      typename boost::graph_traits<Graph>::vertex_descriptor * vertices, // this is actually a plain array
      size_t n
      )
  {
    incremental_induced_subgraph<Graph> subg(g);

    for (size_t i = 0; i < n; ++i)
      subg.add_vertex(vertices[i]);

    return subg.template make_graph<Graph>();
  }

} // conan

#endif //INDUCED_SUBGRAPH_HPP
