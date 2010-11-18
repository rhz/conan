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

#ifndef EDGE_SET_HPP
#define EDGE_SET_HPP
#include <conan/graphs.hpp>
#include <utility> // for std::pair and std::make_pair

namespace conan {

  template <class Graph>
  class incremental_subgraph_from_edge_set
  {
  public:
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::edge_descriptor edge_descriptor;
    typedef typename GraphTraits::vertex_descriptor vertex_descriptor;

    typedef typename std::set<vertex_descriptor> VertexList;
    typedef typename std::set<edge_descriptor> EdgeList;

  protected:
    Graph & _parent_graph;
    mutable VertexList _tmp_vertex_set;
    EdgeList _edge_set;

  public:
    // Ctor.
    incremental_subgraph_from_edge_set(Graph & g)
      : _parent_graph(g) { }

  protected:
    void construct_vertex_set();

  public:
    inline Graph parent_graph()
    { return _parent_graph; }

    inline Graph parent_graph() const
    { return _parent_graph; }

    std::pair<edge_descriptor, bool> add_edge(edge_descriptor e)
    { _edge_set.insert( e ); return std::make_pair( e, true ); }

    void remove_edge(edge_descriptor e)
    { _edge_set.remove( e ); return; }

    template <class OutGraph>
    OutGraph make_graph();

    size_t num_edges() const
    { return _edge_set.size(); }
  };

} // conan

namespace boost {

  template <typename G>
  typename conan::incremental_subgraph_from_edge_set<G>::vertex_descriptor
  add_edge(typename conan::incremental_subgraph_from_edge_set<G>::edge_descriptor e,
             conan::incremental_subgraph_from_edge_set<G> & g)
  { return g.add_edge( e ); }

  template <typename G>
  typename conan::incremental_subgraph_from_edge_set<G>::vertex_descriptor
  remove_edge(typename conan::incremental_subgraph_from_edge_set<G>::edge_descriptor e,
              conan::incremental_subgraph_from_edge_set<G> & g)
  { return g.remove_edge( e ); }

} // boost

namespace conan {

  template <class Graph>
  void incremental_subgraph_from_edge_set<Graph>::construct_vertex_set()
  {
    typedef typename EdgeList::iterator edge_iter;

    _tmp_vertex_set.clear(); // clear vertex set

    for (edge_iter ei = _edge_set.begin(), eiend = _edge_set.end(); ei != eiend; ++ei)
    {
      _tmp_vertex_set.insert( boost::source( *ei, _parent_graph ) );
      _tmp_vertex_set.insert( boost::target( *ei, _parent_graph ) );
    }

    return;
  }


  template <class Graph>
  template <class OutGraph>
  OutGraph incremental_subgraph_from_edge_set<Graph>::make_graph()
  {
    typedef typename EdgeList::iterator edge_iter;
    typedef typename VertexList::iterator vertex_iter;

    construct_vertex_set();

    OutGraph g(_tmp_vertex_set.size());

    std::map<vertex_descriptor, typename OutGraph::vertex_descriptor> vertex_index_map;

    // Copy vertex properties
    size_t vertex_index = 0;

    for (vertex_iter vi = _tmp_vertex_set.begin(), viend = _tmp_vertex_set.end(); vi != viend; ++vi)
    {
      g[ vertex_index ].name = _parent_graph[ *vi ].name;
      vertex_index_map[ *vi ] = vertex_index++;
    }

    // Create edges and copy their properties
    for (edge_iter ei = _edge_set.begin(), eiend = _edge_set.end(); ei != eiend; ++ei)
    {
      typename OutGraph::edge_descriptor e =
          boost::add_edge( vertex_index_map[ boost::source( *ei, _parent_graph ) ],
                           vertex_index_map[ boost::target( *ei, _parent_graph ) ],
                           g ).first;
      g[ e ].weight = _parent_graph[ *ei ].weight;
    }

    return g;
  }
  

  template <class Graph, class List>
  Graph make_subgraph_from_edge_set(
      Graph & g,
      List edges // elements of this list must be of type boost::graph_traits<Graph>::edge_descriptor
      )
  {
    typedef typename List::iterator ListIter;

    incremental_subgraph_from_edge_set<Graph> subg(g);

    for (ListIter li = edges.begin(), liend = edges.end(); li != liend; ++li)
      subg.add_edge(*li);

    return subg.template make_graph<Graph>();
  }


  template <class Graph>
  Graph make_subgraph_from_edge_set(
      Graph & g,
      typename boost::graph_traits<Graph>::edge_descriptor * edges, // this is actually a plain array
      size_t n
      )
  {
    incremental_subgraph_from_edge_set<Graph> subg(g);

    for (size_t i = 0; i < n; ++i)
      subg.add_edge(edges[i]);

    return subg.template make_graph<Graph>();
  }

} // conan

#endif // EDGE_SET_HPP
