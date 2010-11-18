#ifndef VERTEX_AND_EDGE_WRAP_HPP
#define VERTEX_AND_EDGE_WRAP_HPP

#include <conan/graphs.hpp>
#include <conan/properties.hpp>

#include <boost/python.hpp>

using namespace boost::python;


template <class Graph>
struct EdgeWrap;


// VertexWrap class
template <class Graph>
struct VertexWrap
{
  typedef typename Graph::vertex_descriptor vertex_descriptor;
  typedef typename Graph::vertex_iterator vertex_iter;

  VertexWrap(vertex_descriptor vd, Graph & g)
    : v(vd), parent_graph(g)
  {
    bool vertex_found = false;
    vertex_iter vi, viend;
    for (tie(vi, viend) = boost::vertices(parent_graph); vi != viend; ++vi)
    {
      if (v == *vi)
        vertex_found = true;
    }

    if (not vertex_found)
    {
      throw std::runtime_error("vertex " + conan::to_string(v) + " not found");
    }
  }

  vertex_descriptor get_id()
  { return v; }

  std::string get_name()
  { return parent_graph[v].name; }

  void set_name(std::string name)
  { parent_graph[v].name = name; }

  size_t get_degree()
  { return boost::out_degree(v, parent_graph); }

  conan::decimal get_clustering()
  { return conan::vertex_clustering(v, parent_graph); }

  conan::decimal get_asp()
  { return conan::vertex_avg_shortest_path(v, parent_graph); }

  conan::decimal get_closeness()
  { return conan::vertex_closeness(v, parent_graph); }

  conan::decimal get_degree_centrality()
  { return conan::vertex_degree_centrality(v, parent_graph); }

  list get_incident_edges()
  {
    typedef typename Graph::out_edge_iterator out_edge_iterator;

    list l;

    out_edge_iterator oei, oei_end;
    for (tie(oei, oei_end) = boost::out_edges(v, parent_graph); oei != oei_end; ++oei)
      l.append( EdgeWrap<Graph>(*oei, parent_graph) );

    return l;
  }

  list get_adjacent_vertices()
  {
    typedef typename Graph::adjacency_iterator adjacency_iterator;

    list l;

    adjacency_iterator ai, aiend;
    for (tie(ai, aiend) = boost::adjacent_vertices(v, parent_graph); ai != aiend; ++ai)
      l.append( VertexWrap<Graph>(*ai, parent_graph) );

    return l;
  }

  str get_str()
  { return str( "(" + str(v) + ")" ); }

  vertex_descriptor v;
  Graph & parent_graph;
};


template <class Graph>
struct EdgeWrap
{
  typedef typename Graph::edge_descriptor edge_descriptor;
  typedef typename Graph::edge_iterator edge_iter;
  typedef typename Graph::vertex_descriptor vertex_descriptor;

  EdgeWrap(
      edge_descriptor ed,
      Graph & g
      )
    : e(ed), parent_graph(g)
  {
    bool edge_found = false;
    edge_iter ei, eiend;
    for (tie(ei, eiend) = boost::edges(g); ei != eiend; ++ei)
      if (e == *ei)
        edge_found = true;

    if (not edge_found)
      throw std::runtime_error("edge not found");
  }

  EdgeWrap(
      vertex_descriptor a,
      vertex_descriptor b,
      Graph & g
      )
    : parent_graph(g)
  {
    if (boost::edge(a, b, g).second)
      e = boost::edge(a, b, g).first;
    else
      throw std::runtime_error("edge between " + conan::to_string(a) + " and " + conan::to_string(b) + " does not exist");
  }

  edge_descriptor get_id()
  { return e; }

  conan::decimal get_weight()
  { return parent_graph[e].weight; }

  void set_weight(conan::decimal w)
  { parent_graph[e].weight = w; }

  std::string get_type()
  { return parent_graph[e].type; }

  void set_type(std::string t)
  { parent_graph[e].type = t; }

  VertexWrap<Graph> get_source_vertex()
  { return VertexWrap<Graph>( boost::source(e, parent_graph), parent_graph ); }

  VertexWrap<Graph> get_target_vertex()
  { return VertexWrap<Graph>( boost::target(e, parent_graph), parent_graph ); }

  str get_str()
  {
    str s = "(";
    s += str( boost::source(e, parent_graph) );

    if (boost::is_directed(parent_graph))
      s += " -> ";
    else
      s += " -- ";

    s += str( boost::target(e, parent_graph) );
    s += ")";

    return s;
  }

  edge_descriptor e;
  Graph & parent_graph;
};

#endif
