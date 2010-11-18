#ifndef MISC_CLASSES_HPP
#define MISC_CLASSES_HPP

#include <conan/properties.hpp>
#include "vertex_and_edge_wrap.hpp"

template <class PythonGraph>
struct degree_distribution
{
  typedef typename PythonGraph::Graph Graph;
  typedef VertexWrap<Graph> Vertex;
  typedef conan::degree_distribution<Graph> CppDegreeDistribution;

  typedef typename CppDegreeDistribution::VertexList VertexList;
  typedef typename VertexList::iterator VertexListIter;

  degree_distribution(PythonGraph & graph)
    : g(graph.g), d(graph.g) { }

  degree_distribution(PythonGraph & graph, bool cumulative)
    : g(graph.g), d(graph.g, cumulative) { }

  str best_fit()
  {
    int best_fit_int = d.best_fit;
    std::string best_fit_str;
    switch(best_fit_int)
    {
    case 0:
      best_fit_str = "linear";
      break;
    case 1:
      best_fit_str = "power";
      break;
    case 2:
      best_fit_str = "exponential";
      break;
    case 3:
      best_fit_str = "gaussian";
      break;
    default:
      best_fit_str = "unknown";
      break;
    }
    return str(best_fit_str);
  }

  list vertices(size_t k)
  {
    list l;

    VertexList vl = d.vertices(k);

    for (VertexListIter vli = vl.begin(); vli != vl.end(); ++vli)
      l.append( Vertex(*vli, g) );

    return l;
  }

  conan::decimal P(size_t k)
  { return d.P(k); }

  list non_zero_degrees()
  { return list(d.non_zero_degrees()); }

  conan::decimal gamma()
  { return d.gamma(); }

  conan::decimal lambda()
  { return d.lambda(); }

  conan::decimal slope()
  { return d.slope(); }

  conan::decimal r(std::string fit_type)
  {
    size_t fit_type_int = 0;

    if (fit_type == "linear")
      fit_type_int = 0;
    else if (fit_type == "power")
      fit_type_int = 1;
    else if (fit_type == "exponential")
      fit_type_int = 2;
    else if (fit_type == "gaussian")
      fit_type_int = 3;
    else
      throw std::runtime_error("Regression type unknown");

    return d.r(fit_type_int);
  }

  conan::decimal rr(std::string fit_type)
  {
    size_t fit_type_int = 0;

    if (fit_type == "linear")
      fit_type_int = 0;
    else if (fit_type == "power")
      fit_type_int = 1;
    else if (fit_type == "exponential")
      fit_type_int = 2;
    else if (fit_type == "gaussian")
      fit_type_int = 3;
    else
      throw std::runtime_error("Regression type unknown");

    return d.rr(fit_type_int);
  }

  conan::decimal entropy()
  { return d.entropy(); }

  Graph & g;
  CppDegreeDistribution d;
};


template <class PythonGraph>
struct shortest_paths
{
  typedef typename PythonGraph::Graph Graph;
  typedef conan::shortest_paths<Graph> CppShortestPaths;
  typedef VertexWrap<Graph> Vertex;

  typedef boost::graph_traits<Graph> GraphTraits;
  typedef typename GraphTraits::vertex_descriptor vertex;

  shortest_paths(PythonGraph & graph)
    : g(graph.g), sp(graph.g) { }

  conan::decimal vertex_asp(Vertex a)
  { return sp.vertex_asp(a.v); }

  conan::decimal graph_asp()
  { return sp.graph_asp(); }

  conan::decimal vertex_eccentricity(Vertex a)
  { return sp.vertex_eccentricity(a.v); }

  conan::decimal graph_diameter()
  { return sp.graph_diameter(); }

  conan::decimal graph_radius()
  { return sp.graph_radius(); }

  list graph_center()
  {
    typedef typename std::vector<vertex> CppVertexVector;
    typedef typename CppVertexVector::iterator CppVertexVectorIter;

    list l;
    CppVertexVector vv( sp.template graph_center< CppVertexVector >() );

    for (CppVertexVectorIter vi = vv.begin(); vi != vv.end(); ++vi)
      l.append( Vertex(*vi, g) );

    return l;
  }

  Graph & g;
  CppShortestPaths sp;
};


template <class PythonGraph>
struct betweenness_centrality
{
  typedef typename PythonGraph::Graph Graph;
  typedef conan::betweenness_centrality<Graph> CppBetweennessCentrality;
  typedef VertexWrap<Graph> Vertex;
  typedef EdgeWrap<Graph> Edge;

  typedef boost::graph_traits<Graph> GraphTraits;
  typedef typename GraphTraits::vertex_descriptor vertex;

  betweenness_centrality(PythonGraph & graph)
    : bc(graph.g) { }

  conan::decimal vertex_centrality(Vertex a)
  { return bc.vertex_centrality(a.v); }

  conan::decimal vertex_relative_centrality(Vertex a)
  { return bc.vertex_relative_centrality(a.v); }

  conan::decimal edge_centrality(Edge e)
  { return bc.edge_centrality(e.e); }

  conan::decimal central_point_dominance()
  { return bc.central_point_dominance(); }

  CppBetweennessCentrality bc;
};


template <class PythonGraph>
struct eigenvector_centrality
{
  typedef typename PythonGraph::Graph Graph;
  typedef conan::eigenvector_centrality<Graph> CppEigenvectorCentrality;
  typedef VertexWrap<Graph> Vertex;
  typedef EdgeWrap<Graph> Edge;

  typedef boost::graph_traits<Graph> GraphTraits;
  typedef typename GraphTraits::vertex_descriptor vertex;

  eigenvector_centrality(PythonGraph & graph)
    : ec(graph.g) { }

  conan::decimal vertex_eigenvector_centrality(Vertex a)
  { return ec(a.v); }

  CppEigenvectorCentrality ec;
};

#endif // MISC_CLASSES_HPP
