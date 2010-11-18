#include <conan/graphs.hpp>
#define BOOST_TEST_MODULE graphs
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( undirected_graph_test )
{
  using namespace conan;
  typedef undirected_graph<adj_listS> Graph;
  typedef graph_traits< Graph > GraphTraits;
  typedef GraphTraits::vertex_descriptor vertex;
  typedef GraphTraits::edge_descriptor edge;

  // 1 //
  Graph g0;
  BOOST_CHECK_EQUAL( num_vertices(g0), (size_t)0 );
  BOOST_CHECK_EQUAL( num_edges(g0), (size_t)0 );
  BOOST_CHECK_EQUAL( is_directed(g0), false );

  vertex v0 = add_vertex("v0", g0);
  BOOST_CHECK_EQUAL( v0, (vertex)0 );
  BOOST_CHECK_EQUAL( out_degree(v0, g0), 0 );
  BOOST_CHECK_EQUAL( g0[ v0 ].name, "v0" );
  BOOST_CHECK_EQUAL( find_vertex_by_name("v0", g0), v0 );

  BOOST_CHECK_EQUAL( add_vertex("v1", g0), (vertex)1 );
  
  edge e0; bool edge_found;
  tie(e0, edge_found) = add_edge( v0, (vertex)1, 1.0, g0 );
  BOOST_CHECK_EQUAL( edge_found, true );

  BOOST_CHECK_EQUAL( source(e0, g0), v0 );
  BOOST_CHECK_EQUAL( target(e0, g0), (vertex)1 );
  BOOST_CHECK_EQUAL( g0[ e0 ].weight, 1.0 );
  BOOST_CHECK_EQUAL( distance( v0, (vertex)1, g0 ), 1.0 );

  BOOST_CHECK_EQUAL( avg_degree(g0), 1.0 );

  // 2 //
  typedef std::pair<int, int> Edge;
  Edge edge_array[] =
    { Edge(0, 1), Edge(0, 2), Edge(0, 3), Edge(0, 4),
      Edge(1, 2), Edge(1, 5), Edge(1, 3),
      Edge(2, 4), Edge(2, 5),
      Edge(3, 2) };

  size_t num_vertices_g1 = 6;
  size_t num_edges_g1 = sizeof(edge_array) / sizeof(Edge);
  Edge * first = edge_array, * last = edge_array + num_edges_g1;

  Graph g1( first, last, 1.0, num_vertices_g1 );

  BOOST_CHECK_EQUAL( num_vertices(g1), num_vertices_g1 );
  BOOST_CHECK_EQUAL( num_edges(g1), num_edges_g1 );

  size_t cnt = 0;
  BOOST_FOREACH( const vertex & v, vertices(g1) )
  {
    BOOST_CHECK( v < num_vertices_g1 );
    BOOST_CHECK_EQUAL( g1[ v ].name, "" );
    ++cnt;
  }
  BOOST_CHECK_EQUAL( cnt, num_vertices_g1 ); cnt = 0;

  BOOST_FOREACH( const edge & e, edges(g1) )
  {
    BOOST_CHECK_EQUAL( boost::edge( source(e, g1), target(e, g1), g1 ).second, true );
    g1[ e ].weight = 1.0;
    g1[ e ].type = "";
    ++cnt;
  }
  BOOST_CHECK_EQUAL( cnt, num_edges_g1 );

  BOOST_FOREACH( const edge & e, out_edges( (vertex)0, g1 ) )
  {
    BOOST_CHECK_EQUAL( boost::edge( source(e, g1), target(e, g1), g1 ).second, true );
  }

  BOOST_FOREACH( const vertex & v, adjacent_vertices( (vertex)2, g1 ) )
  {
    BOOST_CHECK( v < num_vertices_g1 );
  }

  BOOST_CHECK_EQUAL( distance( (vertex)0, (vertex)5, g1 ), 2.0 );

  // 3 //
  Graph g2(10);
  BOOST_CHECK_EQUAL( num_vertices(g2), (size_t)10 );
  BOOST_CHECK_EQUAL( num_edges(g2), (size_t)0 );
}
