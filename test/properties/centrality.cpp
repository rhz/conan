#include <conan/graphs.hpp>
#include <conan/graph_models/barabasi_albert.hpp>
#include <conan/properties/centrality.hpp>
#define BOOST_TEST_MODULE centrality
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( centrality_test )
{
  using namespace conan;
  typedef undirected_graph<adj_listS> Graph;
  typedef graph_traits< Graph > GraphTraits;
  typedef GraphTraits::vertex_descriptor vertex;
  typedef GraphTraits::edge_descriptor edge;

  // 1 //
  Graph g0( generate_scale_free_network<Graph>(10, 1) );
  BOOST_FOREACH( vertex v, boost::vertices(g0) )
  {
    decimal vdc = vertex_degree_centrality(v, g0); // vdc = vertex degree centrality
    BOOST_CHECK( vdc <= 1.0 and vdc >= 0.0 );
  }
  decimal gdc = graph_degree_centrality(g0); // gdc = graph degree centrality
  BOOST_CHECK( gdc <= 1.0 and gdc >= 0.0 );

  // 2 //
  BOOST_FOREACH( vertex v, boost::vertices(g0) )
  {
    decimal closeness = vertex_closeness(v, g0);
    BOOST_CHECK( closeness >= 0.0 );
  }

  eigenvector_centrality<Graph> ec(g0);
  BOOST_FOREACH( vertex v, boost::vertices(g0) )
    BOOST_CHECK( ec(v) >= 0.0 and ec(v) <= 1.0 );

  // 3 //
  Graph g1(2);
  add_edge(0, 1, 1.0, g1);
  BOOST_FOREACH( vertex v, boost::vertices(g1) )
  {
    BOOST_CHECK_EQUAL( vertex_degree_centrality(v, g1), 1.0 );
    BOOST_CHECK_EQUAL( vertex_closeness(v, g1), 1.0 );
  }
  BOOST_CHECK_EQUAL( graph_degree_centrality(g1), 0.0 );

  // 4 //
  add_vertex(g1);
  add_edge(1, 2, 1.0, g1);
  BOOST_CHECK_EQUAL( vertex_degree_centrality(0, g1), 0.5 );
  BOOST_CHECK_EQUAL( vertex_degree_centrality(1, g1), 1.0 );
  BOOST_CHECK_EQUAL( vertex_degree_centrality(2, g1), 0.5 );
  BOOST_CHECK_EQUAL( graph_degree_centrality(g1), 1.0 );
  BOOST_CHECK_EQUAL( vertex_closeness(0, g1), 2.0 / 3.0 );
  BOOST_CHECK_EQUAL( vertex_closeness(1, g1), 1.0 );
  BOOST_CHECK_EQUAL( vertex_closeness(2, g1), 2.0 / 3.0 );

  eigenvector_centrality<Graph> ec1(g1);
  BOOST_CHECK_EQUAL( ec1(0), 0.70710678118654757 );
  BOOST_CHECK_EQUAL( ec1(1), 1.0 );
  BOOST_CHECK_EQUAL( ec1(2), 0.70710678118654757 );

  // 5 //
  add_edge(2, 0, 1.0, g1);
  eigenvector_centrality<Graph> ec2(g1);
  BOOST_FOREACH( vertex v, boost::vertices(g1) )
  {
    BOOST_CHECK_EQUAL( vertex_degree_centrality(v, g1), 1.0 );
    BOOST_CHECK_EQUAL( vertex_closeness(v, g1), 1.0 );
    BOOST_CHECK_EQUAL( ec2(v), 1.0 );
  }
}
