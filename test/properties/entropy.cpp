#include <conan/graphs.hpp>
#include <conan/io.hpp> // DEBUG
#include <conan/graph_models.hpp>
#include <conan/properties/entropy.hpp>
#define BOOST_TEST_MODULE entropy
#include <boost/test/unit_test.hpp>
#include <conan/utils.hpp> // for get_adj_matrix

BOOST_AUTO_TEST_CASE( entropy_test )
{
  using namespace conan;
  typedef undirected_graph<adj_listS> Graph;
  typedef graph_traits< Graph > GraphTraits;
  typedef GraphTraits::vertex_descriptor vertex;
  typedef GraphTraits::edge_descriptor edge;

  // 1 //
  int N = 500;
  Graph g0( generate_erdos_renyi_graph<Graph>(N, 1) );
  decimal H_ev = graph_entropy_using_eigenvalues(g0),
          H_rw = graph_entropy_using_random_walk(g0),
          H_ss = graph_entropy_using_stationary_states(g0),
          H_th = conan_log2(N);
  BOOST_CHECK_EQUAL( H_ev, H_th );
  BOOST_CHECK_EQUAL( H_ss, H_th );
  BOOST_CHECK_EQUAL( H_rw, H_th );

  Graph g1( generate_scale_free_network<Graph>(3, 1) );
  write_dot(g1, "g.dot");

  H_ev = graph_entropy_using_eigenvalues(g1);
  H_rw = graph_entropy_using_random_walk(g1);
  H_ss = graph_entropy_using_stationary_states(g1);
  H_th = 1.5;
  BOOST_CHECK_EQUAL( H_ev, H_th );
  BOOST_CHECK_EQUAL( H_ss, H_th );
  BOOST_CHECK_EQUAL( H_rw, H_th );

}
