#include <conan/utils/combinatoria.hpp>
#define BOOST_TEST_MODULE combinatoria
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( combinatoria_test )
{
  using namespace conan::detail;
  BOOST_CHECK_EQUAL( binomial_coefficient(0, 0), 1.0 );
}
