#ifndef PYTHON_DYNAMICS_HPP
#define PYTHON_DYNAMICS_HPP

#include "python_undirected_graph.hpp"
#include <conan/dynamics/ising.hpp>

// ising struct
struct ising
{
  typedef conan::undirected_graph<conan::adj_listS> Graph;
  typedef conan::ising<Graph> CppIsing;

  ising(python_undirected_graph g)
    : i(g.g) { }

  void equilibrate_system(
      conan::decimal T,
      size_t num_steps = 50000
      )
  { return i.equilibrate_system(T, num_steps); }

  size_t num_edges(
      size_t index = 0
      )
  { return i.num_edges(index); }

  conan::decimal magnetization(
      size_t index = 0
      )
  { return i.magnetization(index); }

  conan::decimal energy(
      size_t index = 0
      )
  { return i.energy(index); }

  python_undirected_graph run_montecarlo(
      conan::decimal T,
      size_t num_steps = 50000
      )
  { return python_undirected_graph(i.run_montecarlo(T, num_steps)); }

  list simulate(
      conan::decimal T_init,
      conan::decimal T_final,
      conan::decimal deltaT,
      size_t num_steps = 50000,
      size_t num_steps_equilibrium = 50000
      )
  {
    std::list<python_undirected_graph> graph_list;
    std::vector<Graph> graph_vector(i.simulate(T_init, T_final, deltaT, num_steps, num_steps_equilibrium));
    for (std::vector<Graph>::iterator veci = graph_vector.begin(); veci != graph_vector.end(); ++veci)
      graph_list.push_back(python_undirected_graph(*veci));
    return list(graph_list);
  }

  CppIsing i;
};

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(equilibrate_system_overloads, equilibrate_system, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(num_edges_overloads, num_edges, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(energy_overloads, energy, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(magnetization_overloads, magnetization, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(run_montecarlo_overloads, run_montecarlo, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(simulate_overloads, simulate, 3, 5)


BOOST_PYTHON_MODULE(_dynamics)
{
  class_< ising >("ising", init<python_undirected_graph>())
    .def("equilibrate_system", &ising::equilibrate_system, equilibrate_system_overloads())
    .def("num_edges", &ising::num_edges, num_edges_overloads())
    .def("energy", &ising::energy, energy_overloads())
    .def("magnetization", &ising::magnetization, magnetization_overloads())
    .def("run_montecarlo", &ising::run_montecarlo, run_montecarlo_overloads())
    .def("simulate", &ising::simulate, simulate_overloads())
    ;
}

#endif //PYTHON_DYNAMICS_HPP
