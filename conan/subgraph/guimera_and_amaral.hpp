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

#ifndef GUIMERA_AND_AMARAL_HPP
#define GUIMERA_AND_AMARAL_HPP
#include <conan/graphs.hpp>
#include <conan/utils/sa.hpp>
#include <conan/utils/stacktrace.hpp>
#include <conan/subgraph/modularity.hpp>

namespace conan {

  typedef std::vector<size_t> Module;
  typedef std::vector<Module> ModuleVector;

  namespace detail {

    struct ModulesAndGraph
    {
      typedef undirected_graph<> Graph;

      ModulesAndGraph(
          const ModuleVector & m,
          const Graph & g
          )
      : modules(m), graph_ptr(&g)
      { L = boost::num_edges(*graph_ptr); }

      ModulesAndGraph(
          const ModulesAndGraph & other
          )
      : modules(other.modules), graph_ptr(other.graph_ptr), L(other.L) { }

      ModulesAndGraph & operator=(
          const ModulesAndGraph & other
          )
      {
        if (this != &other)
        {
          modules = other.modules;
          graph_ptr = other.graph_ptr;
          L = other.L;
        }
        return *this;
      }

      bool operator==(
          const ModulesAndGraph & other
          )
      {
        if (graph_ptr != other.graph_ptr)
          return false;
        else if (modules != other.modules)
          return false;

        return true;
      }

      const Graph & graph() const
      { return *graph_ptr; }

#if defined CONAN_INFO || !defined CONAN_NON_INTERACTIVE
      std::string tostring()
      {
        std::string s = "||";
        BOOST_FOREACH( Module m, modules )
        {
          BOOST_FOREACH( size_t v, m )
          {
            s += " " + to_string(v);
          }
          s += " ||";
        }
        return s;
      }
#endif

      ModuleVector modules;
      const Graph * graph_ptr;
      size_t L; // for optimization
    };


    double modularity_cost_fcn(const ModulesAndGraph & m)
    { return -modularity( m.graph(), m.modules, m.L ); }


    template <int MaxTriesMod> // = 2>
    ModulesAndGraph individual_node_movement(
        const ModulesAndGraph & m,
        const SA<ModulesAndGraph> & caller_sa
        )
    {
      using namespace detail;
      typedef SA<ModulesAndGraph> SA;

      if ( m.modules.size() == 1 )
      {
#ifdef CONAN_DEBUG
        std::cerr << "conan::detail::individual_node_movement: There is only one module. Nothing to do." << std::endl;
#endif
        return ModulesAndGraph( m );
      }

      ModulesAndGraph nb( m );

      // Choose modules
      SA::RNG & rng = caller_sa.rng;
      boost::uniform_int<size_t> index_range( 0, nb.modules.size() - 1 );
      boost::variate_generator< SA::RNG &, boost::uniform_int<size_t> >
        random_index( rng, index_range );

      size_t i1 = 0, i2 = 0, max_tries = nb.modules.size() * MaxTriesMod;
      while ( i1 == i2 || ( nb.modules[ i1 ].size() == 1 && max_tries > 0 ) )
      { i1 = random_index(); i2 = random_index(); --max_tries; }

      double M_i1 = module_modularity( nb.graph(), nb.modules[ i1 ], nb.L ),
             M_i2 = module_modularity( nb.graph(), nb.modules[ i2 ], nb.L );

      // Choose vertex
      boost::uniform_int<size_t> index_subrange_i1( 0, nb.modules[ i1 ].size() - 1 );
      boost::variate_generator< SA::RNG &, boost::uniform_int<size_t> >
        random_subindex_i1( rng, index_subrange_i1 );

      size_t j1 = random_subindex_i1();
      bool module_removed = false;

      // Move vertex
      nb.modules[ i2 ].push_back( nb.modules[ i1 ][ j1 ] );
      if ( nb.modules[ i1 ].size() == 1 )
      {
        nb.modules.erase( nb.modules.begin() + i1 );
        module_removed = true;
        if ( i2 > i1 )
          --i2;
      }
      else
      {
        nb.modules[ i1 ].erase( nb.modules[ i1 ].begin() + j1 );
      }

#ifdef CONAN_EXTRA_DEBUG
      std::cerr << "conan::detail::individual_node_movement: " << nb.tostring() << std::endl;
#endif

      if ( caller_sa.neighbour_fcn_computes_energy )
      {
        double M_i2_new = module_modularity( nb.graph(), nb.modules[ i2 ], nb.L ),
               M_i1_new = ( module_removed? 0.0 : module_modularity( nb.graph(), nb.modules[ i1 ], nb.L ) );
        caller_sa.neighbour_energy = caller_sa.current_energy + M_i1 + M_i2 - M_i2_new - M_i1_new;
      }

      return nb;
    }


    template <int MaxTriesMod> // = 10>
    ModulesAndGraph split_modules(
        const ModulesAndGraph & m,
        const SA<ModulesAndGraph> & caller_sa
        )
    {
      using namespace detail;
      typedef neighbour_fcn<ModulesAndGraph> NeighbourFcn;
      typedef SA<ModulesAndGraph> SA;

      ModulesAndGraph nb( m );

      // Choose module
      SA::RNG & rng = caller_sa.rng;
      boost::uniform_int<size_t> index_range( 0, nb.modules.size() - 1 );
      boost::variate_generator< SA::RNG &, boost::uniform_int<size_t> >
        random_index( rng, index_range );

      size_t i = 0, max_tries = nb.modules.size() * MaxTriesMod;
      
      do
      {
        i = random_index();

        if ( ! --max_tries ) // if there are so many modules with just one node inside them and we cannot find one suitable for splitting
        { // then we should escape from this function
#ifdef CONAN_DEBUG
          std::cerr << "conan::detail::split_modules: Could not find any suitable module for splitting, ie., a module with more than one node inside it."
                    << std::endl;
#endif
          caller_sa.neighbour_energy = caller_sa.current_energy;
          return nb;
        }
      }
      while ( nb.modules[ i ].size() == 1 );

      size_t module_size = nb.modules[ i ].size();

      Module m1, m2;

      if ( module_size == 2 )
      { // Shortcut
        m1.push_back( nb.modules[ i ][ 0 ] );
        m2.push_back( nb.modules[ i ][ 1 ] );
      
        nb.modules.push_back( m1 );
        nb.modules.push_back( m2 );
        nb.modules.erase( nb.modules.begin() + i );

#ifdef CONAN_EXTRA_DEBUG
        std::cerr << "conan::detail::split_modules: " << nb.tostring() << std::endl;
#endif

        if ( caller_sa.neighbour_fcn_computes_energy )
        { caller_sa.neighbour_energy = caller_sa.E( nb ); }

        return nb;
      }

      size_t cnt = 0;
      BOOST_FOREACH( size_t v, nb.modules[ i ] )
      {
        if ( ++cnt % 2 )
          m1.push_back( v );
        else
          m2.push_back( v );
      }
      
      ModuleVector new_module_vector;
      new_module_vector.push_back( m1 );
      new_module_vector.push_back( m2 );

      ModulesAndGraph new_modules( new_module_vector, nb.graph() );

#ifdef CONAN_EXTRA_DEBUG
      std::cerr << "conan::detail::split_modules: Beginning nested SA:"
                << " initial temp = " << caller_sa.get_initial_temperature() << ","
                << " final temp = " << caller_sa.get_current_temperature() << std::endl;
#endif

      SA nested_sa(
        temperature_schedule( caller_sa.get_initial_temperature(), caller_sa.get_current_temperature() ),
        caller_sa.E,
        caller_sa.P,
        new_modules,
        rng
        );

      nested_sa.add_neighbour_fcn( NeighbourFcn( &individual_node_movement<2>, 2.0 * module_size ) );

#ifdef CONAN_EXTRA_DEBUG
      std::cerr << "conan::detail::split_modules: Running nested SA." << std::endl
                << "conan::detail::split_modules: initial state: " << new_modules.tostring() << std::endl;
#endif

      nested_sa.neighbour_fcn_computes_energy = caller_sa.neighbour_fcn_computes_energy;

      nested_sa.run();

#ifdef CONAN_EXTRA_DEBUG
      std::cerr << "conan::detail::split_modules: Nested SA finished: best state: " << nested_sa.best_state.tostring() << std::endl;
#endif

      // Remove the old module and add the two new ones
      nb.modules.erase( nb.modules.begin() + i );

      nb.modules.push_back( nested_sa.best_state.modules[ 0 ] );
      if ( nested_sa.best_state.modules.size() > 1 )
        nb.modules.push_back( nested_sa.best_state.modules[ 1 ] );

#ifdef CONAN_EXTRA_DEBUG
      std::cerr << "conan::detail::split_modules: " << nb.tostring() << std::endl;
#endif

      // Compute new energy
      if ( caller_sa.neighbour_fcn_computes_energy )
      {
        if ( nested_sa.best_state.modules.size() == 2 )
          caller_sa.neighbour_energy = caller_sa.E( nb );
        else
          caller_sa.neighbour_energy = caller_sa.current_energy;
      }

      return nb;
    }


    ModulesAndGraph merge_modules(
        const ModulesAndGraph & m,
        const SA<ModulesAndGraph> & caller_sa
        )
    {
      using namespace detail;
      typedef SA<ModulesAndGraph> SA;

      if ( m.modules.size() == 1 )
      {
#ifdef CONAN_DEBUG
        std::cerr << "conan::detail::merge_modules: There is only one module. Nothing to do." << std::endl;
#endif
        return ModulesAndGraph( m );
      }

      ModulesAndGraph nb( m );

      // Choose modules
      SA::RNG & rng = caller_sa.rng;
      boost::uniform_int<size_t> index_range( 0, nb.modules.size() - 1 );
      boost::variate_generator< SA::RNG &, boost::uniform_int<size_t> >
        random_index( rng, index_range );

      size_t i1 = 0, i2 = 0;
      while (i1 == i2)
      {
        i1 = random_index();
        i2 = random_index();
      }

      Module & m1 = nb.modules[ i1 ],
             & m2 = nb.modules[ i2 ];

      Module new_module( m1 );
      BOOST_FOREACH( size_t v, m2 )
      {
        new_module.push_back( v );
      }

      // Compute new energy
      if ( caller_sa.neighbour_fcn_computes_energy )
      {
        caller_sa.neighbour_energy = caller_sa.current_energy + module_modularity( nb.graph(), m1, nb.L ) + module_modularity( nb.graph(), m2, nb.L ) -
          module_modularity( nb.graph(), new_module, nb.L );
      }

      // Remove old modules
      nb.modules.erase( nb.modules.begin() + i1 );
      if ( i1 < i2 )
        nb.modules.erase( nb.modules.begin() + i2 - 1 );
      else
        nb.modules.erase( nb.modules.begin() + i2 );

      // Add the new one
      nb.modules.push_back( new_module );

#ifdef CONAN_EXTRA_DEBUG
      std::cerr << "conan::detail::merge_modules: " << nb.tostring() << std::endl;
#endif

      return nb;
    }


    ModulesAndGraph collective_movement(
        const ModulesAndGraph & m,
        const SA<ModulesAndGraph> & caller_sa
        )
    {
      typedef SA<ModulesAndGraph> SA;

      // Choose method
      boost::uniform_int<size_t> range( 0, 1 );
      boost::variate_generator< SA::RNG &, boost::uniform_int<size_t> >
        random( caller_sa.rng, range );

      size_t method = random();

      if ( method == 0 )
        return split_modules<10>( m, caller_sa );
      else
        return merge_modules( m, caller_sa );
    }

  } // detail

  namespace globals {

    conan::SA<conan::detail::ModulesAndGraph> * running_sa_ptr;

  } // globals


#ifdef __linux__
  void print_current_stack_and_sa_state( int signum )
  {
    std::cerr << "Simulated annealing interrupted..." << std::endl;
#ifdef CONAN_INFO
    // Print SA state
    std::cerr << "[sa] Last state: " << conan::globals::running_sa_ptr->current_state.tostring() << std::endl;
#endif
#ifdef PRINT_STACKTRACE
    // Print stack trace
    std::cerr << "[bt] Execution path:" << std::endl;
    show_stackframe();
#endif
    conan::globals::running_sa_ptr->stop();
    return;
  }
#endif


  template <class Graph>
  void sa_communities(
      Graph & g,
      std::vector< std::vector<size_t> > & modules,
      SA<detail::ModulesAndGraph>::EnergyFcn E,
      bool nb_fcn_computes_energy = false,
      double initial_temp = .0001,
      double final_temp = .00006,
      double f = 1.0
      )
  {
    using namespace detail;
    typedef neighbour_fcn<ModulesAndGraph> NeighbourFcn;
    typedef SA<ModulesAndGraph> SA;

    size_t S = boost::num_vertices( g );

    if ( modules.size() == 0 )
    {
      size_t cnt = 1, module_size = ceil( log(S) );
      Module m;

      BOOST_FOREACH( size_t vertex, boost::vertices(g) )
      {
        if ( cnt % module_size )
        {
          m.push_back( vertex );
        }
        else
        {
          m.push_back( vertex );
#ifdef CONAN_INFO
          std::cout << "vertices ";
          BOOST_FOREACH( size_t v, m )
            std::cout << v << " ";
          std::cout << "-> module " << cnt / module_size << std::endl;
#endif
          modules.push_back( m );
          m.clear();
        }

        ++cnt;
      }

      if ( ! m.empty() )
      {
#ifdef CONAN_INFO
        std::cout << "vertices ";
        BOOST_FOREACH( size_t v, m )
          std::cout << v << " ";
        std::cout << "-> module " << cnt / module_size << std::endl;
#endif
        modules.push_back( m );
      }
    }

    ModulesAndGraph m( modules, g );

    SA sa(
      temperature_schedule(initial_temp, final_temp),
      E,
      &default_acceptance_probability,
      m,
      size_t( time(NULL) ),
      true // DEBUG
      );

    sa.add_neighbour_fcn( NeighbourFcn( &individual_node_movement<2>, S * S * f ) );
    sa.add_neighbour_fcn( NeighbourFcn( &collective_movement, S * f ) );

    sa.neighbour_fcn_computes_energy = nb_fcn_computes_energy;

#ifdef __linux__
    struct sigaction sigact, prev_sigact;
    sigact.sa_handler = print_current_stack_and_sa_state;
    sigaction( SIGINT, &sigact, &prev_sigact );
    conan::globals::running_sa_ptr = &sa;
#endif

    sa.run(); // print summary and statistics

#ifdef __linux__
    sigaction( SIGINT, &prev_sigact, NULL );
#endif

    modules = sa.best_state.modules;

    return;
  }


  template <class Graph>
  inline void guimera_amaral_communities(
      Graph & g,
      std::vector< std::vector<size_t> > & modules,
      double initial_temp = .0001,
      double final_temp = .00006,
      double f = 1.0
      )
  { return sa_communities(g, modules, &detail::modularity_cost_fcn, true, initial_temp, final_temp, f); }

} // conan

#endif // GUIMERA_AND_AMARAL_HPP
