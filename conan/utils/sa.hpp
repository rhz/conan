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

#ifndef SA_HPP
#define SA_HPP
#include <conan/config.hpp>
#include <conan/utils/choice.hpp>
#ifdef HAVE_TERMIOS
  #include <conan/utils/kbhit.hpp>
#endif
#include <cmath>
#include <ctime>

namespace conan {

  template <class State, class _RNG>
  struct SA;

  template <class State, class RNG = boost::mt19937>
  struct neighbour_fcn
  {
    typedef State ( * NeighbourGenerator ) ( const State &, const SA<State, RNG> & );

    size_t num_calls_per_iteration;
    size_t remaining_calls_cnt;
    NeighbourGenerator neighbour_generator;

    neighbour_fcn(
        NeighbourGenerator neighbour_gen,
        size_t n
        )
      : num_calls_per_iteration(n), remaining_calls_cnt(n), neighbour_generator(neighbour_gen)
    { }

    State operator()(
        const State & s,
        const SA<State, RNG> & caller_sa
        )
    {
      --remaining_calls_cnt;
      return neighbour_generator(s, caller_sa);
    }

    size_t remaining_calls()
    { return remaining_calls_cnt; }

    void reset_counter()
    {
      remaining_calls_cnt = num_calls_per_iteration;
      return;
    }
  };

  namespace detail {

    double default_next_temperature(
        double T
        )
    { return 0.995 * T; }

  } // detail

  struct temperature_schedule
  {
    double initial_T;
    double ( * next_T_fcn ) ( double T );
    double final_T;
    double current_T;
    double next_T;

    temperature_schedule(
        double Tinitial,
        double Tfinal,
        double ( * Tnext ) ( double T ) = &detail::default_next_temperature
        )
      : initial_T(Tinitial), next_T_fcn(Tnext), final_T(Tfinal), current_T(Tinitial), next_T(Tinitial)
    { }

    double next_temperature()
    {
      current_T = next_T;
      next_T = next_T_fcn( current_T );
      return current_T;
    }

    void reset_schedule()
    {
      current_T = initial_T;
      return;
    }

    bool finished()
    {
      return (current_T < final_T);
    }
  };

  template <class State, class _RNG = boost::mt19937>
  struct SA
  {
    // Typedefs
    typedef _RNG RNG;
    typedef neighbour_fcn<State, RNG>
      NeighbourFcn;
    typedef typename std::vector<NeighbourFcn>
      NeighbourFcnSet;
    typedef double ( * EnergyFcn ) ( const State & );
    typedef double ( * ProbFcn ) ( double, double, double );

    // Data members
    NeighbourFcnSet neighbour_fcn_set;
    temperature_schedule temp_sched;
    EnergyFcn E;
    ProbFcn P;
    State current_state;
    double current_energy;
    State best_state;
    double best_energy;
    bool neighbour_fcn_computes_energy;
    mutable double neighbour_energy;

    // For statistics
    size_t accepted;
    size_t total;
    size_t non_changed;

    // For summary
    size_t total_sa_iterations;
    time_t start, end;

    // RNG
    mutable RNG rng;
    boost::uniform_real<double> prob_range;
    mutable boost::variate_generator< RNG &, boost::uniform_real<double> > random;

    bool paused;
    bool debug;
    size_t responsiveness_level;

    // Ctors.
    SA(
        const temperature_schedule & ts,
        const EnergyFcn & e,
        const ProbFcn & p,
        const State & initial_state,
#ifdef CONAN_DEBUG
        bool _debug = true,
#else
        bool _debug = false,
#endif
        size_t _responsiveness_level = 100
      )
      : temp_sched(ts), E(e), P(p),
        current_state(initial_state), best_state(initial_state),
        neighbour_fcn_computes_energy(false), neighbour_energy(0.0),
        accepted(0), total(0), non_changed(0),
        rng(), prob_range(0, 1), random(rng, prob_range),
        paused(false), debug(_debug),
        responsiveness_level(_responsiveness_level)
    {
      current_energy = best_energy = E( initial_state );
    }

    template <typename RNGSeedType> // = typename RNG::result_type>
    SA(
        const temperature_schedule & ts,
        const EnergyFcn & e,
        const ProbFcn & p,
        const State & initial_state,
        RNGSeedType seed,
#ifdef CONAN_DEBUG
        bool _debug = true,
#else
        bool _debug = false,
#endif
        size_t _responsiveness_level = 100
      )
      : temp_sched(ts), E(e), P(p),
        current_state(initial_state), best_state(initial_state),
        neighbour_fcn_computes_energy(false), neighbour_energy(0.0),
        accepted(0), total(0), non_changed(0),
        rng(seed), prob_range(0, 1), random(rng, prob_range),
        paused(false), debug(_debug),
        responsiveness_level(_responsiveness_level)
    {
      current_energy = best_energy = E( initial_state );
#ifdef CONAN_DEBUG
      if ( debug )
      {
        std::cerr << "conan::SA::SA(): initial state: " << current_state.tostring() << std::endl;
        std::cerr << "conan::SA::SA(): initial energy: " << current_energy << std::endl;
      }
#endif
    }

    // member functions
    void add_neighbour_fcn(
        const neighbour_fcn<State, RNG> & neighbour
        )
    {
#ifdef CONAN_INFO
      if ( debug )
        std::cerr << "conan::SA::add_neighbour_fcn(): Neighbour-generator function added. It will be called "
                  << neighbour.num_calls_per_iteration << " times per SA iteration." << std::endl;
#endif
      neighbour_fcn_set.push_back( neighbour );
      return;
    }

    void run()
    {
      using namespace detail;

      if (neighbour_fcn_set.size() < 1)
        throw std::runtime_error( "No functions found for generating a neighbour state." );

      total_sa_iterations = 0; // for summary
      time( &start ); // for summary

      while ( ! temp_sched.finished() )
      {
        decimal current_temp = temp_sched.next_temperature();

#ifdef CONAN_INFO
        if ( debug )
          std::cerr << "conan::SA::run(): current temperature = " << current_temp << std::endl;
#endif

        std::vector<size_t> num_calls;
        size_t total_num_calls = 0;

        BOOST_FOREACH( NeighbourFcn & nb, neighbour_fcn_set )
        {
          num_calls.push_back( nb.num_calls_per_iteration );
          total_num_calls += nb.num_calls_per_iteration;
        }

        int last_nb_fcn_index = -1;

        while ( ! paused )
        {
#if defined(HAVE_TERMIOS) && !defined(CONAN_NON_INTERACTIVE)
          if ( ! (total % responsiveness_level) ) // to not waste so much time with kbhit()
          {
            if ( kbhit() )
            {
              int c = getchar();
              if ( char(c) == 'p' )
              {
                std::cerr << current_state.tostring() << std::endl;
              }
              else if ( char(c) == 's' )
              {
                print_statistics();
              }
              else if ( char(c) == 'b' )
              {
                std::cerr << best_energy << std::endl;
              }
            }
          }
#endif

          size_t i = 0;

          if ( last_nb_fcn_index >= 0 )
            i = last_nb_fcn_index;
          else
            i = weighted_choice(num_calls, total_num_calls, rng);

          NeighbourFcn & nb = neighbour_fcn_set[ i ];

          if ( nb.remaining_calls() )
          {
            State neighbour( nb( current_state, *this ) );
            if ( ! neighbour_fcn_computes_energy )
              neighbour_energy = E( neighbour );

#ifdef CONAN_INFO
            std::cout << "current energy = " << current_energy << ", neighbour energy = " << neighbour_energy << std::endl;
            if ( current_energy == neighbour_energy )
            {
              std::cout << "current energy = neighbour energy = " << neighbour_energy << std::endl;
              std::cout << "current state: " << current_state.tostring() << std::endl;
              std::cout << "neighbour state: " << neighbour.tostring() << std::endl;
            }
#endif

            if ( neighbour_energy < best_energy )
            {
              best_state = current_state = neighbour;
              best_energy = current_energy = neighbour_energy;
#ifdef CONAN_INFO
              if ( debug )
                std::cerr << "conan::SA::run(): best energy = " << best_energy << " found!" << std::endl;
#endif
              ++accepted; // for statistics
            }
            else if ( P( current_energy, neighbour_energy, current_temp ) > random() )
            {
              if (current_energy == neighbour_energy) ++non_changed; // for statistics

              current_state = neighbour;
              current_energy = neighbour_energy;
              ++accepted; // for statistics
            }

            ++total; // for statistics
          }
          else
          {
            int remaining_nb_fcns = neighbour_fcn_set.size(),
                possible_last_nb_fcn_index = -1, index = 0;

            BOOST_FOREACH( NeighbourFcn & nb_check, neighbour_fcn_set )
            {
              if ( nb_check.remaining_calls() )
                possible_last_nb_fcn_index = index;
              else
                --remaining_nb_fcns;

              ++index;
            }

            if ( remaining_nb_fcns == 0 )
              break;
            else if ( remaining_nb_fcns == 1 && possible_last_nb_fcn_index != -1 )
              last_nb_fcn_index = possible_last_nb_fcn_index;

          } // if nb.remaining_calls()

        } // while true

        if ( paused ) return;

        BOOST_FOREACH( NeighbourFcn & nb, neighbour_fcn_set )
          nb.reset_counter();

        ++total_sa_iterations; // for summary

      } // while !temp_sched.finished()

      time( &end ); // for summary

#ifdef CONAN_INFO
      if (debug)
      {
        print_summary();
        print_statistics();
      }
#endif

      return;
    } // void run()

    double get_current_temperature() const
    { return temp_sched.current_T; }

    double get_final_temperature() const
    { return temp_sched.final_T; }

    double get_initial_temperature() const
    { return temp_sched.initial_T; }

    State & get_current_state()
    { return current_state; }

    const State & get_current_state() const
    { return current_state; }

    State & get_best_state()
    { return best_state; }

    const State & get_best_state() const
    { return best_state; }

    void print_summary()
    {
      std::string method_name = (debug? "conan::SA::print_summary(): ": "");
      std::cerr << method_name << "Completed." << std::endl
                << method_name << "Summary:" << std::endl
                << method_name << "  Total number of SA iterations = " << total_sa_iterations << std::endl
                << method_name << "  Time elapsed = " << difftime( end, start ) << " segs" << std::endl
                << method_name << "  Best energy found = " << best_energy << std::endl;
      return;
    }

    void print_statistics()
    {
      std::cerr << (debug? "conan::SA::print_statistics(): " : "")
                << "Statistics: accepted = " << accepted
                << " (non changed = " << non_changed << "),"
                << " rejected = " << total - accepted << ","
                << " total = " << total << std::endl;
      return;
    }

    void stop()
    { paused = true; return; }

    void resume()
    {
      BOOST_FOREACH( NeighbourFcn & nb, neighbour_fcn_set )
        nb.reset_counter();
      run();
    }
  };

  double default_acceptance_probability(
      double e1,
      double e2,
      double T
      )
  {
    if (e2 <= e1)
      return 1;
    else
      return exp( (e1 - e2) / T );
  }

} // conan

#endif // SA_HPP
