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

#ifndef ROLES_HPP
#define ROLES_HPP
#include <conan/graphs.hpp>


namespace conan {

  namespace detail {

    int index_of(
        size_t vertex,
        std::vector< std::vector<size_t> > & modules
        )
    {
      typedef std::vector<size_t> Module;

      int i = 0;

      // Search vertex in modules
      BOOST_FOREACH( Module m, modules )
      {
        BOOST_FOREACH( size_t v, m )
        {
          if (v == vertex)
            return i;
        }

        ++i;
      }

      return -1; // not found
    }

  } // detail


  template <class Graph>
  size_t within_module_connectivity(
      Graph & g,
      std::vector<size_t> & module,
      size_t vertex
      )
  {
    size_t k_i = 0;

    BOOST_FOREACH( size_t v, module )
    {
      if ( boost::edge(vertex, v, g).second )
        ++k_i;
    }

    return k_i;
  }


  template <class Graph>
  decimal z_score(
      Graph & g,
      std::vector< std::vector<size_t> > & modules,
      size_t vertex
      )
  {
    using namespace detail;
    typedef typename std::vector<size_t> Module;
    typedef typename std::vector<Module> ModuleVector;

    static Graph * prev_g_ptr = NULL;
    static ModuleVector * prev_module_vector_ptr = NULL;
    static std::vector<decimal> mean_k_s_i( modules.size(), -1.0 ),
                                sd_k_s_i( modules.size(), -1.0 );
    decimal k_i = -1.0;

    int index = index_of( vertex, modules );
    Module s_i( modules[ index ] );

#ifdef CONAN_DEBUG
    std::cerr << "vertex = " << vertex << std::endl
              << "index = " << index << std::endl;
#endif

    if ( prev_g_ptr != & g ||
         prev_module_vector_ptr != & modules )
    {
      mean_k_s_i = std::vector<decimal>( modules.size(), -1.0 );
      sd_k_s_i = std::vector<decimal>( modules.size(), -1.0 );
    }

    if ( mean_k_s_i[ index ] != -1.0 &&
         sd_k_s_i[ index ] != -1.0 )
    {
      k_i = within_module_connectivity( g, s_i, vertex );
    }
    else
    {
      size_t size = s_i.size();
      double vec[ size ];
      size_t i = 0;

      BOOST_FOREACH( size_t v, s_i )
      {
        vec[ i ] = within_module_connectivity( g, s_i, v );

        if ( v == vertex )
          k_i = vec[ i ];

        ++i;
      }

      mean_k_s_i[ index ] = gsl_stats_mean( vec, 1, size );
      sd_k_s_i[ index ] = gsl_stats_sd_m( vec, 1, size, mean_k_s_i[ index ] );
    }

#ifdef CONAN_DEBUG
    std::cerr << "(k_i - mean_k_s_i) / sd_k_s_i = (" << k_i << " - " << mean_k_s_i[ index ] << ") / " << sd_k_s_i[ index ]
              << " = " << (k_i - mean_k_s_i[ index ]) / sd_k_s_i[ index ] << std::endl;
#endif

    return (k_i - mean_k_s_i[ index ]) / sd_k_s_i[ index ];
  }


  template <class Graph>
  decimal participation_coefficient(
      Graph & g,
      std::vector< std::vector<size_t> > & modules,
      size_t vertex
      )
  {
    typedef typename std::vector<size_t> Module;

    size_t vertex_degree = boost::out_degree(vertex, g);
    decimal sum = 0.0;

    BOOST_FOREACH( Module m, modules )
    {
      size_t k_is = 0;

      BOOST_FOREACH( size_t v, m )
      {
        if ( boost::edge( vertex, v, g ).second )
          ++k_is;
      }

      sum += conan_pow( k_is / decimal(vertex_degree), 2 );
    }

    return 1 - sum;
  }

} // conan

#endif // MODULARITY_HPP
