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

#ifndef ROSVALL_AND_BERGSTROM_HPP
#define ROSVALL_AND_BERGSTROM_HPP
#include <conan/graphs.hpp>
#include <conan/subgraph/guimera_and_amaral.hpp>
#include <conan/utils/combinatoria.hpp>

namespace conan {

  namespace detail {

    template <class Graph>
    inline bool l_ij(
        const Module & source_module,
        const Module & target_module,
        Graph & g
        )
    {
      size_t l = 0;
      BOOST_FOREACH( const size_t & source_vertex, source_module )
      {
        BOOST_FOREACH( const size_t & target_vertex, target_module )
        {
          if ( boost::edge(source_vertex, target_vertex, g).second )
            l += 1;
        }
      }
      return l;
    }


    double mutual_information_cost_fcn(
        const ModulesAndGraph & m
        )
    {
      if ( m.modules.size() < 1 )
        return 0.0;

      size_t n = boost::num_vertices( m.graph() ),
             l = boost::num_edges( m.graph() ),
             m_size = m.modules.size(),
             n_0 = m.modules[ 0 ].size();

      double HZ = binomial_coefficient( n_0 * (n_0 - 1) / 2, l_ij( m.modules[ 0 ], m.modules[ 0 ], m.graph() ) );

      for (size_t i = 1; i < m_size; ++i)
      {
        size_t n_i = m.modules[ i ].size();
        double l_ii = l_ij( m.modules[ i ], m.modules[ i ], m.graph() );

        double prod2 = binomial_coefficient( n_i * m.modules[ 0 ].size(), l_ij( m.modules[ i ], m.modules[ 0 ], m.graph() ) );
        for (size_t j = 1; j < m_size; ++j)
        {
          if ( i > j )
          {
            size_t n_j = m.modules[ j ].size();
            prod2 *= binomial_coefficient( n_i * n_j, l_ij( m.modules[ i ], m.modules[ j ], m.graph() ) );
          }
        }

        HZ *= binomial_coefficient( n_i * (n_i - 1) / 2, l_ii ) * prod2;
      }

      HZ = conan_log2(HZ);

      return (n * conan_log2(m_size) + m_size * (m_size + 1) / 2 * conan_log2(l) + HZ );
    }

  } // detail

  template <class Graph>
  inline void rosvall_bergstrom_communities(
      Graph & g,
      std::vector< std::vector<size_t> > & modules,
      double initial_temp = .0001,
      double final_temp = .00006,
      double f = 1.0
      )
  { return sa_communities(g, modules, &detail::mutual_information_cost_fcn, false, initial_temp, final_temp, f); }

} // conan

#endif // ROSVALL_AND_BERGSTROM_HPP
