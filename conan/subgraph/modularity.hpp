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

#ifndef MODULARITY_HPP
#define MODULARITY_HPP
#include <conan/graphs.hpp>

#define HAVE_INLINE
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>


namespace conan {

  /// @brief Computes the modularity matrix as defined by Newman (ref.)
  template <class Graph>
  void get_modularity_matrix(
      Graph & g,
      gsl_matrix * B,
      std::vector<size_t> & vertex_subset
      )
  {
    size_t vertex_subset_size = vertex_subset.size();

    if (B->size1 != B->size2 || B->size1 != vertex_subset_size)
      throw std::runtime_error("invalid size for modularity matrix B");

    size_t num_vertices = boost::num_vertices(g);
    size_t num_edges = boost::num_edges(g);
    gsl_matrix * A = gsl_matrix_alloc(num_vertices, num_vertices);

    get_adj_matrix(g, A);

    if (num_vertices == vertex_subset_size)
    {
      // Analize the whole graph
      for (size_t r = 0; r < vertex_subset_size; ++r)
      {
        for (size_t c = 0; c < vertex_subset_size; ++c)
        {
          decimal avg_random_edges = boost::out_degree(vertex_subset[r], g) * boost::out_degree(vertex_subset[c], g) / (2.0 * num_edges);
          gsl_matrix_set(B, r, c, gsl_matrix_get(A, vertex_subset[r], vertex_subset[c]) - avg_random_edges);
        }
      }
    }
    else
    {
      // Analize only a subgraph
      std::vector<size_t> local_degrees( vertex_subset_size );
      decimal total_local_degrees = 0;

      for (size_t r = 0; r < vertex_subset_size; ++r)
      {
        for (size_t c = 0; c < vertex_subset_size; ++c)
        {
          local_degrees[r] += gsl_matrix_get(A, vertex_subset[r], vertex_subset[c]);
        }
        total_local_degrees += local_degrees[r];
      }

      for (size_t r = 0; r < vertex_subset_size; ++r)
      {
        for (size_t c = 0; c < vertex_subset_size; ++c)
        {
          decimal avg_random_edges = boost::out_degree(vertex_subset[r], g) * boost::out_degree(vertex_subset[c], g) / (2.0 * num_edges);
          gsl_matrix_set(B, r, c, gsl_matrix_get(A, vertex_subset[r], vertex_subset[c]) - avg_random_edges);
        }
        gsl_matrix_set(B, r, r, gsl_matrix_get(B, r, r) - (local_degrees[r] - (boost::out_degree(vertex_subset[r], g) * total_local_degrees) /(2*num_edges)));
      }
    }

    gsl_matrix_free(A);

    return;
  }


  template <class Graph>
  void get_modularity_matrix(
      Graph & g,
      gsl_matrix * B
      )
  {
    typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iter;

    std::vector<size_t> all_vertices( boost::num_vertices(g) );

    vertex_iter vi, viend;
    tie(vi, viend) = boost::vertices(g);
    std::copy(vi, viend, all_vertices.begin());

    return get_modularity_matrix(g, B, all_vertices);
  }


  namespace detail {

    template <class Graph>
    inline decimal l_s(
        const Graph & g,
        const std::vector<size_t> & module
        )
    {
      typedef typename std::vector<size_t> Module;
      typedef typename Module::const_iterator ModuleIter;

      bool is_directed = boost::is_directed(g);
      size_t edges_cnt = 0, source_cnt = 1;

      ModuleIter source_vvi_end = module.end() - (is_directed? 0 : 1),
                 target_vvi_end = module.end();
      for (ModuleIter source_vvi = module.begin(); source_vvi != source_vvi_end; ++source_vvi, ++source_cnt)
      {
        for (ModuleIter target_vvi = module.begin() + (is_directed? 0 : source_cnt); target_vvi != target_vvi_end; ++target_vvi)
        {
          if (boost::edge(*source_vvi, *target_vvi, g).second) // FIXME: this could be optmized with an unordered_map
            ++edges_cnt;
        }
      }

      return edges_cnt;
    }


    template <class Graph>
    inline decimal d_s(
        const Graph & g,
        const std::vector<size_t> & module
        )
    {
      size_t total_degree = 0;

      BOOST_FOREACH( size_t v, module )
        total_degree += boost::out_degree( v, g );

      return total_degree;
    }


    template <class Graph>
    inline decimal module_modularity(
        const Graph & g,
        const std::vector<size_t> & module
        )
    {
      size_t L = boost::num_edges(g);
      return ( l_s(g, module) / L - conan_pow( d_s(g, module) / (2 * L), 2 ) );
    }

    template <class Graph>
    inline decimal module_modularity(
        const Graph & g,
        const std::vector<size_t> & module,
        size_t L
        )
    {
      return ( l_s(g, module) / L - conan_pow( d_s(g, module) / (2 * L), 2 ) );
    }

  } // detail


  template <class Graph>
  decimal modularity(
      const Graph & g,
      const std::vector< std::vector<size_t> > & modules
      )
  {
    using namespace detail;
    typedef typename std::vector<size_t> Module;

    size_t L = boost::num_edges( g );
    decimal Q = 0.0;

    BOOST_FOREACH( Module m, modules )
      Q += module_modularity( g, m, L );

    return Q;
  }


  template <class Graph>
  decimal modularity(
      const Graph & g,
      const std::vector< std::vector<size_t> > & modules,
      size_t L
      )
  {
    using namespace detail;
    typedef typename std::vector<size_t> Module;

    decimal Q = 0.0;

    BOOST_FOREACH( Module m, modules )
      Q += module_modularity( g, m, L );

    return Q;
  }


  template <class Graph>
  double delta_modularity(
      const Graph & g,
      const std::vector< std::vector<size_t> > & two_modules
      )
  {
    using namespace detail;
    typedef typename std::vector<size_t> Module;

    if ( two_modules.size() != 2 )
      throw std::runtime_error("conan::delta_modularity: two_modules must have 2 elements.");

    Module prev_module( two_modules[0].size() + two_modules[1].size() );

    BOOST_FOREACH( size_t v, two_modules[0] )
      prev_module.push_back( v );
    BOOST_FOREACH( size_t v, two_modules[1] )
      prev_module.push_back( v );

    return ( module_modularity( g, two_modules[0] ) + module_modularity( g, two_modules[1] ) - module_modularity( g, prev_module ) ); // Q_2 - Q_1
  }

} // conan

#endif // MODULARITY_HPP
