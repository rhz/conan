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

#ifndef SPECTRAL_HPP
#define SPECTRAL_HPP
#include <conan/config.hpp>
#include <cmath> // for sqrt

// GSL
#define HAVE_INLINE
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>


namespace conan {

  // *** Declarations ***

  template <class Graph>
  void graph_laplacian(
      gsl_matrix * L,
      const Graph & g);

  template <class Graph>
  std::vector<double> spectrum(
      const Graph & g);

  template <class Graph>
  void graph_weighted_laplacian(
      gsl_matrix * L,
      const Graph & g);

  template <class Graph>
  size_t graph_volume(
      const Graph & g);

  template <class Graph>
  decimal weighted_degree(
      typename boost::graph_traits<Graph>::vertex_descriptor v,
      const Graph & g);

  template <class Graph>
  decimal graph_weighted_volume(
      const Graph & g);


  // *** Definitions ***

  template <class Graph>
  void graph_laplacian(
      gsl_matrix * L,
      const Graph & g
      )
  {
    size_t size = L->size1;

    if (size != L->size2)
      throw std::runtime_error("gsl_matrix *L must be a square matrix (size1 == size2)");

    for (size_t i = 0; i < size; ++i)
    {
      for (size_t j = 0; j < size; ++j)
      {
        if (i == j)
          gsl_matrix_set(L, i, j, 1.0);
        else if (boost::edge(i, j, g).second) // if i and j are adjacent
          gsl_matrix_set(L, i, j, - 1.0 / sqrt(boost::degree(i, g) * boost::degree(j, g)));
        else
          gsl_matrix_set(L, i, j, 0.0);
      }
    }

    return;
  }


  template <class Graph>
  std::vector<double> spectrum(
      const Graph & g
      )
  {
    size_t V = boost::num_vertices(g);

    gsl_matrix * L = gsl_matrix_alloc(V, V);

    graph_laplacian(L, g);

    // compute eigenvalues
    gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(V);
    gsl_vector *eval = gsl_vector_calloc(V);
    gsl_eigen_symm(A, eval, w);

    std::vector<double> spectrum;
    for (size_t i = 0; i < eval->size; ++i)
      spectrum.push_back(gsl_vector_get(eval, i));

    gsl_vector_free(eval);
    gsl_eigen_symm_free(w);
    gsl_matrix_free(L);

    return spectrum;
  }


  template <class Graph>
  void graph_weighted_laplacian(
      gsl_matrix * L,
      const Graph & g
      )
  {
    size_t size = L->size1;

    if (size != L->size2)
      throw std::runtime_error("gsl_matrix *L must be a square matrix (size1 == size2)");

    for (size_t i = 0; i < size; ++i)
    {
      for (size_t j = 0; j < size; ++j)
      {
        if (i == j)
        {
          if (boost::degree(i, g) != 0)
          {
            double w = (boost::edge(i, i, g).second? g[boost::edge(i, i, g).first].weight : 0.0);
            gsl_matrix_set(L, i, j, 1.0 - (w / boost::degree(i, g)));
          }
          else
            gsl_matrix_set(L, i, j, 0.0);
        }
        else if (boost::edge(i, j, g).second) // if i and j are adjacent
        {
          double w = g[boost::edge(i, i, g).first].weight;
          gsl_matrix_set(L, i, j, - w / sqrt(boost::degree(i, g) * boost::degree(j, g)));
        }
        else
          gsl_matrix_set(L, i, j, 0.0);
      }
    }

    return;
  }


  template <class Graph>
  size_t graph_volume(
      const Graph & g
      )
  { // Do this function have any sense for directed graphs?
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_iterator vertex_iter;

    size_t sum_degree = 0;

    vertex_iter vi, viend;
    for (tie(vi, viend) = boost::vertices(g); vi != viend; ++vi)
      sum_degree += boost::degree(*vi, g);

    return sum_degree;
  }

  template <class Graph>
  decimal weighted_degree(
      typename boost::graph_traits<Graph>::vertex_descriptor v,
      const Graph & g
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::out_edge_iterator out_edge_iter;

    decimal sum_weight = 0.0;

    out_edge_iter oei, oei_end;
    for (tie(oei, oei_end) = boost::out_edges(v, g); oei != oei_end; ++oei)
      sum_weight += g[*oei].weight;

    return sum_weight;
  }

  template <class Graph>
  decimal graph_weighted_volume(
      const Graph & g
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_iterator vertex_iter;

    decimal sum_degree = 0;

    vertex_iter vi, viend;
    for (tie(vi, viend) = boost::vertices(g); vi != viend; ++vi)
      sum_degree += weighted_degree(*vi, g);

    return sum_degree;
  }

} // conan

#endif // SPECTRAL_HPP
