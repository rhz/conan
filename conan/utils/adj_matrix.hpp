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

#ifndef ADJ_MATRIX_HPP
#define ADJ_MATRIX_HPP
#include <conan/graphs.hpp>
#include <fstream>

#define HAVE_INLINE
#include <gsl/gsl_matrix.h>

namespace conan {

  // ***** Declarations *****

  /**
   * @brief Get the adjacency matrix from a graph.
   *
   * @param g a graph
   * @return a matrix containing the adjacency matrix of g.
   */
  template <class Graph, class Matrix>
  Matrix get_adj_matrix(
      const Graph &
      );

  template <class Graph>
  void get_adj_matrix(
      const Graph &,
      gsl_matrix *
      );

  /**
   * @brief Construct a graph based on the adjacency matrix.
   *
   * @param m an adjacency matrix
   * @return the graph
   */
  template <class Graph, class Matrix>
  Graph make_graph_from_adj_matrix(
      const Matrix & m
      );

  template <class Graph>
  Graph make_graph_from_adj_matrix(
      double ** m,
      size_t size
      );

  template <class Graph>
  Graph make_graph_from_adj_matrix(
      gsl_matrix * M
      );


  // ***** Definitions *****


  template <class Graph, class Matrix>
  Matrix get_adj_matrix(
      const Graph & g
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::edge_iterator edge_iter;

    size_t V = boost::num_vertices(g);

    if (V == 0)
      throw std::runtime_error("g must not be a null graph");

    Matrix A(V, V);

    for (size_t i = 0; i < A.size1(); ++i)
      for (size_t j = 0; j < A.size2(); ++j)
        A(i, j) = 0;

    bool is_directed = boost::is_directed(g);

    edge_iter ei, ei_end;
    for (tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
    {
      A( boost::source(*ei, g), boost::target(*ei, g) ) = 1.0;

      if (!is_directed)
        A( boost::target(*ei, g), boost::source(*ei, g) ) = 1.0;
    }

    return A;
  }


  template <class Graph>
  void get_adj_matrix(
      const Graph & g,
      gsl_matrix * A
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::edge_iterator edge_iter;

    size_t V = boost::num_vertices(g);

    if (A->size1 != A->size2 || A->size1 != V)
      throw std::runtime_error("invalid size for adjacency matrix");

    bool is_undirected = boost::is_undirected(g);

    for (size_t r = 0; r < A->size1; ++r)
      for (size_t c = 0; c < A->size2; ++c)
        gsl_matrix_set(A, r, c, 0.0);

    edge_iter ei, ei_end;
    for (tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
    {
      gsl_matrix_set(A, boost::source(*ei, g), boost::target(*ei, g), 1.0);

      if (is_undirected)
        gsl_matrix_set(A, boost::target(*ei, g), boost::source(*ei, g), 1.0);
    }

    return;
  }


  template <class Graph, class Matrix>
  Graph make_graph_from_adj_matrix(
      const Matrix & m
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::edge_descriptor edge;

    if (m.size1() != m.size2())
    {
      std::cerr << "m.size1() != m.size2()" << std::endl;
      throw FormatError();
    }

    size_t V = m.size1();

    Graph g(V);

    bool undirected = boost::is_undirected(g);

    for (size_t i = 0; i < V; ++i)
    {
      for (size_t j = (undirected? i + 1 : 0); j < V; ++j)
      {
        if (m(i, j) != 0)
          boost::add_edge(i, j, m(i, j), g);
      }
    }

    return g;
  }


  template <class Graph>
  Graph make_graph_from_adj_matrix(
      double ** m,
      size_t size
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::edge_descriptor edge;

    Graph g(size);

    bool undirected = boost::is_undirected(g);

    for (size_t i = 0; i < size; ++i)
    {
      for (size_t j = (undirected? i + 1 : 0); j < size; ++j)
      {
        if (m[i][j] != 0.0)
          boost::add_edge(i, j, m[i][j], g);
      }
    }

    return g;
  }


  template <class Graph>
  Graph make_graph_from_adj_matrix(
      gsl_matrix * M
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::edge_descriptor edge;

    size_t size = M->size1;

    if (size != M->size2)
      throw std::runtime_error("invalid lenght of input matrix.");

    Graph g(size);

    bool undirected = boost::is_undirected(g);

    for (size_t i = 0; i < size; ++i)
    {
      for (size_t j = (undirected? i + 1 : 0); j < size; ++j)
      {
        if (gsl_matrix_get(M, i, j) != 0.0)
          boost::add_edge( i, j, gsl_matrix_get(M, i, j), g );
      }
    }

    return g;
  }

} //conan

#endif // ADJ_MATRIX_HPP
