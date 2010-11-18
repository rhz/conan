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

#ifndef COLOR_HPP
#define COLOR_HPP
#include <limits>
#include <conan/config.hpp>
#include <conan/utils.hpp>
#include <boost/graph/sequential_vertex_coloring.hpp>
#include <boost/numeric/ublas/matrix.hpp>


namespace conan {

  /**
   * @brief Compute the chromatic number of the graph.
   * The chromatic number is the minimum number of colors needed to draw each vertex with
   * a color different from those of his neighbors.
   * @warning This function does not compute the minimum.
   *
   * @param g A graph.
   * @return The chromatic number.
   */
  template <class Graph>
  size_t chromatic_number(
      Graph& g);


  template <class Graph>
  typename std::vector<size_t> smallest_last_ordering(
      Graph& g)
  {
    using namespace boost::numeric;
    typedef typename ublas::matrix<decimal> cmatrix;
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex;

    std::vector<size_t> out_vec;

    cmatrix adj_matrix = get_adj_matrix<Graph, cmatrix>(g);
    size_t num_vertices = boost::num_vertices(g);
    for (size_t step = 0; step < num_vertices; ++step)
    {
      size_t min_degree = std::numeric_limits<size_t>::max(),
           min_degree_index = 0;
      for (size_t i = 0; i < num_vertices; ++i)
      {
        // check if i has been added to out_vec
        for (std::vector<size_t>::iterator veci = out_vec.begin(); veci != out_vec.end(); ++veci)
          if (*veci == i)
            continue;

        size_t i_degree = 0;
        for (size_t j = 0; j < num_vertices; j++)
        {
          // check if j has been added to out_vec
          for (std::vector<size_t>::iterator veci = out_vec.begin(); veci != out_vec.end(); ++veci)
            if (*veci == j)
              continue;

          i_degree += adj_matrix(i, j);
        }

        if (i_degree < min_degree)
        {
          min_degree = i_degree;
          min_degree_index = i;
        }
      } // for i < num_vertices
      out_vec.push_back(min_degree_index);
    } // for step < num_vertices

    return out_vec;
  }
  
  template <class Graph>
  size_t chromatic_number(
      Graph& g
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertices_size_type vertices_size_type;
    typedef typename GraphTraits::vertex_iterator vertex_iter;
    typedef typename boost::property_map<Graph, vertex_index_t>::const_type vertex_index_map;

    size_t V = boost::num_vertices(g);
    std::vector<vertices_size_type> color_vec(V);
    iterator_property_map<vertices_size_type*, vertex_index_map>
      color(&color_vec.front(), get(vertex_index, g));
    std::vector<size_t> order_vec(V);
    vertex_iter vi = boost::vertices(g).first;
    for (size_t i = 0; i < V && vi != boost::vertices(g).second; ++i, ++vi)
      order_vec[i] = *vi;
    iterator_property_map<size_t*, vertex_index_map>
      order(&order_vec.front(), get(vertex_index, g));
    return sequential_vertex_coloring(g, order, color);
  }


  template <class Graph>
  size_t chromatic_number(
      Graph & g,
      std::vector<size_t> order_vec
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertices_size_type vertices_size_type;
    typedef typename GraphTraits::vertex_iterator vertex_iter;
    typedef typename boost::property_map<Graph, vertex_index_t>::const_type vertex_index_map;

    size_t V = boost::num_vertices(g);
    if (order_vec.size() != V)
      throw std::runtime_error("chromatic_number: size of vector order_vec must be equal to the number of vertices of graph g");

    std::vector<vertices_size_type> color_vec(V);
    iterator_property_map<vertices_size_type*, vertex_index_map>
      color(&color_vec.front(), get(vertex_index, g));
    iterator_property_map<size_t*, vertex_index_map>
      order(&order_vec.front(), get(vertex_index, g));
    return sequential_vertex_coloring(g, order, color);
  }

} // conan

#endif // COLOR_HPP
