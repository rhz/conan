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

#ifndef RANDOM_PROPS_HPP
#define RANDOM_PROPS_HPP

#include <conan/config.hpp>
#include <cmath>

namespace conan {

  /**
   * @brief Computes the average degree for a random graph analytically.
   * @param num_vertices Number of vertices.
   * @param num_edges Number of edges.
   */
  decimal random_graph_average_degree(size_t num_vertices,
                                      size_t num_edges);
  
  /**
   * @brief Computes the average shortest path for a random graph analytically.
   * @param num_vertices Number of vertices.
   * @param num_edges Number of edges.
   */
  decimal random_graph_avg_shortest_path (size_t num_vertices,
                                          size_t num_edges);

  /**
   * @brief Computes the average clustering coefficient for a random graph analytically.
   * @param num_vertices Number of vertices.
   * @param num_edges Number of edges.
   */
  decimal random_graph_avg_clustering (size_t num_vertices,
                                       size_t num_edges);


  inline
  decimal random_graph_average_degree(
      size_t num_vertices,
      size_t num_edges
      )
  {
    return (2.0 * static_cast<decimal>(num_edges) / num_vertices);
  }


  inline
  decimal random_graph_avg_shortest_path(
      size_t num_vertices,
      size_t num_edges
      )
  {
    if ((num_vertices == 2 && num_edges == 1) || num_vertices == 1)
      return 1;
    else
      return conan_log(num_vertices) / conan_log(random_graph_average_degree(num_vertices, num_edges));
  }


  inline
  decimal random_graph_avg_clustering(
      size_t num_vertices,
      size_t num_edges
      )
  {
    if (num_vertices == 1)
      return 1;
    else
      return random_graph_average_degree(num_vertices, num_edges) / num_vertices;
  }
} // conan

#endif //RANDOM_PROPS_HPP
