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

#ifndef STROGRATZ_WATTS_HPP
#define STROGRATZ_WATTS_HPP
#include <conan/config.hpp>
#include <conan/graph_models/canonical_graphs.hpp>
#include <conan/transformations.hpp>

namespace conan {

  /**
   * @brief Return a graph created based on the Strogratz-Watts model.
   *
   * @param V number of vertices.
   * @param half_degree the same as in generate_regular_graph.
   * @param percent_replaced_edges the number of edges from a regular graph generated with V and k which will be replaced by random edges.
   * @return a small-world graph
   */
  template <class Graph> Graph generate_small_world_network (size_t V,
                                                             size_t degree,
                                                             decimal percent_replaced_edges);


  template <class Graph>
  Graph generate_small_world_network(
      size_t V,
      size_t degree,
      decimal percent_replaced_edges
      )
  { // O ( ? )
    // if k < V this function returns an empty graph (ie. num_vertices(returned_graph) = 0)

    if (V == 0)
      return Graph(0); // shortcut

    Graph g = generate_regular_lattice<Graph>(V, degree);
    return randomize_graph_edges(g, percent_replaced_edges, true);
  }

} //conan

#endif //STROGRATZ_WATTS_HPP
