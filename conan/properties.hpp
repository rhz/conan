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

#ifndef PROPERTIES_HPP
#define PROPERTIES_HPP
#include <conan/graphs.hpp>
#include <conan/properties/random_props.hpp>
#include <conan/properties/clustering.hpp>
#include <conan/properties/connectance.hpp>
#include <conan/properties/entropy.hpp>
#include <conan/properties/fractal_dim.hpp>
#include <conan/properties/distributions.hpp>
#include <conan/properties/shortest_paths.hpp>
#include <conan/properties/centrality.hpp>

namespace conan {

  /**
   * @brief Computes the mu factor (also known as small-world coefficient) for a graph.
   *
   * @param g A graph.
   * @return The value of the mu factor.
   */
  template <class Graph>
  decimal mu_factor (const Graph& g);

  
  template <class Graph>
  decimal mu_factor(
      const Graph & g
      )
  {
    int V = boost::num_vertices(g),
        E = boost::num_edges(g);
    decimal normalized_asp = graph_avg_shortest_path(g) / random_graph_avg_shortest_path(V, E),
            normalized_clustering = graph_avg_clustering(g) / random_graph_avg_clustering(V, E);
    return normalized_clustering / normalized_asp;
  }

} //namespace conan

#endif //PROPERTIES_HPP
