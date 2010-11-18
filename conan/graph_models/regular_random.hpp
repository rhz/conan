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

#ifndef REGULAR_RANDOM_HPP
#define REGULAR_RANDOM_HPP
#include <conan/config.hpp>
#include <conan/properties/distributions.hpp>

namespace conan {

  /**
   * generate_random_quasi_regular_graph
   * Generates a random regular graph according to the pairing model, but enforces it to be simple.
   * This causes that the function can not always create a perfect regular graph.
   * 
   * To overcome this issue exists the function generate_random_regular_graph, which calls the function
   * generate_random_quasi_regular_graph until it has returned a perfect regular graph (with a delta
   * function degree distribution).
   *
   * @param V Number of vertices.
   * @param k Degree.
   */
  template <class Graph>
  Graph generate_random_quasi_regular_graph(
      size_t V,
      size_t k,
      size_t max_num_iterations = 0);

  /**
   * @brief Generate a random regular graph according to the pairing model, but enforces it to be simple.
   *
   * @param V Number of vertices.
   * @param k Degree.
   */
  template <class Graph>
  Graph generate_random_regular_graph(
      size_t V,
      size_t k);


  /***** Definitions *****/

  template <class Graph>
  Graph generate_random_regular_graph(
      size_t V,
      size_t k)
  {
    Graph g;
    while (true)
    {
      g = generate_random_quasi_regular_graph<Graph>(V, k);
      degree_distribution<Graph> d(g);
      if (d.P(k) == 1)
        break;
    }
    return g;
  }


  template <class Graph>
  Graph generate_random_quasi_regular_graph(
      size_t V,
      size_t k,
      size_t max_num_iterations // = 0
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex;
    typedef typename GraphTraits::vertex_iterator vertex_iter;
    typedef typename GraphTraits::edge_descriptor edge;

    typename boost::uniform_real<decimal> percent(0, 1);
    typename boost::uniform_int<size_t> int_dist(0, V - 1);
#ifdef DETERMINISTIC_RNG
    typename boost::mt19937 rng(time(NULL) - V - k - max_num_iterations);
    typename boost::variate_generator< boost::mt19937&, boost::uniform_real<decimal> >
      probability(rng, percent);
    typename boost::variate_generator< boost::mt19937&, boost::uniform_int<size_t> >
      random_vertex(rng, int_dist);
#else
    typename boost::random_device rng;
    typename boost::variate_generator< boost::random_device&, boost::uniform_real<decimal> >
      probability(rng, percent);
    typename boost::variate_generator< boost::random_device&, boost::uniform_int<size_t> >
      random_vertex(rng, int_dist);
#endif

    if (max_num_iterations == 0)
      max_num_iterations = V;

    Graph g(V);
    vertex_iter vi, viend;
    for (tie(vi, viend) = boost::vertices(g); vi != viend; ++vi)
    {
      size_t counter = max_num_iterations;
      while (boost::degree(*vi, g) < k)
      {
        if (counter-- == 0)
        {
#ifdef CONAN_INFO
          std::cerr << "generate_regular_random_graph: vertex " << *vi << " could not find a partner" << std::endl;
#endif
          break;
        }

        vertex target_vertex;
        // To select a target vertex, it must not be the same as the source one (*vi),
        // it must not have degree k and it must not be linked with the source vertex.
        // If all these conditions are not satisfied within the first max_num_iterations,
        // we give up and inform the user.
        target_vertex = vertex(random_vertex());
#if 0
        while ((target_vertex = vertex(random_vertex())) == *vi
                || boost::degree(target_vertex, g) == k
                || boost::edge(*vi, target_vertex, g).second)
          { }
#endif

        if (*vi == target_vertex)
        {
#ifdef CONAN_INFO
          std::cerr << "generate_regular_random_graph: source and target vertices are the same (" << *vi
                    << " and " << target_vertex << ")" << std::endl;
#endif
          continue;
        }

        if (boost::edge(*vi, target_vertex, g).second)
        {
#ifdef CONAN_INFO
          std::cerr << "generate_regular_random_graph: edge between vertices " << *vi << " and "
                    << target_vertex << " already exists" << std::endl;
#endif
          continue;
        }

        if (boost::degree(target_vertex, g) == k)
        {
#ifdef CONAN_INFO
          std::cerr << "generate_regular_random_graph: vertex " << target_vertex << " have degree " << k
                    << " (*vi = " << *vi << ")" << std::endl;
#endif
          continue;
        }

        edge e = boost::add_edge(*vi, target_vertex, g).first;
        g[e].weight = 1.0;

      } // while (boost::degree(*vi, g) < k)
    }

    return g;
  }


} // conan

#endif //REGULAR_RANDOM_HPP
