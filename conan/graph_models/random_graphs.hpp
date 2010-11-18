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

#ifndef RANDOM_GRAPHS_HPP
#define RANDOM_GRAPHS_HPP
#include <conan/config.hpp>
#include <boost/graph/random.hpp>

namespace conan {

  // ***** Declarations *****
  
  /**
   * @brief Return a graph with V vertices and E randomly-generated edges between those vertices.
   *
   * @param V Number of vertices.
   * @param E Number of edges.
   * @return A random graph.
   */
  template <class Graph, class RNG>
  Graph generate_random_graph(
      size_t V,
      size_t E,
      RNG gen
      );
  
  /**
   * @brief Return a graph in which there is a probability 'p' of connecting any two vertices.
   *
   * @param V Number of vertices.
   * @param p Probability of any two vertices to be connected.
   * @return A random graph.
   */
  template <class Graph>
  Graph generate_erdos_renyi_graph(
      size_t V,
      decimal p
      );
  

  // ***** Definitions ******

  template <class Graph, class RNG>
  Graph generate_random_graph(
      size_t V,
      size_t E,
      RNG gen
      )
  {
    Graph g;
    boost::generate_random_graph(g, V, E, gen);
    return booleanize_graph(g);
  }


  template <class Graph>
  Graph generate_random_graph(
      size_t V,
      size_t E
      )
  {
    typename boost::mt19937 rng(time(NULL) - V - E);
    return generate_random_graph<Graph>(V, E, rng);
  }


  template <class Graph>
  Graph generate_erdos_renyi_graph(
      size_t V,
      decimal p
      )
  {
    typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iter;
    typedef typename boost::graph_traits<Graph>::edge_descriptor edge;
    typename boost::uniform_real<decimal> percent_dist(0, 1);
#ifdef DETERMINISTIC_RNG
    typename boost::mt19937 rng(time(NULL) - V - int(p * 100));
    typename boost::variate_generator< boost::mt19937&, boost::uniform_real<decimal> >
      percent(rng, percent_dist);
#else
    typename boost::random_device dev;
    typename boost::variate_generator< boost::random_device&, boost::uniform_real<decimal> >
      percent(dev, percent_dist);
#endif

    if (p < 0 || p > 1)
      throw std::runtime_error("p must be a real number between 0 and 1");

    Graph g(V);
    bool undirected = boost::is_undirected(g);
    size_t counter = 1;
    vertex_iter source_vi, viend;
    for (tie(source_vi, viend) = boost::vertices(g); source_vi != viend; ++source_vi, ++counter)
      for (vertex_iter target_vi = boost::vertices(g).first + (undirected? counter : 0);
           target_vi != viend; ++target_vi)
      {
        if (source_vi == target_vi)
          continue;
        if (percent() < p and not boost::edge(*source_vi, *target_vi, g).second)
        {
          edge e = boost::add_edge(*source_vi, *target_vi, g).first;
          g[e].weight = 1.0;
        }
      }
    
    return g;
  }

} // conan

#endif //RANDOM_GRAPHS_HPP
