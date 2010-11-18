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

#ifndef TRANSFORMATIONS_HPP
#define TRANSFORMATIONS_HPP
#include <conan/graphs.hpp>
#include <conan/utils.hpp>

#include <boost/graph/connected_components.hpp>


namespace conan {

  /**
   * @brief Convert any weighted graph in a boolean graph, ie. a graph in which each edge weights 1.
   *
   * @param g a graph
   * @return the transformed graph
   */
  template <class Graph>
  Graph booleanize_graph(
      Graph g);


  /**
   * @brief Replace edges with a certain probability.
   *
   * @param g a graph.
   * @param p the probability of replacing any edge by a random one.
   * @param maintain_connected Whether to force the graph to remain connected after the edge replacement. False by default.
   * @return a graph.
   */
  template <class Graph>
  Graph randomize_graph_edges(
      Graph g,
      decimal p,
      bool maintain_connected = false);


  /**
   * @brief Swap a certain number of edges by random edges.
   *
   * @param g a graph.
   * @param number_edges the number of edges to be replaced by random edges (or the percentage of them).
   * @param maintain_connected Whether to force the graph to remain connected after the edge replacement. False by default.
   * @return a graph where part of its edges has been replaced by random edges.
   */
  template <class Graph>
  Graph randomize_graph_edges2(
      Graph g,
      decimal number_edges,
      bool maintain_connected = false);


  /**
   * @brief Swap the directionality of all edges.
   *
   * @param g A graph.
   * @return The transformed graph.
   */
  template <class Graph>
  Graph transpose(
      const Graph & g);


  /**
   * @brief Compute the complement of an %undirected graph.
   * The complement of a graph is defined as the graph which has connected all vertices that are not connected in the
   * original graph and viceversa.
   *
   * @param g A graph.
   * @return The complement graph.
   */
  template <class Graph>
  Graph complement(
      const Graph & g);


  /**
   * @brief Remove any number of vertices from the graph or a percentage of them.
   *
   * @param g A graph.
   * @param num_vertices The number of vertices to be removed. It could be also a percentage.
   * @return The resulting graph.
   */
  template <class Graph>
  Graph random_attack(
      Graph g,
      decimal num_vertices);


  /***** DEFINITIONS *****/

  template <class Graph>
  Graph booleanize_graph(
      Graph g
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::edge_iterator edge_iter;

    edge_iter ei, ei_end;
    for (tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
      g[*ei].weight = 1.0;

    return g;
  }


  template <class Graph>
  Graph randomize_graph_edges(
      Graph out_g,
      decimal p,
      bool maintain_connected // = false
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex;
    typedef typename GraphTraits::edge_descriptor edge;
    typedef typename GraphTraits::edge_iterator edge_iter;

    typedef typename std::vector< std::pair<vertex, vertex> > EdgeVector;
    typedef typename EdgeVector::iterator EdgeVectorIter;

#ifdef CONAN_DEBUG
    std::cerr << "conan::randomize_graph_edges: p = " << p << std::endl
              << "conan::randomize_graph_edges: maintain_connected = " << maintain_connected << std::endl;
#endif

    if (p < 0)
      throw std::runtime_error("p must be a probability, ie. it must be equal or greater than 0 and equal or less than 1");
    else if (p > 1)
      return randomize_graph_edges2(out_g, p, maintain_connected);

    typename boost::uniform_real<decimal> percent(0, 1);
    typename boost::uniform_int<size_t> int_dist(0, boost::num_vertices(out_g) - 1);
#ifdef DETERMINISTIC_RNG
    typename boost::mt19937 rng(time(NULL) - int(p * 100));
    typename boost::variate_generator< boost::mt19937&, boost::uniform_real<decimal> >
      random_prob(rng, percent);
    typename boost::variate_generator< boost::mt19937&, boost::uniform_int<size_t> >
      random_vertex(rng, int_dist);
#else
    typename boost::random_device rng;
    typename boost::variate_generator< boost::random_device&, boost::uniform_real<decimal> >
      random_prob(rng, percent);
    typename boost::variate_generator< boost::random_device&, boost::uniform_int<size_t> >
      random_vertex(rng, int_dist);
#endif

    size_t num_edges_to_remove = 0;
    EdgeVector edges_to_remove;

    edge_iter ei, eiend;
    for (tie(ei, eiend) = boost::edges(out_g); ei != eiend; ++ei)
    {
      decimal random_probability = random_prob();

#ifdef CONAN_DEBUG
      std::cerr << "conan::randomize_graph_edges: current_source_vertex = " << boost::source(*ei, out_g)
                << ", current_target_vertex = " << boost::target(*ei, out_g) << std::endl
                << "conan::randomize_graph_edges: random_probability = " << random_probability << (random_probability < p? "... OK" : "") << std::endl;
#endif

      if (random_probability < p)
      {
        edges_to_remove.push_back(
            std::make_pair( boost::source(*ei, out_g), boost::target(*ei, out_g) ) );
        ++num_edges_to_remove;
      }
    }

    // actually replace edges
    for (EdgeVectorIter ei = edges_to_remove.begin(); ei != edges_to_remove.end(); ++ei)
    {
      vertex old_source_vertex = ei->first,
             old_target_vertex = ei->second;

      edge e = boost::edge( old_source_vertex, old_target_vertex, out_g ).first;

      decimal old_weight = out_g[ e ].weight;

      boost::remove_edge(e, out_g);

      vertex source_vertex, target_vertex;

      do
      {
        source_vertex = random_vertex();
        target_vertex = random_vertex();
      }
      while ( source_vertex == target_vertex || boost::edge( source_vertex, target_vertex, out_g ).second );

#ifdef CONAN_DEBUG
      std::cerr << "conan::randomize_graph_edges: source_vertex = " << source_vertex << ", target_vertex = " << target_vertex << std::endl;
#endif

      e = boost::add_edge( source_vertex, target_vertex, out_g ).first;

      if (maintain_connected)
      {
        size_t max_num_tries = 31;
        std::vector<int> component(boost::num_vertices(out_g));

        while (boost::connected_components(out_g, &component[0]) != 1)
        {
          if (--max_num_tries == 0)
          {
            boost::add_edge(old_source_vertex, old_target_vertex, old_weight, out_g);
            continue;
          }

          boost::remove_edge( e, out_g );
          source_vertex = random_vertex();
          target_vertex = random_vertex();

#ifdef CONAN_DEBUG
          std::cerr << "conan::randomize_graph_edges: source_vertex = " << source_vertex << ", target_vertex = " << target_vertex << std::endl;
#endif

          e = boost::add_edge( source_vertex, target_vertex, out_g ).first;
        }
      }

      out_g[ e ].weight = 1.0;
    }

    return out_g;
  }

  
  namespace detail {

    template <typename T1>
    bool is_in_edge_vector(
        T1 & value,
        std::vector<T1> & vec
        )
    {
      typedef typename std::vector<T1>::iterator VectorIter;

      for (VectorIter vi = vec.begin(), viend = vec.end(); vi != viend; ++vi)
      {
        if ( vi->m_source == value.m_source )
        {
          if ( vi->m_target == value.m_target )
          {
            return true;
          }
        }
        else if ( vi->m_source == value.m_target )
        {
          if ( vi->m_target == value.m_source )
          {
            return true;
          }
        }
      }

      return false;
    }

  } // detail

    
  template <class Graph>
  Graph randomize_graph_edges2(
      Graph out_g,
      decimal number_edges,
      bool maintain_connected // = false
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex;
    typedef typename GraphTraits::edge_descriptor edge;
    typedef typename GraphTraits::edge_iterator edge_iter;
    
#ifdef CONAN_DEBUG
    std::cerr << "conan::randomize_graph_edges2: number_edges = " << number_edges << std::endl
              << "conan::randomize_graph_edges2: maintain_connected = " << maintain_connected << std::endl;
#endif

    std::vector<edge> random_edges; // to check whether an edge has been created by the algorithm (and thus not remove it) or it is from the original graph

    size_t V = boost::num_vertices(out_g);

    if (V <= 2)
      throw std::runtime_error("conan::randomize_graph_edges2: There are no edges to randomize");
    if (number_edges < 0)
      throw std::runtime_error("conan::randomize_graph_edges2: number_edges must be a positive non-zero real number");
    else if (number_edges <= 1)
      number_edges *= V;

    // Replace some edges by random edges
    typename boost::uniform_int<size_t> int_dist(0, V - 1);
#ifdef DETERMINISTIC_RNG
    typename boost::mt19937 rng(time(NULL) - number_edges);
    typename boost::variate_generator< boost::mt19937&, boost::uniform_int<size_t> >
      random_vertex(rng, int_dist);
#else
    typename boost::random_device rng;
    typename boost::variate_generator< boost::random_device&, boost::uniform_int<size_t> >
      random_vertex(rng, int_dist);
#endif

    int num_replaced_edges = 0;
    while (num_replaced_edges < number_edges)
    {
      vertex current_vertex;
      edge edge_to_remove;

      do
      {
        current_vertex = random_vertex();

        size_t current_vertex_out_degree = boost::out_degree(current_vertex, out_g);

        int offset = 0;
        if (current_vertex_out_degree == 0)
        {
          if (maintain_connected)
            std::cerr << "conan::randomize_graph_edges2: Warning: a vertex have been disconnected" << std::endl;

          continue;
        }
        else if (current_vertex_out_degree != 1)
        {
          std::vector<decimal> p(current_vertex_out_degree, 1.0);

          offset = detail::weighted_choice(p, decimal(current_vertex_out_degree), rng);
        }

        edge_to_remove = *( boost::out_edges(current_vertex, out_g).first + offset );
      }
      while ( detail::is_in_edge_vector( edge_to_remove, random_edges ) );

#ifdef CONAN_DEBUG
      std::cerr << "conan::randomize_graph_edges2: source_vertex = " << boost::source( edge_to_remove, out_g )
                << ", target_vertex = " << boost::target( edge_to_remove, out_g ) << std::endl;
#endif

      vertex old_source_vertex = boost::source(edge_to_remove, out_g),
             old_target_vertex = boost::target(edge_to_remove, out_g);

      if (old_target_vertex == current_vertex)
      { // swap vertices
        vertex tmp = old_source_vertex;
        old_source_vertex = old_target_vertex;
        old_target_vertex = tmp;
      }

      if (maintain_connected)
      {
        decimal weight = out_g[ edge_to_remove ].weight;

        boost::remove_edge(edge_to_remove, out_g);

        std::vector<int> component(boost::num_vertices(out_g));

        size_t num_components = boost::connected_components(out_g, &component[0]);

        if (num_components > 1)
        { // Undo the removement
          boost::add_edge(old_source_vertex, old_target_vertex, weight, out_g);

          continue;
        }
      }
      else // maintain_connected == false
      {
        boost::remove_edge(edge_to_remove, out_g);
      }

      vertex target_random_vertex;

      do
      {
        target_random_vertex = random_vertex();
      }
      while ( boost::edge( current_vertex, target_random_vertex, out_g ).second ||
              target_random_vertex == current_vertex ||
              target_random_vertex == old_target_vertex );

      edge new_edge = boost::add_edge( current_vertex, target_random_vertex, 1.0, out_g ).first; // FIXME: why weight = 1.0?
      random_edges.push_back( new_edge );

#ifdef CONAN_DEBUG
      std::cerr << "conan::randomize_graph_edges2: Added new random edge between " << boost::source( new_edge, out_g )
                << " and " << boost::target( new_edge, out_g ) << std::endl;
#endif

      ++num_replaced_edges;

    } // while num_replaced_edges < number_edges

    return out_g;
  }


  template <class Graph>
  Graph transpose(
      const Graph & g
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_iterator vertex_iter;
    typedef typename GraphTraits::edge_descriptor edge;
    typedef typename GraphTraits::directed_category directed_category;

    BOOST_STATIC_ASSERT(!(boost::is_same<directed_category, boost::undirected_tag>::value));

    size_t V = boost::num_vertices(g);

    Graph out_g(V);

    vertex_iter source_vi, viend;
    for (tie(source_vi, viend) = boost::vertices(g); source_vi != viend; ++source_vi)
    {
      for (vertex_iter target_vi = boost::vertices(g).first; target_vi != viend; ++target_vi)
      {
        if (source_vi == target_vi)
          continue;
        if (boost::edge(*source_vi, *target_vi, g).second)
        {
          edge e = boost::add_edge(*target_vi, *source_vi, out_g).first;
          out_g[e].weight = 1.0;
        }
      }
    }

    return out_g;
  }


  template <class Graph>
  Graph complement(
      const Graph & g
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_iterator vertex_iter;
    typedef typename GraphTraits::edge_descriptor edge;
    
    size_t V = boost::num_vertices(g);

    Graph out_g(V);

    bool undirected = boost::is_undirected(g);

    size_t counter = 1;
    vertex_iter source_vi, viend;
    for (tie(source_vi, viend) = boost::vertices(g); source_vi != viend; ++source_vi, ++counter)
    {
      for (vertex_iter target_vi = boost::vertices(g).first + (undirected? counter : 0);
           target_vi != viend; ++target_vi)
      {
        if (source_vi == target_vi) // for directed graphs
          continue;

        if (not boost::edge(*source_vi, *target_vi, g).second)
        {
          if (boost::edge(*source_vi, *target_vi, out_g).second)
            std::cerr << "conan::complement: Warning!! complement edge has already been added to output graph." << std::endl;

          boost::add_edge(*source_vi, *target_vi, 1.0, out_g);
        }
      }
    }

    return out_g;
  }


  template <class Graph>
  Graph random_attack(
      Graph g,
      decimal num_vertices
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex;

    size_t num_vertices_uint = size_t(num_vertices);

    if (num_vertices < 0)
      throw std::runtime_error("random_attack: num_vertices cannot be negative");
    else if (num_vertices <= 1)
        num_vertices_uint = size_t(num_vertices * boost::num_vertices(g)); // FIXME: el redondeo es medio extranho

#ifdef DETERMINISTIC_RNG
    typename boost::mt19937 rng(time(NULL) - num_vertices);
#else
    typename boost::random_device rng;
#endif

    while (num_vertices_uint--)
    {
      if (boost::num_vertices(g) <= 1)
        return Graph();

      typename boost::uniform_int<size_t> int_dist(0, boost::num_vertices(g) - 1);
#ifdef DETERMINISTIC_RNG
      typename boost::variate_generator< boost::mt19937&, boost::uniform_int<size_t> >
        random_vertex(rng, int_dist);
#else
      typename boost::variate_generator< boost::random_device&, boost::uniform_int<size_t> >
        random_vertex(rng, int_dist);
#endif

      vertex v = vertex(random_vertex());

      boost::clear_vertex(v, g);
      boost::remove_vertex(v, g);
    }

    return g;
  }

} // conan

#endif //TRANSFORMATIONS_HPP
