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

#ifndef DEMETRIUS_MANKE_HPP
#define DEMETRIUS_MANKE_HPP
#include <conan/graphs.hpp>
#include <conan/properties/entropy.hpp>
#include <conan/utils.hpp>
#include <conan/utils/combinatoria.hpp>

#include <limits>
#include <iomanip>


namespace conan {

  // ***** Declarations *****
  
  /**
   * @brief Return a graph generated based on the Demetrius-Manke model.
   *
   * @param V Number of edges.
   * @param T
   * @param edges_by_vertex
   * @param each_edge_as_separete_step
   * @return Demetrius-Manke graph.
   */
  template <class Graph>
  Graph generate_Demetrius_Manke_network(
      size_t V,
      decimal T,
      size_t edges_by_vertex,
      bool each_edge_as_separete_step = false
      );
  
  /**
   * @brief Expand a graph based on Demetrius-Manke model's rules
   * (see conan::generate_Demetrius_Manke_graph function).
   *
   * @param g Base graph.
   * @param T
   * @param edges_by_vertex The number of edges to add to the base graph.
   * @param each_edge_as_separete_step
   */
  template <class Graph>
  void add_vertex_as_Demetrius_Manke(
      Graph & g,
      decimal T,
      size_t edges_by_vertex,
      bool each_edge_as_separete_step = false
      );


  // ***** Definitions ******
  
  template <class Graph>
  struct ChildGraph
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex;

    Graph g;
    std::vector<vertex> linked_vertices;
    decimal H;
    decimal probability;

    ChildGraph(const Graph & child_g, const std::vector<vertex> & v)
      : g(child_g), linked_vertices(v), H(0.0), probability(0.0) { }
  };

  namespace detail {

    template <class Graph>
    inline decimal get_probability(
        const ChildGraph<Graph> & cg
        )
    { return cg.probability; }

  } // detail


  template <class Graph>
  size_t generate_all_combinations_vertices(
      Graph & g,
      size_t edges_by_vertex,
      std::vector< std::vector<size_t> > & all_combinations
      )
  { // This function assumes that graph's vertex_descriptors are positive integers
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex;

    typedef typename std::vector<vertex> VertexSubset;
    typedef typename std::vector< VertexSubset > VertexSubsetList;

    size_t V = boost::num_vertices(g);

    decimal num_possible_combinations_float = detail::binomial_coefficient(V, edges_by_vertex);
    if (num_possible_combinations_float > std::numeric_limits<size_t>::max())
      throw std::runtime_error("Number of possible combinations is too big");

    size_t num_possible_combinations = size_t(num_possible_combinations_float);
    size_t real_num_possible_combinations = 0;

    all_combinations.resize(num_possible_combinations);

    VertexSubset aux(edges_by_vertex);

    for (size_t i = 0; i < edges_by_vertex; ++i)
      aux[i] = i;

    vertex last_vertex = *(boost::vertices(g).second - 1);

    for (typename VertexSubsetList::iterator li = all_combinations.begin(); li != all_combinations.end(); ++li)
    {
#ifdef CONAN_DEBUG
      std::cout << "aux = ";
      for (typename VertexSubset::iterator vi = aux.begin(); vi != aux.end(); ++vi)
        std::cout << *vi << ' ';
      std::cout << std::endl;
#endif

      *li = aux;
      ++real_num_possible_combinations;

      for (int pos = edges_by_vertex - 1; pos >= 0; --pos)
      {
        ++aux[pos];

        if (aux[pos] > last_vertex && edges_by_vertex == 1)
        { // little workaround for edges_by_vertex == 1
          all_combinations.resize(real_num_possible_combinations);
          return real_num_possible_combinations;
        }

        if (aux[pos] <= last_vertex)
        {
          for (++pos; pos < int(edges_by_vertex); ++pos)
          {
            aux[pos] = aux[pos - 1] + 1;

            if (aux[pos] > last_vertex)
            {
              all_combinations.resize(real_num_possible_combinations);
              return real_num_possible_combinations;
            }
          }
          break;
        }
      }
    }

    return real_num_possible_combinations;
  }


  template <class Graph>
  void generate_all_child_graphs(
      Graph & g,
      size_t edges_by_vertex,
      std::vector< ChildGraph<Graph> > & child_graphs,
      size_t v = 0
      )
  { // If v == 0, then a new vertex will be added to the graph, otherwise this vertex will be used.
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_iterator vertex_iter;
    typedef typename GraphTraits::vertex_descriptor vertex;

    typedef typename std::vector<vertex> VertexSubset;
    typedef typename std::vector< VertexSubset > VertexSubsetList;

    if (edges_by_vertex > boost::num_vertices(g))
      throw std::runtime_error("edges_by_vertex parameter must be less than the number of vertices of the graph");

    child_graphs.clear();

    // generate all possible combinations of vertices
    VertexSubsetList all_combinations;
#ifdef CONAN_DEBUG
    size_t n = generate_all_combinations_vertices(g, edges_by_vertex, all_combinations);
    std::cout << "n = " << n << std::endl;
#else
    generate_all_combinations_vertices(g, edges_by_vertex, all_combinations);
#endif

    // generate all child graphs
    for (typename VertexSubsetList::iterator li = all_combinations.begin(); li != all_combinations.end(); ++li)
    {
      Graph child_g(g);
      vertex new_vertex;
      if (v == 0)
        new_vertex = boost::add_vertex(child_g);
      else
        new_vertex = v;

      bool vertex_subset_contains_new_vertex = false;
      for (typename VertexSubset::iterator si = li->begin(); si != li->end(); ++si)
      {
        if (*si == new_vertex)
        {
          vertex_subset_contains_new_vertex = true;
          break;
        }
        boost::add_edge(new_vertex, *si, 1.0, child_g);
      }

      if (vertex_subset_contains_new_vertex)
        continue;

      ChildGraph<Graph> cg(child_g, *li);

      // Compute entropy... maybe this should stay in another function
      decimal entropy = graph_entropy(child_g);

      cg.H = entropy;

      child_graphs.push_back(cg);
    }

    return;
  }


  template <class Graph>
  decimal compute_child_graph_probabilities(
      std::vector< ChildGraph<Graph> > & child_graphs,
      decimal T,
      decimal tolerance = 1E-8
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex;

#if 0
    decimal H_min = std::numeric_limits<decimal>::max(),
            H_max = 0; // these values will be useful for normalizing H

    for (typename std::vector< ChildGraph<Graph> >::iterator child = child_graphs.begin(); child != child_graphs.end(); ++child)
    {
      if (H_max < child->H)
        H_max = child->H;
      if (H_min > child->H)
        H_min = child->H;
    }
#endif

    decimal prob_sum = 0.0,
            H_sum = 0.0;

    for (typename std::vector< ChildGraph<Graph> >::iterator child = child_graphs.begin(); child != child_graphs.end(); ++child)
      H_sum += child->H;

    for (typename std::vector< ChildGraph<Graph> >::iterator child = child_graphs.begin(); child != child_graphs.end(); ++child)
    {
      // Normalize H... H_i_normalized = (H_i - H_min) / (H_max - H_min)
      /*
      decimal H_i_normalized = 0.0;

      if (T < 0 && H_min > H_max - (.001 * H_max)) // workaround
        H_i_normalized = 0;
      else if (T >= 0 && H_min > H_max - (.001 * H_max)) // workaround
        H_i_normalized = 1;
      else
        H_i_normalized = (child->H - H_min) / (H_max - H_min);
      */
      decimal H_i_normalized = child->H / H_sum; // my own normalization

      // Compute the probability of selecting each new network
      // Note: It's not truly a probability, since sum(P_i) is not always 1 (maybe never indeed).
      //
      // P_i = H_i_normalized ^ T, for T >= 0
      // P_i = (1 - H_i_normalized) ^ T, for T < 0
      /*
      if (T >= 0)
        prob_sum += (child->probability = conan_pow(H_i_normalized, T));
      else
        prob_sum += (child->probability =
            (1.0 - H_i_normalized < tolerance)? 0 : conan_pow(1.0 - H_i_normalized, T));
      */
      prob_sum += (child->probability = conan_pow(H_i_normalized, T));

#ifdef CONAN_INFO
      std::string vertices_str;
      for (typename std::vector<vertex>::iterator vi = child->linked_vertices.begin(); vi != child->linked_vertices.end(); ++vi)
        vertices_str += to_string(*vi) + ',';
      vertices_str.erase( vertices_str.end() - 1 );

      std::cout << "H(" << vertices_str << ") = " << child->H
                << ", Hnorm(" << vertices_str << ") = " << H_i_normalized
                << ", P(" << vertices_str << ") = " << child->probability
                << std::endl;
#endif
    }

    return prob_sum;
  }


  template <class Graph, class RNG>
  size_t select_child_graph(
      std::vector< ChildGraph<Graph> > & child_graphs,
      decimal prob_sum,
      RNG & rng
      )
  {
    using namespace detail;
    return weighted_choice(child_graphs, prob_sum, rng, &get_probability<Graph>);
  }


  template <class Graph>
  void add_vertex_as_Demetrius_Manke(
      Graph & g,
      decimal T,
      size_t edges_by_vertex,
      bool each_edge_as_separete_step // = false
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_iterator vertex_iter;
    typedef typename GraphTraits::vertex_descriptor vertex;
    typedef typename GraphTraits::edge_descriptor edge;

    size_t V = boost::num_vertices(g);

    // check for special cases
    if (edges_by_vertex == 0)
    {
      throw std::runtime_error("It's meaningless to define edges_by_vertex as 0");
    }
    if (V == 0)
    {
      boost::add_vertex(g);
      return;
    }
    else if (V <= edges_by_vertex)
    { // connect new vertex to all vertices in the graph
      vertex new_vertex = boost::add_vertex(g);

      vertex_iter vi, viend;
      for (tie(vi, viend) = boost::vertices(g); vi != viend; ++vi)
        boost::add_edge(new_vertex, *vi, 1.0, g);
    }
    else if (each_edge_as_separete_step)
    {
      vertex new_vertex = 0;

#ifdef DETERMINISTIC_RNG
      boost::mt19937 rng(time(NULL) - int(T * 100) - edges_by_vertex);
#else
      boost::random_device rng;
#endif

      while (edges_by_vertex--)
      {
        std::vector< ChildGraph<Graph> > child_graphs;
        generate_all_child_graphs(g, 1, child_graphs, new_vertex);

        decimal prob_sum = compute_child_graph_probabilities(child_graphs, T);

        size_t selected_index = select_child_graph(child_graphs, prob_sum, rng);

        // Construct the selected child graph
        new_vertex = boost::add_vertex(g);

        for (typename std::vector<vertex>::iterator vi = child_graphs[selected_index].linked_vertices.begin();
             vi != child_graphs[selected_index].linked_vertices.end(); ++ vi)
        {
          vertex selected_vertex = *vi;

#ifdef CONAN_INFO
          std::cout << "selected_vertex = " << selected_vertex << std::endl;
#endif

          boost::add_edge(new_vertex, selected_vertex, 1.0, g);
        }
      } // while edges_by_vertex--
    }
    else
    {
      std::vector< ChildGraph<Graph> > child_graphs;
      generate_all_child_graphs(g, edges_by_vertex, child_graphs);

      decimal prob_sum = compute_child_graph_probabilities(child_graphs, T);

#ifdef DETERMINISTIC_RNG
      boost::mt19937 rng(time(NULL) - int(T * 100) - edges_by_vertex);
#else
      boost::random_device rng;
#endif

      size_t selected_index = select_child_graph(child_graphs, prob_sum, rng);

      // Construct the selected child graph
      vertex new_vertex = boost::add_vertex(g);

      for (typename std::vector<vertex>::iterator vi = child_graphs[selected_index].linked_vertices.begin();
           vi != child_graphs[selected_index].linked_vertices.end(); ++ vi)
      {
        vertex selected_vertex = *vi;

#ifdef CONAN_INFO
        std::cout << "selected_vertex = " << selected_vertex << std::endl;
#endif

        boost::add_edge(new_vertex, selected_vertex, 1.0, g);
      }
    }

    return;
  }


  template <class Graph>
  Graph generate_Demetrius_Manke_network(
      size_t V,
      decimal T,
      size_t edges_by_vertex,
      bool each_edge_as_separete_step // = false
      )
  {
    if (V == 0)
      throw std::runtime_error("The number of vertex must be a positive non-zero integer value");

    Graph g(1);
    --V;

    while (V--)
      add_vertex_as_Demetrius_Manke(g, T, edges_by_vertex, each_edge_as_separete_step);
    
    return g;
  }

} //conan

#endif //DEMETRIUS_MANKE_HPP
