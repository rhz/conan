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

#ifndef ENTROPY_HPP
#define ENTROPY_HPP
#include <conan/config.hpp>
#include <conan/subgraph.hpp>
#include <conan/utils.hpp> // for get_adj_matrix

#define HAVE_INLINE
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>


namespace conan {
  namespace detail {
    bool are_vectors_different(
        gsl_vector * v1,
        gsl_vector * v2,
        decimal tolerance
        )
    {
      size_t v1_size = v1->size,
             v2_size = v2->size;
             
      if (v1_size != v2_size)
      {
        throw std::runtime_error("Invalid size for input vectors");
      }

      for (size_t i = 0; i < v1_size; ++i)
      {
        double v1_val = gsl_vector_get(v1, i),
               v2_val = gsl_vector_get(v2, i);

        if ( (v1_val < v2_val && v1_val < v2_val - tolerance) || (v1_val > v2_val && v1_val > v2_val + tolerance) )
          return true;
      }
      
      return false;
    }


    template <class Graph>
    bool is_graph_not_bipartite_and_connected(
        const Graph & g
        )
    { 
      typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iter;
      std::pair<vertex_iter, vertex_iter> vp = vertices(g);

      // Could we use V to store the number of colored vertexes and test for it being 0?
      size_t V = boost::num_vertices(g);
      size_t visited = 0;

      std::vector<bool> colors(V, false);
      std::vector<bool> was_painted(V, false);
      std::vector<bool> was_visited(V, false);
      if ( visit_vertex(g, colors, was_painted, was_visited, *vp.first, visited) )
        return false;
        
      return V == visited;
    }


    // FIXME: Depth first would be better but memory usage raises concerns
    template <class Graph>
    bool visit_vertex(
        const Graph & g,
        std::vector<bool> & colors,
        std::vector<bool> & was_painted,
        std::vector<bool> & was_visited,
        typename boost::graph_traits<Graph>::vertex_descriptor v,
        size_t & visited
      )
    {
      typedef typename boost::graph_traits<Graph> GraphTraits;
      typedef typename GraphTraits::vertex_descriptor vertex;
      typedef typename GraphTraits::adjacency_iterator adj_vertex_iter;
      typename boost::property_map<Graph, boost::vertex_index_t>::type index =
          boost::get(boost::vertex_index, g);
      
      int vi = index[v];
      int ni;
      visited++;
      was_visited[vi] = true;
      
      adj_vertex_iter ai, ai_end;
      for (tie(ai, ai_end) = boost::adjacent_vertices(v, g); ai != ai_end; ++ai)
      {
         ni = index[*ai];
         if (was_painted[ni])
         {
           if (colors[ni] == colors[vi])
           {
             return false;
           }
         } else
         {
           was_painted[ni] = true;
           colors[ni] = ! colors[vi];
         }
      }
      
      for (tie(ai, ai_end) = boost::adjacent_vertices(v, g); ai != ai_end; ++ai)
      {
        ni = index[*ai];
        if (not was_visited[ni])
          if (not visit_vertex(g, colors, was_painted, was_visited, *ai, visited) )
            return false;
      }
      
      return true;
    }


    template <class Graph>
    void compute_weighted_transition_matrix(
        const Graph & g,
        gsl_matrix * P
        )
    {
      typedef typename boost::graph_traits<Graph> GraphTraits;
      typedef typename GraphTraits::out_edge_iterator out_edge_iter;

      size_t V = boost::num_vertices(g);

      if (P->size1 != V || P->size2 != V)
        throw std::runtime_error("Invalid size for input matrix P");

      // Compute matrix P (transitions), where p_ij = 0 iff a_ij = 0, otherwise p_ij = w_ij / sum (w_ij, j = 0...N)
      for (size_t r = 0; r < V; ++r)
      {
        decimal sum_w_ij = 0;
        out_edge_iter oie, oie_end;
        
        for (tie(oie, oie_end) = boost::out_edges(r, g); oie != oie_end; ++oie)
          sum_w_ij += g[*oie].weight;

        for (size_t c = 0; c < V; ++c)
        {
          if (boost::edge(r, c, g).second)
          {
            decimal w_ij = g[ boost::edge(r, c, g).first ].weight;

            gsl_matrix_set(P, r, c, w_ij / sum_w_ij);
          }
          else
          {
            gsl_matrix_set(P, r, c, 0.0);
          }
        }
      }

#ifdef CONAN_DEBUG
      for (size_t r = 0; r < V; ++r)
      {
        double sum = 0.0;
        for (size_t c = 0; c < V; ++c)
          sum += gsl_matrix_get(P, r, c);
        std::cout << "sum = " << sum << std::endl;
      }
#endif

      return;
    }


    void compute_stationary_states(
        gsl_matrix * P, // transition matrix
        gsl_vector * pi, // row vector
        decimal tolerance,
        int max_steps = 1000000
        )
    {
      size_t V = P->size1;

      if (P->size2 != V || pi->size != V)
        throw std::runtime_error("Invalid size for input matrix or vector");

      // Compute matrix pi (stationary states)
      gsl_vector * last_pi = gsl_vector_calloc(V); // row vector

      //gsl_matrix_set(pi, 0, 0, 1.0);
      /*
      boost::mt19937 rng(time(NULL));
      boost::uniform_int<int> dist(1, V);
      boost::variate_generator< boost::mt19937&, boost::uniform_int<int> >
        random_value(rng, dist);
      double sum = 0.0;
      for (size_t i = 0; i < V; ++i)
      {
        double v = random_value();
        gsl_matrix_set(pi, 0, i, v);
        sum += v;
      }

      for (size_t i = 0; i < V; ++i)
        gsl_matrix_set( pi, 0, i, gsl_matrix_get(pi, 0, i) / sum );
      */
      for (size_t i = 0; i < V; ++i)
        gsl_vector_set( pi, i, double(i + 1) * 2.0 / double(V * (V + 1)) );

#ifdef CONAN_DEBUG
      double sum = 0.0;
      for (size_t i = 0; i < V; ++i)
        sum += gsl_vector_get(pi, i);
      if (sum != 1.0)
        std::cout << "compute_stationary_states: Warning! sum [" << sum << "] != 1.0" << std::endl;
#endif

      do
      {
#ifdef CONAN_DEBUG
        double sum = 0.0;
        for (size_t i = 0; i < V; ++i)
        {
          std::cout << gsl_vector_get(pi, i) << " ";
          sum += gsl_vector_get(pi, i);
        }
        std::cout << "= " << sum << std::endl;
#endif
        gsl_vector_memcpy(last_pi, pi);

        gsl_blas_dgemv(CblasNoTrans, 1.0, P, last_pi, 0.0, pi);
      }
      while (are_vectors_different(pi, last_pi, tolerance) && max_steps-- > 0);

      gsl_vector_free(last_pi);

      return;
    }

  } // detail


  template <class Graph>
  decimal graph_entropy_using_stationary_states(
      const Graph & g,
      decimal tolerance = 1E-8
      )
  {
    using namespace detail;
    size_t V = boost::num_vertices(g);

    std::vector<int> components(V);
    
    if ( is_graph_not_bipartite_and_connected(g) )
    {
      throw std::runtime_error("graph_entropy_using_stationary_states: input graph is bipartite or reducible.");
    }

    gsl_matrix * P = gsl_matrix_calloc(V, V);
    compute_weighted_transition_matrix(g, P);

    gsl_vector * pi = gsl_vector_calloc(V); // row vector
    compute_stationary_states(P, pi, tolerance);

    // Return H(P) = sum(pi_i * H_i, i = 0...N), where H_i = - sum(p_ij * log(p_ij), j = 0...N)
    decimal H_P = 0.0;

    for (size_t i = 0; i < V; ++i)
    {
      decimal pi_i = gsl_vector_get(pi, i);
      H_P -= (pi_i * conan_log2(pi_i));
    }

    gsl_matrix_free(P);
    gsl_vector_free(pi);

    return H_P;
  }

  // FIXME: Why not use the adjacency iterator?
  template <class Graph>
  decimal graph_entropy_using_random_walk(
      Graph & g,
      long long int num_steps = 0,
      long long int jump = 1000
      )
  {
    using namespace detail;
    typedef long long int int64;
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex;
    typedef typename GraphTraits::out_edge_iterator out_edge_iter;

    int V = boost::num_vertices(g);
    
    std::vector<int> components(V);
    
    //FIXME: an ad-hoc method isReducible would be a lot faster than calculating the components
    if (boost::connected_components(g, &components[0]) != 1)
    {
      throw std::runtime_error("graph_entropy_using_random_walk: input graph is reducible.");
    }
    
    std::vector<int64> counter(V, 0);

    gsl_matrix * P = gsl_matrix_calloc(V, V);
    compute_weighted_transition_matrix(g, P);

    boost::mt19937 rng(time(NULL));
    boost::uniform_int<int> dist(0, V-1);
    boost::variate_generator< boost::mt19937&, boost::uniform_int<int> >
      random_vertex(rng, dist);

    if (num_steps == 0)
      num_steps = V * 1E3;

    std::vector< std::vector<decimal> > weights(V);
    std::vector< std::vector<vertex> > adj_vertices(V);
    std::vector< decimal > sums(V);
    BOOST_FOREACH( vertex v, boost::vertices(g) )
    {
      decimal sum = 0.0;
      for (int i = 0; i < V; ++i)
      {
        double w = gsl_matrix_get(P, v, i);
        if (w != 0.0)
        {
          weights[v].push_back(w); sum += w;
          adj_vertices[v].push_back(i);
        }
      }
      sums[v] = sum;
    }

    vertex v;
    for (int64 i = 0; i < num_steps; ++i)
    {
      if (i % jump == 0)
      { // jump to a random vertex, this allows correct sampling even in the presence of strong nodes
        v = vertex(random_vertex());
      }
      else
      { // choose neighbour
        int index = detail::weighted_choice(weights[v], sums[v]);
        v = adj_vertices[v][index];
      }

      counter[v] += 1;
    }

    decimal H = 0.0; double sum = double(num_steps);
    for (int i = 0; i < V; ++i)
    {
      decimal p = double(counter[i]) / sum; // FIXME numerical instabilities?
      H += p * conan_log2(p);
    }

    return -H;
  }


  template <class Graph>
  decimal graph_entropy_using_eigenvalues(
      Graph& g
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::directed_category directed_category;
    typedef typename GraphTraits::edge_iterator edge_iter;

    BOOST_STATIC_ASSERT((boost::is_same<directed_category, boost::undirected_tag>::value));

    size_t V = boost::num_vertices(g);

    // Create adjacency matrix
    gsl_matrix * A = gsl_matrix_calloc(V, V);

    edge_iter ei, ei_end;
    for (tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
    {
      gsl_matrix_set(A, boost::source(*ei, g), boost::target(*ei, g), 1.0);
      gsl_matrix_set(A, boost::target(*ei, g), boost::source(*ei, g), 1.0);
    }

    // Compute eigenvalues
    gsl_eigen_symm_workspace * w = gsl_eigen_symm_alloc(V);
    gsl_vector * eval = gsl_vector_calloc(V);

    gsl_eigen_symm(A, eval, w);

    gsl_eigen_symm_free(w);
    gsl_matrix_free(A);

    decimal max_eigenvalue = 0.0;
    for (size_t i = 0; i < eval->size; ++i)
      if (fabs(gsl_vector_get(eval, i)) > max_eigenvalue)
        max_eigenvalue = fabs(gsl_vector_get(eval, i));

    gsl_vector_free(eval);
    
    return conan_log2(max_eigenvalue);
  }


#if 0
  struct eigenvalues { };

  struct stationary_states
  {
    decimal tolerance;

    stationary_states(decimal tol = 1E-8)
      : tolerance(tol) { }
  };

  struct random_walkers
  {
    size_t num_steps;
    size_t num_walkers;

    random_walkers(size_t s, size_t w)
      : num_steps(s), num_walkers(w) { }
  };


  template <class Graph, typename Method>
  decimal graph_entropy(
      Graph & g,
      const Method & m = eigenvalues()
      )
  {
    return graph_entropy_using_eigenvalues(g);
  }


  template <class Graph>
  decimal graph_entropy<Graph, stationary_states>(
      Graph & g,
      const stationary_states & m
      )
  {
    return graph_entropy_using_stationary_states(g, m.tolerance);
  }


  template <class Graph>
  decimal graph_entropy<Graph, random_walkers>(
      Graph & g,
      const random_walkers & m
      )
  {
    return graph_entropy_using_random_walkers(g, m.num_steps, m.num_walkers);
  }
#endif


  template <class Graph>
  decimal graph_entropy(
      Graph & g
      )
  {
    return graph_entropy_using_eigenvalues(g);
  }

  template <class Graph>
  decimal graph_entropy(
      Graph & g,
      decimal tolerance
      )
  {
    return graph_entropy_using_stationary_states(g, tolerance);
  }

  template <class Graph>
  decimal graph_entropy(
      Graph & g,
      size_t num_steps,
      size_t num_walkers
      )
  {
    return graph_entropy_using_eigenvalues(g, num_steps, num_walkers);
  }

} // conan

#endif //ENTROPY_HPP
