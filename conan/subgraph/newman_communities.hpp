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
 * Important!!
 * Thought I copyright Conan and all its files,
 * the REAL AUTHOR of this file is Bryan Reynaert.
 * Thanks Bryan for making these functions work =)
 */

/* Modularity Algorithm due Mark Newman published as
 * "Modularity and community structure in networks" in PNAS
 * www.pnas.org/cgi/doi/10.1073/pnas.0601602103
*/

#ifndef NEWMAN_COMMUNITIES_HPP
#define NEWMAN_COMMUNITIES_HPP
#include <conan/graphs.hpp>
#include <conan/subgraph/modularity.hpp>

#define HAVE_INLINE
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>


namespace conan {

  /**
   * @brief Search for modules or communities within the graph.
   * It computes these modules with the algorithm published in
   * "Modularity and community structure in networks", Mark Newman, PNAS, 2006.
   *
   * @param g Input graph.
   * @param module Vector to store the module number to which each vertex belongs.
   *   This number is like a module ID. Each vertex will take part of only one module.
   * @return Number of modules found.
   */
  template <class Graph>
  void newman_communities(
      Graph &,
      std::vector< std::vector<size_t> > &,
      bool check_quality_function = true
      );

  template <class Graph>
  void newman_two_communities(
      Graph &,
      std::vector< std::vector<size_t> > &
      );

  template <class Graph>
  void get_submodules(
      Graph & g,
      std::vector<size_t> & vertex_subset,
      std::vector< std::vector<size_t> > & all_modules,
      bool check_quality_function
      );


  template <class Graph>
  void newman_two_communities(
      Graph & g,
      std::vector< std::vector<size_t> > & modules
      )
  {
    size_t V = boost::num_vertices(g);

    modules.clear();

    gsl_matrix * B = gsl_matrix_alloc(V, V);

    get_modularity_matrix(g, B);

    // compute the eigenvalues and eigenvectors of B
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(V);
    gsl_matrix * evec = gsl_matrix_alloc(V, V);
    gsl_vector * eval = gsl_vector_alloc(V);

    gsl_eigen_symmv(B, eval, evec, w);

    gsl_eigen_symmv_free(w);

    size_t max_eigenvalue_index = 0;
    double max_eigenvalue = 0;
    for (size_t i = 0; i < eval->size; ++i)
    {
      if (gsl_vector_get(eval, i) > max_eigenvalue)
      {
        max_eigenvalue = gsl_vector_get(eval, i);
        max_eigenvalue_index = i;
      }
    }

    gsl_vector_view max_eval_eigenvector = gsl_matrix_column(evec, max_eigenvalue_index);

    std::vector<size_t> community1, community2;

    for (size_t i = 0; i < max_eval_eigenvector.vector.size; ++i)
    {
      double val = gsl_vector_get(&max_eval_eigenvector.vector, i);

      if (val > 0)
        community1.push_back(i);
      else
        community2.push_back(i);
    }

    gsl_matrix_free(B);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);

    modules.push_back(community1);
    modules.push_back(community2);

    return;
  }

  template <class Graph>
  void get_submodules(
      Graph & g,
      std::vector<size_t> & vertex_subset,
      std::vector< std::vector<size_t> > & all_modules,
      bool check_quality_function
      )
  {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex;
    typedef typename std::vector<vertex> Module;
    typedef typename std::vector<Module> ModuleVector;

    Module module1;
    Module module2;

    double delta_modularity_value;

    size_t subgraph_size = vertex_subset.size();

    gsl_matrix * B = gsl_matrix_alloc(subgraph_size, subgraph_size);
    get_modularity_matrix(g, B, vertex_subset);

    // compute the eigenvalues and eigenvectors of B
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(subgraph_size);
    gsl_matrix * evec = gsl_matrix_alloc(subgraph_size, subgraph_size);
    gsl_vector * eval = gsl_vector_alloc(subgraph_size);

    gsl_eigen_symmv(B, eval, evec, w);

    gsl_eigen_symmv_free(w);

    int max_eigenvalue_index = -1;
    double max_eigenvalue = 0.0;
    for (size_t i = 0; i < eval->size; ++i)
    {
      if (gsl_vector_get(eval, i) > max_eigenvalue)
      {
        max_eigenvalue = gsl_vector_get(eval, i);
        max_eigenvalue_index = i;
      }
    }

#ifdef CONAN_DEBUG
    std::cerr << "conan::get_submodules: max_eigenvalue = " << max_eigenvalue << std::endl;
#endif

    if (max_eigenvalue_index == -1)
    {
      // Indivisible graph
      all_modules.push_back(vertex_subset);
      gsl_vector_free(eval);
      gsl_matrix_free(evec);
      gsl_matrix_free(B);
      return;
    }

    gsl_vector_view max_eval_eigenvector = gsl_matrix_column(evec, max_eigenvalue_index);

    for (size_t i = 0; i < max_eval_eigenvector.vector.size; ++i)
    {
      double val = gsl_vector_get(&max_eval_eigenvector.vector, i);
      if (val > 0)
        module1.push_back(vertex_subset[i]);
      else
        module2.push_back(vertex_subset[i]);
    }

    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_matrix_free(B);

    if ( check_quality_function && (module1.size() != 0 && module2.size() != 0) )
    {
      ModuleVector new_modules;
      new_modules.push_back(module1);
      new_modules.push_back(module2);
      delta_modularity_value = delta_modularity(g, new_modules);

#ifdef CONAN_DEBUG
      std::cerr << "conan::get_submodules: delta_modularity_value = " << delta_modularity_value << std::endl;
#endif

      if (delta_modularity_value <= 0.0)
      { // No good division found
        all_modules.push_back(vertex_subset);
        return;
      }
    }

    if ( module1.size() != 0 && module2.size() != 0)
    {
      get_submodules(g, module1, all_modules, check_quality_function);
      get_submodules(g, module2, all_modules, check_quality_function);
      return;
    }
    else if (module1.size() == 0)
    {
      all_modules.push_back(module2);
    }
    else
    {
      all_modules.push_back(module1);
    }

    return;
  }


  template <class Graph>
  void newman_communities(
      Graph & g,
      std::vector< std::vector<size_t> > & all_modules,
      bool check_quality_function // = true
      )
  {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex;
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_iterator vertex_iter;

    size_t V = boost::num_vertices(g);

    std::vector<vertex> all_vertices(V);

    vertex_iter vi, viend;
    tie(vi, viend) = boost::vertices(g);
    std::copy(vi, viend, all_vertices.begin());

    get_submodules(g, all_vertices, all_modules, check_quality_function);

    return;
  }

} // conan

#endif // NEWMAN_COMMUNITIES_HPP
