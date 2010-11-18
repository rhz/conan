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

#ifndef MAXIMUM_ENTROPY_HPP
#define MAXIMUM_ENTROPY_HPP
#include <conan/utils/pinv.hpp>
#include <conan/inference/base.hpp>
#include <conan/inference/bootstrap.hpp>
#include <ctime>


namespace conan { namespace inference {

  /**
   * @brief Calculate the interaction matrix from a data matrix using the covariance.
   * @warning You should free the output matrix after using it.
   *
   * @param A Data matrix.
   * @param covariance If it's true, the covariance matrix will be used. Otherwise, the correlation matrix
   *   (default = true).
   * @return Pointer to the interaction matrix (gsl_matrix *). You should free it.
   */
  template <class CoexpressionMeasure>
  void maximum_entropy(
      gsl_matrix * data_matrix,
      gsl_matrix * interaction_matrix,
      CoexpressionMeasure);

  /**
   * @brief Infer a network with the maximum entropy method using a given threshold for the interaction
   * strenght. If abs_values = true, the values for the interaction strenght will be taken absolute.
   * If you use threshold = 0 and abs_values = true, then no threshold will be applied over the interaction
   * matrix.
   * The data matrix must have one row for each gene/bird/whatever and one column
   * for each time-point/observation.
   *
   * @param A Data matrix.
   * @param threshold Threshold.
   * @param abs_values Whether the interaction strength values should be taken absolute (default = true).
   * @param normalize Whether to normalize the data matrix at first (default = true).
   */
  template <class Graph>
  Graph maximum_entropy_with_fixed_threshold(
      gsl_matrix * A,
      decimal threshold,
      bool abs_values = true,
      bool normalize = true);

  /**
   * @brief Infer a network with the maximum entropy method and select only the more powerful interactions.
   * If abs_values = true, the values for the interaction strenght will be taken absolute.
   * The data matrix must have one row for each gene/bird/whatever and one column
   * for each time-point/observation.
   *
   * @param A Data matrix.
   * @param num_edges Number of interactions to be selected.
   * @param abs_values Whether the interaction strength values should be taken absolute (default = true).
   * @param normalize Whether to normalize the data matrix at first (default = true).
   */
  template <class Graph>
  Graph maximum_entropy_with_fixed_num_edges(
      gsl_matrix * A,
      size_t num_edges,
      bool abs_values = true,
      bool normalize = true);

  /**
   * @brief Infer a network with the maximum entropy method using a threshold for the interaction strenght
   * calculated from the mean and standard deviation of the values.
   *
   * @param A Data matrix.
   * @param threshold All the values within the range [mean - threshold * sd; mean + threshold * sd] will be
   *   considered as non-interacting. The values out of this range will maintain their original values for the
   *   edge weight.
   * @param normalize Whether to normalize the data matrix at first (default = true).
   */
  template <class Graph>
  Graph maximum_entropy_with_gaussian_threshold(
      gsl_matrix * A,
      decimal threshold,
      bool normalize = true);

  /**
   * @brief Infer a network with the maximum entropy method using a threshold for the interaction strenght
   * calculated from the mean and standard deviation of the bootstrap values. In contrast with the function
   * maximum_entropy_with_gaussian_threshold, a different threshold will be computed for every pair of
   * genes/birds/whatever.
   *
   * @param A Data matrix.
   * @param num_steps Number of bootstrap steps.
   * @param threshold All the values within the range [mean - threshold * sd; mean + threshold * sd] will be
   *   considered as non-interacting. The values out of this range will maintain their original values for the
   *   edge weight.
   * @param normalize Whether to normalize the data matrix at first (default = true).
   */
  template <class Graph, class CoexpressionMeasure>
  Graph maximum_entropy_with_bootstrapping(
      gsl_matrix * A,
      size_t num_steps,
      decimal threshold,
      CoexpressionMeasure filter,
      bool normalize = true,
      bool use_tmp_file = false);



  template <class CoexpressionMeasure>
  void maximum_entropy(
      gsl_matrix * data_matrix,
      gsl_matrix * interaction_matrix,
      CoexpressionMeasure filter
      )
  {
    size_t num_elems = data_matrix->size1;

    gsl_matrix * coexpr = gsl_matrix_calloc(num_elems, num_elems);

    coexpr_matrix_generator(data_matrix, coexpr, filter);

#ifdef CONAN_DEBUG
    // DEBUG
    static bool first_time = true; // DEBUG

    if (first_time)
    {
      FILE * coexpr_matrix_outfile = fopen("coexpr.out", "wb");
      if (coexpr_matrix_outfile)
      {
        std::cout << "Writing coexpression matrix to file coexpr.out (binary)" << std::endl;
        gsl_matrix_fwrite(coexpr_matrix_outfile, coexpr);
        fclose(coexpr_matrix_outfile);
      }

      std::ofstream coexpr_matrix_ofstream("coexpr.txt", std::ios_base::app);
      if (coexpr_matrix_ofstream.is_open())
      {
        std::cout << "Writing coexpression matrix to file coexpr.tmp (plain-text)" << std::endl;
        for (size_t r = 0; r < num_elems; ++r)
        {
          coexpr_matrix_ofstream << gsl_matrix_get(coexpr, r, 0);
          for (size_t c = 1; c < num_elems; ++c)
            coexpr_matrix_ofstream << '\t' << gsl_matrix_get(coexpr, r, c);
          coexpr_matrix_ofstream << std::endl;
        }
        coexpr_matrix_ofstream << std::endl;
        coexpr_matrix_ofstream.close();
      }
    }
    // END DEBUG
#endif

    conan::linalg::pinv(coexpr, interaction_matrix); // compute the interaction matrix

#ifdef CONAN_DEBUG
    // DEBUG
    if (first_time)
    {
      FILE * interaction_matrix_outfile = fopen("interaction.out", "wb");
      if (interaction_matrix_outfile)
      {
        std::cout << "Writing interaction matrix to file interaction.out (binary)" << std::endl;
        gsl_matrix_fwrite(interaction_matrix_outfile, interaction_matrix);
        fclose(interaction_matrix_outfile);
      }

      std::ofstream interaction_matrix_ofstream("interaction_matrix.txt", std::ios_base::app);
      if (interaction_matrix_ofstream.is_open())
      {
        std::cout << "Writing interaction matrix to file interaction_matrix.txt (plain-text)" << std::endl;
        for (size_t r = 0; r < num_elems; ++r)
        {
          interaction_matrix_ofstream << gsl_matrix_get(interaction_matrix, r, 0);
          for (size_t c = 1; c < num_elems; ++c)
            interaction_matrix_ofstream << '\t' << gsl_matrix_get(interaction_matrix, r, c);
          interaction_matrix_ofstream << std::endl;
        }
        interaction_matrix_ofstream << std::endl;
        interaction_matrix_ofstream.close();
      }
    }

    first_time = false;
    // END DEBUG
#endif

    gsl_matrix_free(coexpr);

    return;
  }


  template <class Graph>
  Graph maximum_entropy_with_fixed_threshold(
      gsl_matrix * A,
      decimal threshold,
      bool abs_values, // = true
      bool normalize // = true
      )
  {
    size_t num_elems = A->size1;

    if (normalize)
      normalize_data_matrix(A);

    gsl_matrix * interaction = gsl_matrix_calloc(num_elems, num_elems);
    maximum_entropy(A, interaction, covarianceS());

    gsl_matrix * adj_matrix = gsl_matrix_alloc(num_elems, num_elems);
    calc_adj_matrix_with_fixed_threshold(interaction, adj_matrix, threshold, abs_values);

    Graph out_g( conan::make_graph_from_adj_matrix<Graph>(adj_matrix) );

    gsl_matrix_free(interaction);
    gsl_matrix_free(adj_matrix);

    return out_g;
  }


  template <class Graph>
  Graph maximum_entropy_with_fixed_num_edges(
      gsl_matrix * A,
      size_t num_edges,
      bool abs_values, // = true
      bool normalize // = true
      )
  {
    size_t num_elems = A->size1;

    if (normalize)
      normalize_data_matrix(A);

    gsl_matrix * interaction = gsl_matrix_calloc(num_elems, num_elems);
    maximum_entropy(A, interaction, covarianceS());

    // Compute the threshold
    double threshold = 0;

    if (interaction->tda != interaction->size2)
      throw std::runtime_error("interaction->tda != interaction->size2");

    {
      size_t n_square = num_elems * num_elems;

      gsl_vector_view values = gsl_vector_view_array(interaction->data, n_square);
      gsl_vector * val_vec = &values.vector;

      std::vector<double> values_vector(n_square);
      std::copy(val_vec->data, val_vec->data + val_vec->size + 1, values_vector.begin());

      if (abs_values)
        for (std::vector<double>::iterator vi = values_vector.begin(); vi != values_vector.end(); ++vi)
          *vi = fabs(*vi);

#ifdef CONAN_DEBUG
      std::cout << "values =";
      for (std::vector<double>::iterator vi = values_vector.begin(); vi != values_vector.end(); ++vi)
        std::cout << ' ' << *vi;
      std::cout << std::endl << std::endl;
#endif
      
      std::sort(values_vector.begin(), values_vector.end());
      threshold = values_vector[n_square - num_edges];
    }

    // Compute the adjacency matrix
    gsl_matrix * adj_matrix = gsl_matrix_alloc(num_elems, num_elems);
    calc_adj_matrix_with_fixed_threshold(interaction, adj_matrix, threshold, abs_values);

#ifdef CONAN_DEBUG
    std::cout << "adjacency matrix:" << std::endl;
    for (size_t r = 0; r < num_elems; ++r)
    {
      std::cout << gsl_matrix_get(adj_matrix, r, 0);
      for (size_t c = 1; c < num_elems; ++c)
        std::cout << '\t' << gsl_matrix_get(adj_matrix, r, c);
      std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;
#endif

    Graph out_g( conan::make_graph_from_adj_matrix<Graph>(adj_matrix) );

    gsl_matrix_free(interaction);
    gsl_matrix_free(adj_matrix);

    return out_g;
  }


  template <class Graph>
  Graph maximum_entropy_with_gaussian_threshold(
      gsl_matrix * A,
      decimal threshold,
      bool normalize // = true
      )
  {
    size_t num_elems = A->size1;

    if (normalize)
      normalize_data_matrix(A);

    gsl_matrix * interaction = gsl_matrix_calloc(num_elems, num_elems);
    maximum_entropy(A, interaction, covarianceS());

#ifdef CONAN_DEBUG
    std::cout << "interaction matrix:" << std::endl;
    for (size_t r = 0; r < num_elems; ++r)
    {
      std::cout << gsl_matrix_get(interaction, r, 0);
      for (size_t c = 1; c < num_elems; ++c)
        std::cout << '\t' << gsl_matrix_get(interaction, r, c);
      std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;
#endif

    gsl_matrix * adj_matrix = gsl_matrix_alloc(num_elems, num_elems);
    calc_adj_matrix_with_gaussian_threshold(interaction, adj_matrix, threshold);

    gsl_matrix_free(interaction);

    Graph out_g( conan::make_graph_from_adj_matrix<Graph>(adj_matrix) );

    gsl_matrix_free(adj_matrix);

    return out_g;
  }


  template <class Graph, class CoexpressionMeasure>
  Graph maximum_entropy_with_bootstrapping(
      gsl_matrix * A,
      size_t num_steps,
      decimal threshold,
      CoexpressionMeasure filter,
      bool normalize, // = true
      bool use_tmp_file // = false
      )
  {
    if (normalize)
      normalize_data_matrix(A);

    gsl_matrix * adj_matrix = bootstrap_and_calc_adj_matrix(A, num_steps, threshold, filter,
        &maximum_entropy<CoexpressionMeasure>, use_tmp_file);

    // Create graph
    Graph out_g( conan::make_graph_from_adj_matrix<Graph>(adj_matrix) );

    gsl_matrix_free(adj_matrix);

    return out_g;
  }

}} // conan::inference

#endif // MAXIMUM_ENTROPY_HPP
