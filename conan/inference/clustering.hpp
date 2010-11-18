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

#ifndef INFERENCE_CLUSTERING_HPP
#define INFERENCE_CLUSTERING_HPP
#include <conan/inference/base.hpp>
#include <conan/inference/bootstrap.hpp>


namespace conan { namespace inference {

  /**
   * @brief Infer a network using the co-expression clustering method.
   *
   * @param A Data matrix.
   * @param threshold Threshold.
   * @param abs_values Whether the interaction strength values should be taken absolute (default = true).
   * @param covariance If it's true, the covariance matrix will be used. Otherwise, the correlation matrix
   *   (default = true).
   * @param normalize Whether to normalize the data matrix at first (default = true).
   * @return The resulting graph.
   */
  template <class Graph>
  Graph clustering_with_fixed_threshold(
      gsl_matrix * A,
      decimal threshold,
      bool abs_values = true,
      bool covariance = true,
      bool normalize = true);

  /**
   * @brief Infer a network using the co-expression clustering method with a threshold for the interaction
   * strenght (covariance matrix) calculated from the mean and standard deviation of the bootstrap values.
   *
   * @param A Data matrix.
   * @param num_steps Number of bootstrap steps.
   * @param threshold All the values within the range [mean - threshold * sd; mean + threshold * sd] will be
   *   considered as non-interacting. The values out of this range will maintain their original values for the
   *   edge weight.
   * @param normalize Whether to normalize the data matrix at first (default = true).
   */
  template <class Graph, class CoexpressionMeasure>
  Graph clustering_with_bootstrapping(
      gsl_matrix * A,
      size_t num_steps,
      decimal threshold,
      CoexpressionMeasure filter,
      bool normalize = true);


  template <class Graph>
  Graph clustering_with_fixed_threshold(
      gsl_matrix * A,
      decimal threshold,
      bool abs_values, // = true
      bool covariance, // = true
      bool normalize // = true
      )
  {
    size_t num_elems = A->size1;

    if (normalize)
      normalize_data_matrix(A);

    gsl_matrix * covar = gsl_matrix_calloc(num_elems, num_elems);
    if (covariance)
      calc_covariance_matrix(A, covar);
    else
      calc_correlation_matrix(A, covar); // for simplicity the correlation matrix will be temporarily named covar.

    gsl_matrix * adj_matrix = gsl_matrix_alloc(num_elems, num_elems);
    calc_adj_matrix_with_fixed_threshold(covar, adj_matrix, threshold, abs_values);

    gsl_matrix_free(covar);

    Graph out_g( conan::make_graph_from_adj_matrix<Graph>(adj_matrix) );

    gsl_matrix_free(adj_matrix);

    return out_g;
  }

#if 0
  namespace detail {

    template <class CoexpressionMeasure>
    inline void covariance(
        gsl_matrix * A,
        gsl_matrix * covar,
        CoexpressionMeasure
        )
    {
      BOOST_STATIC_ASSERT((boost::is_same<CoexpressionMeasure, covarianceS>::value));
      calc_covariance_matrix(A, covar);
      return covar;
    }

    template <class CoexpressionMeasure>
    inline void correlation(
        gsl_matrix * A,
        gsl_matrix * corr,
        CoexpressionMeasure
        )
    {
      BOOST_STATIC_ASSERT((boost::is_same<CoexpressionMeasure, correlationS>::value));
      calc_correlation_matrix(A, corr);
      return corr;
    }

  } // detail
#endif

  template <class Graph, class CoexpressionMeasure>
  Graph clustering_with_bootstrapping(
      gsl_matrix * A,
      size_t num_steps,
      decimal threshold,
      CoexpressionMeasure filter,
      bool normalize // = true
      )
  {
    if (normalize)
      normalize_data_matrix(A);

    gsl_matrix * adj_matrix;
    adj_matrix = bootstrap_and_calc_adj_matrix(A, num_steps, threshold, filter,
        &coexpr_matrix_generator<CoexpressionMeasure>);

#if 0
    if (covariance)
      adj_matrix = bootstrap_and_calc_adj_matrix(A, num_steps, threshold, covarianceS(), &detail::covariance);
    else
      adj_matrix = bootstrap_and_calc_adj_matrix(A, num_steps, threshold, correlationS(), &detail::correlation);
#endif

    // Create graph
    Graph out_g( conan::make_graph_from_adj_matrix<Graph>(adj_matrix) );

    gsl_matrix_free(adj_matrix);

    return out_g;
  }

}} // conan::inference

#endif // INFERENCE_CLUSTERING_HPP
