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

#ifndef BASE_HPP
#define BASE_HPP
#include <conan/graphs.hpp>
#include <conan/utils/matrix_mask.hpp>

#define HAVE_INLINE
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics_double.h>


namespace conan { namespace inference {

  using namespace conan::linalg;

#if 0
  enum CoexpressionMeasure {
    covariance = 0,
    covariance_pairwise_complete_obs = 1,
    correlation = 2
  };
#endif

  struct covarianceS
  {
    bool shuffle_by_column;

    covarianceS(
        bool s = false
        )
      : shuffle_by_column(s) { }
  };

  struct correlationS
  {
    bool shuffle_by_column;

    correlationS(
        bool s = false
        )
      : shuffle_by_column(s) { }
  };

  struct covariance_pairwise_complete_obs_S
  {
    matrix_mask & mask;
    size_t min_num_complete_obs;
    bool shuffle_by_column;

    covariance_pairwise_complete_obs_S(
        matrix_mask & m,
        size_t min,
        bool s = false
        )
      : mask(m), min_num_complete_obs(min), shuffle_by_column(s) { }
  };


  /**
   * @brief Normalize input data matrix.
   * Normalize each row of a matrix, ie., mean = 0, std. dev. = 1.
   *
   * @param A Data matrix.
   */
  void normalize_data_matrix(
      gsl_matrix * A);

  /**
   * @brief Calculate the covariance matrix from a data matrix.
   * The data matrix must have one row for each gene/bird/whatever and one column for each time point.
   *
   * @param A Data matrix.
   * @param covar Matrix for storing the computed covariance matrix.
   */
  void calc_covariance_matrix(
      gsl_matrix * A,
      gsl_matrix * covar);

  /**
   * @brief Calculate the Pearson correlation matrix for a data matrix.
   * The data matrix must have one row for each gene/bird/whatever and one column
   * for each time-point/observation.
   *
   * @param A Data matrix.
   * @param covar Matrix for storing the computed covariance matrix.
   */
  void calc_correlation_matrix(
      gsl_matrix * A,
      gsl_matrix * corr);

  /**
   *
   */
  void calc_covariance_matrix_from_pairwise_complete_obs(
      gsl_matrix * A,
      matrix_mask & Amask,
      size_t min_num_complete_obs,
      gsl_matrix * covar);

  /**
   * @brief Calculate adjacency matrix from an interaction matrix.
   * If you use threshold = 0 and abs_values = true, then no threshold will be applied over the interaction
   * matrix.
   *
   * @param interaction Interaction matrix.
   * @param adj_matrix Adjacency matrix.
   * @param threshold Threshold (default = 0).
   * @param abs_values Whether the interaction strength values should be taken absolute (default = true).
   */
  void calc_adj_matrix_with_fixed_threshold(
      gsl_matrix * interaction,
      gsl_matrix * adj_matrix,
      decimal threshold = 0,
      bool abs_values = true);

  void calc_adj_matrix_with_gaussian_threshold(
      gsl_matrix * interaction,
      gsl_matrix * adj_matrix,
      decimal threshold);


  void normalize_data_matrix(
      gsl_matrix * A
      )
  {
    size_t nrows = A->size1,
           ncols = A->size2;

    for (size_t r = 0; r < nrows; ++r)
    {
      gsl_vector_view row = gsl_matrix_row(A, r);
      gsl_vector * row_vector = &row.vector;

      double mean = gsl_stats_mean(row_vector->data, row_vector->stride, ncols),
             stddev = gsl_stats_sd(row_vector->data, row_vector->stride, ncols);

      for (size_t c = 0; c < ncols; ++c)
        gsl_matrix_set(A, r, c, (gsl_matrix_get(A, r, c) - mean) / stddev);
    }

    return;
  }



  template <class CoexpressionMeasure>
  void coexpr_matrix_generator(
      gsl_matrix * A,
      gsl_matrix * coexpr,
      CoexpressionMeasure
      )
  { // Invalid coexpression measure selected
    BOOST_STATIC_ASSERT((boost::is_same<CoexpressionMeasure, covarianceS>::value ||
                         boost::is_same<CoexpressionMeasure, correlationS>::value ||
                         boost::is_same<CoexpressionMeasure, covariance_pairwise_complete_obs_S>::value));
  }


  template <>
  void coexpr_matrix_generator(
      gsl_matrix * A,
      gsl_matrix * coexpr,
      covarianceS
      )
  { return calc_covariance_matrix(A, coexpr); }


  template <>
  void coexpr_matrix_generator(
      gsl_matrix * A,
      gsl_matrix * coexpr,
      correlationS
      )
  { return calc_correlation_matrix(A, coexpr); }


  template <>
  void coexpr_matrix_generator(
      gsl_matrix * A,
      gsl_matrix * coexpr,
      covariance_pairwise_complete_obs_S filter
      )
  { return calc_covariance_matrix_from_pairwise_complete_obs(
      A, filter.mask, filter.min_num_complete_obs, coexpr); }



  void calc_covariance_matrix(
      gsl_matrix * A,
      gsl_matrix * covar
      )
  {
    size_t num_obs = A->size2;
    size_t num_elems = A->size1;

    if (covar->size1 != num_elems || covar->size2 != num_elems)
      throw std::runtime_error("invalid lenght for covariance matrix");

    for (size_t r1 = 0; r1 < A->size1; ++r1) // for each row
    {
      gsl_vector_view row1 = gsl_matrix_row(A, r1);
      gsl_vector * row1_vector = &row1.vector;

      for (size_t r2 = 0; r2 < A->size1; ++r2) // for each row
      {
        gsl_vector_view row2 = gsl_matrix_row(A, r2);
        gsl_vector * row2_vector = &row2.vector;

        double covar1_2 = gsl_stats_covariance(row1_vector->data, row1_vector->stride,
                                               row2_vector->data, row2_vector->stride,
                                               num_obs);
        gsl_matrix_set(covar, r1, r2, covar1_2);
      }
    }

    return;
  }


  void calc_correlation_matrix(
      gsl_matrix * A,
      gsl_matrix * corr
      )
  { // the same as before, but calls gsl_stats_correlation instead of gsl_stats_covariance
    size_t num_obs = A->size2;
    size_t num_elems = A->size1;

    if (corr->size1 != num_elems || corr->size2 != num_elems)
      throw std::runtime_error("invalid lenght for correlation matrix");

    for (size_t r1 = 0; r1 < A->size1; ++r1) // for each row
    {
      gsl_vector_view row1 = gsl_matrix_row(A, r1);
      gsl_vector * row1_vector = &row1.vector;

      for (size_t r2 = 0; r2 < A->size1; ++r2) // for each row
      {
        gsl_vector_view row2 = gsl_matrix_row(A, r2);
        gsl_vector * row2_vector = &row2.vector;

        double corr1_2 = gsl_stats_correlation(row1_vector->data, row1_vector->stride,
                                               row2_vector->data, row2_vector->stride,
                                               num_obs);
        gsl_matrix_set(corr, r1, r2, corr1_2);
      }
    }

    return;
  }


  void calc_covariance_matrix_from_pairwise_complete_obs(
      gsl_matrix * A,
      matrix_mask & Amask,
      size_t min_num_complete_obs,
      gsl_matrix * covar)
  {
    size_t num_obs = A->size2;
    size_t num_elems = A->size1;

    if (covar->size1 != num_elems || covar->size2 != num_elems)
      throw std::runtime_error("Invalid lenght for covariance matrix");

    for (size_t r1 = 0; r1 < A->size1; ++r1) // for each row
    {
      for (size_t r2 = 0; r2 < A->size1; ++r2) // for each row
      {
        // count the number of pairwise-complete observations
        size_t n = 0;
        for (size_t c = 0; c < num_obs; ++c)
        {
          if (Amask(r1, c) == false && Amask(r2, c) == false)
            ++n;
        }

        if (n < min_num_complete_obs)
        {
          gsl_matrix_set(covar, r1, r2, 0.0);
        }
        else
        {
          gsl_vector * data1 = gsl_vector_calloc(n);
          gsl_vector * data2 = gsl_vector_calloc(n);

          // fill the vectors with the pairwise-complete observations
          for (size_t c = 0, i = 0; c < num_obs; ++c)
          {
            if (Amask(r1, c) == false && Amask(r2, c) == false)
            {
              gsl_vector_set(data1, i, gsl_matrix_get(A, r1, c));
              gsl_vector_set(data2, i, gsl_matrix_get(A, r2, c));
              ++i;
            }
          }

          double covar1_2 = gsl_stats_covariance(data1->data, data1->stride,
                                                 data2->data, data2->stride,
                                                 n);

          gsl_matrix_set(covar, r1, r2, covar1_2);

          gsl_vector_free(data1);
          gsl_vector_free(data2);
        }
      }
    }

    return;
  }


  void calc_adj_matrix_with_fixed_threshold(
      gsl_matrix * interaction,
      gsl_matrix * adj_matrix,
      decimal threshold,
      bool abs_values
      )
  {
    size_t n = interaction->size1;

    if (n != interaction->size2)
      throw std::runtime_error("invalid lenght for interaction matrix.");
    else if (n != adj_matrix->size1 || n != adj_matrix->size2)
      throw std::runtime_error("invalid lenght for adjacency matrix.");

    for (size_t r = 0; r < n; ++r)
    {
      for (size_t c = 0; c < n; ++c)
      {
        double value = gsl_matrix_get(interaction, r, c);
        if (abs_values)
          value = fabs(value);

        gsl_matrix_set(adj_matrix, r, c, (value >= threshold? value : 0.0));
      }
    }

    return;
  }


  void calc_adj_matrix_with_gaussian_threshold(
      gsl_matrix * interaction,
      gsl_matrix * adj_matrix,
      decimal threshold
      )
  {
    size_t n = interaction->size1;

    if (n != interaction->size2)
      throw std::runtime_error("invalid lenght for interaction matrix.");
    else if (n != adj_matrix->size1 || n != adj_matrix->size2)
      throw std::runtime_error("invalid lenght for adjacency matrix.");

    // Calculate mean and standard deviation (sd) of all values.
    if (interaction->tda != interaction->size2)
      throw std::runtime_error("interaction->tda != interaction->size2");

    double mean, sd;
    {
      size_t n_square = n * n;
      gsl_vector_view values = gsl_vector_view_array(interaction->data, n_square);
      gsl_vector * val_vec = &values.vector;

#ifdef CONAN_DEBUG
      std::cout << "values =";
      for (size_t i = 0; i < n_square; ++i)
        std::cout << ' ' << gsl_vector_get(val_vec, i);
      std::cout << std::endl << std::endl;
#endif

      mean = gsl_stats_mean(val_vec->data, val_vec->stride, n_square);
      sd = gsl_stats_sd(val_vec->data, val_vec->stride, n_square);
    }

#ifdef CONAN_DEBUG
    std::cout << "mean = " << mean << std::endl
              << "sd = " << sd << std::endl << std::endl;
#endif

    // Compute adjacency matrix
    for (size_t r = 0; r < n; ++r)
    {
      for (size_t c = 0; c < n; ++c)
      {
        double value = gsl_matrix_get(interaction, r, c);

        if (value < (mean - threshold * sd) || value > (mean + threshold * sd))
          gsl_matrix_set(adj_matrix, r, c, value);
        else
          gsl_matrix_set(adj_matrix, r, c, 0.0);
      }
    }

    return;
  }

}} // conan::inference

#endif // BASE_HPP
