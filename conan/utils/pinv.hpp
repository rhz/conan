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

#ifndef PINV_HPP
#define PINV_HPP
#include <conan/config.hpp>

#define HAVE_INLINE
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

namespace conan { namespace linalg {

  /**
   * @brief Compute the generalized inverse for the input matrix using Singular Value Decomposition (SVD).
   * Moore-Penrose method.
   * @warning The input matrix will not be the same after using this function.
   * @warning You should free the output matrix after using it.
   *
   * @param U Input matrix of m x n dimensions.
   * @return Pointer to the gsl_matrix containing the pseudo-inverse (generalized inverse).
   *   It has n x m dimensions. You should free the output matrix after using it.
   */
  gsl_matrix * pinv(
      gsl_matrix * U,
      gsl_matrix * ginv,
      double tolerance = 1E-8
      )
  {
    size_t m = U->size1,
           n = U->size2;

    if (ginv->size1 != n || ginv->size2 != m)
      throw std::runtime_error("Invalid size of ginv");

    gsl_vector * Sdiag = gsl_vector_calloc(n);
    gsl_vector * work = gsl_vector_calloc(n);
    gsl_matrix * V = gsl_matrix_calloc(n, n);

    // Singular Value Decomposition (SVD) of U
    gsl_linalg_SV_decomp(U, V, Sdiag, work);
    gsl_vector_free(work);

    // Editing Sdiag by choosing a suitable tolerance
    for (size_t i = 0; i < n; ++i)
      if (fabs(gsl_vector_get(Sdiag, i)) < tolerance)
        gsl_vector_set(Sdiag, i, 0.0);

    // V = V * Spinv(n, n) ( ifelse(row == col && Sdiag(col) != 0, 1.0 / Sdiag(col), 0.0) )
    for (size_t r = 0; r < n; ++r)
    {
      for (size_t c = 0; c < n; ++c)
      {
        double Spinv_c_c = 0.0;
        if (gsl_vector_get(Sdiag, c) != 0.0)
          Spinv_c_c = 1.0 / gsl_vector_get(Sdiag, c);

        gsl_matrix_set(V, r, c, Spinv_c_c * gsl_matrix_get(V, r, c));
      }
    }

    gsl_vector_free(Sdiag);

    // return V * U'
    gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, V, U, 0.0, ginv);

    gsl_matrix_free(V);
    return ginv;
  }

}} // conan::linalg

#endif // PINV_HPP
