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

#ifndef SHUFFLE_HPP
#define SHUFFLE_HPP
#include <conan/utils/matrix_mask.hpp>

namespace conan { namespace linalg {

  struct no_matrices { };

  struct two_matrix_mask
  {
    matrix_mask & orig;
    matrix_mask & shuffled;

    two_matrix_mask(
        matrix_mask & o,
        matrix_mask & s
        )
      : orig(o), shuffled(s)
    { }
  };


  template <class MatrixMasks>
  void matrix_mask_copy_value(
      MatrixMasks masks,
      size_t source_row,
      size_t source_col,
      size_t target_row,
      size_t target_col
      )
  { }


  template <>
  void matrix_mask_copy_value<two_matrix_mask>(
      two_matrix_mask masks,
      size_t source_row,
      size_t source_col,
      size_t target_row,
      size_t target_col
      )
  {
    masks.shuffled(target_row, target_col) = masks.orig(source_row, source_col);
  }


  template <class RNG, class MatrixMasks>
  void shuffle(
      gsl_matrix * A,
      gsl_matrix * A_shuffled,
      MatrixMasks masks,
      RNG & gen
      )
  {
    size_t num_elems = A->size1,
           num_obs = A->size2;

    size_t nrows = A_shuffled->size1,
           ncols = A_shuffled->size2;

    if (num_elems * num_obs != nrows * ncols)
      throw std::runtime_error("the total number of elements of input matrices must be the same");

    size_t remaining_entries = nrows * ncols;

    char ready[remaining_entries];
    for (size_t i = 0; i < remaining_entries; ++i)
      ready[i] = 0;

    --remaining_entries;
    for (size_t r = 0; r < num_elems; ++r)
    {
      for (size_t c = 0; c < num_obs; ++c, --remaining_entries)
      {
        boost::uniform_int<size_t> dist(0, remaining_entries);
        boost::variate_generator< RNG&, boost::uniform_int<size_t> >
          random_entry(gen, dist);

        size_t target_pos = random_entry();

        while (ready[target_pos] == 1)
          ++target_pos;

        size_t target_row = target_pos / ncols,
               target_col = target_pos % ncols;

        gsl_matrix_set(A_shuffled, target_row, target_col, gsl_matrix_get(A, r, c));

        matrix_mask_copy_value(masks, r, c, target_row, target_col);

        ready[target_pos] = 1;
      }
    }

    return;
  }


  template <class RNG>
  void shuffle(
      gsl_matrix * A,
      gsl_matrix * A_shuffled,
      RNG& gen
      )
  {
    no_matrices dumb;
    shuffle(A, A_shuffled, dumb, gen);
    return;
  }


  template <class RNG, class MatrixMasks>
  void shuffle_by_column(
      gsl_matrix * A,
      gsl_matrix * A_shuffled,
      MatrixMasks masks,
      RNG & gen
      )
  {
    size_t num_elems = A->size1,
           num_obs = A->size2;

    size_t nrows = A_shuffled->size1,
           ncols = A_shuffled->size2;

    if (num_elems != nrows || num_obs != ncols)
      throw std::runtime_error("invalid size for matrices A and A_shuffled");

    for (size_t c = 0; c < num_obs; ++c)
    {
      size_t remaining_entries = nrows;

      char ready[remaining_entries];
      for (size_t i = 0; i < remaining_entries; ++i)
        ready[i] = 0;

      --remaining_entries;

      gsl_vector_view vec = gsl_matrix_column(A, c);

      for (size_t i = 0; i < num_elems; ++i, --remaining_entries)
      {
        boost::uniform_int<size_t> dist(0, remaining_entries);
        boost::variate_generator< RNG&, boost::uniform_int<size_t> >
          random_entry(gen, dist);

        size_t target_row = random_entry();

        while (ready[target_row] == 1)
          ++target_row;

        gsl_matrix_set(A_shuffled, target_row, c, gsl_vector_get(&vec.vector, i));

        matrix_mask_copy_value(masks, i, c, target_row, c);

        ready[target_row] = 1;
      }
    }

    return;
  }



  template <class RNG>
  void shuffle_by_column(
      gsl_matrix * A,
      gsl_matrix * A_shuffled,
      RNG& gen
      )
  {
    no_matrices dumb;
    shuffle_by_column(A, A_shuffled, dumb, gen);
    return;
  }

}} // conan::linalg

#endif // SHUFFLE_HPP
