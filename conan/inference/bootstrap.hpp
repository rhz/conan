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

#ifndef BOOTSTRAP_HPP
#define BOOTSTRAP_HPP
#include <conan/inference/base.hpp>
#include <conan/inference/matrix_file.hpp>
#include <conan/utils/shuffle.hpp>

namespace conan { namespace inference {

  template <class CoexpressionMeasure>
  void bootstrap_and_calc_adj_matrix_helper(
      gsl_matrix * A,
      std::list<gsl_matrix *> & boot_matrix,
      detail::write_matrix_3d & write_boot_matrix,
      CoexpressionMeasure filter,
      boost::mt19937 & rng,
      void (* interaction_matrix_generator_fnc)(gsl_matrix *, gsl_matrix *, CoexpressionMeasure),
      bool use_tmp_file
      )
  {
    size_t num_elems = A->size1,
           num_obs = A->size2;

    gsl_matrix * A_shuffled = gsl_matrix_calloc(num_elems, num_obs);

    if (filter.shuffle_by_column)
      shuffle_by_column(A, A_shuffled, rng);
    else
      shuffle(A, A_shuffled, rng);

    gsl_matrix * shuffled_interaction_matrix = gsl_matrix_calloc(num_elems, num_elems);
    interaction_matrix_generator_fnc(A_shuffled, shuffled_interaction_matrix, filter);

    gsl_matrix_free(A_shuffled);

    if (use_tmp_file)
    {
      write_boot_matrix.write(shuffled_interaction_matrix);
      gsl_matrix_free(shuffled_interaction_matrix);
    }
    else
    {
      boot_matrix.push_back(shuffled_interaction_matrix);
    }

    return;
  }


  template <>
  void bootstrap_and_calc_adj_matrix_helper<covariance_pairwise_complete_obs_S>(
      gsl_matrix * A,
      std::list<gsl_matrix *> & boot_matrix,
      detail::write_matrix_3d & write_boot_matrix,
      covariance_pairwise_complete_obs_S filter,
      boost::mt19937 & rng,
      void (* interaction_matrix_generator_fnc)(gsl_matrix *, gsl_matrix *, covariance_pairwise_complete_obs_S),
      bool use_tmp_file
      )
  {
    size_t num_elems = A->size1,
           num_obs = A->size2;

    gsl_matrix * A_shuffled = gsl_matrix_calloc(num_elems, num_obs);
    matrix_mask A_shuffled_mask(num_elems, num_obs);
    covariance_pairwise_complete_obs_S A_shuffled_filter(A_shuffled_mask, filter.min_num_complete_obs);

    two_matrix_mask masks(filter.mask, A_shuffled_filter.mask);

    if (filter.shuffle_by_column)
      shuffle_by_column(A, A_shuffled, rng);
    else
      shuffle(A, A_shuffled, masks, rng);

    gsl_matrix * shuffled_interaction_matrix = gsl_matrix_calloc(num_elems, num_elems);
    interaction_matrix_generator_fnc(A_shuffled, shuffled_interaction_matrix, A_shuffled_filter);

    gsl_matrix_free(A_shuffled);

    if (use_tmp_file)
    {
      write_boot_matrix.write(shuffled_interaction_matrix);
      gsl_matrix_free(shuffled_interaction_matrix);
    }
    else
    {
      boot_matrix.push_back(shuffled_interaction_matrix);
    }

    return;
  }


  void write_data_matrix_with_matrix_mask_to_file(
      gsl_matrix * A,
      matrix_mask & m,
      std::string output_filename
      )
  {
    std::ofstream data_matrix_outfile(output_filename.c_str());
    for (size_t r = 0; r < A->size1; ++r)
    {
      if (m(r, 0))
        data_matrix_outfile << "NA";
      else
        data_matrix_outfile << gsl_matrix_get(A, r, 0);

      for (size_t c = 1; c < A->size2; ++c)
      {
        if (m(r, c))
          data_matrix_outfile << '\t' << "NA";
        else
          data_matrix_outfile << '\t' << gsl_matrix_get(A, r, c);
      }

      data_matrix_outfile << std::endl;
    }

    return;
  }


  template <class CoexpressionMeasure>
  gsl_matrix * bootstrap_and_calc_adj_matrix(
      gsl_matrix * A,
      size_t num_steps,
      decimal threshold,
      CoexpressionMeasure filter,
      void (* interaction_matrix_generator_fnc)(gsl_matrix *, gsl_matrix *, CoexpressionMeasure),
      bool use_tmp_file = false
      )
  {
    size_t num_elems = A->size1;
           //num_obs = A->size2;

    gsl_matrix * interaction = gsl_matrix_calloc(num_elems, num_elems);
    interaction_matrix_generator_fnc(A, interaction, filter);

    // Bootstrapping
    std::list<gsl_matrix *> boot_matrix;
    {
      detail::write_matrix_3d write_boot_matrix;

      if (use_tmp_file)
        write_boot_matrix.set_filename("tmp_file.out");

      boost::mt19937 rng(time(NULL) - num_steps - int(threshold * 100));

      for (size_t i = 0; i < num_steps; ++i)
      {
        std::cout << "step = " << i << std::endl; // DEBUG
        bootstrap_and_calc_adj_matrix_helper(A, boot_matrix, write_boot_matrix, filter, rng, interaction_matrix_generator_fnc, use_tmp_file);
      }
    }

    // Compute adjacency matrix
    detail::read_matrix_3d read_boot_matrix(num_elems, num_elems, num_steps);

    if (use_tmp_file)
      read_boot_matrix.set_filename("tmp_file.out");

    gsl_matrix * adj_matrix = gsl_matrix_alloc(num_elems, num_elems);

    for (size_t r = 0; r < num_elems; ++r)
    {
      for (size_t c = 0; c < num_elems; ++c)
      {
        double value = gsl_matrix_get(interaction, r, c);

        gsl_vector * boot_z_vec = gsl_vector_calloc(num_steps);

        if (use_tmp_file)
        {
          read_boot_matrix.read(r, c, boot_z_vec);
        }
        else
        {
          size_t i = 0;
          for (std::list<gsl_matrix *>::iterator li = boot_matrix.begin(); li != boot_matrix.end(); ++li, ++i)
            gsl_vector_set(boot_z_vec, i, gsl_matrix_get(*li, r, c));
        }

        double mean = gsl_stats_mean(boot_z_vec->data, boot_z_vec->stride, num_steps),
               sd = gsl_stats_sd(boot_z_vec->data, boot_z_vec->stride, num_steps);

        // DEBUG
        std::string filename = "boot_" + to_string(r) + "_" + to_string(c) + ".out";
        FILE * outfile = fopen(filename.c_str(), "wb");
        gsl_vector_fwrite(outfile, boot_z_vec);
        fclose(outfile);
        // END DEBUG
        
        if (value < (mean - threshold * sd) || value > (mean + threshold * sd))
          gsl_matrix_set(adj_matrix, r, c, value);
        else
          gsl_matrix_set(adj_matrix, r, c, 0.0);

        gsl_vector_free(boot_z_vec);
      }
    }

    // Free interaction matrices
    gsl_matrix_free(interaction);
    if (!use_tmp_file)
      for (std::list<gsl_matrix *>::iterator li = boot_matrix.begin(); li != boot_matrix.end(); ++li)
        gsl_matrix_free(*li);

    return adj_matrix;
  }
    

}} // conan::inference

#endif // BOOTSTRAP_HPP
