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

#ifndef CONAN_MATRIX_HPP
#define CONAN_MATRIX_HPP
#include <fstream>

#define HAVE_INLINE
#include <gsl/gsl_matrix.h>

namespace conan { namespace linalg {

  gsl_matrix * alloc_and_read_matrix(
      std::string input_filename,
      char sep_char = ' '
      )
  {
    std::ifstream infile(input_filename.c_str());

    size_t num_cols = 0,
           num_rows = 0,
           i = 0,
           j = 0;

    std::string line;

    bool last_char_is_sep_char = false;

    // Compute number of columns and number of rows
    getline(infile, line);
    ++num_rows;

    for (std::string::iterator chr = line.begin(); chr != line.end(); ++chr)
      if (*chr == sep_char)
        ++num_cols;

    if (*(line.end() - 1) != sep_char)
      ++num_cols; // the last one
    else
      last_char_is_sep_char = true;

    while (not std::getline(infile, line).eof())
      ++num_rows;

    // Allocate space for the matrix
    gsl_matrix * M = gsl_matrix_calloc(num_rows, num_cols);

    // Come back to the begining of the file
    infile.clear();
    infile.seekg(0, std::ios_base::beg);

    // Read the file again and fill in the matrix
    while (not std::getline(infile, line).eof())
    {
      std::string aux;

      for (std::string::iterator chr = line.begin(); chr != line.end(); ++chr)
      {
        if (*chr != sep_char)
        {
          aux += *chr;
        }
        else
        {
          gsl_matrix_set(M, i, j++, from_string<decimal>(aux));
          aux.clear();
        }

        if (j == num_cols)
          break;
      }

      if (!last_char_is_sep_char)
      {
        gsl_matrix_set(M, i, j++, from_string<decimal>(aux));
      }

      ++i;
      j = 0;

      if (i == num_rows) // there's no more space in adj_matrix
        break;
    }

    infile.close();
    
    return M;
  }

}} // conan::linalg

#endif // CONAN_MATRIX_HPP
