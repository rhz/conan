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

#ifndef MATRIX_MASK_HPP
#define MATRIX_MASK_HPP
#include <fstream>

namespace conan { namespace linalg {

  struct matrix_mask
  {
    bool * masked;
    size_t size1;
    size_t size2;

    matrix_mask(
        size_t s1,
        size_t s2
        )
      : size1(s1), size2(s2)
    {
      masked = new bool[size1 * size2];
      for (size_t i = 0; i < size1 * size2; ++i)
        masked[i] = false;
    }

    matrix_mask(
        const matrix_mask & m
        )
      : size1(m.size1), size2(m.size2)
    {
      masked = new bool[size1 * size2];
      for (size_t i = 0; i < size1 * size2; ++i)
        masked[i] = m.masked[i];
    }

    ~matrix_mask()
    {
      delete[] masked;
    }

    bool& operator()(
        size_t row,
        size_t col
        )
    {
      if (row >= size1 || col >= size2)
        throw std::runtime_error("matrix_mask subscript out of bounds");

      return masked[size2 * row + col];
    }

    bool operator()(
        size_t row,
        size_t col
        )
      const
    {
      if (row >= size1 || col >= size2)
        throw std::runtime_error("const matrix_mask subscript out of bounds");

      return masked[size2 * row + col];
    }

    bool operator==(
        const matrix_mask & other
        ) const
    {
      if (size1 != other.size1 || size2 != other.size2)
        return false;

      for (size_t i = 0; i < size1 * size2; ++i)
        if (masked[i] != other.masked[i])
          return false;

      return true;
    }

    void write(
        std::string output_filename
        )
    {
      std::ofstream outfile(output_filename.c_str());

      outfile << size1 << ' ' << size2 << std::endl;

      for (size_t i = 0; i < size1 * size2; ++i)
        outfile << masked[i] << std::endl;

      outfile.close();
      return;
    }

    void read(
        std::string input_filename,
        bool tolerant_mode = true
        )
    {
      std::ifstream infile(input_filename.c_str());
      size_t s1, s2;

      infile >> s1 >> s2;
      if (!tolerant_mode && (size1 != s1 || size2 != s2))
        throw std::runtime_error("input matrix mask file doesn't have the same size as this matrix_mask object");

      size_t i = 0;
      while(infile >> masked[i++]) { }

      infile.close();
      return;
    }

  };

}} // conan::linalg

#endif // MATRIX_MASK_HPP
