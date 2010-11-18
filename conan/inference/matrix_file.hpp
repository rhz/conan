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

#ifndef MATRIX_3D_DB_HPP
#define MATRIX_3D_DB_HPP

namespace conan { namespace inference { namespace detail {

  struct read_matrix_3d
  {
    FILE * infile;
    size_t size1,
           size2,
           size3;

    read_matrix_3d(
        size_t s1,
        size_t s2,
        size_t s3
        )
    {
      infile = NULL;
      size1 = s1;
      size2 = s2;
      size3 = s3;
      return;
    }

    ~read_matrix_3d()
    {
      if (infile != NULL)
        fclose(infile);
    }

    void set_filename(
        std::string input_filename
        )
    {
      infile = fopen(input_filename.c_str(), "rb");

      if (infile == NULL)
        throw std::runtime_error(to_string("error opening file ") + input_filename);

      return;
    }

    void read(
        size_t r,
        size_t c,
        gsl_vector * vec
        )
    {
      if (infile == NULL)
        throw std::runtime_error("Filename was not set yet");

      fseek(infile, ((r * size2) + c) * sizeof(double), SEEK_SET); // set the initial position

      for (size_t i = 0; i < size3; ++i)
      {
        double tmp = 0;
        fread(&tmp, sizeof(double), 1, infile);

        gsl_vector_set(vec, i, tmp);

        fseek(infile, ((size1 * size2) - 1) * sizeof(double), SEEK_CUR); // jump to the next pos

        if (feof(infile))
          throw std::runtime_error("EOF found");
      }

      return;
    }

  }; // struct read_matrix_3d

  struct write_matrix_3d
  {
    FILE * outfile;
    size_t size1,
           size2;

    write_matrix_3d()
    {
      outfile = NULL;
      size1 = 0;
      size2 = 0;
    }

    ~write_matrix_3d()
    {
      if (outfile != NULL)
        fclose(outfile);
    }

    void set_filename(
        std::string output_filename
        )
    {
      outfile = fopen(output_filename.c_str(), "wb");

      if (outfile == NULL)
        throw std::runtime_error("Error opening file " + output_filename);

      return;
    }

    void write(
        gsl_matrix * m
        )
    {
      if (size1 == 0 && size2 == 0)
      {
        size1 = m->size1;
        size2 = m->size2;
      }
      else if (size1 != m->size1 || size2 != m->size2)
        throw std::runtime_error("Invalid matrix size");

      if (outfile == NULL)
        throw std::runtime_error("Filename was not set yet");

      int status = gsl_matrix_fwrite(outfile, m);

      if (status)
        throw std::runtime_error("There was a problem writing the matrix to a file");

      return;
    }

  }; // struct write_matrix_3d

}}} // conan::inference::detail

#endif // MATRIX_3D_DB_HPP
