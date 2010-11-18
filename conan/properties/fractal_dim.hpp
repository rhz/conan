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

#ifndef FRACTAL_DIM_HPP
#define FRACTAL_DIM_HPP

#include <conan/graphs.hpp>
#include <conan/utils/boxcounting.hpp>

#include <gsl/gsl_fit.h>
#include <gsl/gsl_statistics_double.h>

namespace conan {

  /**
   * @brief Returns the graph fractal dimension computed with the boxcounting method described in the publication:
   * Self-similarity of complex networks, Chaoming Song, Shlomo Havlin & Hernan A. Makse, Nature, 2005.
   * It actually returns an average of performing multiple times the boxcounting algorithm, as it is an stochastic procedure.
   * The start, end and step arguments can be fractionary.
   *
   * @param g A graph.
   * @param start Minimum box size used in the boxcounting procedure, 2 by default.
   * @param end Maximum box size used in the boxcounting procedure, 8 by default.
   * @param step Difference between the box sizes used, 1 by default.
   * @param num_iterations Number of boxcounting procedures to average.
   * @return Value of the graph fractal dimension.
   */
  template <class Graph>
  decimal graph_fractal_dimension(
      Graph&,
      decimal start = 2.0,
      decimal end = 8.0,
      decimal step = 1.0,
      size_t num_iterations = 10);


  template <class Graph>
  decimal graph_fractal_dimension_one_iteration(
      Graph & g,
      decimal start, //= 2.0
      decimal end, //= 8.0
      decimal step //= 1.0
      )
  {
    int num_iterations = int((end - start) / step) + 1;
    double num_vertices[num_iterations];
    double boxsize_array[num_iterations];
    int count = 0;
    for (decimal boxsize = start; boxsize <= end; boxsize += step)
    {
      Graph renormalized_graph = detail::boxcounting(g, boxsize);
      boxsize_array[count] = conan_log(boxsize);
      num_vertices[count++] = conan_log(boost::num_vertices(renormalized_graph));
    }

    double c0, c1, cov00, cov01, cov11, sumsq;
    gsl_fit_linear(boxsize_array, 1, // x = boxsize
                   num_vertices, 1,  // y = num_vertices
                   num_iterations, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);

    /*
    std::cout << "y = " << c1 << " * x + " << c0 << std::endl;
    std::cout << "cov00 = " << cov00 << std::endl;
    std::cout << "cov01 = " << cov01 << std::endl;
    std::cout << "cov11 = " << cov11 << std::endl;
    std::cout << "sumsq = " << sumsq << std::endl;
    */

    /*
    double r =
      gsl_stats_correlation(boxsize_array, 1, num_vertices, 1, num_iterations);
    std::cout << "r = " << r << std::endl;
    std::cout << "r^2 = " << r*r << std::endl;
    */

    return (decimal) -c1;
  }


  template <class Graph>
  inline decimal graph_fractal_dimension(
      Graph & g,
      decimal start, //= 2.0
      decimal end, //= 8.0
      decimal step, //= 1.0
      size_t num_iterations //=10
      )
  {
    decimal sum = 0, n = num_iterations;
    while (--num_iterations)
      sum += graph_fractal_dimension_one_iteration(g, start, end, step);

    return sum / n;
  }

} // conan

#endif //FRACTAL_DIM_HPP
