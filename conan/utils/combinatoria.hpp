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

/*********************************************
 * The author of this file is Hector Urbina. *
 * Thanks for the code, Hector.              *
 *********************************************/

#ifndef COMBINATORIA_HPP
#define COMBINATORIA_HPP
#include <conan/config.hpp>

namespace conan { namespace detail {

  /// @brief Returns the greater comun denominator between a and b
  int gcd(
      int a,
      int b
      )
  {
    if (a < b)
    {
      a = a + b;
      b = a - b;
      a = a - b;
    }

    int a_tmp;

    while (b)
    {
      a_tmp = a;
      a = b;
      b = a_tmp % b;
    }

    return a;
  }


  void rationalize(
      int num_factors[],
      int den_factors[],
      size_t n
      )
  {
    for (size_t i = 0; i < n; ++i) // iterating over denominator factors
    {
      for (size_t j = 0; j < n && den_factors[i] > 1; ++j) // iterating over numerator factors
      {
        int divisor = gcd(num_factors[j], den_factors[i]);
        num_factors[j] = num_factors[j] / divisor;
        den_factors[i] = den_factors[i] / divisor;
      }
    }

    return;
  }


  decimal binomial_coefficient(
      int a,
      int b
      )
  {
    // doing b > a - b
    if (a - b > b)
      b = a - b;

    int c = a - b;
    int numerator_factors[c];
    int denominator_factors[c];
    
    for (int i = 0; i < c; ++i)
    {
      numerator_factors[i] = b + 1 + i;
      denominator_factors[i] = i + 1;
    }
    
    rationalize(numerator_factors, denominator_factors, c);

    decimal result = 1.0;
    for (int i = 0; i < c; ++i)
    {
      result *= numerator_factors[i];
    }

    return result;
  }

}} // conan::detail

#endif // COMBINATORIA_HPP
