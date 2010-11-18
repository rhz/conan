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

#ifndef CHOICE_HPP
#define CHOICE_HPP
#include <ctime>

namespace conan { namespace detail {

    /**
     * @brief Compute a random choice based on each element weight.
     *
     * Used for Preferential Attachment model.
     *
     * @param vec a vector containing the weight of each element
     * @return the index of the chosen element.
     * 
     * TODO: en que rango esta el resultado? [0,sum) ? [0,sum] ? [1,sum] ?
     */
    template <typename T>
    inline size_t weighted_choice(
        const T &
        );

    template <typename T>
    inline size_t weighted_choice(
        const T &,
        typename T::value_type sum
        );

    template <class RNG, typename T>
    size_t weighted_choice(
        const T &,
        typename T::value_type sum,
        RNG &
        );


    template <typename T>
    inline
    size_t weighted_choice(
        const T & vec
        )
    {
      typename T::value_type sum = 0;
      BOOST_FOREACH( typename T::value_type val, vec )
        sum += val;
      return weighted_choice(vec, sum);
    }


    template <typename T>
    inline
    size_t weighted_choice(
        const T & vec,
        typename T::value_type sum
        )
    {
#ifdef DETERMINISTIC_RNG
      boost::mt19937 rng(time(NULL) - sum);
#else
      boost::random_device rng;
#endif
      return weighted_choice(vec, sum, rng);
    }


    template <class RNG, typename T>
    size_t weighted_choice(
        const T & vec,
        typename T::value_type sum,
        RNG & gen
        )
    {
      boost::uniform_real<decimal> dist(0, decimal(sum));
      boost::variate_generator< RNG &, boost::uniform_real<decimal> >
        dice(gen, dist);

      decimal random_number = dice(),
              accumulated_value = 0;
#ifdef CONAN_DEBUG
      std::cerr << "conan::detail::weighted_choice: sum = " << sum << ", random_number = " << random_number << std::endl;
#endif

      size_t vector_index = 0;
      BOOST_FOREACH( typename T::value_type val, vec )
      {
        accumulated_value += decimal(val);
        if (random_number <= accumulated_value)
          return vector_index;
        ++vector_index;
      }
      return vector_index;
    }


    template <class RNG, typename T>
    size_t weighted_choice(
        const std::vector<T> & vec,
        decimal sum,
        RNG & gen,
        decimal (retrieve_weight) (const T &)
        )
    {
      boost::uniform_real<decimal> dist(0, sum); // [0,sum] uniform distribution. CHECK THIS !!
      boost::variate_generator< RNG &, boost::uniform_real<decimal> > dice(gen, dist);

      decimal random_number = dice(), current_pos = 0;

#ifdef CONAN_DEBUG
      std::cerr << "conan::detail::weighted_choice: sum = " << sum << ", random_number = " << random_number << std::endl;
#endif

      size_t vector_index = 0;

      BOOST_FOREACH( T val, vec )
      {
        current_pos += retrieve_weight(val);

        if (random_number <= current_pos)
          return vector_index;

        ++vector_index;
      }

      // assert vectorIndex < vec.size()
      return vector_index;
    }

}} // conan::detail

#endif // CHOICE_HPP
