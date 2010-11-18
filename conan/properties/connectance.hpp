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

#ifndef CONNECTANCE_HPP
#define CONNECTANCE_HPP

#include <conan/config.hpp>

namespace conan {

  template <class Graph>
  decimal graph_connectance(
      const Graph & g
      )
  {
    if (boost::is_directed(g))
      return boost::num_edges(g) / conan_pow( boost::num_vertices(g), 2.0 );
    else
      return 2.0 * boost::num_edges(g) / ( boost::num_vertices(g) * (boost::num_vertices(g) - 1.0) );
  }

} // conan

#endif // CONNECTANCE_HPP
