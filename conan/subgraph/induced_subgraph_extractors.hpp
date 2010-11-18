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

#ifndef INDUCED_SUBGRAPH_EXTRACTORS_HPP
#define INDUCED_SUBGRAPH_EXTRACTORS_HPP
#include <conan/graphs.hpp>
#include <conan/subgraph/newman_communities.hpp>
#include <conan/subgraph/guimera_and_amaral.hpp>

#include <boost/graph/connected_components.hpp>

namespace conan {

  /// Wrapper for boost::connected_components function. Returns a vector of vectors of vertex_descriptors
  template <class Graph>
  void components(
      Graph & g,
      std::vector< std::vector<size_t> > & all_components
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex;
    typedef typename GraphTraits::vertex_iterator vertex_iter;
    typedef typename std::vector<vertex> Module;
    typedef typename std::vector<Module> ModuleVector;

    all_components.clear();

    std::vector<int> component(boost::num_vertices(g));
    size_t num_components = boost::connected_components(g, &component[0]);

    all_components.resize(num_components);
    for (size_t i = 0; i < num_components; ++i) // initialize all_components elements
      all_components[ i ] = Module();

    vertex_iter vi, viend;
    for (tie(vi, viend) = boost::vertices(g); vi != viend; ++vi)
      all_components[ component[ *vi ] ].push_back( *vi );

    return;
  }

} // conan

#endif // INDUCED_SUBGRAPH_EXTRACTORS_HPP
