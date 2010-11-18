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

#ifndef FIBONACCI_HPP
#define FIBONACCI_HPP
#include <conan/config.hpp>

namespace conan {
  
  template <class Graph> Graph generate_fibonacci_graph (int);

  namespace detail {
#if 0
    template <int N>
    struct Fibonacci
    {
      enum { value = Fibonacci<N - 1>::value + Fibonacci<N - 2>::value };
    };

    template <>
    struct Fibonacci<1>
    {
      enum { value = 1 };
    };

    template <>
    struct Fibonacci<0>
    {
      enum { value = 0 };
    };
#endif

    int fibonacci(
        int n
        )
    {
      if (n == 0)
        return 0;
      else if (n == 1)
        return 1;
      else
      {
        int result = 2;
        for (int i = 2; i <= n; ++i)
          result *= i;
        return result;
      }
    }


    template <class Graph>
    typename boost::graph_traits<Graph>::vertex_descriptor expand_one_step_fibonacci_graph(
        Graph& g,
        typename boost::graph_traits<Graph>::vertex_descriptor v,
        int fibonacci_number
        )
    {
      typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iter;
      typedef typename boost::graph_traits<Graph>::edge_descriptor edge;
      typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex;

      int num_edges_to_create = fibonacci_number - boost::out_degree(v, g);
      vertex_iter target_vi, target_viend;
      for (tie(target_vi, target_viend) = boost::vertices(g); target_vi != target_viend; ++target_vi)
      {
        if (v == *target_vi)
          continue;

        if (!boost::edge(v, *target_vi, g).second)
        {
          edge e = boost::add_edge(v, *target_vi, g).first;
          g[e].weight = 1.0;
          --num_edges_to_create;
        }

        if (num_edges_to_create == 0)
          return v;
      }

      vertex new_v = v;
      for (int i = 0; i < num_edges_to_create; ++i)
      {
        new_v = boost::add_vertex(g);
        edge e = boost::add_edge(new_v, v, g).first;
        g[e].weight = 1.0;
      }
      
      return new_v;
    }
  } // detail


  template <class Graph>
  Graph generate_fibonacci_graph(
      int num_iterations
      )
  {
    Graph g;
    typename boost::graph_traits<Graph>::vertex_descriptor v = boost::add_vertex(g);
    for (int i = 1; i <= num_iterations; ++i)
      v = detail::expand_one_step_fibonacci_graph<Graph>(g, v, detail::fibonacci(i));
    return g;
  }

} // conan


#endif //FIBONACCI_HPP
