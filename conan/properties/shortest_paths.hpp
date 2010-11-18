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

#ifndef SHORTEST_PATH_HPP
#define SHORTEST_PATH_HPP
#include <conan/config.hpp>
#include <conan/subgraph.hpp>

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>

#include <limits>
#define MAX_DECIMAL_VALUE std::numeric_limits<decimal>::max()


namespace conan {

  /**
   * @brief Compute the average shortest path value for a single vertex.
   * @param v Vertex descritor.
   * @param g A graph.
   */
  template <class Graph>
  decimal vertex_avg_shortest_path(
      typename boost::graph_traits<Graph>::vertex_descriptor v,
      Graph & g
      );
  
  /**
   * @brief Return the average shortest path value for the graph.
   * @param g a const graph.
   */
  template <class Graph>
  decimal graph_avg_shortest_path(
      Graph & g
      );

  /**
   * @brief Return the average shortest path value for the graph normalized by the value for a
   *   random graph with the same number of vertices and edges (computed analitically).
   * @param g a const graph.
   */
  template <class Graph>
  decimal graph_avg_shortest_path_normalized(
      Graph & g
      );


  template <class Graph>
  inline decimal graph_avg_shortest_path_normalized(
      Graph & g
      )
  {
    return graph_avg_shortest_path(g) / random_graph_avg_shortest_path(boost::num_vertices(g), boost::num_edges(g));
  }


  template <class Graph>
  decimal vertex_avg_shortest_path(
      typename boost::graph_traits<Graph>::vertex_descriptor v,
      Graph & g
      )
  {
    int V = boost::num_vertices(g);
    decimal sum = 0;
    std::vector<decimal> d(V);

    boost::dijkstra_shortest_paths(g, v,
        boost::weight_map(get(&DefaultEdgeProperties::weight, g)).distance_map(&d[0]));

    BOOST_FOREACH( decimal distance, d )
      sum += distance;
    
    return sum / (decimal)(V - 1);
  }


  template <class Graph>
  inline bool all_pairs_shortest_path_dispatch(
      const Graph & g,
      decimal **D,
      conan::adj_listS
      )
  {
    typedef typename Graph::edge_bundled EdgePropertyType;
    return boost::johnson_all_pairs_shortest_paths(g, D,
        boost::weight_map(get(&EdgePropertyType::weight, g)));
  }


  template <class Graph>
  inline bool all_pairs_shortest_path_dispatch(
      const Graph & g,
      decimal **D,
      conan::adj_matrixS
      )
  {
    typedef typename Graph::edge_bundled EdgePropertyType;
    return boost::floyd_warshall_all_pairs_shortest_paths(g, D,
        boost::weight_map(get(&EdgePropertyType::weight, g)));
  }


  template <class Graph>
  decimal graph_avg_shortest_path(
      Graph & g
      )
  {
    typedef typename Graph::implementation_type implementation_type;
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_iterator vertex_iter;

    decimal graph_shortest_path_sum = 0;

    int V = boost::num_vertices(g);

#if 0
    decimal **D;
    D = (double**) malloc(sizeof(double*) * V);
    for (int i = 0; i < V; ++i)
      D[i] = (double*) malloc(sizeof(double) * V);
    
    // FIXME: conan::adj_matrixS should be replaced by implementation_type
    all_pairs_shortest_path_dispatch<Graph>(g, D, conan::adj_matrixS());

    bool is_undirected = boost::is_undirected(g);

    for (int i = 0; i < V; ++i)
    {
      for (int j = 0; j < V; ++j)
      {
        if (i == j && is_undirected)
          continue;

        graph_shortest_path_sum += D[i][j];
      }
    }

    for (int i = 0; i < V; ++i)
      free(D[i]);
    free(D);
#else
    vertex_iter vi, viend;
    for (tie(vi, viend) = boost::vertices(g); vi != viend; ++vi)
    {
      graph_shortest_path_sum += vertex_avg_shortest_path(*vi, g);
    }
#endif

    return graph_shortest_path_sum / decimal(V);
  }


  /**
   * @brief Class for storing and computing diverse graph properties that require the
   * initial calculation of the average shortest path for every vertex to every other vertex.
   * Examples of these properties are the eccentricity of a vertex and the diameter of the graph.
   */
  template <class Graph>
  struct shortest_paths
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename Graph::implementation_type implementation_type;
    typedef typename Graph::edge_bundled EdgePropertyType;
    typedef typename GraphTraits::vertex_descriptor vertex;
    typedef typename GraphTraits::vertex_iterator vertex_iterator;

    shortest_paths(Graph&);
    ~shortest_paths();

    /// Returns the average shortest path for a given vertex descriptor.
    decimal vertex_asp(vertex);
    /// Returns the graph average shortest path.
    decimal graph_asp();
    /**
     * @brief Returns the eccentricity of a given vertex_descriptor.
     * The eccentricity is defined as the maximum distance between a vertex and the others.
     */
    decimal vertex_eccentricity(vertex);
    /**
     * @brief Returns the diameter of the graph.
     * The diameter is defined as the maximum eccentricity in the graph.
     */
    decimal graph_diameter();
    /**
     * @brief Returns the radius of the graph.
     * The radius is defined as the minimum eccentricity in the graph.
     */
    decimal graph_radius();

    template <class List>
    List graph_center();

    double* eccentricity;
    double* avg_shortest_path;
    size_t V;
  };

  template <class Graph>
  shortest_paths<Graph>::shortest_paths(
      Graph& g
      )
  {
    double** distances;
    V = boost::num_vertices(g);
    distances = (double**) malloc(sizeof(double*) * V);
    eccentricity = (double*) malloc(sizeof(double) * V);
    avg_shortest_path = (double*) malloc(sizeof(double) * V);
    for (size_t i = 0; i < V; ++i)
      distances[i] = (double*) malloc(sizeof(double) * V);

    boost::floyd_warshall_all_pairs_shortest_paths(g, distances,
        boost::weight_map(get(&EdgePropertyType::weight, g)));

    for (size_t i = 0; i < V; ++i)
    {
      double max_distance = 0;
      double vertex_shortest_path_sum = 0;
      for (size_t j = 0; j < V; ++j)
      {
        if (i != j)
        {
          vertex_shortest_path_sum += distances[i][j];
          if (distances[i][j] > max_distance)
            max_distance = distances[i][j];
        }
      }
      eccentricity[i] = max_distance;
      avg_shortest_path[i] = vertex_shortest_path_sum / decimal(V - 1);
    }

    for (size_t i = 0; i < V; ++i)
      free(distances[i]);
    free(distances);

    return;
  }

  template <class Graph>
  shortest_paths<Graph>::~shortest_paths()
  {
    free(eccentricity);
    free(avg_shortest_path);
    return;
  }

  template <class Graph>
  decimal shortest_paths<Graph>::vertex_asp(
      vertex v
      )
  {
    return decimal(avg_shortest_path[v]);
  }

  template <class Graph>
  decimal shortest_paths<Graph>::graph_asp()
  {
    double graph_shortest_path_sum = 0;
    for (size_t i = 0; i < V; ++i)
      graph_shortest_path_sum += avg_shortest_path[i];
    return decimal(graph_shortest_path_sum / decimal(V));
  }

  template <class Graph>
  decimal shortest_paths<Graph>::vertex_eccentricity(
      vertex v
      )
  {
    return decimal(eccentricity[v]);
  }

  template <class Graph>
  decimal shortest_paths<Graph>::graph_diameter()
  {
    double max_eccentricity = 0;
    for (size_t i = 0; i < V; ++i)
      if (max_eccentricity < eccentricity[i])
        max_eccentricity = eccentricity[i];
    return decimal(max_eccentricity);
  }

  template <class Graph>
  decimal shortest_paths<Graph>::graph_radius()
  {
    double min_eccentricity = std::numeric_limits<double>::max();
    for (size_t i = 0; i < V; ++i)
      if (eccentricity[i] < min_eccentricity)
        min_eccentricity = eccentricity[i];
    return decimal(min_eccentricity);
  }

  template <class Graph>
  template <class List>
  List shortest_paths<Graph>::graph_center()
  {
    double min_eccentricity = graph_radius();
    List out_list;
    for (size_t i = 0; i < V; ++i)
      if (eccentricity[i] == min_eccentricity)
        out_list.push_back(i);
    return out_list;
  }

} // conan

#endif //SHORTEST_PATH_HPP
