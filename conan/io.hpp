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

#ifndef CONAN_IO_HPP
#define CONAN_IO_HPP
#include <conan/graphs.hpp>
#include <conan/utils.hpp>
#include <conan/io/graphml.hpp>
#include <fstream> // for write to a file

#include <boost/numeric/ublas/matrix.hpp>

namespace ublas = boost::numeric::ublas;


namespace conan {

  // ***** Declarations *****

  /**
   * @brief Write the graph in DOT language (www.graphviz.org) to a file.
   * DOT is a format used by GraphViz to visualize graphs.
   *
   * @param g A graph.
   * @param output_filename The name of the output DOT file.
   */
  template <class Graph>
  void write_dot(
      const Graph & g,
      std::string output_filename
      );

  /**
   * @brief Read a file containing a graph description in DOT language (www.graphviz.org).
   * DOT is a format used by GraphViz to visualize graphs.
   *
   * @param input_filename The name of the input DOT file.
   * @return The graph.
   */
  template <class Graph>
  Graph read_dot(
      std::string input_filename
      );

  /**
   * @brief Write the graph in GML format to a file.
   * For reference about the GML format, see http://www.infosun.fim.uni-passau.de/Graphlet/GML/
   *
   * @param g A graph.
   * @param output_filename The name of the output GML file.
   */
  template <class Graph>
  void write_gml(
      const Graph & g,
      std::string output_filename
      );

  /**
   * @brief Read a file containing a graph description in GML format.
   * For reference about the GML format, see http://www.infosun.fim.uni-passau.de/Graphlet/GML/
   *
   * @param input_filename The name of the input GML file.
   * @return The graph.
   */
  template <class Graph>
  Graph read_gml(
      std::string input_filename
      );

  /**
   * @brief Write the graph in the format used by the Pajek program to a file.
   * For reference, see the Pajek website: http://vlado.fmf.uni-lj.si/pub/networks/pajek/
   *
   * @param g A graph.
   * @param output_filename The name of the output Pajek file.
   */
  template <class Graph>
  void write_pajek(
      const Graph & g,
      std::string output_filename
      );

  /**
   * @brief Read a file containing a graph description in the format used by the Pajek program.
   * For reference, see the Pajek website: http://vlado.fmf.uni-lj.si/pub/networks/pajek/
   *
   * @param input_filename The name of the input Pajek file.
   * @return The graph.
   */
  template <class Graph>
  Graph read_pajek(
      std::string input_filename
      );

  /**
   * @brief Read a file containing a graph description in a CSV-like format.
   * In the input tabular data file each row represents an edge and the different fields/columns can be:
   * source vertex name (s), target vertex name (t), the weight of the edge (w), and the type of the edge (y).
   * You can define the order in which the fields appear in the input file with the second argument
   * to this function, fields_order.
   * For instance, "stw" means that the first field/column will be interpreted as the source vertex name,
   * the second field as the target vertex name and the last field as the weight of the edge.
   * The weight and type fields are optional.
   *
   * @param input_filename The name of the input CSV-like file.
   * @param fields_order The order of appearance of fields/columns (default = "stw").
   * @param sep_char The separator character used between fields (default = '\t').
   * @return The graph.
   */
  template <class Graph>
  Graph read_csv(
      std::string input_filename,
      std::string fields_order = "stw",
      char sep_char = '\t'
      );

  /**
   * @brief Write the graph in a CSV-like format to a file.
   * For reference about this format, read the description of function conan::read_csv.
   *
   * @param g A graph.
   * @param output_filename The name of the output CSV-like file.
   * @param fields_order The order of appearance of fields/columns (default = "stw").
   * @param sep_char The separator character used between fields (default = '\t').
   */
  template <class Graph>
  void write_csv(
      const Graph & g,
      std::string output_filename,
      std::string fields_order = "stw",
      char sep_char = '\t'
      );

  /**
   * @brief Write the adjacency matrix of a graph to a plain text file.
   *
   * @param g A graph.
   * @param output_filename The name of the output matrix file.
   * @param sep_char The character used to separate columns (default = ' ').
   *   Rows are always separated by newlines.
   */
  template <class Graph>
  void write_adj_matrix_file(
      const Graph & g,
      std::string output_filename,
      char sep_char = ' '
      );

  /**
   * @brief Read a file containing a graph description in matrix format.
   * Each row/column is a row/column of the adjacency matrix of the graph.
   *
   * @param input_filename The name of the input matrix file.
   * @param sep_char The character used to separate columns (default = ' ').
   *   Rows are always separated by newlines.
   * @return The graph.
   */
  template <class Graph>
  Graph read_adj_matrix_file(
      std::string input_filename,
      char sep_char = ' '
      );


  /// @cond

  CONAN_HAS_NAMED_DATA_MEMBER(name);

  namespace detail {
    struct has_name_tag { };
    struct does_not_have_name_tag { };

    template <bool>
    struct has_name : public has_name_tag { };

    template <>
    struct has_name<false> : public does_not_have_name_tag { };


    template <class Graph>
    inline std::string get_name_helper(
        const Graph & g,
        typename boost::graph_traits<Graph>::vertex_descriptor v,
        has_name_tag
        )
    { return (g[v].name != ""? g[v].name : to_string(v)); }
    
    template <class Graph>
    inline std::string get_name_helper(
        const Graph & g,
        typename boost::graph_traits<Graph>::vertex_descriptor v,
        does_not_have_name_tag
        )
    { return to_string(v); }


    template <class Graph>
    inline std::string get_vertex_name(
        const Graph & g,
        typename boost::graph_traits<Graph>::vertex_descriptor v
        )
    {
      return get_name_helper<Graph>(g, v,
          detail::has_name< has_data_member_named_name<std::string, typename Graph::vertex_bundled>::value >());
    }

  } // detail

  /// @endcond
  
  // *************** DOT files ***************
  template <class Graph>
  void write_dot(
      const Graph & g,
      std::string output_filename
      )
  {
    using namespace detail;
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex;
    typedef typename GraphTraits::edge_descriptor edge;

    std::ofstream dotfile(output_filename.c_str());

    bool is_directed = boost::is_directed(g);

    if (is_directed)
      dotfile << "digraph G {" << std::endl;
    else
      dotfile << "graph G {" << std::endl;

    // Write some convenient graph, node and edge properties
    dotfile << "  graph [overlap=scale];" << std::endl  //size=\"10,2\"
            << "  node [shape=\"circle\",style=filled,fillcolor=\"crimson\","
            << "  fontsize=15,fontcolor=\"white\",fixedsize=\"true\"];" << std::endl
            << "  edge [color=\"black\",fontcolor=\"brown\"];" << std::endl;
    
    // Write all nodes
    BOOST_FOREACH( vertex v, boost::vertices(g) )
    {
      dotfile << "  \"" << get_vertex_name(g, v) << "\";" << std::endl;
    }

    // Write all edges
    BOOST_FOREACH( edge e, boost::edges(g) )
    {
      vertex source_vertex = boost::source(e, g);
      vertex target_vertex = boost::target(e, g);

      dotfile << "  \"" << get_vertex_name(g, source_vertex);

      if (is_directed)
        dotfile << "\" -> \"";
      else
        dotfile << "\" -- \"";

      dotfile << get_vertex_name(g, target_vertex) << "\";" << std::endl;
    }

    dotfile << "}" << std::endl;
    return;
  }


  template <class Graph>
  Graph read_dot(
      std::string input_filename
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex;
    typedef typename GraphTraits::edge_descriptor edge;

    std::ifstream infile(input_filename.c_str()); // open input file

    if (!infile.is_open())
      throw std::runtime_error(input_filename + " could not be opened");

    char chr;
    std::string aux;
    Graph g;

    std::string v1_name,
                v2_name;

    bool is_node = true,
         in_quotes = false,
         in_subgraph = false;

    int skip = 0;
    
    // strip the beginning of the file
    {
      std::string buf;
      infile >> buf;

      if (buf != "graph" && buf != "digraph")
        throw std::runtime_error("mal-formatted input file.");

      infile >> buf;

      if (*(buf.end() - 1) != '{')
        infile >> buf;
    }

    infile.unsetf(std::ios_base::skipws);

    while (infile >> chr)
    {
      if (chr == '}' && not in_subgraph)
        break;
      else if (chr == '{')
        in_subgraph = true;
      else if (chr == '[')
        ++skip;
      else if (chr == ']')
        --skip;
      else if (skip)
        continue;
      else if (chr == '"')
        in_quotes = in_quotes? false : true;
      else if (chr == ' ' && not in_quotes)
        continue;
      else if (chr == '-' && not in_quotes)
      {
        v1_name = aux;
        is_node = false;
        // we don't take into account whether the link is defined as directed (->) or undirected (--)
        // in the input file, since this is defined in the type of the graph. So, you can read a graph
        // as directed though it could be written as undirected.
        infile >> chr;
        aux.clear();
      }
      else if ( (chr == ';' || chr == '\n') && not in_quotes )
      {
        if (aux.empty() || aux == "graph" || aux == "node" || aux == "edge")
        {
          aux.clear();
          continue;
        }

        if (is_node)
        {
          vertex v = boost::add_vertex(g);
          g[v].name = aux;
#ifdef CONAN_DEBUG
          std::cerr << "conan::read_dot: Added node " << v << " with name '" << aux << "'" << std::endl;
#endif
        }
        else
        { // is an edge
          v2_name = aux;

          vertex v1 = find_vertex_by_name(v1_name, g);

          if (v1 == GraphTraits::null_vertex())
          {
            v1 = boost::add_vertex(g);
            g[v1].name = v1_name;
          }

          vertex v2 = find_vertex_by_name(v2_name, g);

          if (v2 == GraphTraits::null_vertex())
          {
            v2 = boost::add_vertex(g);
            g[v2].name = v2_name;
          }

          boost::add_edge(v1, v2, 1.0, g);
#ifdef CONAN_DEBUG
          std::cerr << "conan::read_dot: Added edge between vertices " << v1_name << " (" << v1
                    << ") and " << v2_name << " (" << v2 << ")" << std::endl;
#endif

          v1_name.clear();
          v2_name.clear();
        }

        is_node = true;

        aux.clear();
      }
      else
      { // chr is part of the name of a vertex
        aux += chr;
      }
    } // while infile >> chr

    infile.close();

    return g;
  }


  // *************** CSV files ***************
  template <class Graph>
  Graph read_csv(
      std::string input_filename,
      std::string fields_order, //= "stw"
      char sep_char //= '\t'
      )
  {
    std::string fields[4];
    short int source_vertex_index = -1,  // to check if fields_order have source vertex's position info
        target_vertex_index = -1,        // idem
        type_index = -1,                 // to check if edges will have a type
        weight_index = -1,               // idem
        index = 0;

    // Parse fields_order
    for (std::string::iterator chr = fields_order.begin(); chr != fields_order.end(); ++chr)
    {
      switch (*chr)
      {
      case 's':  // source_vertex
        source_vertex_index = index++;
        break;
      case 't':  // target_vertex
        target_vertex_index = index++;
        break;
      case 'y':  // type
        type_index = index++;
        break;
      case 'w':  // weight
        weight_index = index++;
        break;
      default:
        throw FormatError();
      }
    }

    // source_vertex_index and target_vertex_index must be defined in fields_order
    if (source_vertex_index == -1 || target_vertex_index == -1)
      throw FormatError();

    std::ifstream csv_file(input_filename.c_str()); // open input file

    if (!csv_file.is_open())
      throw std::runtime_error("CSV file " + input_filename + " could not be opened");

    csv_file.exceptions ( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );
    
    Graph g;

#ifdef CONAN_DEBUG
    size_t line_num = 1;
#endif

    // Parse input file
    while (true)
    {
      typedef typename boost::graph_traits<Graph> GraphTraits;
      typedef typename GraphTraits::vertex_descriptor vertex;
      typedef typename GraphTraits::edge_descriptor edge;

      std::string line;

#ifdef CONAN_DEBUG
      std::cerr << "conan::read_csv: line_num = " << line_num++ << std::endl;
#endif

      try
      {
        std::getline(csv_file, line);
        line += sep_char;

        std::istringstream iss(line);
        iss.exceptions ( std::istringstream::eofbit | std::istringstream::failbit | std::istringstream::badbit );

        for (int i = 0; i < index; ++i)
          std::getline(iss, fields[i], sep_char);
      }
      catch (std::ios_base::failure & e)
      {
        if (csv_file.good()) // if the problem wasn't with csv_file, it was with iss
        {
          std::cerr << "conan::read_csv: Exception raised: " << e.what() << std::endl;
#ifdef CONAN_DEBUG
          for (int i = 0; i < index; ++i)
            std::cerr << "conan::read_csv: fields[" << i << "] = " << fields[i] << std::endl;
#endif
        }

        break;
      }

      if (fields[source_vertex_index] == fields[target_vertex_index] && boost::is_undirected(g))
      {
#ifdef CONAN_INFO
        std::cerr << "conan::read_csv: Edge between " << fields[source_vertex_index] << " and " << fields[target_vertex_index]
                  << " could not be added to the graph, because undirected graphs does not support it." << std::endl;
#endif
        continue;
      }

      // check if the vertices already exist
      vertex v1 = find_vertex_by_name(fields[source_vertex_index], g);

      if (v1 == GraphTraits::null_vertex())
      { // add the vertex if it doesn't exist yet
        v1 = boost::add_vertex(g);
        g[v1].name = fields[source_vertex_index];
      }

      // the same for the target vertex
      vertex v2 = find_vertex_by_name(fields[target_vertex_index], g);

      if (v2 == GraphTraits::null_vertex())
      {
        v2 = boost::add_vertex(g);
        g[v2].name = fields[target_vertex_index];
      }

      edge e;

      if (!boost::edge(v1, v2, g).second)
      {
        e = boost::add_edge(v1, v2, 1.0, g).first;
      }
      else
      {
#ifdef CONAN_INFO
        std::cerr << "conan::read_csv: Edge between " << g[v1].name << " and " << g[v2].name
                  << " has been already added to the graph. The former weight ("
                  << g[ boost::edge(v1, v2, g).first ].weight << ") will be kept." << std::endl;
#endif
        continue;
      }

      try
      {
        if (weight_index != -1)
          g[e].weight = from_string<decimal>( fields[weight_index] );
      }
      catch (ConversionFailure &)
      {
#ifdef CONAN_INFO
        std::cerr << "conan::read_csv: " << fields[weight_index] << " could not be traslated into a number (double)." << std::endl;
#endif
      }

      if (type_index != -1)
      {
        g[e].type = fields[type_index];
      }
    }

    csv_file.close();

    return g;
  }


  template <class Graph>
  void write_csv(
      const Graph & g,
      std::string output_filename,
      std::string fields_order, //= "stw"
      char sep_char //= '\t'
      )
  {
    using namespace detail;
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::edge_descriptor edge;

    std::ofstream outfile(output_filename.c_str());

    char last_field = fields_order[ fields_order.size() - 1 ];
    BOOST_FOREACH( edge e, boost::edges(g) )
    {
      BOOST_FOREACH( char c, fields_order )
      {
        if (c == 's')
          outfile << get_vertex_name( g, boost::source(e, g) );
        else if (c == 't')
          outfile << get_vertex_name( g, boost::target(e, g) );
        else if (c == 'w')
          outfile << g[ e ].weight;
        else if (c == 'y')
        {
          if ( not g[ e ].type.empty() )
            outfile << '\t' << g[ e ].type;
          else
            outfile << '\t' << "None" << std::endl;
        }

        if (c != last_field)
          outfile << sep_char;
      }
      outfile << std::endl;
    }

    return;
  }


  // *************** Adjacency matrix files ***************
  template <class Matrix>
  void write_matrix(
      Matrix & m,
      std::ostream & out,
      char sep_char = ' '
      )
  {
    for (size_t i = 0; i < m.size1(); ++i)
    {
      out << m(i, 0);

      for (size_t j = 1; j < m.size2(); ++j)
        out << sep_char << m(i, j);

      out << std::endl;
    }

    return;
  }


  void write_matrix(
      gsl_matrix * m,
      std::ostream & out,
      char sep_char = ' '
      )
  {
    for (size_t i = 0; i < m->size1; ++i)
    {
      out << gsl_matrix_get(m, i, 0);

      for (size_t j = 1; j < m->size2; ++j)
        out << sep_char << gsl_matrix_get(m, i, j);

      out << std::endl;
    }

    return;
  }


  template <typename T>
  void write_matrix(
      T ** m,
      size_t nrows,
      size_t ncols,
      std::ostream & out,
      char sep_char = ' '
      )
  {
    for (size_t i = 0; i < nrows; ++i)
    {
      out << m[i][0];

      for (size_t j = 1; j < ncols; ++j)
        out << sep_char << m[i][j];

      out << std::endl;
    }

    return;
  }


  template <class Graph>
  void write_adj_matrix_file(
      const Graph & g,
      std::string output_filename,
      char sep_char // = ' '
      )
  {
    typedef typename ublas::matrix<decimal> cmatrix;

    std::ofstream outfile(output_filename.c_str());

    cmatrix adj_matrix = get_adj_matrix<Graph, cmatrix>(g);

    write_matrix(adj_matrix, outfile, sep_char);

    outfile.close();

    return;
  }


  template <class Graph>
  Graph read_adj_matrix_file(
      std::string input_filename,
      char sep_char // = ' '
      )
  {
    typedef typename ublas::matrix<decimal> matrix;

    std::ifstream infile(input_filename.c_str()); // open input file

    int num_cols = 0,
        i = 0,
        j = 0;

    std::string line;

    bool last_char_is_sep_char = false;

    // Compute number of columns
    std::getline(infile, line);

    for (std::string::iterator chr = line.begin(); chr != line.end(); ++chr)
      if (*chr == sep_char)
        ++num_cols;

    if (*(line.end() - 1) != sep_char)
      ++num_cols; // the last one
    else
      last_char_is_sep_char = true;

    matrix adj_matrix = matrix(num_cols, num_cols);

    infile.seekg(0, std::ios_base::beg);

    while (!std::getline(infile, line).eof())
    {
      std::string aux;
      for (std::string::iterator chr = line.begin(); chr != line.end(); ++chr)
      {
        if (*chr != sep_char)
          aux += *chr;
        else
        {
          adj_matrix(i, j++) = from_string<decimal>(aux);
          aux.clear();
        }
      }

      if (!last_char_is_sep_char)
        adj_matrix(i, j++) = from_string<decimal>(aux);

      ++i;
      j = 0;

      if (i == num_cols) // there's no more space in adj_matrix
        break;
    }

    infile.close();

    return make_graph_from_adj_matrix<Graph, matrix>(adj_matrix);
  }


  // *************** GML ***************
  template <class Graph>
  void write_gml(
      const Graph & g,
      std::string output_filename)
  {
    using namespace detail;
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex;
    typedef typename GraphTraits::edge_descriptor edge;

    bool is_directed = boost::is_directed(g);

    std::ofstream outfile(output_filename.c_str());

    outfile << "graph [" << std::endl
            << "\tdirected " << (is_directed? 1 : 0) << std::endl;
    
    BOOST_FOREACH( vertex v, boost::vertices(g) )
    {
      outfile << "\tnode [" << std::endl
              << "\t\tid " << v << std::endl
              << "\t\tlabel \"" << get_vertex_name( g, v ) << "\"" << std::endl
              << "\t]" << std::endl;
    }
    
    BOOST_FOREACH( edge e, boost::edges(g) )
    {
      outfile << "\tedge [" << std::endl
              << "\t\tsource " << boost::source(e, g) << std::endl
              << "\t\ttarget " << boost::target(e, g) << std::endl
              << "\t\tweight " << g[ e ].weight << std::endl
              << "\t]" << std::endl;
    }

    outfile << "]" << std::endl;
    outfile.close();

    return;
  }


  template <class Graph>
  Graph read_gml(
      std::string input_filename
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex;
    typedef typename GraphTraits::edge_descriptor edge;

    // Mapping between the indices of the input file and the indices within the graph
#if defined(__GNUC__) && (__GNUC__ < 4 || (__GNUC__ == 4 && __GNUC_MINOR__ < 3))
    __gnu_cxx::hash_map<int, vertex> index_map;
#else
    std::unordered_map<int, vertex> index_map;
#endif

    std::string buf;
    std::ifstream infile(input_filename.c_str()); // open input file

    Graph out_g;

    bool in_graph = false,
         in_node = false,
         in_edge = false;

    int id = 0,
        source_id = -1,
        target_id = -1; // source and target id not set

    std::string node_name;
    decimal weight = 1.0;

    while (infile >> buf)
    {
      if (buf[0] == '#')
      {
        getline(infile, buf);
      }
      else if (buf == "graph")
      {
        in_graph = true;
      }
      else if (not in_graph)
      {
        continue;
      }
      else if (in_node and in_edge)
      {
        throw std::runtime_error("read_gml: input file ill-formated");
      }
      else if (buf == "]")
      {
        if (in_node)
        {
          vertex v = add_vertex(out_g);
          index_map[id] = v;

          if (node_name.empty())
            out_g[v].name = to_string(id);
          else
            out_g[v].name = node_name;

          in_node = false;
#ifdef CONAN_DEBUG
          std::cerr << "conan::read_gml: vertex " << v << " with name \"" << out_g[v].name << "\" added" << std::endl;
#endif
        }
        else if (in_edge)
        {
          if (source_id < 0 || target_id < 0)
            throw std::runtime_error("read_gml: input file ill-formated");

          boost::add_edge(index_map[source_id], index_map[target_id], weight, out_g);

#ifdef CONAN_DEBUG
          std::cerr << "conan::read_gml: edge between vertices " << source_id << " and " << target_id
                    << " with weight " << weight << " added" << std::endl;
#endif

          in_edge = false;
          source_id = -1;
          target_id = -1;
          weight = 1.0;
        }
        else if (in_graph)
        {
          in_graph = false;
        }
      }
      else if (buf == "node")
      {
        infile >> buf;

        if (buf != "[")
          throw std::runtime_error("read_gml: input file ill-formated");

        in_node = true;
      }
      else if (buf == "edge")
      {
        infile >> buf;

        if (buf != "[")
          throw std::runtime_error("read_gml: input file ill-formated");

        in_edge = true;
      }
      else if (in_node && buf == "id")
      {
        infile >> id;
      }
      else if (in_node && buf == "label")
      {
        infile >> node_name;

        if (node_name[0] == '"')
        {
          node_name.erase(0, 1); // removes the first character (")
          
          size_t last_pos = node_name.length() - 1;

          if (node_name[last_pos] == '"')
          {
            node_name.erase(last_pos);
          }
          else
          {
            std::string aux;
            getline(infile, aux, '"'); // get the string to the next '"'
            node_name += aux; // append it to the label's part we have already taken into node_name
          }
        }
      }
      else if (in_edge && buf == "source")
      {
        infile >> source_id;
      }
      else if (in_edge && buf == "target")
      {
        infile >> target_id;
      }
      else if (in_edge && buf == "label")
      {
        infile >> weight;

        if (infile.rdstate() & std::ifstream::failbit)
        { // if the string next to "label" could not be converted to double, then discard this attribute
          std::string aux;
          infile >> aux;
          if (aux[0] == '"' && aux[aux.length() - 1] != '"')
            getline(infile, aux, '"');
        }
      }
      else if (in_edge && buf == "weight")
      {
        infile >> weight;
      }
    } // while infile >> buf

    return out_g;
  }

  // *************** Pajek ***************
  template <class Graph>
  void write_pajek(
      const Graph & g,
      std::string output_filename
      )
  {
    using namespace detail;
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex;
    typedef typename GraphTraits::edge_descriptor edge;

    std::ofstream outfile(output_filename.c_str());

    outfile << "*Vertices " << boost::num_vertices(g) << std::endl;

    BOOST_FOREACH( vertex v, boost::vertices(g) )
    {
      outfile << v+1 << " \"" << get_vertex_name( g, v ) << "\"" << std::endl;
    }

    outfile << "*Edges" << std::endl;

    BOOST_FOREACH( edge e, boost::edges(g) )
    {
      outfile << boost::source(e, g)+1 << " "
              << boost::target(e, g)+1 << " "
              << g[ e ].weight << std::endl;
    }

    return;
  }


  inline std::string to_lower(
      std::string s
      )
  {
    for (std::string::iterator si = s.begin(); si != s.end(); ++si)
      *si = tolower(*si);
    return s;
  }

  template <class Graph>
  Graph read_pajek(
      std::string input_filename
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex;

    std::ifstream infile(input_filename.c_str()); // open input file

    Graph out_g;

    bool in_vertices = false,
         in_edges = false,
         in_arcslist = false;

    while (not infile.eof())
    {
      std::string buf, line;

      getline(infile, line);

      if (to_lower(line.substr(0, 9)) == "*vertices")
      {
        in_vertices = true;
        in_edges = false;
        in_arcslist = false;
      }
      else if (to_lower(line) == "*edges" || to_lower(line) == "*arcs")
      {
        in_vertices = false;
        in_edges = true;
        in_arcslist = false;
      }
      else if (to_lower(line) == "*arcslist")
      {
        in_vertices = false;
        in_edges = false;
        in_arcslist = true;
      }
      else if (in_vertices)
      {
        size_t first_space = line.find(' ');
        size_t id;
        if (first_space != std::string::npos)
        {
          id = from_string<size_t>(line.substr(0, first_space));
          vertex v = boost::add_vertex(out_g);
          if (v != id)
            throw std::runtime_error("read_pajek: v != id");

          char sep;
          if (line[first_space + 1] == '"')
            sep = '"';
          else
            sep = ' ';

          size_t next_sep_pos = line.find(sep, first_space + 2);
          out_g[v].name = line.substr(first_space + 2, next_sep_pos - first_space - 2);
        }
        else
        {
          id = from_string<size_t>(line);
          vertex v = boost::add_vertex(out_g);
          if (v != id)
            throw std::runtime_error("read_pajek: v != id");
        }
      }
      else if (in_edges)
      {
        size_t first_space = line.find(' ');
        size_t second_space = line.find(' ', first_space + 1);
        if (second_space != std::string::npos)
        {
          size_t source_id = from_string<size_t>(line.substr(0, first_space - 1)),
                 target_id = from_string<size_t>(line.substr(first_space + 1, second_space - 1));
          decimal weight = from_string<decimal>(line.substr(second_space + 1));
          boost::add_edge(source_id, target_id, weight, out_g);
        }
        else
        {
          size_t source_id = from_string<size_t>(line.substr(0, first_space - 1)),
                 target_id = from_string<size_t>(line.substr(first_space + 1));
          boost::add_edge(source_id, target_id, out_g);
        }
      }
      else if (in_arcslist)
      {
        throw std::runtime_error("read_pajek: implementation missing");
      }
    } // while (!infile.eof())

    return out_g;
  } // read_pajek

} // conan

#endif //CONAN_IO_HPP
