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

#ifndef GRAPHML_HPP
#define GRAPHML_HPP
#include <conan/config.hpp>

#define TIXML_USE_STL
#include <tinyxml/tinyxml.h>

#include <cstring>
#include <fstream> // for write to a file

#if defined(__GNUC__) && (__GNUC__ < 4 || (__GNUC__ == 4 && __GNUC_MINOR__ < 3))
namespace __gnu_cxx
{
  template<>
  struct hash< std::string >
  {
    size_t operator()( const std::string& x ) const
    { return hash< const char* >()( x.c_str() ); }
  };                                                                                          
}
#endif


namespace conan {

  /**
   * @brief Write the graph in GraphML format (XML dialect) to a file.
   * For reference about the GraphML format, see http://graphml.graphdrawing.org/
   *
   * @param g A graph.
   * @param output_filename The name of the output GraphML file.
   */
  template <class Graph>
  void write_graphml(
      Graph & g,
      std::string output_filename
      )
  {
    typedef typename boost::graph_traits<Graph> GraphTraits;
    typedef typename GraphTraits::vertex_iterator vertex_iter;
    typedef typename GraphTraits::edge_iterator edge_iter;

    TiXmlDocument doc;
    TiXmlDeclaration *decl = new TiXmlDeclaration("1.0", "", "");
    doc.LinkEndChild(decl);

    TiXmlElement *graphml_element = new TiXmlElement("graphml");
    graphml_element->SetAttribute("xmlns", "http://graphml.graphdrawing.org/xmlns");
    graphml_element->SetAttribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
    graphml_element->SetAttribute("xsi:schemaLocation", "http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd");
    doc.LinkEndChild(graphml_element);

    TiXmlElement *graph_element = new TiXmlElement("graph");
    graph_element->SetAttribute("id", "G");
    if (boost::is_directed(g))
      graph_element->SetAttribute("edgedefault", "directed");
    else
      graph_element->SetAttribute("edgedefault", "undirected");
    graphml_element->LinkEndChild(graph_element);

    vertex_iter vi, viend;
    for (tie(vi, viend) = boost::vertices(g); vi != viend; ++vi)
    {
      TiXmlElement *node_element = new TiXmlElement("node");
      node_element->SetAttribute("id", to_string(*vi));
      if (!g[*vi].name.empty())
        node_element->SetAttribute("label", g[*vi].name); // FIXME: not standard
      graph_element->LinkEndChild(node_element);
    }

    edge_iter ei, eiend;
    size_t edge_counter = 0;
    for (tie(ei, eiend) = boost::edges(g); ei != eiend; ++ei, ++edge_counter)
    {
      TiXmlElement *edge_element = new TiXmlElement("edge");
      edge_element->SetAttribute("id", to_string(edge_counter));
      edge_element->SetAttribute("source", to_string(boost::source(*ei, g)));
      edge_element->SetAttribute("target", to_string(boost::target(*ei, g)));
      if (g[*ei].weight != 1.0)
        edge_element->SetAttribute("weight", to_string(g[*ei].weight)); // FIXME: not standard
      graph_element->LinkEndChild(edge_element);
    }

    doc.SaveFile(output_filename.c_str());
    return;
  }


  /**
   * @brief Read a file containing a graph description in GraphML format (XML dialect).
   * For reference about the GraphML format, see http://graphml.graphdrawing.org/
   *
   * @param input_filename The name of the input GraphML file.
   * @return The graph.
   */
  template <class Graph>
  Graph read_graphml(
      std::string input_filename
      )
  {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex;
#if defined(__GNUC__) && (__GNUC__ < 4 || (__GNUC__ == 4 && __GNUC_MINOR__ < 3))
    typedef typename __gnu_cxx::hash_map<std::string, vertex> HashMap;
#else
    typedef typename std::unordered_map<std::string, vertex> HashMap;
#endif

    TiXmlDocument doc(input_filename.c_str());
    if (!doc.LoadFile())
      throw std::runtime_error("read_graphml: cannot open input file");

    TiXmlHandle hDoc(&doc);
    TiXmlElement *pElem;

    pElem = hDoc.FirstChildElement().Element();
    if (!pElem || strcmp(pElem->Value(), "graphml"))
      throw std::runtime_error("read_graphml: malformed input graphml file");

    pElem = hDoc.FirstChildElement().FirstChild().Element();
    while (pElem)
    {
      if ( ! strcmp(pElem->Value(), "graph") )
        break;
      pElem = pElem->NextSiblingElement();
    }

    if ( ! pElem )
      throw std::runtime_error("read_graphml: malformed input graphml file");

    TiXmlHandle hRoot(pElem);
    Graph out_g;
    pElem = hRoot.FirstChild().Element();
    HashMap id_hash;

    while (pElem)
    {
      std::string id, source_id, target_id, name;
      double weight = 1.0;

      TiXmlAttribute *pAttrib = pElem->FirstAttribute();
#ifdef CONAN_DEBUG
      std::cerr << "  " << pElem->Value() << ":" << std::endl;
#endif
      while (pAttrib)
      {
        if (!strcmp(pAttrib->Name(), "id"))
        {
          const char *pId = pAttrib->Value();
          if (pId)
            id = pId;
        }
        else if (!strcmp(pAttrib->Name(), "source"))
        {
          const char *pSource = pAttrib->Value();
          if (pSource)
            source_id = pSource;
        }
        else if (!strcmp(pAttrib->Name(), "target"))
        {
          const char *pTarget = pAttrib->Value();
          if (pTarget)
            target_id = pTarget;
        }
        else if (!strcmp(pAttrib->Name(), "weight"))
        {
          if (pAttrib->QueryDoubleValue(&weight) != TIXML_SUCCESS)
            std::cerr << "weight isn't a double presicion number" << std::endl;
        }
        else if (!strcmp(pAttrib->Name(), "label"))
        {
          const char *pName = pAttrib->Value();
          if (pName)
            name = pName;
        }
#ifdef CONAN_DEBUG
        std::cerr << "    " << pAttrib->Name() << " = " << pAttrib->Value() << std::endl;
#endif
        pAttrib = pAttrib->Next();
      }

      if (!strcmp(pElem->Value(), "node"))
      {
        vertex v = boost::add_vertex(out_g);
        id_hash[id] = v;
        if (!name.empty())
          out_g[v].name = name;
        else
          out_g[v].name = id;
#ifdef CONAN_DEBUG
        std::cerr << "Added node with id = " << id << " (remapped to " << v << ") and name = '" << out_g[v].name << "'" << std::endl;
#endif
      }
      else if (!strcmp(pElem->Value(), "edge"))
      {
        boost::add_edge(id_hash[source_id], id_hash[target_id], weight, out_g);
#ifdef CONAN_DEBUG
        std::cerr << "Added edge with id = " << id << ", source_id = " << source_id << " ("
                  << id_hash[source_id] << "), target_id = " << target_id << " ("
                  << id_hash[target_id] << ") and weight = " << weight << std::endl;
#endif
      }
      pElem = pElem->NextSiblingElement();
    }

    return out_g;
  }

} // conan

#endif //GRAPHML_HPP
