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

#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <iostream>
#include <cstdlib>
#include <string>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <cctype>

// STL
#include <vector>

#ifdef __GNUC__
  #if __GNUC__ >= 4 && __GNUC_MINOR__ >= 3
    #include <unordered_map>
  #else
    #include <ext/hash_map>
  #endif
#else
  #include <unordered_map>
#endif

// Boost
#include <boost/foreach.hpp>

// Boost.Random
#include <boost/random.hpp>
#ifndef DETERMINISTIC_RNG
  #ifdef __linux__
    #include <conan/utils/random_device.hpp>
  #else
    #define DETERMINISTIC_RNG
  #endif
#endif

// Macros
#if defined(CONAN_DEBUG) && !defined(CONAN_INFO)
  #define CONAN_INFO
#endif

#define CONAN_HAS_NAMED_DATA_MEMBER(name) \
  template <typename T, typename U> \
  class has_data_member_named_##name \
  { \
  public: \
    typedef char NotFound; \
    struct NonStaticDataMemberFound { char x[2]; }; \
    \
    template <typename Klass, T Klass::*> \
    struct TestNonStatic { }; \
    \
    template <typename Klass> \
    static NotFound TestDataMember(...) ; \
    \
    template <typename Klass> \
    static NonStaticDataMemberFound TestDataMember( TestNonStatic< Klass, &Klass::name >* ); \
    \
    enum { value = sizeof(TestDataMember<U>( NULL )) == 2 }; \
  }

#ifndef DECIMAL
  #define DECIMAL double
#endif

// Magic tricks to convert a macro into a string
#define _QUOTEME(x) #x
#define QUOTEME(x) _QUOTEME(x)


/**
 * @brief Main namespace
 * 
 * @exception ConversionFailure
 * This exception is raised when conversion to or from strings isn't possible.
 *
 * @exception FormatError
 * This exception is raised when a function requires a special format for one or more arguments
 * and the given data doesn't conform with the specified format.
 */
namespace conan {

  // Typedefs
  /// Data type 'decimal'. It can be double, float, and long double. By default is double.
  typedef DECIMAL decimal;


  // Exceptions
  /// @cond
  struct ConversionFailure { }; // : public std::runtime_error {};
  struct FormatError { }; // : public std::runtime_error {};
  /// @endcond


  // Structs
  //
  /// @brief Default set of vertex/vertex properties.
  struct DefaultVertexProperties
  {
    /// vertex name
    std::string name;
  };

  /// @brief Default set of edge properties.
  struct DefaultEdgeProperties
  {
    /// edge weight
    decimal weight;
    /// relation type
    std::string type;
  };

  struct VertexPropertiesViewableGraph :
    public DefaultVertexProperties
  {
    std::string shape;
    std::string fillcolor;
    std::string fontcolor;
    std::string fontname;
    size_t fontsize;
    decimal height;
    decimal width;
    std::string URL;
  };

  struct EdgePropertiesViewableGraph :
    public DefaultEdgeProperties
  {
    std::string arrowhead;
    std::string arrowsize;
    std::string arrowtail;
    std::string color;
    std::string label;
    std::string fontcolor;
    std::string fontname;
    size_t fontsize;
    std::string style;
    size_t weight;
  };


  // Small free functions
  /// @brief Tries to convert anything to string
  template <typename T>
  inline std::string to_string(T value)
  {
    std::stringstream ss;
    if (ss << value)
      return ss.str();
    throw ConversionFailure();
  }

  /// @brief Tries to convert a string to a given data type
  template <typename T>
  inline T from_string(const std::string& s)
  {
    T result;
    std::istringstream stream(s);
    if (stream >> result)
      return result;
    throw ConversionFailure();
  }

  /// @brief Base-2 logarithm for data type 'decimal'
  inline decimal conan_log2(decimal value)
  { // Ideally this should be done with preprocessor directives, but it could not. Anyways, the compiler can strip the "dead code"
    if (!strcmp(QUOTEME(DECIMAL), "float"))
      return log2f(value);
    else if (!strcmp(QUOTEME(DECIMAL), "double"))
      return log2(value);
    else if (!strcmp(QUOTEME(DECIMAL), "long double"))
      return log2l(value);
    else
      throw std::runtime_error(QUOTEME(DECIMAL) " not suppported by conan_log2");
  }

  /// @brief Natural logarithm for data type 'decimal'
  inline decimal conan_log(decimal value)
  { // Ideally this should be done with preprocessor directives, but it could not. Anyways, the compiler can strip the "dead code"
    if (!strcmp(QUOTEME(DECIMAL), "float"))
      return logf(value);
    else if (!strcmp(QUOTEME(DECIMAL), "double"))
      return log(value);
    else if (!strcmp(QUOTEME(DECIMAL), "long double"))
      return logl(value);
    else
      throw std::runtime_error(QUOTEME(DECIMAL) " not suppported by conan_log");
  }

  /// @brief Pow for data type 'decimal'
  inline decimal conan_pow(decimal x, decimal y)
  {
    if (x == 0.0)
      return 0.0;
    else if (x == 1.0 || y == 0.0)
      return 1.0;
    else
    { // Ideally this should be done with preprocessor directives, but it could not. Anyways, the compiler can strip the "dead code"
      if (!strcmp(QUOTEME(DECIMAL), "float"))
        return powf(x, y);
      else if (!strcmp(QUOTEME(DECIMAL), "double"))
        return pow(x, y);
      else if (!strcmp(QUOTEME(DECIMAL), "long double"))
        return powl(x, y);
      else
        throw std::runtime_error(QUOTEME(DECIMAL) " not suppported by conan_pow");
    }
  }
  
} // namespace conan

#endif //CONFIG_HPP
