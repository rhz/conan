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

#ifndef DISTRIBUTIONS_HPP
#define DISTRIBUTIONS_HPP

#include <conan/config.hpp>
#include <conan/utils/chisq.hpp>

#include <cmath>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>

#include <limits>


namespace conan {

  // Distribution Types
  enum DistributionType {
    linear = 0,       // y = m * x + n
    power = 1,        // y = A * x^(-gamma) + b
    exponential = 2,  // y = A * exp(-lambda * x) + b
    logaritmic = 3,   // y = A * log(-lambda * x) + b
    gaussian = 4,     // y = y0 + A exp(-(x - x0)^2 / (2 var))
    gaussian_2peaks,  // y = y0 + A1 exp(-(x - x1)^2 / (2 var1)) + A2 exp(-(x - x2)^2 / (2 var2))
    double_exp,       // y = y0 + A1 exp(-b1 * x) + A2 exp(-b2 * x)
    sinusoidal,       // y = y0 + A sin(f * x + phi)
    lorentzian_peak,  // y = y0 + A / ((x - x0)^2 + B)
    hill,             // y = y0 + (m - y0) / (1 + (xhalf / x)^r)
    sigmoid,          // y = y0 + m / (1 + exp((x0 - x) / r))
    lognormal         // y = y0 + A exp[ -(log(x / x0) / width)^2 ]
  };

  // Degree Types
  enum DegreeType {
    out_degree_type = 0,
    in_degree_type = 1,
    in_out_degree_type = 2
  };

  struct invalid_degree_tag { };
  struct out_degree_tag { };
  struct in_degree_tag { };
  struct in_out_degree_tag { };

  template <int DegreeType>
  struct degree_selector
  {
    typedef invalid_degree_tag degree_type_tag;
  };

  template <>
  struct degree_selector<out_degree_type>
  {
    typedef out_degree_tag degree_type_tag;
  };

  template <>
  struct degree_selector<in_degree_type>
  {
    typedef in_degree_tag degree_type_tag;
  };

  template <>
  struct degree_selector<in_out_degree_type>
  {
    typedef in_out_degree_tag degree_type_tag;
  };


  // Degree distribution class
  /**
   * @brief Class used to store and compute the degree distribution of a graph and
   * its fit parameters.
   */
  template <class Graph,
            int DegreeType = out_degree_type>
  struct degree_distribution
  {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex;
    typedef typename boost::graph_traits<Graph>::degree_size_type degree_size_type;
    typedef typename degree_selector<DegreeType>::degree_type_tag degree_type_tag;

    typedef typename std::vector<vertex> VertexList;
    typedef typename std::vector<degree_size_type> DegreeList;

    typedef typename std::pair<degree_size_type, VertexList> k_vertices_pair;
    typedef typename std::vector<k_vertices_pair> DistVector;

    // Data members
    size_t V;
    bool cum;
    int best_fit;
    DistVector dist_vector;
    double linear_fit_parameters[3];  // r, m(slope), n(constant)... y = m * x + n
    double power_fit_parameters[3];   // r, gamma, A... y = A * x^gamma + b
    double exp_fit_parameters[3];     // r, lambda, A... y = A * exp(lambda * x) + b
    double log_fit_parameters[3];     // r, lambda, A... y = A * log(lambda * x) + b

    // Member functions
    degree_distribution(Graph &, bool cumulative = false);

    /// Returns the probability of taking a random vertex and it to have the given degree.
    decimal P(degree_size_type);
    /// Returns a STL vector with the vertices of the graph which have the given degree.
    VertexList vertices(degree_size_type);
    /// Returns a STL vector with all non-zero degrees present in the degree distribution. They are not sorted.
    DegreeList non_zero_degrees();
    /// Returns the minimum non-zero degree present in the degree distribution.
    degree_size_type min_non_zero_degree();
    /// Returns the maximum non-zero degree present in the degree distribution.
    degree_size_type max_non_zero_degree();
    /// Returns the exponent for the power-law regression y = A * x^(-gamma) + b. Note the minus sign for gamma.
    decimal gamma();
    /// Returns the exponent for the exponential regression y = A * exp(-lambda * x) + b. Note the minus sign for lambda.
    decimal lambda();
    /// Returns the linear slope for the regression y = m * x + n.
    decimal slope();

    /**
     * @brief Returns the Pearson's correlation coefficient for the given regression.
     * @param fit_type The fit type, by default a power-law. See the conan::DistributionType enum for reference.
     */
    decimal r(int fit_type = power);
    /**
     * @ brief Returns the square of the Pearson's correlation coefficient for the given regression.
     * @param fit_type The fit type, by default a power-law. See the conan::DistributionType enum for reference.
     */
    decimal rr(int fit_type = power);
    /**
     * @ brief Returns a three elements C array containing the regression parameters for a given fit type.
     * @param fit_type The fit type, by default a power-law. See the conan::DistributionType enum for reference.
     */
    double* _fit_parameters(int fit_type);
    decimal entropy();

    degree_size_type degree_fcn(vertex v, Graph& g, out_degree_tag)
    { return boost::out_degree(v, g); }

    degree_size_type degree_fcn(vertex v, Graph& g, in_degree_tag)
    { return boost::in_degree(v, g); }

    degree_size_type degree_fcn(vertex v, Graph& g, in_out_degree_tag)
    { return boost::degree(v, g); }

  };

  namespace detail {

    bool k_vertices_pair_less_than(
        const std::pair< size_t, std::vector<size_t> > & left_k_vertices_pair,
        const std::pair< size_t, std::vector<size_t> > & right_k_vertices_pair
        )
    {
      return (left_k_vertices_pair.first > right_k_vertices_pair.first); // inverse order
    }

  } // detail

  template <class Graph, int DegreeType>
  degree_distribution<Graph, DegreeType>::degree_distribution(
      Graph & g,
      bool cumulative // = false
      )
  {
    typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iter;
    typedef typename boost::graph_traits<Graph>::degree_size_type degree_size_type;

    // set V and cumulative
    V = boost::num_vertices(g);
    cum = cumulative;

    // add vertices' degree to dist_vector
    vertex_iter vi, viend;
    for (tie(vi, viend) = boost::vertices(g); vi != viend; ++vi)
    {
      bool vertex_added = false;
      degree_size_type degree = degree_fcn(*vi, g, degree_type_tag());

      // check if already exists a k_vertices_pair with degree = degree
      for (typename DistVector::iterator veci = dist_vector.begin(); veci != dist_vector.end(); ++veci)
      {
        if (degree == (*veci).first)
        {
          // if it exists, then put this vertex in the list
          (*veci).second.push_back(*vi);
          vertex_added = true;
          break;
        }
      }

      if (!vertex_added)
      {
        // if it does not, then create a new k_vertices_pair and put it in dist_vector
        std::vector<vertex> new_vertices_vector;
        new_vertices_vector.push_back(*vi);
        dist_vector.push_back(std::make_pair(degree, new_vertices_vector));
      }
    }

    // sort dist_vector
    std::sort(dist_vector.begin(), dist_vector.end(), detail::k_vertices_pair_less_than);

    // set input data for fitting
    double cov00, cov01, cov11, sumsq;

    DegreeList degrees = non_zero_degrees();
    int num_non_zero_degrees = degrees.size();

    double p_k[num_non_zero_degrees],
           k[num_non_zero_degrees],
           log_k[num_non_zero_degrees],
           log_p_k[num_non_zero_degrees],
           exp_p_k[num_non_zero_degrees];

    int count = 0;
    for (typename DegreeList::iterator li = degrees.begin(); li != degrees.end(); ++li, ++count)
    {
      k[count] = *li;
      log_k[count] = log(*li);

      p_k[count] = P(*li);
      log_p_k[count] = log(p_k[count]);
      exp_p_k[count] = exp(p_k[count]);
    }

    // linear fitting
    {
      double linear_c0, linear_c1;
      gsl_fit_linear(k, 1,    // x = k
                     p_k, 1,  // y = p(k)
                     num_non_zero_degrees, &linear_c0, &linear_c1,
                     &cov00, &cov01, &cov11, &sumsq);
      double r = gsl_stats_correlation(k, 1, p_k, 1, num_non_zero_degrees);
#ifdef DISTRIBUTIONS_DEBUG
      std::cout << "linear regression:" << std::endl
                << '\t' << "r = " << r << std::endl
                << '\t' << "slope = " << linear_c1 << std::endl;
#endif

      linear_fit_parameters[0] = r;
      linear_fit_parameters[1] = linear_c1;
      linear_fit_parameters[2] = linear_c0;
    }

    // power fitting
    {
      double power_c0, power_c1;
      gsl_fit_linear(log_k, 1,    // x = log(k)
                     log_p_k, 1,  // y = log(p(k))
                     num_non_zero_degrees, &power_c0, &power_c1,
                     &cov00, &cov01, &cov11, &sumsq);
      double r = gsl_stats_correlation(log_k, 1, log_p_k, 1, num_non_zero_degrees);
#ifdef DISTRIBUTIONS_DEBUG
      std::cout << "power regression:" << std::endl
                << '\t' << "r = " << r << std::endl
                << '\t' << "gamma = " << -power_c1 << std::endl;
#endif

      power_fit_parameters[0] = r;
      power_fit_parameters[1] = power_c1;
      power_fit_parameters[2] = exp(power_c0);
    }

    // exponential fitting
    {
      double exp_c0, exp_c1;
      gsl_fit_linear(k, 1,        // x = k
                     log_p_k, 1,  // y = log(p(k))
                     num_non_zero_degrees, &exp_c0, &exp_c1,
                     &cov00, &cov01, &cov11, &sumsq);
      double r = gsl_stats_correlation(k, 1, log_p_k, 1, num_non_zero_degrees);
#ifdef DISTRIBUTIONS_DEBUG
      std::cout << "exponential regression:" << std::endl
                << '\t' << "r = " << r << std::endl
                << '\t' << "lambda = " << exp_c1 << std::endl
                << '\t' << "A = " << exp(exp_c0) << std::endl;
#endif

      exp_fit_parameters[0] = r;
      exp_fit_parameters[1] = exp_c1;
      exp_fit_parameters[2] = exp(exp_c0);
    }

    // logaritmic fitting
    {
      double log_c0, log_c1;
      gsl_fit_linear(k, 1,        // x = k
                     exp_p_k, 1,  // y = exp(p(k))
                     num_non_zero_degrees, &log_c0, &log_c1,
                     &cov00, &cov01, &cov11, &sumsq);
      double r = gsl_stats_correlation(k, 1, exp_p_k, 1, num_non_zero_degrees);

      log_fit_parameters[0] = r;
      log_fit_parameters[1] = log_c1;
      log_fit_parameters[2] = log(log_c0);
    }

    // select best fit by comparing theirs squared pearson coefficient
    best_fit = linear;
    double max_rr = linear_fit_parameters[0] * linear_fit_parameters[0];
    if (power_fit_parameters[0] * power_fit_parameters[0] > max_rr)
    {
      best_fit = power;
      max_rr = power_fit_parameters[0] * power_fit_parameters[0];
    }
    if (exp_fit_parameters[0] * exp_fit_parameters[0] > max_rr)
    {
      best_fit = exponential;
      max_rr = exp_fit_parameters[0] * exp_fit_parameters[0];
    }
    if (log_fit_parameters[0] * log_fit_parameters[0] > max_rr)
      best_fit = logaritmic;

#ifdef GOODNESS_OF_FIT
    // goodness-of-fit test
    { // FIXME: this must be reviewed
      size_t df = max_non_zero_degree() - min_non_zero_degree(); // degrees of freedom
#ifdef CONAN_DEBUG
      std::cerr << "degrees of freedom = " << df << std::endl;
#endif
      double linear_chisq = 0,
             power_chisq = 0,
             exp_chisq = 0,
             log_chisq = 0;
      for (size_t i = min_non_zero_degree(); i <= max_non_zero_degree(); ++i)
      {
        double p_k = P(i),
               linear_expected_p_k = linear_fit_parameters[1] * i + linear_fit_parameters[2],
               power_expected_p_k = power_fit_parameters[2] * pow(i, power_fit_parameters[1]),
               exp_expected_p_k = exp_fit_parameters[2] * exp(exp_fit_parameters[1] * i),
               log_expected_p_k = log_fit_parameters[2] * log(log_fit_parameters[1] * i);
#ifdef DISTRIBUTIONS_DEBUG
        std::cerr << "P(" << i << ") = " << p_k << std::endl
                  << "linear_expected_p_k = " << linear_expected_p_k << std::endl
                  << "power_expected_p_k = " << power_expected_p_k << std::endl
                  << "exp_expected_p_k = " << exp_expected_p_k << std::endl
                  << "log_expected_p_k = " << log_expected_p_k << std::endl;
#endif

        linear_chisq += pow(p_k - linear_expected_p_k, 2) / linear_expected_p_k;
        power_chisq += pow(p_k - power_expected_p_k, 2) / power_expected_p_k;
        exp_chisq += pow(p_k - exp_expected_p_k, 2) / exp_expected_p_k;
        log_chisq += pow(p_k - log_expected_p_k, 2) / log_expected_p_k;
      }
#ifdef CONAN_DEBUG
      std::cerr << "linear_chisq = " << linear_chisq << std::endl
                << "power_chisq = " << power_chisq << std::endl
                << "exp_chisq = " << exp_chisq << std::endl
                << "log_chisq = " << log_chisq << std::endl;

      double theoretical_chisq = detail::critchi(0.05, df); // alpha = 0.05
      std::cerr << "theoretical_chisq = " << theoretical_chisq << std::endl;

      if (linear_chisq < theoretical_chisq)
        std::cerr << "The degree distribution is a linear distribution" << std::endl;
      if (power_chisq < theoretical_chisq)
        std::cerr << "The degree distribution is a power distribution" << std::endl;
      if (exp_chisq < theoretical_chisq)
        std::cerr << "The degree distribution is a exponential distribution" << std::endl;
      if (log_chisq < theoretical_chisq)
        std::cerr << "The degree distribution is a logaritmic distribution" << std::endl;
#endif
    }
#endif

    return;
  }

  template <class Graph, int DegreeType>
  inline decimal degree_distribution<Graph, DegreeType>::P(
      degree_size_type k
      )
  {
    decimal prob = 0.0;

    for (typename DistVector::iterator veci = dist_vector.begin(); veci != dist_vector.end(); ++veci)
    {
      if (k == (*veci).first)
        return prob + (decimal)(*veci).second.size() / V;
      else if (k > (*veci).first)
        break;
      else if (cum)
        prob += (decimal)(*veci).second.size() / V;
    }

    return prob;
  }

  template <class Graph, int DegreeType>
  typename degree_distribution<Graph, DegreeType>::VertexList
  degree_distribution<Graph, DegreeType>::vertices(
      degree_size_type k
      )
  {
    VertexList tmp_list;

    for (typename DistVector::iterator veci = dist_vector.begin(); veci != dist_vector.end(); ++veci)
    {
      if (k == (*veci).first)
      {
        tmp_list = (*veci).second;
        return tmp_list;
      }
    }

    return tmp_list;
  }

  template <class Graph, int DegreeType>
  inline typename degree_distribution<Graph, DegreeType>::DegreeList
  degree_distribution<Graph, DegreeType>::non_zero_degrees()
  {
    DegreeList tmp_list;

    for (typename DistVector::iterator veci = dist_vector.begin(); veci != dist_vector.end(); ++veci)
      tmp_list.push_back((*veci).first);

    return tmp_list;
  }

  template <class Graph, int DegreeType>
  inline typename degree_distribution<Graph, DegreeType>::degree_size_type
  degree_distribution<Graph, DegreeType>::min_non_zero_degree()
  {
    degree_size_type min = V;

    for (typename DistVector::iterator veci = dist_vector.begin(); veci != dist_vector.end(); ++veci)
      if ((*veci).first < min)
        min = (*veci).first;

    return min;
  }

  template <class Graph, int DegreeType>
  inline typename degree_distribution<Graph, DegreeType>::degree_size_type
  degree_distribution<Graph, DegreeType>::max_non_zero_degree()
  {
    degree_size_type max = 0;

    for (typename DistVector::iterator veci = dist_vector.begin(); veci != dist_vector.end(); ++veci)
      if ((*veci).first > max)
        max = (*veci).first;

    return max;
  }

  template <class Graph, int DegreeType>
  inline decimal degree_distribution<Graph, DegreeType>::gamma()
  {
    return -power_fit_parameters[1];
  }

  template <class Graph, int DegreeType>
  inline decimal degree_distribution<Graph, DegreeType>::lambda()
  {
    return -exp_fit_parameters[1];
  }

  template <class Graph, int DegreeType>
  inline decimal degree_distribution<Graph, DegreeType>::slope()
  {
    return linear_fit_parameters[1];
  }

  template <class Graph, int DegreeType>
  inline decimal degree_distribution<Graph, DegreeType>::r(
      int fit_type // = power
      )
  {
    if (fit_type == linear)
      return linear_fit_parameters[0];
    else if (fit_type == power)
      return power_fit_parameters[0];
    else if (fit_type == exponential)
      return exp_fit_parameters[0];
    else if (fit_type == logaritmic)
      return log_fit_parameters[0];
    else
      return 0;
  }

  template <class Graph, int DegreeType>
  inline decimal degree_distribution<Graph, DegreeType>::rr(
      int fit_type // = power
      )
  {
    if (fit_type == linear)
      return linear_fit_parameters[0] * linear_fit_parameters[0];
    else if (fit_type == power)
      return power_fit_parameters[0] * power_fit_parameters[0];
    else if (fit_type == exponential)
      return exp_fit_parameters[0] * exp_fit_parameters[0];
    else if (fit_type == logaritmic)
      return log_fit_parameters[0] * log_fit_parameters[0];
    else
      return 0;
  }

  template <class Graph, int DegreeType>
  inline double* degree_distribution<Graph, DegreeType>::_fit_parameters(
      int fit_type
      )
  {
    if (fit_type == linear)
      return linear_fit_parameters;
    else if (fit_type == power)
      return power_fit_parameters;
    else if (fit_type == exponential)
      return exp_fit_parameters;
    else if (fit_type == logaritmic)
      return log_fit_parameters;
    else
      throw std::runtime_error("fit_type unknown");
  }

  template <class Graph, int DegreeType>
  inline decimal degree_distribution<Graph, DegreeType>::entropy()
  {
    DegreeList dl = non_zero_degrees();

    decimal S = 0;

    for (typename DegreeList::iterator dli = dl.begin(); dli != dl.end(); ++dli)
    {
      S += P(*dli) * conan_log2( P(*dli) );
    }
    
    return -S;
  }

  /// @cond
  // General distribution class
  template <class Graph>
  struct distribution
  {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex;
    typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iter;

    typedef typename std::vector<vertex> VertexList;
    typedef typename std::vector<decimal> ValueList;

    struct Interval
    {
      decimal min_value;
      decimal max_value;
      VertexList vertices;

      Interval(decimal min, decimal max)
        : min_value(min), max_value(max), vertices()
      { }

      Interval(decimal min, decimal max, vertex v)
        : min_value(min), max_value(max), vertices()
      {
        vertices.push_back(v);
      }

      Interval(decimal min, decimal max, VertexList vl)
        : min_value(min), max_value(max), vertices(vl)
      { }

      size_t num_vertices()
      { return vertices.size(); }

      decimal mean()
      { return (max_value - min_value) / 2; }
    };
    typedef typename std::vector<Interval> DistVector;

    distribution(
        Graph& g,
        size_t num_bins,
        decimal (*fcn)(vertex, Graph&)
        )
      : dist_vector()
    {
      // set V and num_bins
      this->V = boost::num_vertices(g);
      this->num_bins = num_bins;

      // add vertices' degree to dist_vector
#if defined(__GNUC__) && (__GNUC__ < 4 || (__GNUC__ == 4 && __GNUC_MINOR__ < 3))
      __gnu_cxx::hash_map<vertex, decimal> value_hash;
#else
      std::unordered_map<vertex, decimal> value_hash;
#endif
      min_value = std::numeric_limits<decimal>::max();
      max_value = 0;
      vertex_iter vi, viend;
      for (tie(vi, viend) = boost::vertices(g); vi != viend; ++vi)
      {
        decimal value = fcn(*vi, g);
        value_hash[*vi] = value;
        if (min_value > value)
          min_value = value;
        if (max_value < value)
          max_value = value;
      }

      step = (max_value - min_value) / num_bins;
      for (decimal i = min_value; i < max_value; i += step)
        dist_vector.push_back( Interval(i, i + step) );

      for (tie(vi, viend) = boost::vertices(g); vi != viend; ++vi)
      {
        bool vertex_added = false;
        for (typename DistVector::iterator veci = dist_vector.begin(); veci != dist_vector.end(); ++veci)
          if (value_hash[*vi] >= (*veci).min_value && value_hash[*vi] <= (*veci).max_value)
          {
            (*veci).vertices.push_back(*vi);
            vertex_added = true;
            break;
          }

        if (!vertex_added)
          std::cerr << "distribution::distribution: problem has occur adding vertex " << *vi << std::endl;
      }
      return;
    }

    decimal P(
        size_t bin
        )
    { return decimal(dist_vector[bin].num_vertices()) / V; }

    decimal P(
        decimal value
        )
    {
      for (typename DistVector::iterator veci = dist_vector.begin(); veci != dist_vector.end(); ++veci)
      {
        if (value >= (*veci).min_value && value <= (*veci).max_value)
          return decimal( (*veci).num_vertices() ) / V;
      }
      return 0.0;
    }

    VertexList vertices(
        size_t bin
        )
    { return dist_vector[bin].vertices; }

    VertexList vertices(
        decimal value
        )
    {
      for (typename DistVector::iterator veci = dist_vector.begin(); veci != dist_vector.end(); ++veci)
      {
        if (value >= (*veci).min_value && value <= (*veci).max_value)
          return (*veci).vertices;
      }
      return VertexList();
    }

    // Data members
    size_t V;
    size_t num_bins;
    decimal min_value, max_value, step;
    DistVector dist_vector;

#if 0
    int best_fit;
    double linear_fit_parameters[3];  // r, m(slope), n(constant)... y = m * x + n
    double power_fit_parameters[3];   // r, gamma, A... y = A * x^gamma + b
    double exp_fit_parameters[3];     // r, lambda, A... y = A * exp(lambda * x) + b
    double log_fit_parameters[3];     // r, lambda, A... y = A * log(lambda * x) + b
#endif
  };

  /// @endcond

}

#endif //DISTRIBUTIONS_HPP
