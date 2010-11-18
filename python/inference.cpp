#ifndef PYTHON_INFERENCE_HPP
#define PYTHON_INFERENCE_HPP

#include "python_undirected_graph.hpp"
#include <conan/inference/maximum_entropy.hpp>
#include <conan/inference/clustering.hpp>
#include <conan/utils/matrix.hpp>

#include <boost/python/tuple.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/raw_function.hpp>


using conan::decimal;
using conan::inference::matrix_mask;

gsl_matrix * list_to_matrix(
    list input_matrix
    )
{
  // Convert the input list in a gsl_matrix
  size_t num_elems = extract<size_t>(input_matrix.attr("__len__")());

  size_t num_obs = 0;
  {
    list first_row = extract<list>(input_matrix[0]);
    num_obs = extract<size_t>(first_row.attr("__len__")());
  }

  gsl_matrix * A = gsl_matrix_calloc(num_elems, num_obs);

  for (size_t r = 0; r < num_elems; ++r)
  {
    list row = extract<list>(input_matrix[r]);

    for (size_t c = 0; c < num_obs; ++c)
      gsl_matrix_set(A, r, c, extract<double>(row[c]));
  }

  return A;
}


void list_to_matrix_and_matrix_mask(
    list input_matrix,
    gsl_matrix * A,
    matrix_mask & Amask
    )
{
  // Convert the input list in a gsl_matrix
  size_t num_elems = extract<size_t>(input_matrix.attr("__len__")());

  size_t num_obs = 0;
  {
    list first_row = extract<list>(input_matrix[0]);
    num_obs = extract<size_t>(first_row.attr("__len__")());
  }

  for (size_t r = 0; r < num_elems; ++r)
  {
    list row = extract<list>(input_matrix[r]);

    for (size_t c = 0; c < num_obs; ++c)
    {
      object result = extract<object>(row[c]);
      if (result.ptr() == Py_None)
      { // None case
        Amask(r, c) = true;
        gsl_matrix_set(A, r, c, 0.0);
      }
      else
      {
        gsl_matrix_set(A, r, c, extract<double>(row[c]));
      }
    }
  }

  return;
}


list read_matrix(
    std::string input_filename,
    char sep_char = ' '
    )
{
  gsl_matrix * m = conan::linalg::alloc_and_read_matrix(input_filename, sep_char);
  list out;
  for (int r = 0; r < m->size1; ++r)
  {
    list row;
    for (int c = 0; c < m->size2; ++c)
      row.append( gsl_matrix_get(m, r, c) );
    out.append(row);
  }
  gsl_matrix_free(m);
  return out;
}


list maximum_entropy(
    tuple args,
    dict kw
    )
{
  using namespace conan::inference;

  list data_matrix = extract<list>(args[0]);

  size_t num_elems = extract<size_t>(data_matrix.attr("__len__")());

  size_t num_obs = 0;
  {
    list first_row = extract<list>(data_matrix[0]);
    num_obs = extract<size_t>(first_row.attr("__len__")());
  }

  gsl_matrix * A = gsl_matrix_calloc(num_elems, num_obs);

  matrix_mask Amask(num_elems, num_obs);

  list_to_matrix_and_matrix_mask(data_matrix, A, Amask);

  gsl_matrix * interaction = gsl_matrix_calloc(num_elems, num_elems);

  std::string coexpr_measure = "covariance";
  try
  {
    extract<std::string> get_coexpr_measure(kw["coexpr_measure"]);
    if (get_coexpr_measure.check())
      coexpr_measure = get_coexpr_measure();
  }
  catch (...) { }

  if (coexpr_measure == "covariance")
  {
    maximum_entropy(A, interaction, covarianceS());
  }
  else if (coexpr_measure == "correlation")
  {
    maximum_entropy(A, interaction, correlationS());
  }
  else if (coexpr_measure == "covariance_pairwise_complete_obs")
  {
    size_t min_num_complete_obs = 5;
    try
    {
      extract<size_t> get_min_num_complete_obs(kw["min_num_complete_obs"]);
      if (get_min_num_complete_obs.check())
        min_num_complete_obs = get_min_num_complete_obs();
    }
    catch (...) { }

    covariance_pairwise_complete_obs_S filter(Amask, min_num_complete_obs);

    maximum_entropy(A, interaction, filter);
  }
  else
  {
    throw std::runtime_error("Invalid coexpression measure selected");
  }

  gsl_matrix_free(A);

  // Convert interaction matrix to a list of lists
  list interaction_matrix;

  for (size_t r = 0; r < interaction->size1; ++r)
  {
    list row;

    for (size_t c = 0; c < interaction->size2; ++c)
      row.append( gsl_matrix_get(interaction, r, c) );

    interaction_matrix.append(row);
  }

  gsl_matrix_free(interaction);

  return interaction_matrix;
}


python_undirected_graph maximum_entropy_with_fixed_threshold(
    list data_matrix,
    decimal threshold,
    bool abs_values = true,
    bool normalize = true
    )
{
  typedef python_undirected_graph::Graph Graph;

  gsl_matrix * A = list_to_matrix(data_matrix);

  Graph g(conan::inference::maximum_entropy_with_fixed_threshold<Graph>(
        A, threshold, abs_values, normalize));

  gsl_matrix_free(A);

  return python_undirected_graph(g);
}


python_undirected_graph maximum_entropy_with_fixed_num_edges(
    list data_matrix,
    size_t num_edges,
    bool abs_values = true,
    bool normalize = true
    )
{
  typedef python_undirected_graph::Graph Graph;

  gsl_matrix * A = list_to_matrix(data_matrix);

  Graph g(conan::inference::maximum_entropy_with_fixed_num_edges<Graph>(
        A, num_edges, abs_values, normalize));

  gsl_matrix_free(A);

  return python_undirected_graph(g);
}


python_undirected_graph maximum_entropy_with_gaussian_threshold(
    list data_matrix,
    decimal threshold,
    bool normalize = true
    )
{
  typedef python_undirected_graph::Graph Graph;

  gsl_matrix * A = list_to_matrix(data_matrix);

  Graph g(conan::inference::maximum_entropy_with_gaussian_threshold<Graph>(
        A, threshold, normalize));

  gsl_matrix_free(A);

  return python_undirected_graph(g);
}


python_undirected_graph maximum_entropy_with_bootstrapping(
    tuple args,
    dict kw
    )
{
  using namespace conan::inference;
  typedef python_undirected_graph::Graph Graph;

  // Required parameters
  list data_matrix = extract<list>(args[0]);
  size_t num_steps = extract<size_t>(args[1]);
  decimal threshold = extract<decimal>(args[2]);
  // Optional parameters
  // normalize: can be True or False
  // coexpr_measure: can be "covariance", "correlation" or "covariance_pairwise_complete_obs"
  // min_num_complete_obs (only used if coexpr_measure == "covariance_pairwise_complete_obs"): integer
  // shuffle_by_column (only used if coexpr_measure == "covariance_pairwise_complete_obs"): bool

  bool normalize = false;
  try
  {
    extract<bool> get_normalize(kw["normalize"]);
    if (get_normalize.check())
      normalize = get_normalize();
  }
  catch (...) { }

  size_t num_elems = extract<size_t>(data_matrix.attr("__len__")());

  size_t num_obs = 0;
  {
    list first_row = extract<list>(data_matrix[0]);
    num_obs = extract<size_t>(first_row.attr("__len__")());
  }

  gsl_matrix * A = gsl_matrix_calloc(num_elems, num_obs);

  matrix_mask Amask(num_elems, num_obs);

  list_to_matrix_and_matrix_mask(data_matrix, A, Amask);

  // Get coexpression measure to use
  std::string coexpr_measure = "covariance";
  try
  {
    extract<std::string> get_coexpr_measure(kw["coexpr_measure"]);
    if (get_coexpr_measure.check())
      coexpr_measure = get_coexpr_measure();
  }
  catch (...) { }

  // Compute the graph
  Graph g;

  if (coexpr_measure == "covariance")
  {
    g = maximum_entropy_with_bootstrapping<Graph>(A, num_steps, threshold, covarianceS(), normalize);
  }
  else if (coexpr_measure == "correlation")
  {
    g = maximum_entropy_with_bootstrapping<Graph>(A, num_steps, threshold, correlationS(), normalize);
  }
  else if (coexpr_measure == "covariance_pairwise_complete_obs")
  {
    size_t min_num_complete_obs = 5;
    try
    {
      extract<size_t> get_min_num_complete_obs(kw["min_num_complete_obs"]);
      if (get_min_num_complete_obs.check())
        min_num_complete_obs = get_min_num_complete_obs();
    }
    catch (...) { }

    bool shuffle_by_column = false;
    try
    {
      extract<size_t> get_shuffle_by_column(kw["shuffle_by_column"]);
      if (get_shuffle_by_column.check())
        shuffle_by_column = get_shuffle_by_column();
    }
    catch (...) { }

    conan::inference::covariance_pairwise_complete_obs_S filter(Amask, min_num_complete_obs, shuffle_by_column);

    g = maximum_entropy_with_bootstrapping<Graph>(A, num_steps, threshold, filter, normalize);
  }
  else
  {
    throw std::runtime_error("Invalid coexpression measure selected");
  }

  gsl_matrix_free(A);

  return python_undirected_graph(g);
}


python_undirected_graph clustering_with_fixed_threshold(
    list data_matrix,
    decimal threshold,
    bool abs_values = true,
    bool covariance = true,
    bool normalize = true
    )
{
  typedef python_undirected_graph::Graph Graph;

  gsl_matrix * A = list_to_matrix(data_matrix);

  Graph g(conan::inference::clustering_with_fixed_threshold<Graph>(
        A, threshold, abs_values, covariance, normalize));

  gsl_matrix_free(A);

  return python_undirected_graph(g);
}


python_undirected_graph clustering_with_bootstrapping(
    list data_matrix,
    size_t num_steps,
    decimal threshold,
    bool normalize = true,
    std::string coexpr_measure = "covariance"
    )
{
  using namespace conan::inference;
  typedef python_undirected_graph::Graph Graph;

  size_t num_elems = extract<size_t>(data_matrix.attr("__len__")());

  size_t num_obs = 0;
  {
    list first_row = extract<list>(data_matrix[0]);
    num_obs = extract<size_t>(first_row.attr("__len__")());
  }

  gsl_matrix * A = gsl_matrix_calloc(num_elems, num_obs);

  matrix_mask Amask(num_elems, num_obs);

  list_to_matrix_and_matrix_mask(data_matrix, A, Amask);

  Graph g;
  if (coexpr_measure == "covariance")
    g = clustering_with_bootstrapping<Graph>(A, num_steps, threshold, covarianceS(), normalize);
  else if (coexpr_measure == "correlation")
    g = clustering_with_bootstrapping<Graph>(A, num_steps, threshold, correlationS(), normalize);
  else
  {
    size_t min_num_complete_obs = 5;
    conan::inference::covariance_pairwise_complete_obs_S filter(Amask, min_num_complete_obs);
    g = clustering_with_bootstrapping<Graph>(A, num_steps, threshold, filter, normalize);
  }

  gsl_matrix_free(A);

  return python_undirected_graph(g);
}


BOOST_PYTHON_FUNCTION_OVERLOADS(read_matrix_overloads, read_matrix, 1, 2)
BOOST_PYTHON_FUNCTION_OVERLOADS(maximum_entropy_with_fixed_threshold_overloads, maximum_entropy_with_fixed_threshold, 2, 4)
BOOST_PYTHON_FUNCTION_OVERLOADS(maximum_entropy_with_fixed_num_edges_overloads, maximum_entropy_with_fixed_num_edges, 2, 4)
BOOST_PYTHON_FUNCTION_OVERLOADS(maximum_entropy_with_gaussian_threshold_overloads, maximum_entropy_with_gaussian_threshold, 2, 3)
BOOST_PYTHON_FUNCTION_OVERLOADS(clustering_with_fixed_threshold_overloads, clustering_with_fixed_threshold, 2, 5)
BOOST_PYTHON_FUNCTION_OVERLOADS(clustering_with_bootstrapping_overloads, clustering_with_bootstrapping, 3, 5)


BOOST_PYTHON_MODULE(_inference)
{
  def("read_matrix", read_matrix, read_matrix_overloads());

  def("maximum_entropy", raw_function(maximum_entropy, 1));
  def("maximum_entropy_with_bootstrapping", raw_function(maximum_entropy_with_bootstrapping, 3));

  def("maximum_entropy_with_fixed_threshold", maximum_entropy_with_fixed_threshold, maximum_entropy_with_fixed_threshold_overloads());
  def("maximum_entropy_with_fixed_num_edges", maximum_entropy_with_fixed_num_edges, maximum_entropy_with_fixed_num_edges_overloads());
  def("maximum_entropy_with_gaussian_threshold", maximum_entropy_with_gaussian_threshold, maximum_entropy_with_gaussian_threshold_overloads());

  def("clustering_with_fixed_threshold", clustering_with_fixed_threshold, clustering_with_fixed_threshold_overloads());
  def("clustering_with_bootstrapping", clustering_with_bootstrapping, clustering_with_bootstrapping_overloads());
}

#endif // PYTHON_INFERENCE_HPP
