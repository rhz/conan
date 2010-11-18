#ifndef PYTHON_MI_HPP
#define PYTHON_MI_HPP

#include "python_undirected_graph.hpp"
#include <conan/inference/mutual_information.hpp>

template<class T>
struct std_item
{
  typedef typename T::value_type V;
  static V& get(T & x, int i)
  {
    if( i < 0 ) i += x.size();
    if( i >= 0 && i < (int)x.size() ) return x[i];
    throw std::runtime_error("Index out of range");
  }
  static void set(T & x, int i, const V & v)
  {
    if( i < 0 ) i += x.size();
    if( i >= 0 && i < (int)x.size() ) x[i] = v;
    else throw std::runtime_error("Index out of range");
  }
  static void del(T & x, int i)
  {
    if( i < 0 ) i += x.size();
    if( i >= 0 && i < (int)x.size() )
    {
      typename T::iterator it = x.begin() + i;
      x.erase(it);
    }
    else throw std::runtime_error("Index out of range");
  }
  static void add(T & x, const V & v)
  {
    x.push_back(v);
  }
};

python_undirected_graph mutual_information(
    aracne::Parameter& p
    )
{
  return python_undirected_graph(
      conan::inference::mi::mutual_information<python_undirected_graph::Graph>(p));
}

BOOST_PYTHON_MODULE(_mi)
{
  class_<std::vector<std::string> >("string_list")
    .def("__len__", &std::vector<std::string> ::size)
    .def("clear", &std::vector<std::string> ::clear)
    .def("append", &std_item<std::vector<std::string> >::add,
          with_custodian_and_ward<1,2>()) // to let container keep value
    .def("__getitem__", &std_item<std::vector<std::string> >::get,
          return_value_policy<copy_non_const_reference>())
    .def("__setitem__", &std_item<std::vector<std::string> >::set,
          with_custodian_and_ward<1,2>()) // to let container keep value
    .def("__delitem__", &std_item<std::vector<std::string> >::del)
    .def("__iter__", boost::python::iterator<std::vector<std::string> >())
    ;
  class_< aracne::Parameter >("parameters")
    .def_readwrite("threshold", &aracne::Parameter::threshold)
    .def_readwrite("pvalue", &aracne::Parameter::pvalue)
    .def_readwrite("eps", &aracne::Parameter::eps)
    .def_readwrite("sigma", &aracne::Parameter::sigma)
    .def_readwrite("miSteps", &aracne::Parameter::miSteps)
    .def_readwrite("sample", &aracne::Parameter::sample)
    .def_readwrite("percent", &aracne::Parameter::percent)
    .def_readwrite("mean", &aracne::Parameter::mean)
    .def_readwrite("cv", &aracne::Parameter::cv)
    .def_readwrite("correction", &aracne::Parameter::correction)
    .def_readwrite("verbose", &aracne::Parameter::verbose)
    .def_readwrite("algorithm", &aracne::Parameter::algorithm)
    .def_readwrite("infile", &aracne::Parameter::infile)
    .def_readwrite("outfile", &aracne::Parameter::outfile)
    .def_readwrite("adjfile", &aracne::Parameter::adjfile)
    .def_readwrite("hub", &aracne::Parameter::hub)
    .def_readwrite("controlId", &aracne::Parameter::controlId)
    .def_readwrite("condition", &aracne::Parameter::condition)
    .def_readwrite("subnet", &aracne::Parameter::subnet)
    .def_readwrite("tf_list", &aracne::Parameter::tf_list)
    ;
  def("mutual_information_two_genes", &conan::inference::mi::mutual_information_two_genes);
  def("mutual_information", &mutual_information);
}

#endif // PYTHON_MI_HPP
