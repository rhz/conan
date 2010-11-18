## Makefile.am -- Process this file with automake to produce Makefile.in
INCLUDES = -I../ -I../contrib/
AM_CPPFLAGS = -DBOOST_TEST_DYN_LINK -DBOOST_NO_SLIST -DBOOST_NO_HASH
AM_CXXFLAGS = -std=c++0x
check_PROGRAMS = combinatoria undirected_graph \
    properties/centrality properties/entropy
combinatoria_SOURCES = combinatoria.cpp
undirected_graph_SOURCES = undirected_graph.cpp
properties_centrality_SOURCES = properties/centrality.cpp
properties_entropy_SOURCES = properties/entropy.cpp
TESTS = combinatoria undirected_graph \
    properties/centrality properties/entropy
XFAIL_TESTS = properties/entropy
