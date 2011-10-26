## Makefile.am -- Process this file with automake to produce Makefile.in
INCLUDES = -I../ -I../contrib/
AM_CPPFLAGS = -DBOOST_TEST_DYN_LINK -DBOOST_NO_SLIST -DBOOST_NO_HASH -DDETERMINISTIC_RNG
AM_CXXFLAGS = -std=c++0x
check_PROGRAMS = combinatoria undirected_graph \
    properties/centrality properties/entropy
bin_PROGRAMS = exec/aracne exec/props exec/communities exec/graphml2dot
combinatoria_SOURCES = combinatoria.cpp
undirected_graph_SOURCES = undirected_graph.cpp
properties_centrality_SOURCES = properties/centrality.cpp
properties_entropy_SOURCES = properties/entropy.cpp
exec_aracne_SOURCES = exec/aracne.cpp
exec_aracne_LDFLAGS = @TINYXML_LIB@ -ltinyxml @ARACNE_LIB@ -laracne
exec_props_SOURCES = exec/props.cpp
exec_props_LDFLAGS = @TINYXML_LIB@ -ltinyxml
exec_communities_SOURCES = exec/communities.cpp
exec_communities_LDFLAGS = @TINYXML_LIB@ -ltinyxml
exec_graphml2dot_SOURCES = exec/graphml2dot.cpp
exec_graphml2dot_LDFLAGS = @TINYXML_LIB@ -ltinyxml
TESTS = combinatoria undirected_graph \
    properties/centrality properties/entropy
XFAIL_TESTS = properties/entropy
