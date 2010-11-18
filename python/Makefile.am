## Makefile.am -- Process this file with automake to produce Makefile.in
pkgpython_PYTHON = conan/__init__.py conan/internals.py conan/graph_models.py \
    conan/inference.py conan/dynamics.py
INCLUDES = -I../ -I../contrib/
AM_CXXFLAGS = -std=c++0x
AM_CPPFLAGS = -DBOOST_NO_SLIST -DBOOST_NO_HASH -DBOOST_PYTHON_DYNAMIC_LIB \
    -ftemplate-depth-128 -isystem @PYTHON_INCLUDE_DIR@
AM_LDFLAGS = -avoid-version -module -no-undefined
pkgpyexec_LTLIBRARIES = _conan.la _graph_models.la _inference.la _mi.la _dynamics.la
_conan_la_SOURCES = conan_python.cpp vertex_and_edge_wrap.hpp \
    python_undirected_graph.hpp python_directed_graph.hpp \
    misc_classes.hpp debug.hpp
_conan_la_LDFLAGS = -avoid-version -module -no-undefined @TINYXML_LIB@ -ltinyxml
_graph_models_la_SOURCES = graph_models.cpp
_inference_la_SOURCES = inference.cpp
_mi_la_SOURCES = mi.cpp
_mi_la_LDFLAGS = -avoid-version -module -no-undefined @ARACNE_LIB@ -laracne
_dynamics_la_SOURCES = dynamics.cpp
