## Makefile.am -- Process this file with automake to produce Makefile.in
tinyxmldir = $(includedir)/tinyxml
AM_CPPFLAGS = -DTIXML_USE_STL
lib_LTLIBRARIES = libtinyxml.la
libtinyxml_la_SOURCES = tinyxml.cpp tinystr.cpp tinyxml.h tinystr.h \
    tinyxmlparser.cpp tinyxmlerror.cpp xmltest.cpp
libtinyxml_la_LDFLAGS = -no-undefined
tinyxml_HEADERS = tinyxml.h tinystr.h
