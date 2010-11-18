## Makefile.am -- Process this file with automake to produce Makefile.in
AM_CXXFLAGS = -Wall -Wno-long-long -ansi -pedantic -std=c++0x
borlanddir = ${pkgincludedir}/Borland
referencedir = ${pkgincludedir}/Reference
pkginclude_HEADERS = Matrix.h Matrix_Op.h Microarray_Set.h MutualInfo.h param.h \
    parseargs.h TypeManip.h
borland_HEADERS = Borland/TypeManip.h
reference_HEADERS = Reference/TypeManip.h
pkgdata_DATA = config_kernel.txt config_threshold.txt usage.txt
testdir = ${pkgdatadir}/test
test_DATA = test/arraydata100x336.exp test/arraydata10x336.exp
lib_LTLIBRARIES = libaracne.la
libaracne_la_SOURCES = main.cpp Matrix.cpp Matrix.h Matrix_Op.cpp Matrix_Op.h \
    Microarray_Set.cpp Microarray_Set.h MutualInfo.cpp MutualInfo.h \
    param.cpp param.h parseargs.h TypeManip.h Borland/TypeManip.h Reference/TypeManip.h
libaracne_la_LDFLAGS = -no-undefined
