# test which gcc version we have
gcc_version := $(shell gcc -v 2>&1 |tail -n 1 |cut -d\  -f3 |sed 's/\.//g')
# if gcc version >= 4.3.0 then use the coming standard extensions
gcc43_opts := $(shell [ "$(gcc_version)" -gt 430 ] && echo '-std=gnu++0x')

CXX = g++
INCLUDES = -I$(CONAN_PATH) -I$(CONAN_PATH)/contrib

BOOST_OPT = -DBOOST_NO_SLIST -DBOOST_NO_HASH
FLAGS =
OPTS = $(BOOST_OPT) $(FLAGS) -Wall $(gcc43_opts) $(OPTIMIZATION_FLAGS)

#BLAS_LIBS = -lcblas -latlas
BLAS_LIBS = -lgslcblas
LIBS = -L/usr/lib/sse2/ -lgsl $(BLAS_LIBS) -lm $(EXTRA_LIBS)
CXXFLAGS = $(OPTS) $(INCLUDES) $(LIBS)

TINYXML_PATH = $(CONAN_PATH)/contrib/tinyxml
ARACNE_PATH = $(CONAN_PATH)/contrib/aracne

TINYXML_SRCS = \
	$(TINYXML_PATH)/tinyxml.cpp \
	$(TINYXML_PATH)/tinyxmlparser.cpp \
	$(TINYXML_PATH)/tinyxmlerror.cpp \
	$(TINYXML_PATH)/tinystr.cpp \

TINYXML_OBJS = \
	$(TINYXML_PATH)/tinyxml.o \
	$(TINYXML_PATH)/tinyxmlerror.o \
	$(TINYXML_PATH)/tinyxmlparser.o

ARACNE_SRCS = \
	$(ARACNE_PATH)/Matrix.cpp \
	$(ARACNE_PATH)/Matrix_Op.cpp \
	$(ARACNE_PATH)/Microarray_Set.cpp \
	$(ARACNE_PATH)/MutualInfo.cpp \
	$(ARACNE_PATH)/param.cpp

ARACNE_HEADERS = \
	$(ARACNE_PATH)/Matrix.h \
	$(ARACNE_PATH)/Matrix_Op.h \
	$(ARACNE_PATH)/Microarray_Set.h \
	$(ARACNE_PATH)/MutualInfo.h \
	$(ARACNE_PATH)/param.h

ARACNE_OBJS = \
	$(ARACNE_PATH)/Matrix.o \
	$(ARACNE_PATH)/Matrix_Op.o \
	$(ARACNE_PATH)/Microarray_Set.o \
	$(ARACNE_PATH)/MutualInfo.o \
	$(ARACNE_PATH)/param.o


$(TINYXML_OBJS): $(TINYXML_SRCS)
	cd $(CONAN_PATH)/conan/contrib/tinyxml ; make EXTRA_FLAGS="-fPIC" objs

tinyxml: $(TINYXML_OBJS)

$(ARACNE_OBJS): $(ARACNE_SRCS) $(ARACNE_HEADERS)
	cd $(CONAN_PATH)/conan/contrib/aracne ; make EXTRA_FLAGS="$(gcc43_opts) -fPIC" objs

aracne: $(ARACNE_OBJS)

