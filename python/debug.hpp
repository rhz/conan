#ifndef PYTHON_DEBUG
#define PYTHON_DEBUG
#include <conan/config.hpp>
#ifdef HAVE_ATLAS
  #include <atlas/atlas_buildinfo.h>
#endif
#include <gsl/gsl_version.h>
#include <boost/python.hpp>

#define NL "\n"
#define DATE_INFO "Compilation date: " __DATE__ "."

#ifdef UNAME
  #define SYS_INFO "Conan compiled at " HOSTNAME ", " UNAME "." NL
#else
  #define SYS_INFO
#endif

#if (defined ATL_VERS && defined ATL_ARCH && defined GSL_VERSION)
  #define ATLAS_AND_GSL_INFO NL "Compiled against ATLAS version " ATL_VERS " (optimized for architecture " ATL_ARCH ") and GSL version " GSL_VERSION "."
#else
  #define ATLAS_AND_GSL_INFO
#endif

#ifdef __GNUC__
  #define GCC_VERSION_INFO "GCC version " __VERSION__

#ifdef __OPTIMIZE__
  #define OPT_FLAGS_INFO ", optimization flags enabled"
#else
  #define OPT_FLAGS_INFO
#endif

#ifdef __GXX_EXPERIMENTAL_CXX0X__
  #define CXX0X_INFO ", experimental C++0x standard support enabled"
#else
  #define CXX0X_INFO
#endif

  #define COMPILATION_INFO SYS_INFO GCC_VERSION_INFO OPT_FLAGS_INFO CXX0X_INFO "." ATLAS_AND_GSL_INFO NL DATE_INFO

#else

  #define COMPILATION_INFO SYS_INFO "Unknown compiler used." ATLAS_AND_GSL_INFO NL DATE_INFO

#endif

str version()
{
#ifdef CONAN_VERSION
  return CONAN_VERSION;
#else
  return __DATE__;
#endif
}

str compilation_info()
{
  return COMPILATION_INFO;
}

#endif // PYTHON_DEBUG
