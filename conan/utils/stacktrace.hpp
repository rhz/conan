#ifndef STACKTRACE_HPP
#define STACKTRACE_HPP
#ifdef __linux__
#include <csignal>
#include <execinfo.h>
#include <cxxabi.h>


#ifdef PRINT_STACKTRACE
void show_stackframe()
{
  void * trace[16];
  char ** messages = (char **)NULL;
  int trace_size = 0;
  std::vector<std::string> msgs;

  trace_size = backtrace( trace, 16 );
  messages = backtrace_symbols( trace, trace_size );
  for (int i = 0; i < trace_size; ++i)
    msgs.push_back( messages[ i ] );

  free( messages );

  // Demangle function names
  BOOST_FOREACH( std::string & s, msgs )
  {
    size_t p1 = s.find( '(' );
    if ( p1 == std::string::npos )
    {
      std::cerr << "[bt] " << s << std::endl;
    }
    else
    {
      size_t p2 = s.find( '+', p1 );
      if ( p2 == std::string::npos )
        std::cerr << "[bt] " << s << std::endl;
      else
      {
        char * realname;
        int status;

        realname = abi::__cxa_demangle( s.substr( p1 + 1, p2 - p1 - 1 ).c_str(), 0, 0, &status);

        if ( realname != NULL )
          std::cerr << "[bt] " << s.substr( 0, p1 + 1) << realname << s.substr( p2 ) << std::endl;

        free( realname );
      }
    }
  } // foreach

  return;
}
#endif // PRINT_STACKTRACE
#endif // __linux__
#endif // STACKTRACE_HPP
