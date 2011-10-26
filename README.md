## Conan

Conan is a C++ and Python library for the generation, inference and analysis of complex networks.

#### Documentation
Documentation is available inside the `doc/` directory.
Also <a href="http://www.doxys.dk">DoxyS</a> can be used to generate documentation from the source code.

#### Installation
To install it run the following commands inside the main directory:

    ./bootstrap.sh
    cp /usr/share/libtool/config/ltmain.sh . # this may depend on your Unix distribution
    ./bootstrap.sh
    ./configure
    make
    make install

Note that Conan's compilation requires GCC version 4.3 or higher!
Also, Conan requires Autotools, the GNU Scientific Library (GSL) and the Boost libraries
to be properly set up in your computer.
