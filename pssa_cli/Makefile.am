LINKLIBS = $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) \
  -L$(builddir)/../libpssa/src/.libs -lpssa $(SBML_LDFLAGS)

bin_PROGRAMS = simulator analyzer
simulator_SOURCES = simulator.cpp CmdLineOptions.hpp
analyzer_SOURCES = analyzer.cpp CmdLineOptions.hpp

AM_CFLAGS = -DUNIX -rdynamic -I$(srcdir)/../libpssa/include/ $(BOOST_CPPFLAGS) $(SBML_CPPFLAGS)
AM_CXXFLAGS = -DUNIX -rdynamic -I$(srcdir)/../libpssa/include/ $(BOOST_CPPFLAGS) $(SBML_CPPFLAGS)
LDADD = $(LINKLIBS)
