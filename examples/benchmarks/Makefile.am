LINKLIBS = $(OPT_LIBS) $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) \
  -L$(builddir)/../../libpssa/src/.libs -lpssa $(SBML_LDFLAGS)

noinst_PROGRAMS = benchmarks
benchmarks_SOURCES = main.cpp ColloidalAggregation.hpp  CyclicLinearChain.hpp

benchmarks_CFLAGS = -DUNIX -rdynamic -I$(srcdir)/../../libpssa/include $(GSL_CFLAGS) $(SBML_CPPFLAGS)
benchmarks_CXXFLAGS = -DUNIX -rdynamic  -I$(srcdir)/../../libpssa/include $(GSL_CFLAGS) $(SBML_CPPFLAGS)
benchmarks_LDADD = $(LINKLIBS)
