LINKLIBS = $(OPT_LIBS) $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) \
  -L$(builddir)/../../libpssa/src/.libs -lpssa $(SBML_LDFLAGS) $(GSL_LIBS)

noinst_PROGRAMS = validation
validation_SOURCES = main.cpp Heteroreaction.hpp  Homoreaction.hpp  Hreaction.hpp

validation_CFLAGS = -DUNIX -rdynamic -I$(srcdir)/../../libpssa/include $(GSL_CFLAGS) $(SBML_CPPFLAGS)
validation_CXXFLAGS = -DUNIX -rdynamic  -I$(srcdir)/../../libpssa/include $(GSL_CFLAGS) $(SBML_CPPFLAGS)
validation_LDADD = $(LINKLIBS)
