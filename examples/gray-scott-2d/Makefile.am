LINKLIBS = $(OPT_LIBS) $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) \
  -L$(builddir)/../../libpssa/src/.libs -lpssa $(SBML_LDFLAGS) $(GSL_LIBS)

noinst_PROGRAMS = grayscott
grayscott_SOURCES = main.cpp

grayscott_CFLAGS = -DUNIX -rdynamic -I$(srcdir)/../../libpssa/include $(GSL_CFLAGS) $(BOOST_CPPFLAGS) $(SBML_CPPFLAGS)
grayscott_CXXFLAGS = -DUNIX -rdynamic  -I$(srcdir)/../../libpssa/include $(GSL_CFLAGS) $(BOOST_CPPFLAGS) $(SBML_CPPFLAGS)
grayscott_LDADD = $(LINKLIBS)
