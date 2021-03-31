LINKLIBS = $(OPT_LIBS) $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) \
  -L$(builddir)/../../libpssa/src/.libs -lpssa $(SBML_LDFLAGS)

noinst_PROGRAMS = etlswitch
etlswitch_SOURCES = main.cpp

etlswitch_CFLAGS = -DUNIX -rdynamic -I$(srcdir)/../../libpssa/include $(BOOST_CPPFLAGS) $(SBML_CPPFLAGS)
etlswitch_CXXFLAGS = -DUNIX -rdynamic  -I$(srcdir)/../../libpssa/include $(BOOST_CPPFLAGS) $(SBML_CPPFLAGS)
etlswitch_LDADD = $(LINKLIBS)
