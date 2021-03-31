LINKLIBS = $(OPT_LIBS) -L$(builddir)/../libpssa/src/.libs -lpssa $(SBML_LDFLAGS)

EXTRA_DIST = \
sbml/Diffusion.sbml \
sbml/Multimerization.sbml

noinst_PROGRAMS = pssa_test
pssa_test_SOURCES = \
main.cpp \
TestBase.cpp \
TestBase.h \
TestDiffusion.cpp \
TestDiffusion.h \
TestReaction.cpp \
TestReaction.h \
TestReactionDiffusion.cpp \
TestReactionDiffusion.h

pssa_test_CFLAGS = -DUNIX -rdynamic -I$(srcdir)/../libpssa/include $(GSL_CFLAGS) $(SBML_CPPFLAGS)
pssa_test_CXXFLAGS = -DUNIX -rdynamic  -I$(srcdir)/../libpssa/include $(GSL_CFLAGS) $(SBML_CPPFLAGS)
pssa_test_LDADD = $(LINKLIBS)
