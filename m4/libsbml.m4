#
# Check for availability of libSBML
#
# Oleksandr Ostrenko <oleksandr.ostrenko@tu-dresden.de>
# Pietro Incardona <incardon@mpi-cbg.de>

AC_DEFUN([AX_LIBSBML],
[
  AC_ARG_WITH(libsbml-libdir,[  --with-libsbml-libdir=LIBDIR Path to libSBML shared library (optional)],
              libsbml_libdir="$withval", libsbml_libdir="")
  AC_ARG_WITH(libsbml-include,[  --with-libsbml-include=INCLUDEDIR Path to libSBML header files (optional)],
              libsbml_include="$withval", libsbml_include="")
  AC_ARG_ENABLE(libsbml-test, [  --disable-libsbml-test Do not try to compile and run a test libSBML program],
          , enable_libsbml_test=yes)
  AC_ARG_ENABLE(libsbml-static, [  --enable-libsbml-static Link to libSBML statically],
          , enable_libsbml_static=no)
  AC_REQUIRE([AC_PROG_CXX])
  AC_ARG_VAR(SBML_LIB,[libSBML library file name])
  AC_ARG_VAR(SBML_LDFLAGS,[Linker arguments for libSBML (e.g., dependecies for static linking)])

  if test "x$SBML_LIB" = "x"; then
    SBML_LIB="sbml"
    if test "x$enable_libsbml_static" = "xyes"; then
      SBML_LIB="sbml-static"
    fi
  fi

  SBML_LDFLAGS="-l$SBML_LIB $SBML_LDFLAGS"

  if test "x$libsbml_libdir" != "x"; then
    SBML_LDFLAGS="-L$libsbml_libdir $SBML_LDFLAGS"
  fi

  if test "x$libsbml_include" != "x"; then
    SBML_CPPFLAGS="-I$libsbml_include $SBML_CPPFLAGS"
  fi

  OLD_LIBS=$LIBS
  LIBS="$LIBS $SBML_LDFLAGS"

  OLD_CPPFLAGS=$CPPFLAGS
  CPPFLAGS="$CPPFLAGS $SBML_CPPFLAGS"

  if test "x$enable_libsbml_test" = "xyes"; then
    AC_REQUIRE([AC_PROG_CC])

    AC_CACHE_CHECK([whether SBML library is available],
            ax_cv_libsbml_available,
            [AC_LANG_PUSH(C++)
        AC_LINK_IFELSE([AC_LANG_PROGRAM([[
        #include <iostream>
        #include <sbml/SBMLTypes.h>
                                                          ]],
                                  [[    LIBSBML_CPP_NAMESPACE_QUALIFIER 
    SBMLDocument      *ptrSBMLDocument = LIBSBML_CPP_NAMESPACE_QUALIFIER
      readSBMLFromString("");

    if (ptrSBMLDocument->getNumErrors() > 0)
    {
      ptrSBMLDocument->printErrors(std::cout);
      return -1;
    }
    return 0;]])],
                          ax_cv_libsbml_available=yes, ax_cv_libsbml_available=no)
          AC_LANG_POP([C++])
    ])
  else
    AC_LANG_PUSH(C++)
    AC_CHECK_LIB($SBML_LIB, readSBMLFromString,
                [ax_cv_libsbml_available="yes"],
                [ax_cv_libsbml_available="no"])
    AC_LANG_POP([C++])
  fi

  LIBS=$OLD_LIBS
  CPPFLAGS=$OLD_CPPFLAGS

  AC_SUBST(SBML_CPPFLAGS)
  AC_SUBST(SBML_LDFLAGS)

  if test "x$ax_cv_libsbml_available" != "xyes"; then
    $2
    :
  else
    ifelse([$1],,[AC_DEFINE(HAVE_LIBSBML,1,[Define if you have the SBML library.])],[$1])
    :
  fi
])
