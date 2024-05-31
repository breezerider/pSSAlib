
#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([pSSAlib],[2.0.0],[incardon@mpi-cbg.de])
AC_USE_SYSTEM_EXTENSIONS
AC_CONFIG_MACRO_DIR([m4])
AC_SUBST(CURRENT,[2])
AC_SUBST(REVISION,[0])
AC_SUBST(AGE,[0])
LT_INIT()
AC_CANONICAL_TARGET
AC_CONFIG_SRCDIR([libpssa/include/])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([libpssa/include/config.h])

AC_PROG_INSTALL


m4_ifdef([AX_CHECK_COMPILE_FLAG],,[m4_include([m4/ax_check_compile_flag.m4])])
m4_ifdef([AX_GCC_ARCHFLAG],,[m4_include([m4/ax_gcc_archflag.m4])])
m4_ifdef([AX_GCC_X86_CPUID],,[m4_include([m4/ax_gcc_x86_cpuid.m4])])
m4_ifdef([AX_BOOST_BASE],,[m4_include([m4/ax_boost_base.m4])])
m4_ifdef([AX_BOOST_PROGRAM_OPTIONS],,[m4_include([m4/ax_boost_program_options.m4])])
m4_ifdef([AX_PATH_GSL],,[m4_include([m4/gsl.m4])])
m4_ifdef([ACX_MPI],,[m4_include([m4/acx_mpi.m4])])
m4_ifdef([AX_LIBSBML],,[m4_include([m4/libsbml.m4])])
m4_ifdef([AX_CXX_COMPILE_STDCXX],,[m4_include([m4/ax_cxx_compile_stdcxx.m4])])

# check for doxygen support

DX_DOXYGEN_FEATURE(ON)
DX_DOT_FEATURE(ON)
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(ON)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)
AC_CONFIG_FILES([doc/pSSAlib.doxygen.conf])
DX_INIT_DOXYGEN(${PACKAGE_NAME}, [doc/pSSAlib.doxygen.conf], [doc/doxygen])
AM_CONDITIONAL(ENABLE_DOXYGEN,[test "x${DX_FLAG_doc}" = x1])

CFLAGS="$CFLAGS"
CXXFLAGS="$CXXFLAGS"

# # #
# Check for args
# # #


# Check for debug compilation

AC_MSG_CHECKING(whether to build with debug support)
debug_build=no
AC_ARG_ENABLE(debug,
  AS_HELP_STRING([--enable-debug@<:@=ARG@:>@],[enable debug mode, can be one of: "info" - add "-g3" to CFLAGS (default), "full" - "info" + enable debug output and disable compiler code optimisation, "engine" - "full" + debug the simulation engine
  ]),
  debug_build="$enableval"
)

AC_MSG_RESULT($debug_build)
case "${debug_build}" in
  engine)
    AC_DEFINE([PSSALIB_ENGINE_CHECK],[],[Debug engine])
    ;;&
  engine | full)
    AC_DEFINE([DEBUG],[],[Debug])
    AC_DEFINE([_DEBUG],[],[Debug])
#    CFLAGS=$(echo "$CFLAGS" | sed "s/-O\?[[0-9]]//g")
#    CXXFLAGS=$(echo "$CXXFLAGS" | sed "s/-O\?[[0-9]]//g")
    CFLAGS="$CFLAGS -O0 "
    CXXFLAGS="$CXXFLAGS -O0 "
    ;;&
  engine | full | yes | info )
#    CFLAGS=$(echo "$CFLAGS" | sed "s/-g\?[[0-9]]//g")
#    CXXFLAGS=$(echo "$CXXFLAGS" | sed "s/-g\?[[0-9]]//g")
    CFLAGS="$CFLAGS -g3"
    CXXFLAGS="$CXXFLAGS -g3"
    LDFLAGS="$LDFLAGS -g3"
    ;;
esac


# Check whether to build CLI

AC_MSG_CHECKING(whether to build the CLI)
build_cli=yes
AC_ARG_ENABLE([cli],
  AS_HELP_STRING([--disable-cli],[build without the command line interface]),
  [build_cli="$enableval"],
  [build_cli="yes"]
)
AC_MSG_RESULT($build_cli)

# Check for MPI

AC_ARG_WITH([mpi],
        AS_HELP_STRING([--with-mpi], [build with multiprocessing support using MPI]),
        [use_mpi="$withval"],
        [use_mpi=no]
        )

# check for libSBML

AC_MSG_CHECKING(whether to build with libSBML)
AC_ARG_WITH([libsbml],
  AS_HELP_STRING([--without-libsbml], [build without SBML support and dependency on libSBML]),
  [use_sbml=no],
  [use_sbml=yes]
)
AC_MSG_RESULT($use_sbml)

AM_CONDITIONAL([LIBSBML], [test x$use_sbml = xyes])

# Target information
AC_CANONICAL_TARGET

# disable optimisation

AC_MSG_CHECKING(whether target optimization is disabled during build)
AC_ARG_ENABLE([optimisation],
        AS_HELP_STRING([--disable-optimisation], [disable target optimisation during build]),
        [no_opts=yes],
        [no_opts=no]
        )
AC_MSG_RESULT($no_opts)

# check for headers

AC_PROG_CXX
AC_LANG_PUSH([C++])
AX_CXX_COMPILE_STDCXX(11, ext, optional)
if test "x$use_mpi" = "xyes"; then
  ACX_MPI([AC_DEFINE([HAVE_MPI],[],[Use MPI])],[AC_MSG_ERROR([MPI headers not found])])
fi
AC_CHECK_HEADER([fstream],,[AC_MSG_ERROR([*** This package needs file I/O streams header])])
AC_CHECK_HEADER([iostream],,[AC_MSG_ERROR([*** This package needs standard I/O streams header])])
AC_CHECK_HEADER([sstream],,[AC_MSG_ERROR([*** This package needs standard string streams header])])
AC_CHECK_HEADER([iomanip],,[AC_MSG_ERROR([*** This package needs standard console output manipulation header])])
AC_CHECK_HEADER([stdexcept],,[AC_MSG_ERROR([*** This package needs standard exception header])])
AC_CHECK_HEADER([cerrno],,[AC_MSG_ERROR([*** This package needs C-style errno macro header])])
AC_CHECK_HEADER([string.h],,[AC_MSG_ERROR([*** This package needs STL string header])])
AC_CHECK_HEADER([vector],,[AC_MSG_ERROR([*** This package needs STL vector header])])
AC_CHECK_HEADER([algorithm],,[AC_MSG_ERROR([*** This package needs STL algortihms header])])
AC_CHECK_HEADER([limits],,[AC_MSG_ERROR([*** This package needs maths limits header])])
# platform-specific
case "$target_os" in
  linux*)
    AC_CHECK_HEADER([time.h],,[AC_MSG_ERROR([*** This package needs time header])])
    AC_CHECK_HEADER([sys/stat.h],,[AC_MSG_ERROR([*** This package needs file-system definition header])])
    AC_CHECK_HEADER([sys/types.h],,[AC_MSG_ERROR([*** This package needs common type definitions header])])
  ;;
  macos*|darwin*)
    AC_CHECK_HEADER([mach/mach.h],,[AC_MSG_ERROR([*** This package needs MacOS common header])])
    AC_CHECK_HEADER([mach/mach_time.h],,[AC_MSG_ERROR([*** This package needs MacOS time definitions header])])
    AC_CHECK_HEADER([sys/stat.h],,[AC_MSG_ERROR([*** This package needs file-system definition header])])
    AC_CHECK_HEADER([sys/types.h],,[AC_MSG_ERROR([*** This package needs common type definitions header])])
  ;;
  *mingw*)
    AC_CHECK_HEADER([windows.h],,[AC_MSG_ERROR([*** This package needs Windows common header])])
    AC_CHECK_HEADER([tchar.h],,[AC_MSG_ERROR([*** This package needs generic text mappings header])])
    AC_CHECK_HEADER([direct.h],,[AC_MSG_ERROR([*** This package needs directory definitions header])])
  ;;
esac
AC_LANG_POP([C++])

# Checks for libraries.

if test "x$use_sbml" = "xyes"; then
  AX_LIBSBML([AC_DEFINE([HAVE_LIBSBML],[],[Include libSBML support])],[AC_MSG_ERROR([Could not link against libSBML])])
fi

PKG_CONFIG="pkg-config --static"
AX_PATH_GSL("1.15",
                [],
                [AC_MSG_ERROR("Could not find the GSL development package, use the following command for Red Hat, Fedora or CentOS: "yum install gsl-devel"; for Debian, Ubuntu: "apt-get install gsl-dev", in order to solve the problem")],[],[])
# Make GSL link statically by default
#GSL_LIBS="-Wl,-Bstatic $GSL_LIBS -Wl,-Bdynamic"
PKG_CONFIG="pkg-config"

AX_BOOST_BASE([1.56],, AC_MSG_ERROR([*** Boost library not found. Please install it.]))
AM_CONDITIONAL([BUILD_CLI], [test x"$build_cli" = x"yes"])
if test "x$build_cli" = "xyes"; then
  AX_BOOST_PROGRAM_OPTIONS
fi

# Checks g++ flags

if test x"$is_custom_arch" = x"no"; then 
  case $target_cpu in
  i[[3456]]86*) CXXFLAGS="$CXXFLAGS"
                no_64=yes;;
  x86_64*) CXXFLAGS="$CXXFLAGS -m64"
                   no_64=no;;
  powerpc64*) CXXFLAGS="$CXXFLAGS"
                     no_64=no;;
  powerpc*) CXXFLAGS="$CXXFLAGS"
            LDFLAGS="$LDFLAGS  -qnostaticlink -qnostaticlink=libgcc"
                     use_bgp=yes
                     no_64=yes;;
  *) AC_MSG_ERROR([*** $target_cpu is not a supported target cpu])
  esac
else
  no_64=yes
fi

# Check target architetture

is_custom_arch=no
if test x"$is_custom_arch" = x"no"; then
    AX_GCC_ARCHFLAG([], [CXXFLAGS="$CXXFLAGS"], [])
else
    ax_cv_gcc_archflag="$custom_arch"
fi

no_avx=no
no_sse42=no
no_sse41=no
no_sse3=no
no_sse2=no
no_sse=no
no_mmx=no

if test x"$no_opts" = x"no"; then
  case "${target}" in
    *-*-darwin*)
      AC_MSG_CHECKING([for an avx processor])
      if [ test `sysctl -n hw.optional.avx1_0` == 1]
      then
      AC_MSG_RESULT([yes])
      no_avx=no
      else
      AC_MSG_RESULT([no])
      no_avx=yes
      fi

      AC_MSG_CHECKING([for a sse4_2 processor])
      if [ test `sysctl -n hw.optional.sse4_2` == 1]
      then
        AC_MSG_RESULT([yes])
        no_sse42=no
      else
        AC_MSG_RESULT([no])
        no_sse42=yes
      fi

      AC_MSG_CHECKING([for a sse4_1 processor])
      if [ test `sysctl -n hw.optional.sse4_1` == 1]
      then
        AC_MSG_RESULT([yes])
        no_sse41=no
      else
        AC_MSG_RESULT([no])
        no_sse41=yes
      fi

      AC_MSG_CHECKING([for a sse3 processor])
      if [ test `sysctl -n hw.optional.sse3` == 1]
      then
        AC_MSG_RESULT([yes])
        no_sse3=no
      else
        AC_MSG_RESULT([no])
        no_sse3=yes
      fi

      AC_MSG_CHECKING([for a sse2 processor])
      if [ test `sysctl -n hw.optional.sse2` == 1]
      then
        AC_MSG_RESULT([yes])
        no_sse2=no
      else
        AC_MSG_RESULT([no])
        no_sse2=yes
      fi

      AC_MSG_CHECKING([for a sse processor])
      if [ test `sysctl -n hw.optional.sse` == 1]
      then
        AC_MSG_RESULT([yes])
        no_sse=no
      else
        AC_MSG_RESULT([no])
        no_sse=yes
      fi
    ;;
    *)
      if grep -q "avx" /proc/cpuinfo
      then
        AC_MSG_RESULT([yes])
        no_avx=no
      else
        AC_MSG_RESULT([no])
        no_avx=yes
      fi

      AC_MSG_CHECKING([for a sse4_2 processor])
      if grep -q "sse4_2" /proc/cpuinfo
      then
        AC_MSG_RESULT([yes])
        no_sse42=no
      else
        AC_MSG_RESULT([no])
        no_sse42=yes
      fi

      AC_MSG_CHECKING([for a sse4_1 processor])
      if grep -q "sse4_1" /proc/cpuinfo
      then
        AC_MSG_RESULT([yes])
        no_sse41=no
      else
        AC_MSG_RESULT([no])
        no_sse41=yes
      fi

      AC_MSG_CHECKING([for a sse3 processor])
      if grep -q "sse3" /proc/cpuinfo
      then
        AC_MSG_RESULT([yes])
        no_sse3=no
      else
        AC_MSG_RESULT([no])
        no_sse3=yes
      fi

      AC_MSG_CHECKING([for a sse2 processor])
      if grep -q "sse2" /proc/cpuinfo
      then
        AC_MSG_RESULT([yes])
        no_sse2=no
      else
        AC_MSG_RESULT([no])
        no_sse2=yes
      fi

      AC_MSG_CHECKING([for a sse processor])
      if grep -q "sse" /proc/cpuinfo
      then
        AC_MSG_RESULT([yes])
        no_sse=no
      else
        AC_MSG_RESULT([no])
        no_sse=yes
      fi
    ;;
  esac

  if test x"$no_avx" = x"no"; then
      AX_CHECK_COMPILE_FLAG([-mavx],[CXXFLAGS="$CXXFLAGS -mavx"],[])
  fi

  if test x"$no_sse42" = x"no"; then
      AX_CHECK_COMPILE_FLAG([-msse4.2],[CXXFLAGS="$CXXFLAGS -msse4.2"],[])
  fi

  if test x"$no_avx" = x"no"; then
      AX_CHECK_COMPILE_FLAG([-msse4.1],[CXXFLAGS="$CXXFLAGS -msse4.1"],[])
  fi

  if test x"$no_avx" = x"no"; then
      AX_CHECK_COMPILE_FLAG([-msse3],[CXXFLAGS="$CXXFLAGS -msse3"],[])
  fi

  if test x"$no_avx" = x"no"; then
      AX_CHECK_COMPILE_FLAG([-msse2],[CXXFLAGS="$CXXFLAGS -msse2"],[])
  fi

  if test x"$no_avx" = x"no"; then
      AX_CHECK_COMPILE_FLAG([-msse],[CXXFLAGS="$CXXFLAGS -msse"],[])
  fi

  if test x"$no_64" = x"no"; then
    AX_CHECK_COMPILE_FLAG([-m64],[CXXFLAGS="$CXXFLAGS -m64"],[no_64=yes])
  fi

  AX_CHECK_COMPILE_FLAG([-mmmx],[CXXFLAGS="$CXXFLAGS -mmmx"],[no_mmx=yes])
#AX_CHECK_COMPILE_FLAG([-Wno-unused-but-set-variable],[CXXFLAGS="$CXXFLAGS -Wno-unused-but-set-variable"],[])

#AX_CHECK_COMPILE_FLAG([-frounding-math],[CXXFLAGS="$CXXFLAGS -frounding-math"],[])

# Define avx or sse

  if test x"$no_avx" = x"no"; then
    AC_DEFINE([HAVE_AVX],[],[Enable avx])
  elif test x"$no_sse41" = x"no"; then
    AC_DEFINE([HAVE_SSE],[],[Enable sse])
  fi

fi

# Checks for typedefs, structures, and compiler characteristics.

# Checks for programs.
if [ test x"$use_mpi" = x"yes" ]; then
  CXX="$MPICXX"
fi

# Checks for library functions.

case "${target_os}" in
  darwin*)
    # Nothing to check
    ;;
  *)
    AC_CHECK_LIB([rt],[clock_gettime],[],AC_MSG_ERROR("Install librt"),[])
    OPT_LIBS="$OPT_LIBS -lrt"
    ;;
esac

AC_CONFIG_FILES([Makefile doc/Makefile libpssa/Makefile libpssa/src/Makefile libpssa/include/Makefile pssa_cli/Makefile pssa_test/Makefile examples/Makefile examples/validation/Makefile examples/benchmarks/Makefile examples/etl-endosome-switch-1d/Makefile examples/gray-scott-2d/Makefile])
AC_OUTPUT

echo ""
echo " Compiler flags: "
echo ""
echo " CXX = \"$CXX\" "
echo " CFLAGS = \"$CFLAGS\" "
echo " CXXFLAGS = \"$CXXFLAGS\" "
echo " LDFLAGS = \"$LDFLAGS\" "
echo ""
echo " Software options: "
echo ""
echo   "***********************************"
echo   "*                                 *"
if [ test x"$pic_mode" = x"yes" ]; then
  echo "*    build PIC: yes               *"
else
  echo "*    build PIC: no                *"
fi
if [ test x"$use_mpi" = x"yes" ]; then
  echo "*    use MPI: yes                 *"
else
  echo "*    use MPI: no                  *"
fi
if [ test x"$use_sbml" = x"yes" ]; then
  echo "*    use SBML: yes                *"
else
  echo "*    use SBML: no                 *"
fi
if [ test x"$build_cli" = x"yes" ]; then
  echo "*    build CLI: yes               *"
else
  echo "*    build CLI: no                *"
fi
if [ test x"$debug_build" = x"yes" ]; then
  echo "*    debug build: yes             *"
elif [ test x"$debug_build" = x"info" ]; then
  echo "*    debug build: info            *"
elif [ test x"$debug_build" = x"full" ]; then
  echo "*    debug build: full            *"
elif [ test x"$debug_build" = x"engine" ]; then
  echo "*    debug build: engine          *"
else
  echo "*    debug build: no              *"
fi
echo   "*                                 *"
echo   "***********************************"

echo   ""
echo   "In order to compile:"
echo   "make"
echo   ""
echo   "Getting started documentation: doc/tutorial/tutorial.hml"
echo   ""
