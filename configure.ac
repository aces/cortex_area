dnl Process this file with autoconf to produce a configure script.

# Require autoconf 2.57 or newer.
AC_PREREQ([2.57])

# The arguments are package name, and package version.
AC_INIT([cortex_area],[1.3.0],[Claude Lepage <claude@bic.mni.mcgill.ca>])
AC_CONFIG_SRCDIR([src/cortex_area.c])

AM_INIT_AUTOMAKE

# The argument is the name of the generated header file.
# It is recommended to leave it as "config.h".
AC_CONFIG_HEADER([config.h])

smr_WITH_BUILD_PATH

AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

AC_PROG_CC
AC_PROG_CXX

AC_PATH_PROGS(PERL, perl5 perl)

mni_REQUIRE_BICPL

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT(epm-header)

