dnl Process this file with autoconf to produce a configure script.

AC_INIT
AC_CONFIG_SRCDIR([src/cortex_area.c])
AM_INIT_AUTOMAKE(cortex_area, 1.1)
AM_CONFIG_HEADER(config.h)

smr_WITH_BUILD_PATH

AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

AC_PROG_CC
AC_PROG_CXX

AC_PATH_PROGS(PERL, perl5 perl)

mni_REQUIRE_BICPL

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
