AC_PREREQ(2.59)
AC_INIT([pfit], [0.2], [cngilbreth@gmail.com])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AC_CONFIG_MACRO_DIR([m4])
AC_PROG_FC

dnl Check for BLAS libraries                                                                              
sinclude(ax_blas.m4)
AX_BLAS
if test $ax_blas_ok = no; then
   AC_MSG_ERROR([Cannot find BLAS libraries])
fi

dnl Check for LAPACK libraries
sinclude(ax_lapack.m4)
AX_LAPACK
if test $ax_lapack_ok = no; then
   AC_MSG_ERROR([Cannot find LAPACK libraries])
fi

AC_FC_SRCEXT(F90)
AC_FC_LINE_LENGTH([unlimited])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
