# Edit this file to specify your compiler and LAPACK library

# Fortran compiler
FC = gfortran

# LAPACK and BLAS library

# ATLAS
#LAPACK = -llapack -lcblas -lf77blas -latlas

# Reference LAPACK/BLAS (from netlib.org)
#LAPACK = /usr/local/src/lapack-3.4.2/liblapack.a /usr/local/src/lapack-3.4.2/librefblas.a

# For Mac
LAPACK = -framework accelerate
