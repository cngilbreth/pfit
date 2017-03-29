# Fortran compiler
FC = gfortran -cpp -fbounds-check -g

# LAPACK and BLAS library
# ATLAS
LAPACK = -llapack -lcblas -lf77blas -latlas
# Reference LAPACK/BLAS (from netlib.org)
#LAPACK = /usr/local/src/lapack-3.4.2/liblapack.a /usr/local/src/lapack-3.4.2/librefblas.a
#LAPACK = -framework accelerate


pfit: pfit.o options.o main.o
	$(FC) -o $@ $+ $(LAPACK)
	cp pfit ../bin/

clean:
	rm -f *.o *.mod

%.o: %.f90
	$(FC) $(DEFS) -c $< -o $@

