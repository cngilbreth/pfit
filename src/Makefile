include ../make.inc

pfit: pfit.o options.o main.o
	$(FC) -o $@ $+ $(LAPACK)
	cp pfit ../bin/

clean:
	rm -f *.o *.mod

%.o: %.f90
	$(FC) $(DEFS) -c $< -o $@

