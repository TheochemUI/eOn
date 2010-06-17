POT = lenosky.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $(LIB) $*.f90 

## FC = xlf90
FC = /opt/ibmcmp/xlf/8.1/bin/xlf90

FFLAGS = -qlanglvl=90pure -qcheck -qflttrap -qsuffix=f=f90 -qsource -qextname -O4

## create edip potential in library
lenosky: ${POT}
	${AR} rcs ../../../build/Release/libLenosky_Fortran.a ${POT}
	ranlib ../../../build/Release/libLenosky_Fortran.a

clean:
	rm -f *.o *.a *.lst
