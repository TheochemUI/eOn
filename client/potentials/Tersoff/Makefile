POT = tersoff.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $(LIB) $*.f90 

FC = xlf90

FFLAGS = -qlanglvl=90pure -qcheck -qflttrap -qsuffix=f=f90 -qsource -qextname -O4

## create tersoff potential in library
edip: ${POT}
	${AR} rcs libTersoff_Fortran.a ${POT}
	ranlib libTersoff_Fortran.a

clean:
	rm -f *.o *.a *.lst
