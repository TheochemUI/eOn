
OBJECTS = edipFortran.o EDIP.o

#FC ?= gfortran
#FFLAGS ?= -O2
#AR ?= ar cru

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

LIB := libEDIP.a

$(LIB) : $(OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@
clean :
	rm -f $(OBJECTS) $(LIB)
