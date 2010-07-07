#C compiler settings
CXXFLAGS += -O2 -Wall

#Fortran compiler settings
ifdef NO_FORTRAN
	CXXFLAGS += -DNO_FORTRAN
else
# gfortran
FC = gfortran
FFLAGS += -O2
AR = ar cru
LDFLAGS += -lgfortran
# ifort
#FC = ifort
#FFLAGS += -O2
# Linux linking
#AR = ar cru
#LDFLAGS += -limf -lsvml -lifcore -ldl
# OSX linking
#AR = ifort -staticlib -o
endif

# define NDEBUG to remove debugging lines
ifdef NDEBUG
	DEBUG=No
	CXXFLAGS += -DNDEBUG
else
	DEBUG=Yes
endif

# define BOINC to link against boinc libs
ifdef BOINC
	BOINC=Yes	
	DEBUG=No
	CXXFLAGS += -DNDEBUG -DBOINC
	LDFLAGS += -static -lboinc_api -lboinc -lpthread
else
	BOINC=No
endif

#Potentials
POTDIRS += ./potentials/EMT/
LIBS += ./potentials/EMT/libEMT.a
POTENTIALS += "+EMT "
POTDIRS += ./potentials/Morse/
LIBS += ./potentials/Morse/libMorse.a
POTENTIALS += "+MORSE "
POTDIRS += ./potentials/LennardJones/
LIBS += ./potentials/LennardJones/libLJ.a
POTENTIALS += "+LJ "
#POTDIRS += ./potentials/EAM/
#LIBS += ./potentials/EAM/libEAM.a
#POTENTIALS += "+EAM"

ifdef NO_FORTRAN
	POTENTIALS += "-Aluminum "
else
	FPOTDIRS += ./potentials/Aluminum/
	LIBS += ./potentials/Aluminum/libAL.a
	POTENTIALS += "+Aluminum "
endif

OBJECTS := ClientEON.o SaddlePoint.o Dimer.o\
	Eigenvalues.o EpiCenters.o Prefactors.o ConjugateGradients.o\
	Constants.o debug.o HelperFunctions.o Matter.o Parameters.o\
	PotentialsInterface.o Potentials.o QMBox.o Quickmin.o SaddlePoint.o

all: $(POTDIRS) $(FPOTDIRS) client
	@echo
	@echo "EON Client Compilation Options"
	@echo "BOINC: $(BOINC)" 
	@echo "DEBUG: $(DEBUG)"
	@echo "POTENTIALS: $(POTENTIALS)"

client: $(OBJECTS) $(LIBS)
	g++  -o $@ $^ $(LDFLAGS)

$(LIBS):
	$(MAKE) -C $@

$LIBS: $(POTDIRS) $(FPOTDIRS)

$(POTDIRS):
	$(MAKE) -C $@

$(FPOTDIRS):
	@echo $(AR)
	export FC=$(FC); export FFLAGS=$(FFLAGS); export AR="$(AR)"; $(MAKE) -C $@

clean:
	rm -f $(OBJECTS) client

clean-all: clean
	for pot in $(POTDIRS) $(FPOTDIRS) ; do $(MAKE) -C $$pot clean ; done

.PHONY : all $(POTDIRS) $(FPOTDIRS) clean clean-all
