CXXFLAGS += -O2 -Wall
FC=gfortran
LDFLAGS += -lgfortran

# define NDEBUG to remove debugging lines
ifdef NDEBUG
	NDEBUG=No
	CXXFLAGS += -DNDEBUG
else
	NDEBUG=Yes
endif

# define BOINC to link against boinc libs
ifdef BOINC
	BOINC=Yes	
	NDEBUG=No
	CXXFLAGS += -DNDEBUG -DBOINC
	LDFLAGS += -static -lboinc_api -lboinc -lpthread
else
	BOINC=No
endif

#Potentials
SUBDIRS += ./potentials/EMT/
LIBS += ./potentials/EMT/libEMT.a
POTENTIALS += "+EMT "
SUBDIRS += ./potentials/Morse/
LIBS += ./potentials/Morse/libMorse.a
POTENTIALS += "+MORSE "
SUBDIRS += ./potentials/LennardJones/
LIBS += ./potentials/LennardJones/libLJ.a
POTENTIALS += "+LJ "

ifdef NO_FORTRAN
   CXXFLAGS += -DNO_FORTRAN
   POTENTIALS += "-Aluminum "
else
   SUBDIRS += ./potentials/Aluminum/
   LIBS += ./potentials/Aluminum/libAL.a
   POTENTIALS += "+Aluminum "
endif

OBJECTS := ClientEON.o SaddlePoint.o Dimer.o DimerConstants.o \
	Eigenvalues.o EpiCenters.o Prefactors.o ConjugateGradients.o\
	Constants.o debug.o HelperFunctions.o Matter.o Parameters.o\
	PotentialsInterface.o Potentials.o QMBox.o Quickmin.o SaddlePoint.o

all: $(SUBDIRS) client
	@echo "EON Client Compilation Options"
	@echo "BOINC: $(BOINC)" 
	@echo "DEBUG: $(NDEBUG)"
	@echo "POTENTIALS: $(POTENTIALS)"

client: $(OBJECTS) $(LIBS)
	g++  -o $@ $^ $(LDFLAGS)

$(LIBS):
	$(MAKE) -C $@

$LIBS: $(SUBDIRS)

$(SUBDIRS):
	export FC=$(FC); $(MAKE) -C $@

clean:
	rm -f $(OBJECTS) client

clean-all: clean
	for subdir in $(SUBDIRS) ; do $(MAKE) -C $$subdir clean ; done

.PHONY : all $(SUBDIRS) clean clean-all
