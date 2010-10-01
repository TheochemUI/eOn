#C compiler settings
CXXFLAGS += -O2 -Wall

ifdef COMPRESSION
	CXXFLAGS += -DCOMPRESSION
	LDFLAGS += -larchive
	OBJECTS += Compression.o
endif

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
	CXXFLAGS += -g
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
ifdef WITH_LANCZOS
	CXXFLAGS += -DWITH_LANCZOS
	FPOTDIRS += Lanczos
	LIBS += Lanczos/libLanczos.a
    LDFLAGS += -lgsl -lgslcblas
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
POTDIRS += ./potentials/EAM/
LIBS += ./potentials/EAM/libEAM.a
POTENTIALS += "+EAM"
POTDIRS += ./potentials/QSC/
LIBS += ./potentials/QSC/libQSC.a
POTENTIALS += "+QSC"
POTDIRS += ./potentials/platinum-water/
LIBS += ./potentials/platinum-water/libplatinum-water.a
POTENTIALS += "+Pt-H2O"


ifdef NO_FORTRAN
	POTENTIALS += "-Aluminum -Lenosky -SW -Tersoff -EDIP "
else
	POTENTIALS += "+Aluminum +Lenosky +SW +Tersoff +EDIP "
	FPOTDIRS += ./potentials/Aluminum/ ./potentials/Lenosky/ ./potentials/SW/ ./potentials/Tersoff/ ./potentials/EDIP/
	LIBS += ./potentials/Aluminum/libAL.a ./potentials/Lenosky/libLenosky.a ./potentials/SW/libSW.a ./potentials/Tersoff/libTersoff.a ./potentials/EDIP/libEDIP.a
endif

OBJECTS += ClientEON.o INIFile.o SaddlePoint.o Dimer.o\
	Eigenvalues.o EpiCenters.o Prefactors.o ConjugateGradients.o\
	debug.o HelperFunctions.o Matter.o Parameters.o\
	PotentialsInterface.o Potentials.o QMBox.o Quickmin.o SaddlePoint.o \
	ProcessSearch.o Minimization.o

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
