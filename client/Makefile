#C compiler settings
CC=$(CROSS)gcc$(GCC_SUFFIX)
CXX=$(CROSS)g++$(GCC_SUFFIX)
LD=$(CROSS)ld
AR=$(CROSS)ar cru
RANLIB=$(CROSS)ranlib

CXXFLAGS += -O2 -Wall -Wfatal-errors

ifndef DEBUG
	NDEBUG=1
endif

ifdef EONMPIBGP
	EONMPI=1
endif

ifdef EONMPI
	CXX=mpicxx
	CXXFLAGS += -DEONMPI
endif

#If we aren't doing a cross compile, guess our platform.
ifndef CROSS
#Ensures all code is statically linked on a Linux machine
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
   LINUX=1
   #MPI often doesn't mix well with static linking...
   ifndef EONMPI
	   LDFLAGS += -static
   endif
endif
ifeq ($(UNAME), Darwin)
   OSX=1
endif
ifeq ($(UNAME), CYGWIN_NT-5.1)
   WIN32=1
endif
endif

ifdef LINUX
	CXXFLAGS += -DLINUX
endif
ifdef OSX
	CXXFLAGS += -DOSX
endif
ifdef WIN32
	CXXFLAGS += -DWIN32
endif

#Fortran compiler settings
ifdef NO_FORTRAN
    CXXFLAGS += -DNO_FORTRAN
else
    # gfortran
    FC = $(CROSS)gfortran
    FFLAGS += -O2 -ffpe-trap=invalid,zero,overflow
    FAR = $(CROSS)ar cru
    ifdef BOINC
	ifdef OSX
	    LDFLAGS += /opt/local/lib/gcc45/libgfortran.a
	else
		LDFLAGS += -lgfortran
	endif
    else
        LDFLAGS += -lgfortran
    endif
    # ifort
    #FC = ifort
    #FFLAGS += -O2 -fpe0
    # Linux linking
    #FAR = ar cru
    #LDFLAGS += -limf -lsvml -lifcore -ldl
    # OSX linking
    #FAR = ifort -staticlib -o
endif

#Fortran compiler settings
ifdef BOPFOX
    CXXFLAGS += -DBOPFOX
    LDFLAGS += libbopfox.a -lgfortran
    POTDIRS += ./potentials/bop
    LIBS += ./potentials/bop/libbop.a
    POTENTIALS += "+bop"
endif

ifdef LAMMPS_POT
#   CXXFLAGS += -I./potentials/LAMMPS/src/STUBS
# CAUSING PROBLEMS IN USING -DLAMMPS MUST BE SOMETHING DIFFERENT!
   CXXFLAGS += -DLAMMPS_POT #-I./potentials/filesFromLAMMPS
   POTDIRS += ./potentials/LAMMPS
   LIBS += ./potentials/LAMMPS/liblammps.a ./potentials/LAMMPS/liblmp_serial.a ./potentials/LAMMPS/libmpi.a ./potentials/LAMMPS/libmeam.a 
   POTENTIALS += "+LAMMPS"
endif

# define BOINC to link against boinc libs
ifdef BOINC
    BOINC=Yes
    DEBUG=No
    OBJECTS += Compression.o
    CXXFLAGS += -DNDEBUG -DBOINC -Iboinc/include
    OBJECTS += Compression.o

    ifdef WIN32
        LDFLAGS += -static -Lboinc/lib_win -lboinc_api -lboinc 
	LDFLAGS += -Larchive/lib_win -larchive -lz 
	CXXFLAGS += -Iarchive/include -DLIBARCHIVE_STATIC
    endif

    ifdef OSX
	CC="/usr/bin/gcc"
	CXX="/usr/bin/g++"
        LDFLAGS += -Lboinc/lib_osx -lboinc_api -lboinc -lpthread -larchive -lz
    endif

    ifdef LINUX
	UNAME := $(shell uname -m)
	ifeq ($(UNAME), x86_64)
        LDFLAGS += -static -Lboinc/lib_linux64 -lboinc_api -lboinc -lpthread -larchive -lz
	else
	LDFLAGS += -static -Lboinc/lib_linux32 -lboinc_api -lboinc -lpthread -larchive -lz
	endif
    endif
else
    BOINC=No
	# define NDEBUG to remove debugging lines
	ifdef NDEBUG
		DEBUG=No
		CXXFLAGS += -DNDEBUG
	else
		CXXFLAGS += -g
		DEBUG=Yes
	endif
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
POTENTIALS += "-EAM"
POTDIRS += ./potentials/QSC/
LIBS += ./potentials/QSC/libQSC.a
POTENTIALS += "+QSC"
POTDIRS += ./potentials/platinum-water/
LIBS += ./potentials/platinum-water/libplatinum-water.a
POTENTIALS += "+Pt-H2O"
POTDIRS += ./potentials/IMD/
LIBS += ./potentials/IMD/libIMD.a
POTENTIALS += "+IMD"

ifndef WIN32
	POTDIRS += ./potentials/VASP
	LIBS += ./potentials/VASP/libVASP.a
	POTENTIALS += "+VASP"
else
	POTENTIALS += "-VASP"
endif

POTDIRS += ./potentials/bopfox
LIBS += ./potentials/bopfox/libbopfox.a
POTENTIALS += "+bopfox"

ifdef EONMPI
	POTENTIALS += "+MPIPot"
	POTDIRS += ./potentials/MPIPot
	LIBS += ./potentials/MPIPot/libMPIPot.a
else
	POTENTIALS += "-MPIPot"
endif


ifdef NO_FORTRAN
    POTENTIALS += "-Aluminum -Lenosky -SW -Tersoff -EDIP "
else
    POTENTIALS += "+Aluminum +Lenosky +SW +Tersoff +EDIP "
    FPOTDIRS += ./potentials/Aluminum/ ./potentials/Lenosky/ ./potentials/SW/ ./potentials/Tersoff/ ./potentials/EDIP/
    LIBS += ./potentials/Aluminum/libAL.a ./potentials/Lenosky/libLenosky.a ./potentials/SW/libSW.a ./potentials/Tersoff/libTersoff.a ./potentials/EDIP/libEDIP.a
endif

OBJECTS += ClientEON.o INIFile.o SaddleSearch.o Dimer.o \
    EpiCenters.o Hessian.o ConjugateGradients.o \
    HelperFunctions.o Matter.o Parameters.o \
    Potential.o Quickmin.o QuickminBox.o\
    ProcessSearchJob.o PointJob.o MinimizationJob.o HessianJob.o \
    ParallelReplicaJob.o Dynamics.o  BondBoost.o  DistributedReplicaJob.o \
	FiniteDifferenceJob.o DimerRotationJob.o DisplacementSamplingJob.o NudgedElasticBandJob.o \
	TestJob.o BasinHoppingJob.o SaddleSearchJob.o ImprovedDimer.o \
	NudgedElasticBand.o Lanczos.o ExactMinMode.o Bundling.o Job.o CommandLine.o

ifdef EONMPIBGP
	CXXFLAGS += -I/bgsys/drivers/ppcfloor/gnu-linux/include/python2.6/
	CXXFLAGS += -DEONMPIBGP
	OBJECTS  += EonBGP.o
	LDFLAGS  += -L/bgsys/drivers/ppcfloor/gnu-linux/lib/ -Wl,-dy -lpython2.6 -lm
endif

all: $(POTDIRS) $(FPOTDIRS) client
	@echo
	@echo "EON Client Compilation Options"
	@echo "BOINC: $(BOINC)" 
	@echo "DEBUG: $(DEBUG)"
	@echo "POTENTIALS: $(POTENTIALS)"

client: $(OBJECTS) $(LIBS)
	$(CXX) -o $@ $^ $(LDFLAGS)

libeon: $(filter-out ClientEON.o,$(OBJECTS)) $(POTDIRS) $(FPOTDIRS)
	$(AR) libeon.a $(filter-out ClientEON.o,$(OBJECTS)) potentials/*/*.o potentials/EMT/Asap/*.o

include Depend.mk

$(LIBS):
	$(MAKE) -C $@

$LIBS: $(POTDIRS) $(FPOTDIRS)

$(POTDIRS):
	$(MAKE) -C $@ CC="$(CC)" CXX="$(CXX)" LD="$(LD)" AR="$(AR)" RANLIB="$(RANLIB)"

$(FPOTDIRS):
	$(MAKE) -C $@ CC="$(CC)" CXX="$(CXX)" LD="$(LD)" AR="$(FAR)" FC="$(FC)" FFLAGS="$(FFLAGS)" RANLIB="$(RANLIB)"

clean:
	rm -f $(OBJECTS) client

clean-all: clean
	for pot in $(POTDIRS) $(FPOTDIRS) ; do $(MAKE) -C $$pot clean ; done

.PHONY : all $(POTDIRS) $(FPOTDIRS) clean clean-all
# DO NOT DELETE
