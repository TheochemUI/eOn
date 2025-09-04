#------------------------------------
#Potentials
POTDIRS += ./potentials/EMT/
LIBS += ./potentials/EMT/libEMT.a
POTENTIALS += "+EMT "
POTDIRS += ./potentials/Morse/
LIBS += ./potentials/Morse/libMorse.a
POTENTIALS += "+MORSE "
POTDIRS += ./potentials/LJ/
LIBS += ./potentials/LJ/libLJ.a
POTENTIALS += "+LJ"
POTDIRS += ./potentials/LJCluster/
LIBS += ./potentials/LJCluster/libLJCluster.a
POTENTIALS += "+LJ"
POTDIRS += ./potentials/EAM/
LIBS += ./potentials/EAM/libEAM.a
POTENTIALS += "+EAM"
POTDIRS += ./potentials/QSC/
LIBS += ./potentials/QSC/libQSC.a
POTENTIALS += "+QSC"
POTDIRS += ./potentials/Water/
LIBS += ./potentials/Water/libwater.a
POTENTIALS += "+H2O"
POTDIRS += ./potentials/Water_Pt/
LIBS += ./potentials/Water_Pt/libtip4p_pt.a
POTENTIALS += "+H2O_Pt"
POTDIRS += ./potentials/IMD/
LIBS += ./potentials/IMD/libIMD.a
POTENTIALS += "+IMD"
POTDIRS += ./potentials/ExtPot/
LIBS += ./potentials/ExtPot/libextpot.a
POTENTIALS += "+EXT_POT"
POTDIRS += ./potentials/AMS/
LIBS += ./potentials/AMS/libAMS.a
POTENTIALS += "+AMS"
#POTDIRS += ./potentials/PyAMFF/
#LIBS += ./potentials/PyAMFF/libPyAMFF.a
#POTENTIALS += "+PYAMFF"

#Potentials relying on fortran
ifdef NO_FORTRAN
    POTENTIALS += "-Aluminum -Lenosky -SW -Tersoff -EDIP -H2O_H -FeHe"
    OPOTDIRS += ./potentials/Aluminum/ ./potentials/Lenosky/ ./potentials/SW/ \
                ./potentials/Tersoff/ ./potentials/EDIP/ ./potentials/Water_H/ \
                ./potentials/FeHe/
else
    POTENTIALS += "+Aluminum +Lenosky +SW +Tersoff +EDIP +H2O_H +FeHe"
    FPOTDIRS += ./potentials/Aluminum/ ./potentials/Lenosky/ ./potentials/SW/ \
                ./potentials/Tersoff/ ./potentials/EDIP/ ./potentials/Water_H/ \
                ./potentials/FeHe/
    LIBS += ./potentials/Aluminum/libAL.a ./potentials/Lenosky/libLenosky.a \
            ./potentials/SW/libSW.a ./potentials/Tersoff/libTersoff.a \
            ./potentials/EDIP/libEDIP.a ./potentials/Water_H/libtip4p_h.a \
            ./potentials/FeHe/libFeHe.a

endif

#Optional potentials
ifdef PYAMFF_POT
   CXXFLAGS += -DPYAMFF_POT
   POTDIRS += ./potentials/PyAMFF
   LIBS += ./potentials/PyAMFF/libPyAMFF.a ./potentials/PyAMFF/libAMFF.a
   POTENTIALS += "+PyAMFF"
else
   OPOTDIRS += ./potentials/PyAMFF
   POTENTIALS += "-PyAMFF"
endif

ifdef LAMMPS_POT
    CXXFLAGS += -DLAMMPS_POT
    POTDIRS += ./potentials/LAMMPS
    LIBS += ./potentials/LAMMPS/liblammps.a
    ifdef EONMPI
        LIBS += ./potentials/LAMMPS/liblammps_mpi.a
    else
        LIBS += ./potentials/LAMMPS/liblammps_serial.a
    endif
    ifndef EONMPI
        LIBS += ./potentials/LAMMPS/libmpi_stubs.a
    endif
    POTENTIALS += "+LAMMPS"
    ifdef LAMMPS_MEAM
        LIBS += ./potentials/LAMMPS/libmeam.a
        POTENTIALS += "+LAMMPS_MEAM"
    endif
    ifdef LAMMPS_KIM
        LIBS += /home/graeme/.local/kim-v2-mpi-4/lib64/libkim-api-v2.so
        LDFLAGS += -lcurl
        POTENTIALS += "+LAMMPS_KIM"
    endif

else
   OPOTDIRS += ./potentials/LAMMPS
   POTENTIALS += "-LAMMPS"
endif

ifdef ASE_POT
   # first -I for pybind11, second -I for Python.h
   CXXFLAGS += -DASE_POT -I/path/to/python/rootdir/include -I/path/to/python/rootdir/include/python3.xx  # EDIT
   # for libpython3.xx.so
   LDFLAGS += -L/path/to/python/rootdir/lib -lpython3.xx  # EDIT
   POTDIRS += ./potentials/ASE
   LIBS += ./potentials/ASE/libASE.a
   POTENTIALS += "+ASE"
else
   OPOTDIRS += ./potentials/ASE
   POTENTIALS += "-ASE"
endif

ifdef NEW_POT
   CXXFLAGS += -DNEW_POT
   POTDIRS += ./potentials/New
   LIBS += ./potentials/New/libnew.a
   POTENTIALS += "+NEW_POT"
else
   OPOTDIRS += ./potentials/New
   POTENTIALS += "-NEW_POT"
endif

ifndef WIN32
    POTDIRS += ./potentials/VASP
    LIBS += ./potentials/VASP/libVASP.a
    POTENTIALS += "+VASP"
else
    OPOTDIRS += ./potentials/VASP
    POTENTIALS += "-VASP"
endif

ifdef EONMPI
    POTDIRS += ./potentials/MPIPot
    LIBS += ./potentials/MPIPot/libMPIPot.a
    POTENTIALS += "+MPI_POT"
else
    OPOTDIRS += ./potentials/MPIPot
    POTENTIALS += "-MPI_POT"
endif

#------------------------------------
#MPI settings continued
ifdef EONMPI
    #python_include_path=$(shell python -c "import sys,os; print(os.path.join(sys.prefix, 'include', 'python'+sys.version[:3]))")
    python_include_path=$(shell python -c "from sysconfig import get_paths as gp; print(gp()[\"include\"])")
    #python_lib=$(shell python -c "import sys,os; print('python'+sys.version[:3])")
    python_lib=$(shell python3-config --libs)
    #python_lib_path=$(shell python -c "import sys,os;print(os.path.join(sys.prefix, 'lib'))")
    python_lib_path=$(shell python3-config --prefix)
    ## uncomment for comilation on hopper, comment above definition
    # python_lib_path=$(shell python -c "import sys,os;print os.path.join(sys.prefix, 'lib','python'+sys.version[:3], 'config')")
    #ifneq ($(python_lib_path),/usr/lib)
    #    LDFLAGS += -L${python_lib_path} -l${python_lib} -lm
    #else
    #    LDFLAGS += -l${python_lib} -lm
    #endif
	LDFLAGS += -L$(shell python3-config --prefix) $(shell python3-config --libs)
    CXXFLAGS += -I${python_include_path}
    #CXXFLAGS += ${python_include_path}
endif

#------------------------------------
#client source code
OBJECTS += ClientEON.o INIFile.o MinModeSaddleSearch.o Dimer.o EpiCenters.o \
           Hessian.o ConjugateGradients.o HelperFunctions.o Matter.o \
           Parameters.o Potential.o Quickmin.o ProcessSearchJob.o PointJob.o \
           MinimizationJob.o HessianJob.o ParallelReplicaJob.o \
           SafeHyperJob.o TADJob.o\
           ReplicaExchangeJob.o Dynamics.o BondBoost.o FiniteDifferenceJob.o \
           NudgedElasticBandJob.o TestJob.o BasinHoppingJob.o \
           SaddleSearchJob.o ImprovedDimer.o NudgedElasticBand.o Lanczos.o \
           Bundling.o Job.o CommandLine.o DynamicsJob.o Log.o \
           LBFGS.o LowestEigenmode.o Optimizer.o Prefactor.o \
           DynamicsSaddleSearch.o PrefactorJob.o FIRE.o \
           GlobalOptimizationJob.o GlobalOptimization.o StructureComparisonJob.o \
           MonteCarloJob.o MonteCarlo.o SteepestDescent.o BasinHoppingSaddleSearch.o \
           BiasedGradientSquaredDescent.o ExceptionsEON.o


#ifneq ($(or unitTests,check),)
#CXXFLAGS += -std=c++11
TEMPOBJ := $(filter-out ClientEON.o,$(OBJECTS))
DEPOBJECTS := $(addprefix ../,$(TEMPOBJ))
DEPLIBS := $(addprefix ../,$(LIBS))
export
#endif

#------------------------------------
#Build rules
all: $(POTDIRS) $(FPOTDIRS) eonclient
	@echo
	@echo "EON Client Compilation Options"
	@echo "DEBUG: $(DEBUG)"
	@echo "POTENTIALS: $(POTENTIALS)"

eonclient: $(OBJECTS) $(LIBS)
	$(CXX) -o $(TARGET_NAME) $^ $(LDFLAGS)

libeon: $(filter-out ClientEON.o,$(OBJECTS)) $(POTDIRS) $(FPOTDIRS)
	$(AR) libeon.a $(filter-out ClientEON.o,$(OBJECTS)) potentials/*/*.o potentials/EMT/Asap/*.o

ClientEON.o: version.h
CommandLine.o: version.h

version.h:
	./version.sh > version.h

$(LIBS):
	$(MAKE) -C $@

$LIBS: $(POTDIRS) $(FPOTDIRS)

$(POTDIRS):
	$(MAKE) -C $@ CC="$(CC)" CXX="$(CXX)" LD="$(LD)" AR="$(AR)" RANLIB="$(RANLIB)" CXXFLAGS="$(CXXFLAGS)"

$(FPOTDIRS):
	$(MAKE) -C $@ CC="$(CC)" CXX="$(CXX)" LD="$(LD)" AR="$(FAR)" FC="$(FC)" FFLAGS="$(FFLAGS)" RANLIB="$(RANLIB)" CXXFLAGS="$(CXXFLAGS)"

mkUnitTests: $(filter-out ClientEON.o,$(OBJECTS))
	cd unittests && $(MAKE)

testsClean:
	cd unittests && $(MAKE) clean

testsClobber:
	cd unittests && $(MAKE) clobber

unitTests: mkUnitTests
	cd unittests && sh run_unit_tests.sh

check: unitTests
	@echo "In the future, regression testing will automatically run now."

clean: testsClean
	rm -f $(OBJECTS) $(DEPENDS)

clean-all: clean testsClobber
	for pot in $(POTDIRS) $(FPOTDIRS) $(OPOTDIRS); do $(MAKE) -C $$pot clean ; done

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $<

DEPENDS= $(wildcard *.d)
-include $(DEPENDS)

.PHONY : all $(POTDIRS) $(FPOTDIRS) clean clean-all version.h
# DO NOT DELETE
