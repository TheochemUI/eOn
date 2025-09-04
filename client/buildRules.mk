#------------------------------------
#Potentials
POTDIRS += ./potentials/GPRPotential/
LIBS += ./potentials/GPRPotential/libGPRPot.a
POTENTIALS += "+GPRPot "
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
#POTDIRS += ./potentials/EAM/
#LIBS += ./potentials/EAM/libEAM.a
POTENTIALS += "-EAM"
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
# --- Build
SRC := ./
TEMPDIR := ./.build
FILES :=   ClientEON.cpp INIFile.cpp MinModeSaddleSearch.cpp Dimer.cpp EpiCenters.cpp \
           Hessian.cpp ConjugateGradients.cpp HelperFunctions.cpp Matter.cpp \
           Parameters.cpp Potential.cpp Quickmin.cpp ProcessSearchJob.cpp PointJob.cpp \
           MinimizationJob.cpp HessianJob.cpp ParallelReplicaJob.cpp \
           SafeHyperJob.cpp TADJob.cpp\
           ReplicaExchangeJob.cpp Dynamics.cpp BondBoost.cpp FiniteDifferenceJob.cpp \
           NudgedElasticBandJob.cpp TestJob.cpp BasinHoppingJob.cpp \
           SaddleSearchJob.cpp ImprovedDimer.cpp NudgedElasticBand.cpp Lanczos.cpp \
           Bundling.cpp Job.cpp CommandLine.cpp DynamicsJob.cpp Log.cpp \
           LBFGS.cpp LowestEigenmode.cpp Optimizer.cpp Prefactor.cpp \
           DynamicsSaddleSearch.cpp PrefactorJob.cpp FIRE.cpp \
           GlobalOptimizationJob.cpp GlobalOptimization.cpp StructureComparisonJob.cpp \
           MonteCarloJob.cpp MonteCarlo.cpp SteepestDescent.cpp BasinHoppingSaddleSearch.cpp \
           BiasedGradientSquaredDescent.cpp AtomicGPDimer.cpp \
		   gprdimer/gpr/AtomicDimer.cpp \
		   gprdimer/gpr/AtomicDimerInit.cpp \
		   gprdimer/backend/DistributionFunctions.cpp \
	       gprdimer/gpr/auxiliary/ProblemSetUp.cpp \
	       gprdimer/gpr/auxiliary/Gradient.cpp \
	       gprdimer/gpr/auxiliary/Distance.cpp \
		   gprdimer/gpr/prior/PriorGaussian.cpp \
		   gprdimer/gpr/prior/PriorLogUnif.cpp \
		   gprdimer/gpr/prior/PriorSqrtt.cpp \
		   gprdimer/gpr/prior/PriorT.cpp \
		   gprdimer/gpr/dimer/LBFGS.cpp \
	       gprdimer/gpr/covariance_functions/ConstantCF.cpp \
	       gprdimer/gpr/covariance_functions/SexpatCF.cpp \
	       gprdimer/gpr/observation_models/LikGaussian.cpp \
	       gprdimer/gpr/ml/GaussianProcessRegression.cpp \
	       gprdimer/gpr/ml/SCG.cpp

SOURCES := $(patsubst %, $(SRC)/%, $(FILES))
OBJECTS := $(patsubst %.cpp, $(TEMPDIR)/%.o, $(FILES))
DEPFLAGS = -MT $(TEMPDIR)/$@ -MMD -MP -MF $(TEMPDIR)/$*.d
DEPFILES := $(SRCS:%.c=$(TEMPDIR)/%.d)
$(DEPFILES):
include $(wildcard $(DEPFILES))

$(TEMPDIR)/%.o: $(SRC)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

$(TARGET): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $(TARGET)

TEMPOBJ := $(filter-out ClientEON.o,$(OBJECTS))
DEPOBJECTS := $(addprefix ../,$(TEMPOBJ))
DEPLIBS := $(addprefix ../,$(LIBS))
export

all: $(OBJECTS) $(POTDIRS) $(FPOTDIRS) eonclient
	@echo
	@echo "EON Client Compilation Options"
	@echo "DEBUG: $(DEBUG)"
	@echo "POTENTIALS: $(POTENTIALS)"

eonclient: $(TARGET) $(LIBS)
	$(CXX) -o $(TARGET_NAME) $(OBJECTS) $^ $(LDFLAGS)

libeon: $(filter-out ClientEON.o,$(OBJECTS)) $(POTDIRS) $(FPOTDIRS)
	$(AR) libeon.a $(filter-out ClientEON.o,$(OBJECTS)) potentials/*/*.o potentials/EMT/Asap/*.o

ClientEON.cpp: version.h
CommandLine.cpp: version.h

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
	rm -rf $(TEMPDIR) version.h eonclient

clean-all: clean testsClobber
	for pot in $(POTDIRS) $(FPOTDIRS) $(OPOTDIRS); do $(MAKE) -C $$pot clean ; done

.PHONY : all $(POTDIRS) $(FPOTDIRS) clean clean-all version.h
