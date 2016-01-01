#------------------------------------
#Potentials
POTDIRS += ./potentials/EMT/
LIBS += ./potentials/EMT/libEMT.a
POTENTIALS += "+EMT "
POTDIRS += ./potentials/Morse/
LIBS += ./potentials/Morse/libMorse.a
POTENTIALS += "+MORSE "
POTDIRS += ./potentials/LennardJones/
LIBS += ./potentials/LennardJones/libLJ.a
POTENTIALS += "+LJ"
POTDIRS += ./potentials/LennardJonesCluster/
LIBS += ./potentials/LennardJonesCluster/libLJCluster.a
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
POTDIRS += ./potentials/bopfox/
LIBS += ./potentials/bopfox/libbopfox.a
POTENTIALS += "+bopfox"
POTDIRS += ./potentials/TerminalPotential/
LIBS += ./potentials/TerminalPotential/libterminalpotential.a
POTENTIALS += "+TerminalPotential"


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
ifdef BOPFOX
    CXXFLAGS += -DBOPFOX
    POTDIRS += ./potentials/bop
    LIBS += ./potentials/bop/libbop.a
    ifdef WIN32
        LIBS += bopfox/libbopfox_win32.a
    endif
    ifdef LINUX32
        LIBS += bopfox/libbopfox_lnx32.a
    endif
    ifdef LINUX64
        LIBS += bopfox/libbopfox_lnx64.a
    endif
    ifdef OSX
        LIBS += bopfox/libbopfox_osx.a
    endif
    POTENTIALS += "+bop"
else
    OPOTDIRS += ./potentials/bop
    POTENTIALS += "-bop"
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
else
   OPOTDIRS += ./potentials/LAMMPS
   POTENTIALS += "-LAMMPS"
endif

ifdef NEW_POT
   CXXFLAGS += -DNEW_POT
   POTDIRS += ./potentials/NewPotential
   LIBS += ./potentials/NewPotential/libnewpotential.a
   POTENTIALS += "+NewPotential"
else
   OPOTDIRS += ./potentials/NewPotential
   POTENTIALS += "-NewPotential"
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
    POTENTIALS += "+MPIPot"
    POTDIRS += ./potentials/MPIPot
    LIBS += ./potentials/MPIPot/libMPIPot.a
else
    POTENTIALS += "-MPIPot"
    OPOTDIRS += ./potentials/MPIPot
endif

#------------------------------------
#MPI settings continued
ifdef EONMPI
    ifdef EONMPIBGP
        CXXFLAGS += -I/bgsys/drivers/ppcfloor/gnu-linux/include/python2.7/
        CXXFLAGS += -DEONMPIBGP
        LDFLAGS  += -L/bgsys/drivers/ppcfloor/gnu-linux/lib/ -Wl,-dy -lpython2.7 -lm
    else
        python_include_path=$(shell python -c "import sys,os; print os.path.join(sys.prefix, 'include', 'python'+sys.version[:3])")
        python_lib=$(shell python -c "import sys,os; print 'python'+sys.version[:3]")
python_lib_path=$(shell python -c "import sys,os;print os.path.join(sys.prefix, 'lib')")       
	## uncomment for comilation on hopper, comment above definition	
	# python_lib_path=$(shell python -c "import sys,os;print os.path.join(sys.prefix, 'lib','python'+sys.version[:3], 'config')")
        ifneq ($(python_lib_path),/usr/lib)
        LDFLAGS += -L${python_lib_path} -l${python_lib} -lm
        else
            LDFLAGS += -l${python_lib} -lm
        endif
        CXXFLAGS += -I${python_include_path}
    endif
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
		   BiasedGradientSquaredDescent.o


#ifneq ($(or unitTests,check),)
#CXXFLAGS += -std=c++11
TEMPOBJ := $(filter-out ClientEON.o,$(OBJECTS))
DEPOBJECTS := $(addprefix ../,$(TEMPOBJ))
DEPLIBS := $(addprefix ../,$(LIBS))
export
#endif

#------------------------------------
#Build rules
all: $(POTDIRS) $(FPOTDIRS) eonclient docs
	@echo
	@echo "EON Client Compilation Options"
	@echo "BOINC: $(BOINC)" 
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

docs:
	@echo "Begin making Docs"
	@-doxygen $(DOXYCONFIG) \
	&& ([ $$? -eq 0 ] && echo "Docs have been generated and output to: $(DOXYDIR)")  \
	|| echo "Couldn't compile Docs; doxygen is not installed on this machine" 

docsCheck:
	@echo $(DOXYCONFIG) $(DOXYDIR) $(VERSION) $(BUILDDATE) $(BUILDHOST) $(BUILDUSER)

docsClean:
	rm -rf $(DOXYDIR)/*

clean: testsClean
	rm -f $(OBJECTS) $(DEPENDS) eonclient client

clean-all: clean docsClean testsClobber
	for pot in $(POTDIRS) $(FPOTDIRS) $(OPOTDIRS); do $(MAKE) -C $$pot clean ; done

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $<

DEPENDS= $(wildcard *.d)
-include $(DEPENDS)

.PHONY : all $(POTDIRS) $(FPOTDIRS) clean clean-all version.h
# DO NOT DELETE
