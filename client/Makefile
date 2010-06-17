OBJECTS := ClientEON.o SaddlePoint.o Dimer.o DimerConstants.o \
	Eigenvalues.o EpiCenters.o Prefactors.o ConjugateGradients.o\
	Constants.o debug.o HelperFunctions.o Matter.o Parameters.o\
	PotentialsInterface.o Potentials.o QMBox.o Quickmin.o SaddlePoint.o
#CXXFLAGS += -I..
# Libraries to link the executable 
#LIBS := ./common/libCommon.a
# define NDEBUG to link to boinc and remove debuggin lines
ifdef NDEBUG
   LIBS+= -lboinc_api -lboinc -lpthread
   CXXFLAGS += -DNDEBUG -static
endif
# sub-directories with a makefile and something to build, such as a library
#SUBDIRS :=  ./common/
ifndef NO_POT_EMT
   SUBDIRS += ./potentials/EMT/
   LIBS += ./potentials/EMT/libEMT.a
endif
ifndef NO_MORSE
   SUBDIRS += ./potentials/Morse/
   LIBS += ./potentials/Morse/libMorse.a
endif
ifndef NO_AL
   SUBDIRS += ./potentials/Aluminum/
   LIBS += ./potentials/Aluminum/libAL.a
   CXXFLAGS += -lgfortran
endif

all : client
client : $(OBJECTS) $(LIBS)
# In order to link, libCommon.a has to be inserted both before and after the potentials' libraries. 
	g++  $(CXXFLAGS) -o $@ $^
#$(SUBDIRS) :
#	export CXXFLAGS="$(CXXFLAGS)" ; $(MAKE) -C $@

clean :
	rm -f $(OBJECTS) client

clean-all : clean
	for subdir in $(SUBDIRS) ; do $(MAKE) -C $$subdir clean ; done

.PHONY : all $(SUBDIRS) clean clean-all
