OBJECTS := ClientEON.o SaddlePoint.o Dimer.o DimerConstants.o \
	Eigenvalues.o EpiCenters.o Prefactors.o ConjugateGradients.o\
	Constants.o debug.o HelperFunctions.o Matter.o Parameters.o\
	PotentialsInterface.o Potentials.o QMBox.o Quickmin.o SaddlePoint.o
CXXFLAGS += -O2 -Wall
# Libraries to link the executable 
#LIBS := ./common/libCommon.a
# define NDEBUG to remove debugging lines
ifdef NDEBUG
   CXXFLAGS += -DNDEBUG
endif

# define BOINC to link against boinc libs
ifdef BOINC
   CXXFLAGS += -DNDEBUG -DBOINC
   LDFLAGS += -static -lboinc_api -lboinc -lpthread
endif
# sub-directories with a makefile and something to build, such as a library
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
   LDFLAGS += -lgfortran
else
   CXXFLAGS += -DNO_AL
endif

all : $(SUBDIRS) client

client : $(OBJECTS) $(LIBS)  
	g++  -o $@ $^ $(LDFLAGS)

$(LIBS):
	$(MAKE) -C $@

$LIBS: $(SUBDIRS)

$(SUBDIRS) :
	$(MAKE) -C $@

clean :
	rm -f $(OBJECTS) client

clean-all : clean
	for subdir in $(SUBDIRS) ; do $(MAKE) -C $$subdir clean ; done

.PHONY : all $(SUBDIRS) clean clean-all
