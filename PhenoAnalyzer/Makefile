
# Name of the main file and the target file
MAINFILE=PhenoAnalyzer.cc
#MAINFILE=PhenoAnalyzer.cc
TARGETFILE=PhenoAnalyzer
HEADERS==PhenoAnalyzer.h 
OTHERFILES=ROOTFunctions.h DelphesFunctions.h  #mt2_bisect.h  # Other files
Electron=ElectronA.h
#SOURCEFILES=mt2_bisect.cpp #mt2bl_bisect.cpp mt2w_bisect.cpp 

# Directory where HepMC is installed
HEPMCDIR=/usr/local/HEP/MG5_aMC_v2_6_3_2/HEPTools/hepmc/
HEPMCLIB=$(HEPMCDIR)/lib
HEPMCINC=$(HEPMCDIR)/include/

# Root flags and libs
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)  
ROOTGLIBS=$(shell root-config --glibs)  

# Directory of Exrootanalysis
EXROOTDIR=/usr/local/HEP/ExRootAnalysis
EXROOTLIB=$(EXROOTDIR)/lib
EXROOTINC=$(EXROOTDIR)/ExRootAnalysis


# Directory of Delphes
DELPHESDIR=/usr/local/HEP/Delphes-3.4.1
DELPHESLIB=${DELPHESDIR}
DELPHESINC=${DELPHESDIR}

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

CXXFLAGS += -Wall $(ROOTCFLAGS) -I./ -I$(DELPHESINC) -I$(EXROOTINC)
LDFLAGS += -Wall $(ROOTLIBS) -lGenVector

DIR = .
CXX = g++
LD = g++
CFLAGS = -std=c++11

SRCDIR = $(DIR)/src
OBJDIR = $(DIR)/obj
EXE = PhenoAnalyzer

#OBJECTS:=$(wildcard $(OBJDIR)/*.o)
OBJECT = $(OBJDIR)/PhenoAnalyzer.o $(OBJDIR)/ElectronA.o

all: $(EXE)

$(EXE): $(OBJECT)
	$(CXX) -o $@ $^ -L$(DELPHESLIB) -lDelphes -L$(EXROOTLIB) -lExRootAnalysis $(ROOTLIBS)

%: $(OBJDIR)/%.o
	$(LD) -o $@ $(LDFLAGS) $<  $(LIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cc $(SRCDIR)/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ 


.PHONY: clean

clean:
	@rm $(EXE) $(OBJDIR)/*
