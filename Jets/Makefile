LIBS    =`root-config --libs`
INCS    =`root-config --cflags`
CCOPTS  =-g $(LIBS) $(INCS)

#ROOTPHYLIB=/home/ma.segura10/root/lib
ROOTPHYLIB=/usr/local/root/lib/root

# Directory of Exrootanalysis
# EXROOTLIB=/home/ma.segura10/background/ExRootAnalysis
# EXROOTINC=/home/ma.segura10/background/ExRootAnalysis/ExRootAnalysis

 # Directory of Delphes
 DELPHESDIR=/usr/local/Delphes-3.2.0
 DELPHESDIR1=/home/ma.segura10/background/Delphes-3.3.2
 DELPHESLIB=$(DELPHESDIR1)
 DELPHESINC=$(DELPHESDIR1)

 # Directory of Exrootanalysis
 EXROOTLIB=$(DELPHESDIR1)
 EXROOTINC=$(DELPHESDIR1)/external/ExRootAnalysis

 

compile: jet.cpp
	g++ -o jet jet.cpp $(CCOPTS) \
	-I$(DELPHESINC) -I$(EXROOTINC) \
        -L$(DELPHESLIB) -lDelphes -L$(EXROOTLIB) -lExRootAnalysis -L$(ROOTPHYLIB) -lPhysics
