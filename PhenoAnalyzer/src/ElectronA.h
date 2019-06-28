#ifndef ElectronA_h
#define ElectronA_h

#include "ROOTFunctions.h"
#include "ExRootTreeReader.h"
//#include "HepMCFunctions.h"
#include "DelphesFunctions.h"
#include<iostream>

class ElectronA{

public:
    ElectronA(TClonesArray* branchelectron_, double mass_);
    TLorentzVector p4(uint index) const; 
    TLorentzVector GetMinVectorPt();
    ~ElectronA();

    double GetEnergy(double pt, double eta);

private:

 std::vector<TLorentzVector> ElectronReco;
 std::vector<TLorentzVector> *cur_P;
 double mass; 

};

#endif

