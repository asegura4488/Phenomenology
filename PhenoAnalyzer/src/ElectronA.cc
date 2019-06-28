#include "ElectronA.h"

ElectronA::ElectronA(TClonesArray* branchelectron_, double mass_): mass(mass_){

// branchelectron_ = treeReader->UseBranch("Electron");

 for(Long64_t it = 0; it < branchelectron_->GetEntriesFast(); it++){
 Electron *electron = (Electron*)branchelectron_->At(it);
 TLorentzVector tmp;
 tmp.SetPtEtaPhiE(electron->PT, electron->Eta, electron->Phi, GetEnergy(electron->PT, electron->Eta));
 std::cout << it << " " <<  tmp.Pt() << " "  << tmp.Phi() << tmp.E() << std::endl;
 ElectronReco.push_back(tmp);
 } 

 cur_P = &ElectronReco;
}

TLorentzVector ElectronA::p4(uint index) const {return (cur_P->at(index));}

double ElectronA::GetEnergy(double pt, double eta){

  double theta = TMath::ATan(TMath::Exp(-eta));
  double sin_theta = TMath::Sin(2*theta);
  double p= pt/sin_theta;
  double e = sqrt(pow(p, 2) + pow(mass, 2));
  return e;

}

TLorentzVector ElectronA::GetMinVectorPt(){

  TLorentzVector tmp;  

  double minpt=1e+11;
  std::vector<TLorentzVector>::iterator it;
  for(it = ElectronReco.begin(); it != ElectronReco.end(); it++){
  
  if((*it).Pt() < minpt ){ 
         minpt = (*it).Pt();
         tmp = (*it);
    }
  }

  return tmp;

}

ElectronA::~ElectronA(){
}
