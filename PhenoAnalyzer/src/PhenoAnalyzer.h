////////////////////////////////////////////////////////////////
//                                                            //
// Author: Andrés Flórez, Universidad de los Andes, Colombia  //
//                                                            //
////////////////////////////////////////////////////////////////

#ifndef PHENOANALYZER_H
#define PHENOANALYZER_H

#include "ROOTFunctions.h"
#include "DelphesFunctions.h"
#include "ElectronA.h"

using namespace std;

class PhenoAnalysis {
public :
   PhenoAnalysis(TChain&, TFile*, TDirectory* dir[], int nDir);
   ~PhenoAnalysis();
   TLorentzVector Cuts(ElectronA *Electron_);
   void crateHistoMasps (int);
   bool overlapingObjects(double, double, double, double, double);
   double calculateE(double, double, double);
   double normalizedDphi(double);
   double calculate_deltaR(TLorentzVector,  Track*);

   // For Jets
   std::map<unsigned int, TH1*> _hmap_Nevents;
   std::map<unsigned int, TH1*> _hmap_jet_pt;
   std::map<unsigned int, TH1*> _hmap_Nev_bins;
   //Lina y Fabio 
   std::map<unsigned int, TH1*> _hmap_MET;
   std::map<unsigned int, TH1*> _hmap_jet_pt_cuts;
   std::map<unsigned int, TH1*> _hmap_jet_eta_cuts;
   std::map<unsigned int, TH1*> _hmap_jet_phi_cuts;
   std::map<unsigned int, TH1*> _hmap_muon_pt;
   std::map<unsigned int, TH1*> _hmap_muon_pt_cuts;
   std::map<unsigned int, TH1*> _hmap_muon_eta_cuts;
   std::map<unsigned int, TH1*> _hmap_muon_phi_cuts;
   std::map<unsigned int, TH1*> _hmap_muon_phi;
   std::map<unsigned int, TH1*> _hmap_electron_pt;
   std::map<unsigned int, TH1*> _hmap_electron_eta;
   std::map<unsigned int, TH1*> _hmap_electron_phi;
   std::map<unsigned int, TH1*> _hmap_electron_met_phi;
   std::map<unsigned int, TH1*> _hmap_muon_met_phi;
   std::map<unsigned int, TH1*> _hmap_tau_met_phi;
   std::map<unsigned int, TH1*> _hmap_electron_jet_phi;
   std::map<unsigned int, TH1*> _hmap_muon_jet_phi;
   std::map<unsigned int, TH1*> _hmap_tau_jet_phi;
   std::map<unsigned int, TH1*> _hmap_jet_met_phi;
   std::map<unsigned int, TH1*> _hmap_pt_ratio_electron;
   std::map<unsigned int, TH1*> _hmap_pt_ratio_muon;
   std::map<unsigned int, TH1*> _hmap_pt_ratio_tau;
   std::map<unsigned int, TH1*> _hmap_mt;
   std::map<unsigned int, TH2*> _hmap_ptlead_vs_met;
   
   };

#endif
