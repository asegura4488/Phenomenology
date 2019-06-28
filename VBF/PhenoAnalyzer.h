////////////////////////////////////////////////////////////////
//                                                            //
// Author: Andrés Flórez, Universidad de los Andes, Colombia  //
//                                                            //
////////////////////////////////////////////////////////////////

#ifndef _PHENOANALYZER_H_
#define _PHENOANALYZER_H_

#include "TF2.h"
#include "TH1F.h"
#include "TChain.h"
#include "TEnv.h" 

#include <iostream>
#include <utility>
#include <vector>

#include "TROOT.h"
#include "TSystem.h"

#include "TString.h"

#include "TH2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TApplication.h"
#include "DelphesFunctions.h"
#include "TDirectory.h"
#include "TFile.h"
#include <fstream>

using namespace std;

//std::map<unsigned int, TH1*> _hmap_numero_taus;

class PhenoAnalysis {
public :
   PhenoAnalysis(TChain&, TFile*, TDirectory* dir[], int nDir);
   ~PhenoAnalysis();
   void crateHistoMasps (int);
   double calculateE(double, double, double);
   double calculateM(double, double, double);
   double calculateMt(double, double, double, double);
   double normalizedDphi(double);
   //Numero de eventos
   std::map<unsigned int, TH1*> _hmap_Nevents;

   // For Jets
   std::map<unsigned int, TH1*> _hmap_lead_jet_pT;
   std::map<unsigned int, TH1*> _hmap_lead_jet_eta;
   std::map<unsigned int, TH1*> _hmap_lead_jet_phi;
   std::map<unsigned int, TH1*> _hmap_sublead_jet_pT;
   std::map<unsigned int, TH1*> _hmap_sublead_jet_eta;
   std::map<unsigned int, TH1*> _hmap_sublead_jet_phi;
   std::map<unsigned int, TH1*> _hmap_delta_phi_jj;
   std::map<unsigned int, TH1*> _hmap_MET;
   std::map<unsigned int, TH1*> _hmap_jet_mjj;
   std::map<unsigned int, TH1*> _hmap_jet_mjj_formula;
   std::map<unsigned int, TH1*> _hmap_delta_eta_jj;
      
   // For Taus
   //std::map<unsigned int, TH1*> _hmap_numero_taus;
   std::map<unsigned int, TH1*> _hmap_tau1_pT;
   std::map<unsigned int, TH1*> _hmap_tau1_eta;
   std::map<unsigned int, TH1*> _hmap_tau1_phi;
   std::map<unsigned int, TH1*> _hmap_delta_phi_tau_met;
   std::map<unsigned int, TH1*> _hmap_razon_delta_phi_Pt;
   std::map<unsigned int, TH1*> _hmap_tau2_pT;
   std::map<unsigned int, TH1*> _hmap_tau2_eta;
   std::map<unsigned int, TH1*> _hmap_tau2_phi;
   std::map<unsigned int, TH1*> _hmap_tau3_pT;
   std::map<unsigned int, TH1*> _hmap_tau3_eta;
   std::map<unsigned int, TH1*> _hmap_tau3_phi;
private :

};

#endif
