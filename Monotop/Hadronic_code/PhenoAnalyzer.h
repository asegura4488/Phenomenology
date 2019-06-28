////////////////////////////////////////////////////////////////
//                                                            //
// Author: Andrés Flórez, Universidad de los Andes, Colombia  //
//                                                            //
////////////////////////////////////////////////////////////////

#ifndef _PHENOANALYZER_H_
#define _PHENOANALYZER_H_

#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
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

class PhenoAnalysis {
public :
   PhenoAnalysis(TChain&, TFile*, TDirectory* dir[], int nDir);
   ~PhenoAnalysis();
   void crateHistoMasps (int);
   double calculateE(double, double, double);
   double calculateM(double, double, double);
   double CalculateP(double, double);
   double CalculateP12_square_2particle(double, double, double, double);
   double CalculateP123_square_3particle(double, double, double, double, double, double, double, double, double);
   double normalizedDphi(double);
   double transverseMass(double, double, double);
   // Number of events
   std::map<unsigned int, TH1*> _hmap_Nevents;
   // Jets information
   std::map<unsigned int, TH1*> _hmap_N_jets;
   std::map<unsigned int, TH1*> _hmap_jet1_pT;
   std::map<unsigned int, TH1*> _hmap_jet1_eta;
   std::map<unsigned int, TH1*> _hmap_jet1_phi;

   std::map<unsigned int, TH1*> _hmap_N_jets2;
   std::map<unsigned int, TH1*> _hmap_jet2_pT;
   std::map<unsigned int, TH1*> _hmap_jet2_eta;
   std::map<unsigned int, TH1*> _hmap_jet2_phi;
   // Taus information
   std::map<unsigned int, TH1*> _hmap_N_tau;
   std::map<unsigned int, TH1*> _hmap_tau1_pT;
   std::map<unsigned int, TH1*> _hmap_tau1_eta;
   std::map<unsigned int, TH1*> _hmap_tau1_phi;
   std::map<unsigned int, TH1*> _hmap_tau2_pT;
   std::map<unsigned int, TH1*> _hmap_tau2_eta;
   std::map<unsigned int, TH1*> _hmap_tau2_phi;
   // b-jets information
   std::map<unsigned int, TH1*> _hmap_N_bjet;
   std::map<unsigned int, TH1*> _hmap_bjet_plus_lep;   
   std::map<unsigned int, TH1*> _hmap_bjet_ratio_pT;
   std::map<unsigned int, TH1*> _hmap_Hbjet_ratio_pT;
   std::map<unsigned int, TH1*> _hmap_bjet_pT;
   std::map<unsigned int, TH1*> _hmap_bjet_eta;
   std::map<unsigned int, TH1*> _hmap_bjet_phi;
   std::map<unsigned int, TH1*> _hmap_bjet_energy;
   std::map<unsigned int, TH1*> _hmap_top_energy;
   // electrons information
   std::map<unsigned int, TH1*> _hmap_N_elec;
   std::map<unsigned int, TH1*> _hmap_elec_pT;
   std::map<unsigned int, TH1*> _hmap_elec_eta;
   std::map<unsigned int, TH1*> _hmap_elec_phi;
   // muons information
   std::map<unsigned int, TH1*> _hmap_N_muon;
   std::map<unsigned int, TH1*> _hmap_muon_pT;
   std::map<unsigned int, TH1*> _hmap_muon_eta;
   std::map<unsigned int, TH1*> _hmap_muon_phi;
   // met 
   std::map<unsigned int, TH1*> _hmap_met;
   std::map<unsigned int, TH1*> _hmap_mt;
   std::map<unsigned int, TH1*> _hmap_mt_top;
   std::map<unsigned int, TH1*> _hmap_pt_w;
   std::map<unsigned int, TH1*> _hmap_pt_bl;
   std::map<unsigned int, TH1*> _hmap_pt_t;
   std::map<unsigned int, TH1*> _hmap_pt_top;
   std::map<unsigned int, TH1*> _hmap_pt_W;
   std::map<unsigned int, TH1*> _hmap_e_tb;
   std::map<unsigned int, TH1*> _hmap_e_t;
   std::map<unsigned int, TH1*> _hmap_delta_r;
   std::map<unsigned int, TH1*> _hmap_jet_1_met_dphi;
   std::map<unsigned int, TH1*> _hmap_jet_2_met_dphi;
   std::map<unsigned int, TH1*> _hmap_invariant_mjjb;
   std::map<unsigned int, TH1*> _hmap_invariant_mjj;
   std::map<unsigned int, TH1*> _hmap_invariant_mj1;
   std::map<unsigned int, TH1*> _hmap_invariant_mj2;
   std::map<unsigned int, TH1*> _hmap_invariant_pj1;
   std::map<unsigned int, TH1*> _hmap_invariant_pj2;
   std::map<unsigned int, TH1*> _hmap_invariant_ej1;
   std::map<unsigned int, TH1*> _hmap_invariant_ej2;
   std::map<unsigned int, TH1*> _hmap_invariant_ej1j2;
   std::map<unsigned int, TH1*> _hmap_met_jet1_met_dphi;
   std::map<unsigned int, TH1*> _hmap_mt_mt_top;
   std::map<unsigned int, TH1*> _hmap_jet_consti;
   std::map<unsigned int, TH1*> _hmap_jet_consti2;
   std::map<unsigned int, TH1*> _hmap_top_jet;
private :

};

#endif
