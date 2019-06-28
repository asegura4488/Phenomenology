#include <iostream>
#include "ROOTFunctions.h"
#include "DelphesFunctions.h"
#include "mt2_bisect.h"
#include "mt2w_bisect.h"
#include "mt2bl_bisect.h"

//Global variables
float mt=173.34;
float mw=80.40;
double m_electron=0.511e-3;
double m_muon=105e-3;
double mn    = 0.;
double pi = 4.0*atan(1.0);

using namespace std;

//Functions definitions
//Energy function
double CalculateE(double Eta, double PT, double Mass){
  double theta = 2*TMath::ATan(TMath::Exp(-Eta));
  double sin_theta = TMath::Sin(theta);
  double p = PT/sin_theta;
  double e  = TMath::Sqrt(TMath::Power(Mass,2)+TMath::Power(p,2));
  return e;
}
//P12_square(two particles)
double CalculateP12_square_2particle(double Eta1, double Phi1, double PT1, double Eta2, double Phi2, double PT2){
  double theta1 = 2*TMath::ATan(TMath::Exp(-Eta1));
  double theta2 = 2*TMath::ATan(TMath::Exp(-Eta2));
  double px1 = PT1*TMath::Cos(Phi1+TMath::Pi());
  double py1 = PT1*TMath::Sin(Phi1+TMath::Pi());
  double pz1 = PT1*TMath::Cos(theta1)/TMath::Sin(theta1);
  double px2 = PT2*TMath::Cos(Phi2+TMath::Pi());
  double py2 = PT2*TMath::Sin(Phi2+TMath::Pi());
  double pz2 = PT2*TMath::Cos(theta2)/TMath::Sin(theta2);
  double p2  = TMath::Power((px1+px2),2)+TMath::Power((py1+py2),2)+TMath::Power((pz1+pz2),2);
  return p2;
}

//P123_square(three particles)
double CalculateP123_square_3particle(double Eta1, double Phi1, double PT1, double Eta2, double Phi2, double PT2, double Eta3, double Phi3, double PT3){
  double theta1 = 2*TMath::ATan(TMath::Exp(-Eta1));
  double theta2 = 2*TMath::ATan(TMath::Exp(-Eta2));
  double theta3 = 2*TMath::ATan(TMath::Exp(-Eta3));
  double px1 = PT1*TMath::Cos(Phi1+TMath::Pi());
  double py1 = PT1*TMath::Sin(Phi1+TMath::Pi());
  double pz1 = PT1*TMath::Cos(theta1)/TMath::Sin(theta1);
  double px2 = PT2*TMath::Cos(Phi2+TMath::Pi());
  double py2 = PT2*TMath::Sin(Phi2+TMath::Pi());
  double pz2 = PT2*TMath::Cos(theta2)/TMath::Sin(theta2);
  double px3 = PT3*TMath::Cos(Phi3+TMath::Pi());
  double py3 = PT3*TMath::Sin(Phi3+TMath::Pi());
  double pz3 = PT3*TMath::Cos(theta3)/TMath::Sin(theta3);
  double p2  = TMath::Power((px1+px2+px3),2)+TMath::Power((py1+py2+py3),2)+TMath::Power((pz1+pz2+pz3),2);
  return p2;
}

// Started
int main(int argc, char *argv[]) {
  
  // Load shared library
  gSystem->Load("lib/libExRootAnalysis.so");
  gSystem->Load("libPhysics");

  // Create chains of root trees
  
  // Root ntuplas ttbar semileptonico
  
  TChain chain("Delphes");
  chain.Add(argv[1]);
  
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader_Delphes = new ExRootTreeReader(&chain);
  
  // Number of entries access
  long int numberOfEntries = treeReader_Delphes->GetEntries();
  
  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader_Delphes->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader_Delphes->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader_Delphes->UseBranch("Muon");
  TClonesArray *branchMissingET = treeReader_Delphes->UseBranch("MissingET");
  
  // Book histograms
  
  TH1 *histNumJets = new TH1F("num_jets","numjets",100, 0.0, 10.0);
  TH1 *histNumElectrons = new TH1F("num_electrons","numelectrons",100,0.0,10.0);
  TH1 *histNumMuons = new TH1F("num_muons","nummuons",100,0.0,10.0);
  TH1 *histsigma_123 = new TH1F("SIGMA_123","sigma_123", 100,0.0,300.0);
  TH1 *histNum_btags = new TH1F("Num_btags","num_btags",100,-1.0,4.0);
  TH1 *hist_w_recons_2bjet = new TH1F("W_recons_2bjet","w_recons_2bjet",100,40.0,120.0);
  TH1 *hist_t_recons_2bjet = new TH1F("T_recons_2bjet","t_recons_2bjet",100,110.0,240.0);
  TH1 *hist_w_recons_1bjet = new TH1F("W_recons_1bjet","w_recons_1bjet",100,40.0,120.0);
  TH1 *hist_t_recons_1bjet = new TH1F("T_recons_1bjet","t_recons_1bjet",100,110.0,240.0);
  TH1 *hist_MT_Muon = new TH1F("Transverse_mass_muon_","mt_muon",100,0,200.0);
  TH1 *hist_MT_Electron = new TH1F("Transverse_mass_electron_","mt_electron",100,0,200.0);
  TH1 *hist_MT2 = new TH1F("MT2","mt2",100,150,350.0);
  TH1 *hist_MT2W = new TH1F("MT2W","mt2w",100,0,300.0);
  TH1 *hist_MT2BL = new TH1F("MT2BL","mt2bl",100,0,300.0);
  TH1 *hist_Eta_average = new TH1F("Eta_average","eta_average",100,-7.0,7.0);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TH1 *hist_MT_lepton = new TH1F("MT_lepton","mt_lepton",150,0.0,300.0);
  TH1 *hist_MT_lepton_1btag_corte_200 = new TH1F("MT_lepton_1btag_corte_200","mt_lepton_1btag_corte_200",150,0.0,200.0);
  TH1 *hist_MT_lepton_1btag_corte_300 = new TH1F("MT_lepton_1btag_corte_300","mt_lepton_1btag_corte_300",150,0.0,200.0);
  TH1 *hist_MT_lepton_1btag_corte_400 = new TH1F("MT_lepton_1btag_corte_400","mt_lepton_1btag_corte_400",150,0.0,200.0);


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1F *hist_Cociente_met_raiz_ht_2btag_corte_200 = new TH1F("Cociente_met_raiz_ht_2btag_corte_200","cociente_met_raiz_ht_2btag_corte_200",100,0.0,30.0);
  TH1F *hist_Cociente_met_raiz_ht_2btag_corte_300 = new TH1F("Cociente_met_raiz_ht_2btag_corte_300","cociente_met_raiz_ht_2btag_corte_300",100,0.0,30.0);
  TH1F *hist_Cociente_met_raiz_ht_2btag_corte_400 = new TH1F("Cociente_met_raiz_ht_2btag_corte_400","cociente_met_raiz_ht_2btag_corte_400",100,0.0,30.0);


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1 *hist_Cociente_pt_event_pt_leading_2btag_corte_200 = new TH1F("Cociente_pt_event_pt_leading_2btag_corte_200","cociente_pt_event_pt_leading_2btag_corte_200",100,0.0,2.0);
  TH1 *hist_Cociente_pt_event_pt_leading_2btag_corte_300 = new TH1F("Cociente_pt_event_pt_leading_2btag_corte_300","cociente_pt_event_pt_leading_2btag_corte_300",100,0.0,2.0);
  TH1 *hist_Cociente_pt_event_pt_leading_2btag_corte_400 = new TH1F("Cociente_pt_event_pt_leading_2btag_corte_400","cociente_pt_event_pt_leading_2btag_corte_400",100,0.0,2.0);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TH1 *hist_Delta_phi_w_jets_bjet_hadronic_2btag_corte_200 = new TH1F("Delta_phi_w_jets_bjet_hadronic_2btag_corte_200","delta_phi_w_jets_bjet_hadronic_2btag_corte_200",100,0.0,4.0);
  TH1 *hist_Delta_phi_w_jets_bjet_hadronic_2btag_corte_300 = new TH1F("Delta_phi_w_jets_bjet_hadronic_2btag_corte_300","delta_phi_w_jets_bjet_hadronic_2btag_corte_300",100,0.0,4.0);
  TH1 *hist_Delta_phi_w_jets_bjet_hadronic_2btag_corte_400 = new TH1F("Delta_phi_w_jets_bjet_hadronic_2btag_corte_400","delta_phi_w_jets_bjet_hadronic_2btag_corte_400",100,0.0,4.0);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1 *hist_Delta_phi_hadronic_top_leading_jet_2btag_corte_200 = new TH1F("Delta_phi_hadronic_top_leading_jet_2btag_corte_200","Delta_phi_hadronic_top_leading_jet_2btag_corte_200",100,0.0,4.0);
  TH1 *hist_Delta_phi_hadronic_top_leading_jet_2btag_corte_300 = new TH1F("Delta_phi_hadronic_top_leading_jet_2btag_corte_300","Delta_phi_hadronic_top_leading_jet_2btag_corte_300",100,0.0,4.0);
  TH1 *hist_Delta_phi_hadronic_top_leading_jet_2btag_corte_400 = new TH1F("Delta_phi_hadronic_top_leading_jet_2btag_corte_400","Delta_phi_hadronic_top_leading_jet_2btag_corte_400",100,0.0,4.0);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1 *hist_Delta_phi_leptonic_top_leading_jet_2btag_corte_200 = new TH1F("Delta_phi_leptonic_top_leading_jet_2btag_corte_200","Delta_phi_leptonic_top_leading_jet_2btag_corte_200",100,0.0,4.0);
  TH1 *hist_Delta_phi_leptonic_top_leading_jet_2btag_corte_300 = new TH1F("Delta_phi_leptonic_top_leading_jet_2btag_corte_300","Delta_phi_leptonic_top_leading_jet_2btag_corte_300",100,0.0,4.0);
  TH1 *hist_Delta_phi_leptonic_top_leading_jet_2btag_corte_400 = new TH1F("Delta_phi_leptonic_top_leading_jet_2btag_corte_400","Delta_phi_leptonic_top_leading_jet_2btag_corte_400",100,0.0,4.0);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1 *hist_Delta_Phi_lepton_bl_2btag_corte_200 = new TH1F("Delta_Phi_lepton_bl_2btag_corte_200","delta_Phi_lepton_bl_2btag_corte_200",100,0.0,4.0);
  TH1 *hist_Delta_Phi_lepton_bl_2btag_corte_300 = new TH1F("Delta_Phi_lepton_bl_2btag_corte_300","delta_Phi_lepton_bl_2btag_corte_300",100,0.0,4.0);
  TH1 *hist_Delta_Phi_lepton_bl_2btag_corte_400 = new TH1F("Delta_Phi_lepton_bl_2btag_corte_400","delta_Phi_lepton_bl_2btag_corte_400",100,0.0,4.0);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Delta_phi_leading_jet_muon
  TH1 *hist_Delta_Phi_leading_jet_muon_2btag = new TH1F("DeltaPhi_leading_jet_muon_2btag","deltaphi_leading_jet_muon",100,0.0,7.0);
  TH1 *hist_Delta_Phi_leading_jet_muon_2btag_corte_200 = new TH1F("DeltaPhi_leading_jet_muon_2btag_corte_200","deltaphi_leading_jet_muon_corte_200",100,0.0,7.0);
  TH1 *hist_Delta_Phi_leading_jet_muon_2btag_corte_300 = new TH1F("DeltaPhi_leading_jet_muon_2btag_corte_300","deltaphi_leading_jet_muon_corte_300",100,0.0,7.0);
  TH1 *hist_Delta_Phi_leading_jet_muon_2btag_corte_400 = new TH1F("DeltaPhi_leading_jet_muon_2btag_corte_400","deltaphi_leading_jet_muon_corte_400",100,0.0,7.0);
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1 *hist_Delta_Phi_leading_jet_electron_2btag = new TH1F("DeltaPhi_leading_jet_electron_2btag","deltaphi_leading_jet_electron",100,-7.0,7.0);
  TH1 *hist_cociente_HT_PT_leading_jet = new TH1F("Cociente_HT_PT_leading_jet","cociente_HT_PT_leading_jet", 100, 0.0,10.0);
  TH1 *hist_cociente_HT_nobjet_PT_leading_jet = new TH1F("Cociente_HT_nobjet_PT_leading_jet","cociente_HT_nobjet_PT_leading_jet", 100, 0.0, 10.0);
  TH1 *hist_cociente_HT_bjet_PT_leading_jet = new TH1F("Cociente_HT_bjet_PT_leading_jet","cociente_HT_bjet_PT_leading_jet", 100, 0.0, 10.0);
  TH1 *hist_cociente_PT_Electron_PT_leading_jet = new TH1F("Cociente_PT_electron_PT_leading_jet","cociente_PT_electron_PT_leading_jet",100,0.0,3.0);
  TH1 *hist_cociente_PT_Muon_PT_leading_jet = new TH1F("Cociente_PT_muon_PT_leading_jet","cociente_PT_muon_PT_leading_jet",100,0.0,3.0);
  TH1 *hist_rt_2btag = new TH1F ("RT_2btag","rt_2btag",100,0.0,5.0);
  TH1 *hist_delta_phi_th_leading_jet_2btag = new TH1F("Delta_phi_th_leading_jet_2btag","delta_phi_th_leading_jet_2btag",100,0.0,7.0);
  TH1* hist_delta_eta_th_leading_jet_2btag = new TH1F("Delta_eta_th_leading_jet_2btag","delta_eta_th_leading_jet_2btag",100,0.0,8.0);
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Interesting variables

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1 *hist_Cociente_pt_visible_st_2btag_corte_200 = new TH1F("Cociente_pt_visible_st_2btag_corte_200","cociente_pt_visible_st_2btag_corte_200",100.0,0.0,550.0);
  TH1 *hist_Cociente_pt_visible_st_2btag_corte_300 = new TH1F("Cociente_pt_visible_st_2btag_corte_300","cociente_pt_visible_st_2btag_corte_300",100.0,0.0,550.0);
  TH1 *hist_Cociente_pt_visible_st_2btag_corte_400 = new TH1F("Cociente_pt_visible_st_2btag_corte_400","cociente_pt_visible_st_2btag_corte_400",100.0,0.0,550.0);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1 *hist_Cociente_met_st_1btag_corte_200 = new TH1F("Cociente_met_st_1btag_corte_200","cociente_met_st_1btag_corte_200",100.0,0.0,2.0);
  TH1 *hist_Cociente_met_st_1btag_corte_300 = new TH1F("Cociente_met_st_1btag_corte_300","cociente_met_st_1btag_corte_300",100.0,0.0,2.0);
  TH1 *hist_Cociente_met_st_1btag_corte_400 = new TH1F("Cociente_met_st_1btag_corte_400","cociente_met_st_1btag_corte_400",100.0,0.0,2.0);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1F *hist_Delta_eta_ave_1btag_corte_200 = new TH1F("Delta_eta_ave_1btag_corte_200","delta_eta_ave_1btag_corte_200",100,0.0,7.0);
  TH1F *hist_Delta_eta_ave_1btag_corte_300 = new TH1F("Delta_eta_ave_1btag_corte_300","delta_eta_ave_1btag_corte_300",100,0.0,7.0);
  TH1F *hist_Delta_eta_ave_1btag_corte_400 = new TH1F("Delta_eta_ave_1btag_corte_400","delta_eta_ave_1btag_corte_400",100,0.0,7.0);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Cociente PT_top_lep/PT_leading_jet
  TH1F * hist_Cociente_PT_top_lept_PT_leading_jet_2btag_corte_200 = new TH1F("Cociente_PT_top_lept_PT_leading_jet_2btag_corte_200","cociente_PT_top_lept_PT_leading_jet_2btag_corte_200",100,0.0,3.0);
  TH1F * hist_Cociente_PT_top_lept_PT_leading_jet_2btag_corte_300 = new TH1F("Cociente_PT_top_lept_PT_leading_jet_2btag_corte_300","cociente_PT_top_lept_PT_leading_jet_2btag_corte_300",100,0.0,3.0);
  TH1F * hist_Cociente_PT_top_lept_PT_leading_jet_2btag_corte_400 = new TH1F("Cociente_PT_top_lept_PT_leading_jet_2btag_corte_400","cociente_PT_top_lept_PT_leading_jet_2btag_corte_400",100,0.0,3.0);
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Cociente PT_top_lep/PT_top_hadronico
  TH1F *hist_Cociente_PT_top_lept_PT_top_had_2btag_corte_200 = new TH1F("Cociente_PT_top_lept_PT_top_had_2btag_corte_200","Cociente_PT_top_lept_PT_top_had_2btag_corte_200",100,0.0,3.0);
  TH1F *hist_Cociente_PT_top_lept_PT_top_had_2btag_corte_300 = new TH1F("Cociente_PT_top_lept_PT_top_had_2btag_corte_300","Cociente_PT_top_lept_PT_top_had_2btag_corte_300",100,0.0,3.0);
  TH1F *hist_Cociente_PT_top_lept_PT_top_had_2btag_corte_400 = new TH1F("Cociente_PT_top_lept_PT_top_had_2btag_corte_400","Cociente_PT_top_lept_PT_top_had_2btag_corte_400",100,0.0,3.0);
  /////// //////////////////////////////////////////////////////////////////////////////////////////////////
  //MT_top
  TH1F *hist_MT_top_2btag = new TH1F("MT_top_2btag","mt_top_2btag",100,0.0,300.0);
  TH1F *hist_MT_top_2btag_corte_200 = new TH1F("MT_top_2btag_corte_200","mt_top_2btag_corte_200",100,0.0,600.0);
  TH1F *hist_MT_top_2btag_corte_300 = new TH1F("MT_top_2btag_corte_300","mt_top_2btag_corte_300",100,0.0,600.0);
  TH1F *hist_MT_top_2btag_corte_400 = new TH1F("MT_top_2btag_corte_400","mt_top_2btag_corte_400",100,0.0,600.0);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //MT_lb
  TH1F *hist_MT_lb_2btag =new TH1F("MT_lb_2btag","mt_lb_2btag",100,0.0,600.0);
  TH1F *hist_MT_lb_2btag_corte_200 =new TH1F("MT_lb_2btag_corte_200","mt_lb_2btag_corte_200",100,0.0,600.0);
  TH1F *hist_MT_lb_2btag_corte_300 =new TH1F("MT_lb_2btag_corte_300","mt_lb_2btag_corte_300",100,0.0,600.0);
  TH1F *hist_MT_lb_2btag_corte_400 =new TH1F("MT_lb_2btag_corte_400","mt_lb_2btag_corte_400",100,0.0,600.0);
  TH1F *hist_MT_lb_1btag_corte_200 =new TH1F("MT_lb_1btag_corte_200","mt_lb_1btag_corte_200",100,0.0,600.0);
  TH1F *hist_MT_lb_1btag_corte_300 =new TH1F("MT_lb_1btag_corte_300","mt_lb_1btag_corte_300",100,0.0,600.0);
  TH1F *hist_MT_lb_1btag_corte_400 =new TH1F("MT_lb_1btag_corte_400","mt_lb_1btag_corte_400",100,0.0,600.0);
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  //ST
  TH1 *hist_ST= new TH1F("ST","st",100,0.0,10.0);
  TH1 *hist_ST_2btag_corte_200= new TH1F("ST_2btag_corte_200","st_2btag_corte_200",100,0.0,4.0);
  TH1 *hist_ST_2btag_corte_300= new TH1F("ST_2btag_corte_300","st_2btag_corte_300",100,0.0,4.0);
  TH1 *hist_ST_2btag_corte_400= new TH1F("ST_2btag_corte_400","st_2btag_corte_400",100,0.0,4.0);
  TH1 *hist_ST_1btag_corte_200= new TH1F("ST_1btag_corte_200","st_1btag_corte_200",100,0.0,4.0);
  TH1 *hist_ST_1btag_corte_300= new TH1F("ST_1btag_corte_300","st_1btag_corte_300",100,0.0,4.0);
  TH1 *hist_ST_1btag_corte_400= new TH1F("ST_1btag_corte_400","st_1btag_corte_400",100,0.0,4.0);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Cociente p(visible)/leading_jet_PT
  TH1 *hist_Cociente_p_visible_leading_jet_pt = new TH1F("Cociente_p_visible_leading_jet_pt","cociente_p_visible_leading_jet_pt", 100,0.0,3.0);
  TH1 *hist_Cociente_p_visible_leading_jet_pt_2btag_corte_200 = new TH1F("Cociente_p_visible_leading_jet_pt_2btag_corte_200","cociente_p_visible_leading_jet_pt_2btag_corte_200", 100, 0.0,3.0);
  TH1 *hist_Cociente_p_visible_leading_jet_pt_2btag_corte_300 = new TH1F("Cociente_p_visible_leading_jet_pt_2btag_corte_300","cociente_p_visible_leading_jet_pt_2btag_corte_300", 100, 0.0,3.0);
  TH1 *hist_Cociente_p_visible_leading_jet_pt_2btag_corte_400 = new TH1F("Cociente_p_visible_leading_jet_pt_2btag_corte_400","cociente_p_visible_leading_jet_pt_2btag_corte_400", 100, 0.0,3.0);
  TH1 *hist_Cociente_p_visible_leading_jet_pt_1btag_corte_200 = new TH1F("Cociente_p_visible_leading_jet_pt_1btag_corte_200","cociente_p_visible_leading_jet_pt_1btag_corte_200", 100, 0.0,3.0);
  TH1 *hist_Cociente_p_visible_leading_jet_pt_1btag_corte_300 = new TH1F("Cociente_p_visible_leading_jet_pt_1btag_corte_300","cociente_p_visible_leading_jet_pt_1btag_corte_300", 100, 0.0,3.0);
  TH1 *hist_Cociente_p_visible_leading_jet_pt_1btag_corte_400 = new TH1F("Cociente_p_visible_leading_jet_pt_1btag_corte_400","cociente_p_visible_leading_jet_pt_1btag_corte_400", 100, 0.0,3.0);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Delta_phi_p_visible_leading_jet
  TH1 *hist_Delta_phi_p_visible_leading_jet = new TH1F("Delta_phi_p_visible_leading_jet","delta_phi_p_visible_leading_jet",100,0.0,8.0);
  TH1 *hist_Delta_phi_p_visible_leading_jet_2btag_corte_200 = new TH1F("Delta_phi_p_visible_leading_jet_2btag_corte_200","delta_phi_p_visible_leading_jet_2btag_corte_200",100,0.0,8.0);
  TH1 *hist_Delta_phi_p_visible_leading_jet_2btag_corte_300 = new TH1F("Delta_phi_p_visible_leading_jet_2btag_corte_300","delta_phi_p_visible_leading_jet_2btag_corte_300",100,0.0,8.0);
  TH1 *hist_Delta_phi_p_visible_leading_jet_2btag_corte_400 = new TH1F("Delta_phi_p_visible_leading_jet_2btag_corte_400","delta_phi_p_visible_leading_jet_2btag_corte_400",100,0.0,8.0);
  TH1 *hist_Delta_phi_p_visible_leading_jet_1btag_corte_200 = new TH1F("Delta_phi_p_visible_leading_jet_1btag_corte_200","delta_phi_p_visible_leading_jet_1btag_corte_200",100,0.0,8.0);
  TH1 *hist_Delta_phi_p_visible_leading_jet_1btag_corte_300 = new TH1F("Delta_phi_p_visible_leading_jet_1btag_corte_300","delta_phi_p_visible_leading_jet_1btag_corte_300",100,0.0,8.0);
  TH1 *hist_Delta_phi_p_visible_leading_jet_1btag_corte_400 = new TH1F("Delta_phi_p_visible_leading_jet_1btag_corte_400","delta_phi_p_visible_leading_jet_1btag_corte_400",100,0.0,8.0);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Delta_eta_p_visible_leading_jet
  TH1 *hist_Delta_eta_p_visible_leading_jet = new TH1F("Delta_eta_p_visible_leading_jet","hist_delta_eta_p_visible_leading_jet",100,0.0,8.0);
  TH1 *hist_Delta_eta_p_visible_leading_jet_2btag_corte_200 = new TH1F("Delta_eta_p_visible_leading_jet_2btag_corte_200","hist_delta_eta_p_visible_leading_jet_2btag_corte_200",100,0.0,8.0);
  TH1 *hist_Delta_eta_p_visible_leading_jet_2btag_corte_300 = new TH1F("Delta_eta_p_visible_leading_jet_2btag_corte_300","hist_delta_eta_p_visible_leading_jet_2btag_corte_300",100,0.0,8.0);
  TH1 *hist_Delta_eta_p_visible_leading_jet_2btag_corte_400 = new TH1F("Delta_eta_p_visible_leading_jet_2btag_corte_400","hist_delta_eta_p_visible_leading_jet_2btag_corte_400",100,0.0,8.0);
  TH1 *hist_Delta_eta_p_visible_leading_jet_1btag_corte_200 = new TH1F("Delta_eta_p_visible_leading_jet_1btag_corte_200","hist_delta_eta_p_visible_leading_jet_1btag_corte_200",100,0.0,8.0);
  TH1 *hist_Delta_eta_p_visible_leading_jet_1btag_corte_300 = new TH1F("Delta_eta_p_visible_leading_jet_1btag_corte_300","hist_delta_eta_p_visible_leading_jet_1btag_corte_300",100,0.0,8.0);
  TH1 *hist_Delta_eta_p_visible_leading_jet_1btag_corte_400 = new TH1F("Delta_eta_p_visible_leading_jet_1btag_corte_400","hist_delta_eta_p_visible_leading_jet_1btag_corte_400",100,0.0,8.0);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  //MT
  TH1 *hist_MT_Muon_1btag_corte_200 = new TH1F("MT_Muon_1btag_corte_200","mt_Muon_1btag_corte_200",100,0.0,200.0);
  TH1 *hist_MT_Muon_1btag_corte_300 = new TH1F("MT_Muon_1btag_corte_300","mt_Muon_1btag_corte_300",100,0.0,200.0);
  TH1 *hist_MT_Muon_1btag_corte_400 = new TH1F("MT_Muon_1btag_corte_400","mt_Muon_1btag_corte_400",100,0.0,200.0);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1F *Cociente_PT_top_PT_leading_corte_200 =  new TH1F("Cociente_PT_top_PT_leading_corte_200","cociente_PT_top_PT_leading_corte_200",100,0.0,4.0);
  TH1F *Cociente_PT_top_PT_leading_corte_300 =  new TH1F("Cociente_PT_top_PT_leading_corte_300","cociente_PT_top_PT_leading_corte_300",100,0.0,4.0);
  TH1F *Cociente_PT_top_PT_leading_corte_400 =  new TH1F("Cociente_PT_top_PT_leading_corte_400","cociente_PT_top_PT_leading_corte_400",100,0.0,4.0);
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //corte 200-2 btag
  TH1 *hist_rt_2btag_corte_200 = new TH1F ("RT_2btag_corte_200","rt_2btag_corte_200",100,0.0,4.0);
  TH1 *hist_delta_phi_th_leading_jet_2btag_corte_200 = new TH1F("Delta_phi_th_leading_jet_2btag_corte_200","delta_phi_th_leading_jet_2btag_corte_200",100,0.0,7.0);
  TH1* hist_delta_eta_th_leading_jet_2btag_corte_200 = new TH1F("Delta_eta_th_leading_jet_2btag_corte_200","delta_eta_th_leading_jet_2btag_corte_200",100,0.0,8.0);
  TH1 *hist_rm_2btag_corte_200 = new TH1F("RM_2btag_corte_200","rm_2btag_corte_200",100,0.0,2.0);
  TH1 *hist_Delta_Phi_leading_jet_met_2btag_corte_200 = new TH1F("DeltaPhi_leading_jet_met_2btag_corte_200","deltaphi_leading_jet_met_2btag_corte_200",100,-7.0,7.0);
  TH1 *hist_Delta_Eta_leading_jet_met_2btag_corte_200 = new TH1F("DeltaEta_leading_jet_met_2btag_corte200","deltaeta_leading_jet_met_2btag_corte_200",50,0.0,8.0);
  TH1 *hist_MissingET_2btag_corte_200 = new TH1F("MissingET_2btag_corte_200","missinget_2btag_corte_200", 100,0.0,300.0);
  TH1 *hist_MT2_2btag_corte_200 = new TH1F("MT2_2btag_corte_200","mt2_2btag_corte_200",100,150,350.0);
  TH1 *hist_MT2W_2btag_corte_200 = new TH1F("MT2W_2btag_corte_200","mt2w_2btag_corte_200",100,70.0,400.0);
  TH1 *hist_MT2BL_2btag_corte_200 = new TH1F("MT2BL_2btag_corte_200","mt2bl_2btag_corte_200",100,0,300.0);
  TH1 *hist_Eta_average_2btag_corte_200 = new TH1F("Eta_average_2btag_corte_200","eta_average_2btag_corte_200",100,-7.0,7.0);
  TH2 *h_RM_Delta_Phi_leading_jet_met_2btag_corte_200 = new TH2F("RM_Delta_Phi_leading_jet_met_2btag_corte_200","rm_Delta_phi_leading_met_2btag_corte_200",200, 0.0, 1.5, 200, 0, 6.0);
  TH1 *hist_cociente_HT_nobjet_PT_leading_jet_2btag_corte_200 = new TH1F("Cociente_HT_nobjet_PT_leading_jet_2btag_corte_200","cociente_HT_nobjet_PT_leading_jet_2btag_corte_200", 100, 0.0, 3.0);
  TH1 *hist_cociente_HT_bjet_PT_leading_jet_2btag_corte_200 = new TH1F("Cociente_HT_bjet_PT_leading_jet_2btag_corte_200","cociente_HT_bjet_PT_leading_jet_2btag_corte_200", 100, 0.0, 2.0);
  TH1 *hist_cociente_PT_Electron_PT_leading_jet_2btag_corte_200 = new TH1F("Cociente_PT_electron_PT_leading_jet_2btag_corte_200","cociente_PT_electron_PT_leading_jet_2btag_corte_200",100,0.0,3.0);
  TH1 *hist_cociente_PT_Muon_PT_leading_jet_2btag_corte_200 = new TH1F("Cociente_PT_muon_PT_leading_jet_2btag_corte_200","cociente_PT_muon_PT_leading_jet_2btag_corte_200",100,0.0,3.0);

  //corte 300-2btag
  TH1 *hist_rt_2btag_corte_300 = new TH1F ("RT_2btag_corte_300","rt_2btag_corte_300",100,0.0,4.0);
  TH1 *hist_delta_phi_th_leading_jet_2btag_corte_300 = new TH1F("Delta_phi_th_leading_jet_2btag_corte_300","delta_phi_th_leading_jet_2btag_corte_300",100,0.0,7.0);
  TH1* hist_delta_eta_th_leading_jet_2btag_corte_300 = new TH1F("Delta_eta_th_leading_jet_2btag_corte_300","delta_eta_th_leading_jet_2btag_corte_300",100,0.0,8.0);
  TH1 *hist_rm_2btag_corte_300 = new TH1F("RM_2btag_corte_300","rm_2btag_corte_300",100,0.0,2.0);
  TH1 *hist_Delta_Phi_leading_jet_met_2btag_corte_300 = new TH1F("DeltaPhi_leading_jet_met_2btag_corte_300","deltaphi_leading_jet_met_2btag_corte_300",100,-7.0,7.0);
  TH1 *hist_Delta_Eta_leading_jet_met_2btag_corte_300 = new TH1F("DeltaEta_leading_jet_met_2btag_corte300","deltaeta_leading_jet_met_2btag_corte_300",50,0.0,8.0);
  TH1 *hist_MissingET_2btag_corte_300 = new TH1F("MissingET_2btag_corte_300","missinget_2btag_corte_300", 100,0.0,300.0);
  TH1 *hist_MT2_2btag_corte_300 = new TH1F("MT2_2btag_corte_300","mt2_2btag_corte_300",100,150,350.0);
  TH1 *hist_MT2W_2btag_corte_300 = new TH1F("MT2W_2btag_corte_300","mt2w_2btag_corte_300",100,70.0,400.0);
  TH1 *hist_MT2BL_2btag_corte_300 = new TH1F("MT2BL_2btag_corte_300","mt2bl_2btag_corte_300",100,0,300.0);
  TH1 *hist_Eta_average_2btag_corte_300 = new TH1F("Eta_average_2btag_corte_300","eta_average_2btag_corte_300",100,-7.0,7.0);
  TH2 *h_RM_Delta_Phi_leading_jet_met_2btag_corte_300 = new TH2F("RM_Delta_Phi_leading_jet_met_2btag_corte_300","rm_Delta_phi_leading_met_2btag_corte_300",200, 0.0, 1.5, 200, 0.0, 6.0);
  TH1 *hist_cociente_HT_nobjet_PT_leading_jet_2btag_corte_300 = new TH1F("Cociente_HT_nobjet_PT_leading_jet_2btag_corte_300","cociente_HT_nobjet_PT_leading_jet_2btag_corte_300", 100, 0.0, 3.0);
  TH1 *hist_cociente_HT_bjet_PT_leading_jet_2btag_corte_300 = new TH1F("Cociente_HT_bjet_PT_leading_jet_2btag_corte_300","cociente_HT_bjet_PT_leading_jet_2btag_corte_300", 100, 0.0, 2.0);
  TH1 *hist_cociente_PT_Electron_PT_leading_jet_2btag_corte_300 = new TH1F("Cociente_PT_electron_PT_leading_jet_2btag_corte_300","cociente_PT_electron_PT_leading_jet_2btag_corte_300",100,0.0,3.0);
  TH1 *hist_cociente_PT_Muon_PT_leading_jet_2btag_corte_300 = new TH1F("Cociente_PT_muon_PT_leading_jet_2btag_corte_300","cociente_PT_muon_PT_leading_jet_2btag_corte_300",100,0.0,3.0);

  //corte 400-2btag
  TH1 *hist_rt_2btag_corte_400 = new TH1F ("RT_2btag_corte_400","rt_2btag_corte_400",100,0.0,4.0);
  TH1 *hist_delta_phi_th_leading_jet_2btag_corte_400 = new TH1F("Delta_phi_th_leading_jet_2btag_corte_400","delta_phi_th_leading_jet_2btag_corte_400",100,0.0,7.0);
  TH1* hist_delta_eta_th_leading_jet_2btag_corte_400 = new TH1F("Delta_eta_th_leading_jet_2btag_corte_400","delta_eta_th_leading_jet_2btag_corte_400",100,0.0,8.0);
  TH1 *hist_rm_2btag_corte_400 = new TH1F("RM_2btag_corte_400","rm_2btag_corte_400",100,0.0,2.0);
  TH1 *hist_Delta_Phi_leading_jet_met_2btag_corte_400 = new TH1F("DeltaPhi_leading_jet_met_2btag_corte_400","deltaphi_leading_jet_met_2btag_corte_400",100,-7.0,7.0);
  TH1 *hist_Delta_Eta_leading_jet_met_2btag_corte_400 = new TH1F("DeltaEta_leading_jet_met_2btag_corte400","deltaeta_leading_jet_met_2btag_corte_400",50,0.0,8.0);
  TH1 *hist_MissingET_2btag_corte_400 = new TH1F("MissingET_2btag_corte_400","missinget_2btag_corte_400", 100,0.0,300.0);
  TH1 *hist_MT2_2btag_corte_400 = new TH1F("MT2_2btag_corte_400","mt2_2btag_corte_400",100,150,350.0);
  TH1 *hist_MT2W_2btag_corte_400 = new TH1F("MT2W_2btag_corte_400","mt2w_2btag_corte_400",100,70.0,400.0);
  TH1 *hist_MT2BL_2btag_corte_400 = new TH1F("MT2BL_2btag_corte_400","mt2bl_2btag_corte_400",100,0,300.0);
  TH1 *hist_Eta_average_2btag_corte_400 = new TH1F("Eta_average_2btag_corte_400","eta_average_2btag_corte_400",100,-7.0,7.0);
  TH2 *h_RM_Delta_Phi_leading_jet_met_2btag_corte_400 = new TH2F("RM_Delta_Phi_leading_jet_met_2btag_corte_400","rm_Delta_phi_leading_met_2btag_corte_400",200, 0.0, 1.5, 200, 0.0, 6.0);
  TH1 *hist_cociente_HT_nobjet_PT_leading_jet_2btag_corte_400 = new TH1F("Cociente_HT_nobjet_PT_leading_jet_2btag_corte_400","cociente_HT_nobjet_PT_leading_jet_2btag_corte_400", 100, 0.0, 3.0);
  TH1 *hist_cociente_HT_bjet_PT_leading_jet_2btag_corte_400 = new TH1F("Cociente_HT_bjet_PT_leading_jet_2btag_corte_400","cociente_HT_bjet_PT_leading_jet_2btag_corte_400", 100, 0.0, 2.0);
  TH1 *hist_cociente_PT_Electron_PT_leading_jet_2btag_corte_400 = new TH1F("Cociente_PT_electron_PT_leading_jet_2btag_corte_400","cociente_PT_electron_PT_leading_jet_2btag_corte_400",100,0.0,3.0);
  TH1 *hist_cociente_PT_Muon_PT_leading_jet_2btag_corte_400 = new TH1F("Cociente_PT_muon_PT_leading_jet_2btag_corte_400","cociente_PT_muon_PT_leading_jet_2btag_corte_400",100,0.0,3.0);

 //corte 200-1 btag
  TH1 *hist_rm_1btag_corte_200 = new TH1F("RM_1btag_corte_200","rm_1btag_corte_200",100,0.0,1.5);
  TH1 *hist_Delta_Phi_leading_jet_met_1btag_corte_200 = new TH1F("DeltaPhi_leading_jet_met_1btag_corte_200","deltaphi_leading_jet_met_1btag_corte_200",100,0.0,6.0);
  TH1 *hist_Delta_Eta_leading_jet_met_1btag_corte_200 = new TH1F("DeltaEta_leading_jet_met_1btag_corte200","deltaeta_leading_jet_met_1btag_corte_200",50,0.0,8.0);
  TH1 *hist_MissingET_1btag_corte_200 = new TH1F("MissingET_1btag_corte_200","missinget_1btag_corte_200", 100,0.0,600.0);
  TH1 *hist_Eta_average_1btag_corte_200 = new TH1F("Eta_average_1btag_corte_200","eta_average_1btag_corte_200",100,-7.0,7.0);
  TH1 *hist_cociente_HT_PT_leading_jet_1btag_corte_200 = new TH1F("Cociente_HT_PT_leading_jet_1btag_corte_200","cociente_HT_PT_leading_jet_1btag_corte_200", 100, 0.0,4.0);
  TH2 *h_RM_Delta_Phi_leading_jet_met_1btag_corte_200 = new TH2F("RM_Delta_Phi_leading_jet_met_1btag_corte_200","rm_Delta_phi_leading_met_1btag_corte_200",200, 0.0, 1.5, 200, 0.0, 6.0);
  TH2 *hist_cociente_HT_PT_leading_jet_RM_1btag_corte_200 = new TH2F("Cociente_HT_PT_leading_jet_RM_1btag_corte_200","cociente_HT_PT_leading_jet_RM_1btag_corte_200", 200, 0.0,4.0,200,0.0,1.5);
  TH1 *hist_cociente_HT_nobjet_PT_leading_jet_1btag_corte_200 = new TH1F("Cociente_HT_nobjet_PT_leading_jet_1btag_corte_200","cociente_HT_nobjet_PT_leading_jet_1btag_corte_200", 100, 0.0, 3.0);
  TH1 *hist_cociente_HT_bjet_PT_leading_jet_1btag_corte_200 = new TH1F("Cociente_HT_bjet_PT_leading_jet_1btag_corte_200","cociente_HT_bjet_PT_leading_jet_1btag_corte_200", 100, 0.0, 1.2);
  TH1 *hist_cociente_PT_Electron_PT_leading_jet_1btag_corte_200 = new TH1F("Cociente_PT_electron_PT_leading_jet_1btag_corte_200","cociente_PT_electron_PT_leading_jet_1btag_corte_200",100,0.0,3.0);
  TH1 *hist_cociente_PT_Muon_PT_leading_jet_1btag_corte_200 = new TH1F("Cociente_PT_muon_PT_leading_jet_1btag_corte_200","cociente_PT_muon_PT_leading_jet_1btag_corte_200",100,0.0,3.0);


  //corte 300-1btag
  TH1 *hist_rm_1btag_corte_300 = new TH1F("RM_1btag_corte_300","rm_1btag_corte_300",100,0.0,1.5);
  TH1 *hist_Delta_Phi_leading_jet_met_1btag_corte_300 = new TH1F("DeltaPhi_leading_jet_met_1btag_corte_300","deltaphi_leading_jet_met_1btag_corte_300",100,0.0,6.0);
  TH1 *hist_Delta_Eta_leading_jet_met_1btag_corte_300 = new TH1F("DeltaEta_leading_jet_met_1btag_corte300","deltaeta_leading_jet_met_1btag_corte_300",50,0.0,8.0);
  TH1 *hist_MissingET_1btag_corte_300 = new TH1F("MissingET_1btag_corte_300","missinget_1btag_corte_300", 100,0.0,600.0);
  TH1 *hist_Eta_average_1btag_corte_300 = new TH1F("Eta_average_1btag_corte_300","eta_average_1btag_corte_300",100,-7.0,7.0);
  TH1 *hist_cociente_HT_PT_leading_jet_1btag_corte_300 = new TH1F("Cociente_HT_PT_leading_jet_1btag_corte_300","cociente_HT_PT_leading_jet_1btag_corte_300", 100, 0.0,4.0);
  TH2 *h_RM_Delta_Phi_leading_jet_met_1btag_corte_300 = new TH2F("RM_Delta_Phi_leading_jet_met_1btag_corte_300","rm_Delta_phi_leading_met_1btag_corte_300",200, 0.0, 1.5, 100, 0.0, 6.0);
  TH2 *hist_cociente_HT_PT_leading_jet_RM_1btag_corte_300 = new TH2F("Cociente_HT_PT_leading_jet_RM_1btag_corte_300","cociente_HT_PT_leading_jet_RM_1btag_corte_300", 200, 0.0,4.0,200,0.0,1.5);
  TH1 *hist_cociente_HT_nobjet_PT_leading_jet_1btag_corte_300 = new TH1F("Cociente_HT_nobjet_PT_leading_jet_1btag_corte_300","cociente_HT_nobjet_PT_leading_jet_1btag_corte_300", 100, 0.0, 3.0);
  TH1 *hist_cociente_HT_bjet_PT_leading_jet_1btag_corte_300 = new TH1F("Cociente_HT_bjet_PT_leading_jet_1btag_corte_300","cociente_HT_bjet_PT_leading_jet_1btag_corte_300", 100, 0.0, 1.2);
  TH1 *hist_cociente_PT_Electron_PT_leading_jet_1btag_corte_300 = new TH1F("Cociente_PT_electron_PT_leading_jet_1btag_corte_300","cociente_PT_electron_PT_leading_jet_1btag_corte_300",100,0.0,3.0);
  TH1 *hist_cociente_PT_Muon_PT_leading_jet_1btag_corte_300 = new TH1F("Cociente_PT_muon_PT_leading_jet_1btag_corte_300","cociente_PT_muon_PT_leading_jet_1btag_corte_300",100,0.0,3.0);


  //corte 400-1btag
  TH1 *hist_rm_1btag_corte_400 = new TH1F("RM_1btag_corte_400","rm_1btag_corte_400",100,0.0,1.5);
  TH1 *hist_Delta_Phi_leading_jet_met_1btag_corte_400 = new TH1F("DeltaPhi_leading_jet_met_1btag_corte_400","deltaphi_leading_jet_met_1btag_corte_400",100,0.0,6.0);
  TH1 *hist_Delta_Eta_leading_jet_met_1btag_corte_400 = new TH1F("DeltaEta_leading_jet_met_1btag_corte400","deltaeta_leading_jet_met_1btag_corte_400",50,0.0,8.0);
  TH1 *hist_MissingET_1btag_corte_400 = new TH1F("MissingET_1btag_corte_400","missinget_1btag_corte_400", 100,0.0,600.0);
  TH1 *hist_Eta_average_1btag_corte_400 = new TH1F("Eta_average_1btag_corte_400","eta_average_1btag_corte_400",100,-7.0,7.0);
  TH1 *hist_cociente_HT_PT_leading_jet_1btag_corte_400 = new TH1F("Cociente_HT_PT_leading_jet_1btag_corte_400","cociente_HT_PT_leading_jet_1btag_corte_400", 100, 0.0,4.0);
  TH2 *h_RM_Delta_Phi_leading_jet_met_1btag_corte_400 = new TH2F("RM_Delta_Phi_leading_jet_met_1btag_corte_400","rm_Delta_phi_leading_met_1btag_corte_400",200, 0.0, 1.5, 200, 0.0, 6.0);
  TH2 *hist_cociente_HT_PT_leading_jet_RM_1btag_corte_400 = new TH2F("Cociente_HT_PT_leading_jet_RM_1btag_corte_400","cociente_HT_PT_leading_jet_RM_1btag_corte_400", 100, 0.0,4.0,100,0.0,1.5);
  TH1 *hist_cociente_HT_nobjet_PT_leading_jet_1btag_corte_400 = new TH1F("Cociente_HT_nobjet_PT_leading_jet_1btag_corte_400","cociente_HT_nobjet_PT_leading_jet_1btag_corte_400", 100, 0.0, 3.0);
  TH1 *hist_cociente_HT_bjet_PT_leading_jet_1btag_corte_400 = new TH1F("Cociente_HT_bjet_PT_leading_jet_1btag_corte_400","cociente_HT_bjet_PT_leading_jet_1btag_corte_400", 100, 0.0, 1.2);
  TH1 *hist_cociente_PT_Electron_PT_leading_jet_1btag_corte_400 = new TH1F("Cociente_PT_electron_PT_leading_jet_1btag_corte_400","cociente_PT_electron_PT_leading_jet_1btag_corte_400",100,0.0,3.0);
  TH1 *hist_cociente_PT_Muon_PT_leading_jet_1btag_corte_400 = new TH1F("Cociente_PT_muon_PT_leading_jet_1btag_corte_400","cociente_PT_muon_PT_leading_jet_1btag_corte_400",100,0.0,3.0);


  //Jets counters
  int bjet_cntr=0;
  int wjet_cntr=0;
  int one_electron=0;
  int one_muon=0;
  int more_electron=0;
  int more_muon=0;
  int more_muon_electron=0;
  int bjet_counter=0;
  int contador[50]={50*0};

  //Scan counters
  int nscans=40;
  double MET_min, MET_max, MET_2btag_max,RM_min, RM_max;
  double RT_min, RT_max, MT2_min, MT2_max;
  double MT_min, MT_max;
  double Delta_Phi_leading_jet_met_min,Delta_Phi_leading_jet_met_max;
  double Cociente_p_visible_leading_jet_pt_min,Cociente_p_visible_leading_jet_pt_max;
  double Cociente_HT_bjet_PT_leading_jet_min, Cociente_HT_bjet_PT_leading_jet_max;
  double  PT_leading_jet_min, PT_leading_jet_max;
  double valuev,deltav;
  int MET_counter[100]={100*0};
  int RM_counter[100]={100*0};
  int RT_counter[100]={100*0};
  int MT2_counter[100]={100*0};
  int MT_counter[100]={100*0};
  int PT_leading_jet_counter[100]={100*0};
  int Delta_Phi_leading_jet_met_counter[100]={100*0};
  int Cociente_p_visible_leading_jet_pt_counter[100]={100*0};
  int Cociente_HT_bjet_PT_leading_jet_counter[100]={100*0};
  // Loop over all events
  int nevents=0;
  for ( Int_t entry = 0; entry < numberOfEntries; ++entry ){
    //Primer contador
    contador[0]++;
    // Load selected branches with data from specified event
    treeReader_Delphes->ReadEntry(entry);
  
    //Number of jets, muons and electrons in the event
    
    int NumJets = branchJet->GetEntries();
    histNumJets->Fill(NumJets);
    int NumElectrons = branchElectron->GetEntries();
    histNumElectrons->Fill(NumElectrons);
    int NumMuons = branchMuon->GetEntries();
    histNumMuons->Fill(NumMuons);
    nevents++;
    //Topology Analysis in Delphes
    MissingET *metpointer;
    Muon *Muonpointer;
    Electron *Electronpointer;
    Jet *Jetpointer1;
    Jet *Jetpointer2;
    Jet *Jetpointer3;
    Jet *Jetpointerb;
    Jet *Jetpointerbl;
    Jet *leading_jet;
    Jet *Jet_i;
    Jet *i_Jet;
    int jjet =0;
    int kjet =0;
    int bjet[5]={5*0};
    int wjet[5]={5*0};
    int p1w[5]={5*0};
    int p2w[5]={5*0};
    int p_bjet1[5]={5*0};
    double PT_Muon;
    double PT_Electron;
    double PT_leading_jet;
    double transverse_mass_Muon=0;
    double transverse_mass_Electron=0;
    double transverse_mass_Lepton=0;
    double transverse_mass_square_Muon;
    double transverse_mass_square_Electron;
    double Cociente_PT_Muon_PT_leading_jet;
    double Delta_Eta_Muon_leading_jet;
    double Delta_Phi_Muon_leading_jet;
    double Cociente_PT_Electron_PT_leading_jet;
    double Delta_Eta_Electron_leading_jet;
    double Delta_Phi_Electron_leading_jet;
    double Cociente_MissingET_PT_leading_jet;
    double Delta_Eta_leading_jet_met;
    double Delta_Phi_leading_jet_met;
    double theta_lepton, px_lepton, py_lepton, pz_lepton, E_lepton;
    double Delta_Phi_leading_jet_muon,  Delta_Phi_leading_jet_electron;
    double Delta_Phi_electron_met, Delta_Phi_muon_met;

    Jet *wjet1 = new Jet;
    Jet *wjet2 = new Jet;
    Jet *bjet1 = new Jet;
    Jet *bjet2 = new Jet;
    double max_j_pt = 200.;
    double rm;
    bool  esolo=false;
    bool musolo=false;
    double pl[4];           // lepton cuadrimomentum(MT2W, MT2BL)
    double pb1[4];         // bjet1 cuadrimomentum  (MT2W, MT2BL)
    double pb2[4];        // bjet2 cuadrimomentum   (MT2W, MT2BL)
    double pmiss[3];     // MissingEt cuadrimomentum(MT2W, MT2BL)
    double pa[3]={0};   //                            (MT2)
    double pb[3]={0};  //                             (MT2)
    double pth[4];    // Cuadrimomentum of top hadronico
    double ptl[2];
    double ptw[2];
    //Topology analysis
    
    if ((NumElectrons == 1) && (NumMuons == 0) ) {
      esolo=true;
      one_electron++;
    }
    if ((NumElectrons == 0) && (NumMuons == 1)) {
      musolo=true;
      one_muon++;
    }
    if ((NumElectrons > 1) && (NumMuons == 0) ) {
      more_electron++;
    }
    if ((NumElectrons == 0) && (NumMuons > 1) ) {
      more_muon++;
    }
    if ((NumElectrons > 1) || (NumMuons > 1) ) {
      more_muon_electron++;
    }
    
    /*  
    metpointer = (MissingET*)branchMissingET->At(0);
    leading_jet = (Jet*)branchJet->At(0);
    Muonpointer = (Muon*)branchMuon->At(0);
    Electronpointer = (Electron*)branchElectron->At(0);
    if(NumJets<0) continue;
    contador[1]++;
    if(metpointer<0) continue;
    contador[2]++;
    if((NumElectrons==0) && (NumMuons==0)||(NumElectrons>0) && (NumMuons>0)) continue;
    contador[3]++;
    if ((esolo&&Electronpointer->PT<30)||(musolo&&Muonpointer->PT<30)) continue;
    contador[4]++;
    */
    if((NumJets>= 0) && (branchMissingET->GetEntries()>0)&&((NumElectrons>0 &&NumMuons==0)||(NumElectrons==0 && NumMuons>0)))  {
      contador[1]++;
      metpointer = (MissingET*)branchMissingET->At(0);
      leading_jet = (Jet*)branchJet->At(0);
      Muonpointer = (Muon*)branchMuon->At(0);
      Electronpointer = (Electron*)branchElectron->At(0);
    
      if (esolo&&Electronpointer->PT<30) continue;
      if (musolo&&Muonpointer->PT<30) continue;
      contador[2]++;

      //Eta Average and HT/PT leading
      double eta_visible;
      double phi_visible;
      double theta_visible;
      double pvisible;
      double p_lepton;
      double sum=0.0;
      double HT=0.0;
      double ST=0.0;
      double rST=0.0;
      double HT_nobjet=0.0;
      double HT_bjet=0.0;
      double num_bjets=0.0;
      double delta_phi_visible_leading_jet;
      double delta_eta_visible_leading_jet;
      bjet_cntr=0;
      double p_visible[3]={3*0};
      double p_event[3]={3*0};
      double pt_event;
      double delta_eta_ave=0;
      
      for(int i=0; i<NumJets; i++){
	Jet_i = (Jet*)branchJet->At(i);
	
	// add all visible momenta
	if (i>0) {
	  p_visible[0]+=Jet_i->PT*cos(Jet_i->Phi+pi);
	  p_visible[1]+=Jet_i->PT*sin(Jet_i->Phi+pi);
	  double theta_i=2*atan(exp(-Jet_i->Eta));
	  p_visible[2]+=Jet_i->PT*cos(theta_i)/sin(theta_i);
	}
	

	sum+=fabs((leading_jet->Eta)-(Jet_i->Eta));
	
	//No leading_jet
	if(i>0){
	  HT += Jet_i->PT;
	}
	//No leading_jet and bjet
	if(i>0 && Jet_i->BTag==0){
	  HT_nobjet += Jet_i->PT;
	}
	//No leading_jet and only considering bjet
	if (Jet_i->BTag>0&&i>0){
	  num_bjets++;
	  HT_bjet += Jet_i->PT;
	}
	//if (num_bjets>0) HT_bjet=HT_bjet/(num_bjets);
	
      } //end Jet_i cycle

      for(int i=1; i<NumJets; i++){
	i_Jet = (Jet*)branchJet->At(i);
	delta_eta_ave+=fabs(i_Jet->Eta-leading_jet->Eta);

      }

      delta_eta_ave = delta_eta_ave/(NumJets-1);
      hist_cociente_HT_PT_leading_jet->Fill(HT/leading_jet->PT);
      hist_cociente_HT_nobjet_PT_leading_jet->Fill(HT_nobjet/leading_jet->PT);
      hist_cociente_HT_bjet_PT_leading_jet->Fill(HT_bjet/leading_jet->PT);
      
      sum = sum/(NumJets-1);
      hist_Eta_average->Fill(sum);
      
      //Rm Delta phi Delta Eta
      rm = metpointer->MET/leading_jet->PT;
      Delta_Phi_leading_jet_met = fabs((leading_jet->Phi)-(metpointer->Phi));
      Delta_Eta_leading_jet_met = fabs((leading_jet->Eta)-(metpointer->Eta));
      
      //MT Variables for muon and electron
      if(musolo ){
	transverse_mass_Lepton = TMath::Sqrt(Muonpointer->PT * metpointer->MET * (1 - TMath::Cos(Muonpointer->Phi - metpointer->Phi)));
	Delta_Phi_leading_jet_muon = (leading_jet->Phi+pi)-(Muonpointer->Phi+pi);
	hist_Delta_Phi_leading_jet_muon_2btag->Fill(Delta_Phi_leading_jet_muon);
	hist_cociente_PT_Muon_PT_leading_jet->Fill(Muonpointer->PT/leading_jet->PT);
	
	//Relevant variables to MT2W and MT2BL (Muon)
	theta_lepton = 2*TMath::ATan(TMath::Exp(-Muonpointer->Eta));
	px_lepton = Muonpointer->PT*TMath::Cos(Muonpointer->Phi+pi);
	py_lepton = Muonpointer->PT*TMath::Sin(Muonpointer->Phi+pi);
	pz_lepton = Muonpointer->PT*TMath::Cos(theta_lepton)/TMath::Sin(theta_lepton);
	E_lepton = TMath::Sqrt(TMath::Power((Muonpointer->PT),2) + TMath::Power(pz_lepton,2) + TMath::Power(m_muon,2));
	p_lepton=sqrt(px_lepton*px_lepton+py_lepton*py_lepton+pz_lepton*pz_lepton);
	ST=HT+Muonpointer->PT;
	rST=ST/leading_jet->PT;

	
	//Add p_lepton to p_visible
	p_visible[0]+=px_lepton;
	p_visible[1]+=py_lepton;
	p_visible[2]+=pz_lepton;

	//Add p_lepton to p_event
	p_event[0]+=px_lepton;
	p_event[1]+=py_lepton;
	p_event[2]+=pz_lepton;
	
      }

      if(esolo){
	transverse_mass_Lepton = TMath::Sqrt(Electronpointer->PT * metpointer->MET * (1 - TMath::Cos(Electronpointer->Phi - metpointer->Phi)));
	Delta_Phi_leading_jet_electron = (leading_jet->Phi)-(Electronpointer->Phi);
	hist_Delta_Phi_leading_jet_electron_2btag->Fill(Delta_Phi_leading_jet_electron);
        hist_cociente_PT_Electron_PT_leading_jet->Fill(Electronpointer->PT/leading_jet->PT);
	
	//Relevant variables to MT2W and MT2BL (Electron)
	theta_lepton = 2*TMath::ATan(TMath::Exp(-Electronpointer->Eta));
	px_lepton = Electronpointer->PT*TMath::Cos(Electronpointer->Phi+pi);
	py_lepton = Electronpointer->PT*TMath::Sin(Electronpointer->Phi+pi);
	pz_lepton = Electronpointer->PT*TMath::Cos(theta_lepton)/TMath::Sin(theta_lepton);
	E_lepton = TMath::Sqrt(TMath::Power((Electronpointer->PT),2) + TMath::Power(pz_lepton,2) + TMath::Power( m_electron,2));
	p_lepton=sqrt(px_lepton*px_lepton+py_lepton*py_lepton+pz_lepton*pz_lepton);
	ST=HT+Electronpointer->PT;
	rST=ST/leading_jet->PT;
	//Add p_lepton to p_visible
	p_visible[1]+=px_lepton;
	p_visible[2]+=py_lepton;

	//Add p_lepton to p_event
	p_event[1]+=px_lepton;
	p_event[2]+=py_lepton;
	
	
      }
      hist_MT_lepton->Fill(transverse_mass_Lepton);

      
      if ((esolo||musolo)&&ST>20) hist_ST->Fill(ST);
      //Magnitud of p_visible
      //pvisible = sqrt(p_visible[0]*p_visible[0]+p_visible[1]*p_visible[1]+p_visible[2]*p_visible[2]);
      pvisible = sqrt(p_visible[0]*p_visible[0]+p_visible[1]*p_visible[1]);

      hist_Cociente_p_visible_leading_jet_pt->Fill(pvisible/leading_jet->PT);
      
      //Phi_visible
      phi_visible=acos(p_visible[0]/sqrt(p_visible[0]*p_visible[0]+p_visible[1]*p_visible[1]));
      if (p_visible[1]<0) phi_visible=2*pi-phi_visible;
      delta_phi_visible_leading_jet = fabs(phi_visible-(leading_jet->Phi+pi));
      hist_Delta_phi_p_visible_leading_jet->Fill(delta_phi_visible_leading_jet);

      //Eta_visible
      theta_visible=acos(pth[2]/sqrt(pth[0]*pth[0]+pth[1]*pth[1]+pth[2]*pth[2]));
      eta_visible = -log(tan((theta_visible)/2.0));
      delta_eta_visible_leading_jet = fabs(eta_visible-leading_jet->Eta);
      hist_Delta_eta_p_visible_leading_jet->Fill(delta_eta_visible_leading_jet);

      //Saving information of lepton(MT2W, MT2BL)
      pl[0]=E_lepton;
      pl[1]=px_lepton;
      pl[2]=py_lepton;
      pl[3]=pz_lepton;
      
      //Saving information of MissingET (MT2W, MT2BL)
      pmiss[0]=0;
      pmiss[1]= metpointer->MET*cos(metpointer->Phi);
      pmiss[2]= metpointer->MET*sin(metpointer->Phi);
      
      for(int ijet=0; ijet<NumJets; ijet++){
	// Checking type of jets
	// Bjet identification
	Jetpointer1 = (Jet*)branchJet->At(ijet);
	if (Jetpointer1->BTag>0){
	  bjet[bjet_cntr]= ijet;
	  bjet_cntr++;
	}
      }
      histNum_btags->Fill(bjet_cntr);
      //Definition of some variables
      
      double theta_bjet1,px_bjet1,py_bjet1,pz_bjet1,e_bjet1;
      double theta_bjet2,px_bjet2,py_bjet2,pz_bjet2,e_bjet2;
      double chi_dijet,chi_dijet1, chi_dijet2;
      double  p2_j1_j2_b1, p2_j1_j2_b2,e_j1_j2_b1,e_j1_j2_b2, p2_j1_j2_j3;
      
      double pxj1, pyj1, pzj1, ej1;
      double pxj2, pyj2, pzj2, ej2;
      double pxj3, pyj3, pzj3, ej3;

      double e_j1_j2, p2_j1_j2, M_12, M_12b1, M_12b2, M_12b, e_j1_j2_j3, M_123;
      double theta_Jetpointer1;
      double theta_Jetpointer2;
      double theta_Jetpointer3;

      double theta_Jetpointerb;
      double jw1, jw2, jbh, jbl;
      bool found_hbranch;
      double thetab, pxb, pyb, pzb, ejb, eb, ej12b, p2_j12b;
      double bj;
      float sigma_12b = 20.0;
      double px_top1, py_top1, px_top2, py_top2;
      double rt, delta_phi_th_leading_jet, delta_eta_th_leading_jet;
      double theta1, theta2;
      double phi_th, eta_th;
      double PT_top, mt_top;
      double MT_lb;
      double cociente_PT_top_PT_leading, PT_top_lep;
      double phi_tl, phi_lb;
      double phi_bjethadronic, phi_wjets;
      // Now we find the b from the hadronic branch
      if (bjet_cntr==2){
	//contador[5];
	found_hbranch=false;
	bjet1 = (Jet*)branchJet->At(bjet[0]);
	bjet2 = (Jet*)branchJet->At(bjet[1]);
	e_bjet1 = CalculateE(bjet1->Eta, bjet1->PT, bjet1->Mass);
	e_bjet2 = CalculateE(bjet2->Eta, bjet2->PT, bjet2->Mass);

	if ((bjet1->PT<30)||(bjet2->PT<30)) continue;  // Demand  a threshold on bjet PT
	contador[3]++;

	////////////////////////////////////////////////////////////
	//Saving information of bjet1 and bjet2 (MT2W, MT2BL)
	
	theta_bjet1 = 2*TMath::ATan(TMath::Exp(-bjet1->Eta));
	pb1[1] = bjet1->PT*TMath::Cos(bjet1->Phi);
	pb1[2] = bjet1->PT*TMath::Sin(bjet1->Phi);
	pb1[3] = bjet1->PT*TMath::Cos(theta_bjet1)/sin(theta_bjet1);
	pb1[0] = TMath::Sqrt(bjet1->PT*bjet1->PT+pb1[3]*pb1[3]+bjet1->Mass);
	
	theta_bjet2 = 2*TMath::ATan(TMath::Exp(-bjet2->Eta));
	pb2[1] = bjet2->PT*TMath::Cos(bjet2->Phi);
	pb2[2] = bjet2->PT*TMath::Sin(bjet2->Phi);
	pb2[3] = bjet2->PT-TMath::Cos(theta_bjet2)/sin(theta_bjet2);
	pb2[0] = TMath::Sqrt(bjet2->PT*bjet2->PT+pb2[3]*pb2[3]+bjet2->Mass);
	
	float chi_tw = 9999999999.99;
	float chi_w1 = 9999999999.99;
	for(int ijet=0; ijet<NumJets; ijet++){
	  Jetpointer1 = (Jet*)branchJet->At(ijet);
          if(Jetpointer1->BTag==0) {
	    for (int kjet=0; kjet<NumJets; kjet++){
	      Jetpointer2 = (Jet*)branchJet->At(kjet);
	      if((kjet!=ijet)&& (Jetpointer2->BTag==0)){
		
		//W Invariant mass
		ej1 = CalculateE(Jetpointer1->Eta, Jetpointer1->PT, Jetpointer1->Mass);
	        ej2 = CalculateE(Jetpointer2->Eta, Jetpointer2->PT, Jetpointer2->Mass);
                e_j1_j2 = ej1+ej2;
		p2_j1_j2 = CalculateP12_square_2particle(Jetpointer1->Eta, Jetpointer1->Phi, Jetpointer1->PT, Jetpointer2->Eta, Jetpointer2->Phi, Jetpointer2->PT);
		M_12 = TMath::Sqrt(TMath::Power((e_j1_j2),2)-p2_j1_j2);
		chi_dijet = TMath::Power((M_12-mw),2);
		
		// also reconstruct simultaneously the top mass from hadronic decay
		
		p2_j1_j2_b1 = CalculateP123_square_3particle(Jetpointer1->Eta, Jetpointer1->Phi, Jetpointer1->PT, Jetpointer2->Eta, Jetpointer2->Phi, Jetpointer2->PT, bjet1->Eta, bjet1->Phi, bjet1->PT);
		p2_j1_j2_b2 = CalculateP123_square_3particle(Jetpointer1->Eta, Jetpointer1->Phi, Jetpointer1->PT, Jetpointer2->Eta, Jetpointer2->Phi, Jetpointer2->PT, bjet2->Eta, bjet2->Phi, bjet2->PT);
		
		e_j1_j2_b1 = ej1+ej2+e_bjet1;
		e_j1_j2_b2 = ej1+ej2+e_bjet2;
		M_12b1= TMath::Sqrt(TMath::Power((e_j1_j2_b1),2)-p2_j1_j2_b1);
      		M_12b2= TMath::Sqrt(TMath::Power((e_j1_j2_b2),2)-p2_j1_j2_b2);
		chi_dijet1=chi_dijet+TMath::Power((M_12b1-mt),2);
		chi_dijet2=chi_dijet+TMath::Power((M_12b2-mt),2);
		
		if(chi_dijet1<chi_tw){
		  found_hbranch=true;
		  chi_tw=chi_dijet1;
		  //store jet from w
		  jw1=ijet;
		  jw2=kjet;
		  jbh=bjet[0];
		  jbl=bjet[1];
		}

		else if(chi_dijet2<chi_tw){
		  found_hbranch=true;
		  chi_tw=chi_dijet2;
		  //store jet from w
		  jw1=ijet;
		  jw2=kjet;
		  jbh=bjet[1];
		  jbl=bjet[0];
		}
		
	      }
	    }
	  }
	}


	// here we know the jets coming from w
	//plot mW
	if (found_hbranch){
	  Jetpointer1 = (Jet*)branchJet->At(jw1);
	  Jetpointer2 = (Jet*)branchJet->At(jw2);
	  Jetpointerb = (Jet*)branchJet->At(jbh);
 	  Jetpointerbl = (Jet*)branchJet->At(jbl);
	  
	  // Relevant information to p_event
	  p_event[1]=(Jetpointer1->PT*cos(Jetpointer1->Phi+pi))+(Jetpointer2->PT*cos(Jetpointer2->Phi+pi))+(Jetpointerb->PT*cos(Jetpointerb->Phi+pi))+(Jetpointerbl->PT*cos(Jetpointerbl->Phi+pi));

	  p_event[2]=(Jetpointer1->PT*sin(Jetpointer1->Phi+pi))+(Jetpointer2->PT*sin(Jetpointer2->Phi+pi))+(Jetpointerb->PT*sin(Jetpointerb->Phi+pi))+(Jetpointerbl->PT*sin(Jetpointerbl->Phi+pi));

	  pt_event = sqrt(p_event[1]*p_event[1]+p_event[2]*p_event[2]);
      	  
	  // MT_lb
	  if (musolo){
	    ptl[1]=(Jetpointerbl->PT*cos(Jetpointerbl->Phi+pi))+(Muonpointer->PT*cos(Muonpointer->Phi+pi));
	    ptl[2]=(Jetpointerbl->PT*sin(Jetpointerbl->Phi+pi))+(Muonpointer->PT*sin(Muonpointer->Phi+pi));
	    
	  }
	  if(esolo){
	    ptl[1]=(Jetpointerbl->PT*cos(Jetpointerbl->Phi+pi))+(Electronpointer->PT*cos(Electronpointer->Phi+pi));
	    ptl[2]=(Jetpointerbl->PT*sin(Jetpointerbl->Phi+pi))+(Electronpointer->PT*sin(Electronpointer->Phi+pi));
	    
	  }
	  
	  phi_lb=acos(ptl[1]/sqrt(ptl[1]*ptl[1]+ptl[2]*ptl[2]));
	  double PT_lb = sqrt(ptl[1]*ptl[1]+ptl[2]*ptl[2]);
	  if (ptl[2]<0) phi_lb=2*pi-fabs(phi_lb);
	  MT_lb = sqrt(2*metpointer->MET*PT_lb*(1-cos(phi_lb-(metpointer->Phi+pi))));
	  hist_MT_lb_2btag->Fill(MT_lb);
	  ej1 = CalculateE(Jetpointer1->Eta, Jetpointer1->PT, Jetpointer1->Mass);
	  ej2 = CalculateE(Jetpointer2->Eta, Jetpointer2->PT, Jetpointer2->Mass);
	  e_j1_j2 = ej1+ej2;
	  p2_j1_j2 = CalculateP12_square_2particle(Jetpointer1->Eta, Jetpointer1->Phi, Jetpointer1->PT, Jetpointer2->Eta, Jetpointer2->Phi, Jetpointer2->PT);
	  M_12 = TMath::Sqrt(TMath::Power((e_j1_j2),2)-p2_j1_j2);
	  ejb = CalculateE(Jetpointerb->Eta, Jetpointerb->PT, Jetpointerb->Mass);
	  ej12b = ej1+ej2+ejb;
	  p2_j12b = CalculateP123_square_3particle(Jetpointer1->Eta, Jetpointer1->Phi, Jetpointer1->PT, Jetpointer2->Eta, Jetpointer2->Phi, Jetpointer2->PT, Jetpointerb->Eta, Jetpointerb->Phi, Jetpointerb->PT);
	  
	  M_12b = TMath::Sqrt(TMath::Power((ej12b),2)-p2_j12b);
	  hist_w_recons_2bjet->Fill(M_12);
	  hist_t_recons_2bjet->Fill(M_12b);

	  //Reconstruction of information wjets
	  ptw[1]=(Jetpointer1->PT*cos(Jetpointer1->Phi+pi))+(Jetpointer2->PT*cos(Jetpointer2->Phi+pi));
	  ptw[2]=(Jetpointer1->PT*sin(Jetpointer1->Phi+pi))+(Jetpointer2->PT*sin(Jetpointer2->Phi+pi));
	  phi_wjets = acos(ptw[1]/sqrt(ptw[1]*ptw[1]+ptw[2]*ptw[2]));
          phi_bjethadronic = Jetpointerb->Phi+pi;

	  //Reconstruction of information of top leptonic
	  PT_top_lep = sqrt(ptl[1]*ptl[1]+ptl[2]*ptl[2]);
	  phi_tl=acos(ptl[1]/sqrt(ptl[1]*ptl[1]+ptl[2]*ptl[2]));
	  if (ptl[2]<0) phi_tl=2*pi-fabs(phi_tl);


	  //Reconstruction of information of top hadronic
	  pth[0]=ej1+ej2+ejb;
	  pth[1]=(Jetpointer1->PT*cos(Jetpointer1->Phi+pi))+(Jetpointer2->PT*cos(Jetpointer2->Phi+pi))+(Jetpointerb->PT*cos(Jetpointerb->Phi+pi));
	  pth[2]=(Jetpointer1->PT*sin(Jetpointer1->Phi+pi))+(Jetpointer2->PT*sin(Jetpointer2->Phi+pi))+(Jetpointerb->PT*sin(Jetpointerb->Phi+pi));
	  double theta1=2*atan(exp(-Jetpointer1->Eta));
	  double theta2=2*atan(exp(-Jetpointer2->Eta));
	  double thetab=2*atan(exp(-Jetpointerb->Eta));
	  pth[3]= (Jetpointer1->PT*cos(theta1)/sin(theta1))+(Jetpointer2->PT*cos(theta2)/sin(theta2))+(Jetpointerb->PT*cos(thetab)/sin(thetab));
	  rt=pth[0]/leading_jet->PT;
	  double phi_th=acos(pth[1]/sqrt(pth[1]*pth[1]+pth[2]*pth[2]));
	  
	  if (pth[2]<0) phi_th=2*pi-fabs(phi_th);
	  
	  double theta_th=acos(pth[3]/sqrt(pth[1]*pth[1]+pth[2]*pth[2]+pth[3]*pth[3]));
	  double eta_th = -log(tan((theta_th)/2.0));
	  delta_phi_th_leading_jet = fabs(phi_th-leading_jet->Phi-pi);
	  delta_eta_th_leading_jet = fabs(eta_th-leading_jet->Eta);
	  hist_rt_2btag->Fill(rt);
	  hist_delta_phi_th_leading_jet_2btag->Fill(delta_phi_th_leading_jet);
	  hist_delta_eta_th_leading_jet_2btag->Fill(delta_eta_th_leading_jet);
	  
	  // Top transverse mass
	  PT_top = sqrt(pth[1]*pth[1]+pth[2]*pth[2]);
	  cociente_PT_top_PT_leading = PT_top/leading_jet->PT;
	  mt_top = sqrt(2*metpointer->MET*PT_top*(1-cos(phi_th-(metpointer->Phi+pi))));
	  hist_MT_top_2btag->Fill(mt_top);
	  
	  //Relevant information to MT2
	  px_top1 = bjet1->PT*cos(bjet1->Phi)+(Jetpointer2->PT*cos(Jetpointer2->Phi))+(Jetpointer1->PT*cos(Jetpointer1->Phi));
	  py_top1 = bjet1->PT*sin(bjet1->Phi)+(Jetpointer2->PT*sin(Jetpointer2->Phi))+(Jetpointer1->PT*sin(Jetpointer1->Phi));
	  px_top2 = bjet2->PT*cos(bjet2->Phi)+px_lepton;
	  py_top2 = bjet2->PT*sin(bjet2->Phi)+py_lepton;
	  
	  // Saving information of pa and pb to MT2
	  pa[0] = mt;
	  pa[1] = px_top1;
	  pa[2] = py_top1;
	    
	  pb[0] = mt;
	  pb[1] = px_top2;
	  pb[2] = py_top2;
	  
	}
	
	// Calcula MT2
	mt2_bisect::mt2 mt2_event;
	mt2_event.set_momenta(pa,pb,pmiss);
	mt2_event.set_mn(mn);
	double value_mt2=mt2_event.get_mt2();
	hist_MT2->Fill(value_mt2);
	
	// Calcula MT2W
	mt2w_bisect::mt2w mt2w_event;
	mt2w_event.set_momenta(pl,pb1,pb2,pmiss);
	double value_mt2w=mt2w_event.get_mt2w();
	hist_MT2W->Fill(value_mt2w);
	
	//Calcula MT2BL
	mt2bl_bisect::mt2bl mt2bl_event;
	mt2bl_event.set_momenta(pl,pb1,pb2,pmiss);
	double value_mt2bl=mt2bl_event.get_mt2bl();
	hist_MT2BL->Fill(value_mt2bl);
	
	//now we try to find the b jet from top hadronic decay
	//for this we reconstruct a 3 jet invariant mass
	// Some cuts
	//if((rm>0.5)&& (metpointer->MET>200)&&(Delta_Phi_leading_jet_met<pi)){
	    if (leading_jet->PT>200){
	      if (esolo||musolo) hist_Cociente_pt_visible_st_2btag_corte_200->Fill(pvisible/ST);
	      hist_Cociente_met_raiz_ht_2btag_corte_200->Fill(metpointer->MET/sqrt(HT));
	      hist_MT2BL_2btag_corte_200->Fill(value_mt2bl);
	      if(musolo)hist_Delta_Phi_leading_jet_muon_2btag_corte_200->Fill(Delta_Phi_leading_jet_muon);
	      hist_Delta_phi_hadronic_top_leading_jet_2btag_corte_200->Fill(fabs(phi_th-(leading_jet->Phi+pi)));
	      hist_Delta_phi_leptonic_top_leading_jet_2btag_corte_200->Fill(fabs(phi_tl-(leading_jet->Phi+pi)));
	      if(musolo) hist_Delta_Phi_lepton_bl_2btag_corte_200->Fill(fabs((Muonpointer->Phi+pi)-phi_lb));
	      hist_Delta_phi_w_jets_bjet_hadronic_2btag_corte_200->Fill(fabs(phi_wjets-phi_bjethadronic));
	      hist_Cociente_pt_event_pt_leading_2btag_corte_200->Fill(pt_event/leading_jet->PT);
	      hist_MT_lb_2btag_corte_200->Fill(MT_lb);
	      hist_MT_top_2btag_corte_200->Fill(mt_top);
	      hist_rt_2btag_corte_200->Fill(rt);
	      hist_delta_phi_th_leading_jet_2btag_corte_200->Fill(delta_phi_th_leading_jet);
	      hist_delta_eta_th_leading_jet_2btag_corte_200->Fill(delta_eta_th_leading_jet);
	      hist_ST_2btag_corte_200->Fill(rST);
	      hist_Cociente_p_visible_leading_jet_pt_2btag_corte_200->Fill(pvisible/leading_jet->PT);
	      hist_Delta_phi_p_visible_leading_jet_2btag_corte_200->Fill(delta_phi_visible_leading_jet);
	      hist_Delta_eta_p_visible_leading_jet_2btag_corte_200->Fill(delta_eta_visible_leading_jet);
	      hist_MissingET_2btag_corte_200->Fill(metpointer->MET);
	      hist_Delta_Phi_leading_jet_met_2btag_corte_200->Fill(Delta_Phi_leading_jet_met);
	      hist_Delta_Eta_leading_jet_met_2btag_corte_200->Fill(Delta_Eta_leading_jet_met);
	      hist_rm_2btag_corte_200->Fill(rm);
	      hist_MT2_2btag_corte_200->Fill(value_mt2);
	      hist_MT2BL_2btag_corte_200->Fill(value_mt2bl);
	      hist_MT2W_2btag_corte_200->Fill(value_mt2w);
	      hist_Eta_average_2btag_corte_200->Fill(sum);
	      h_RM_Delta_Phi_leading_jet_met_2btag_corte_200->Fill(rm,Delta_Phi_leading_jet_met);
	      hist_cociente_HT_nobjet_PT_leading_jet_2btag_corte_200->Fill(HT_nobjet/leading_jet->PT);
	      Cociente_PT_top_PT_leading_corte_200->Fill(cociente_PT_top_PT_leading);
	      hist_Cociente_PT_top_lept_PT_leading_jet_2btag_corte_200->Fill(PT_top_lep/leading_jet->PT);
              hist_Cociente_PT_top_lept_PT_top_had_2btag_corte_200->Fill(PT_top_lep/PT_top);
	      if (HT_bjet>0.1) hist_cociente_HT_bjet_PT_leading_jet_2btag_corte_200->Fill(HT_bjet/leading_jet->PT);
	      
	      
	    }
	    
	    if (leading_jet->PT>300){
	      if (esolo||musolo) hist_Cociente_pt_visible_st_2btag_corte_300->Fill(pvisible/ST);
	      hist_Cociente_met_raiz_ht_2btag_corte_300->Fill(metpointer->MET/sqrt(HT));
	      hist_MT2BL_2btag_corte_300->Fill(value_mt2bl);
	      if(musolo)hist_Delta_Phi_leading_jet_muon_2btag_corte_300->Fill(Delta_Phi_leading_jet_muon);
	      hist_Delta_phi_hadronic_top_leading_jet_2btag_corte_300->Fill(fabs(phi_th-(leading_jet->Phi+pi)));
	      hist_Delta_phi_leptonic_top_leading_jet_2btag_corte_300->Fill(fabs(phi_tl-(leading_jet->Phi+pi)));
	      if(musolo) hist_Delta_Phi_lepton_bl_2btag_corte_300->Fill(fabs((Muonpointer->Phi+pi)-phi_lb));
	      hist_Delta_phi_w_jets_bjet_hadronic_2btag_corte_300->Fill(fabs(phi_wjets-phi_bjethadronic));
	      hist_Cociente_pt_event_pt_leading_2btag_corte_300->Fill(pt_event/leading_jet->PT);
	      hist_MT_lb_2btag_corte_300->Fill(MT_lb);
	      hist_MT_top_2btag_corte_300->Fill(mt_top);
	      hist_rt_2btag_corte_300->Fill(rt);
	      hist_delta_phi_th_leading_jet_2btag_corte_300->Fill(delta_phi_th_leading_jet);
	      hist_delta_eta_th_leading_jet_2btag_corte_300->Fill(delta_eta_th_leading_jet);
	      hist_ST_2btag_corte_300->Fill(rST);
	      hist_Cociente_p_visible_leading_jet_pt_2btag_corte_300->Fill(pvisible/leading_jet->PT);
	      hist_Delta_phi_p_visible_leading_jet_2btag_corte_300->Fill(delta_phi_visible_leading_jet);
	      hist_Delta_eta_p_visible_leading_jet_2btag_corte_300->Fill(delta_eta_visible_leading_jet);
	      hist_MissingET_2btag_corte_300->Fill(metpointer->MET);
	      hist_Delta_Phi_leading_jet_met_2btag_corte_300->Fill(Delta_Phi_leading_jet_met);
	      hist_Delta_Eta_leading_jet_met_2btag_corte_300->Fill(Delta_Eta_leading_jet_met);
	      hist_rm_2btag_corte_300->Fill(rm);
	      hist_MT2_2btag_corte_300->Fill(value_mt2);
	      hist_MT2BL_2btag_corte_300->Fill(value_mt2bl);
	      hist_MT2W_2btag_corte_300->Fill(value_mt2w);
	      hist_Eta_average_2btag_corte_300->Fill(sum);
	      h_RM_Delta_Phi_leading_jet_met_2btag_corte_300->Fill(rm,Delta_Phi_leading_jet_met);
	      hist_cociente_HT_nobjet_PT_leading_jet_2btag_corte_300->Fill(HT_nobjet/leading_jet->PT);
	      Cociente_PT_top_PT_leading_corte_300->Fill(cociente_PT_top_PT_leading);
	      hist_Cociente_PT_top_lept_PT_leading_jet_2btag_corte_300->Fill(PT_top_lep/leading_jet->PT);
              hist_Cociente_PT_top_lept_PT_top_had_2btag_corte_300->Fill(PT_top_lep/PT_top);
	      if (HT_bjet>0.1) hist_cociente_HT_bjet_PT_leading_jet_2btag_corte_300->Fill(HT_bjet/leading_jet->PT);
	      
	    }
	    
	    if (leading_jet->PT>400){
	      if (esolo||musolo) hist_Cociente_pt_visible_st_2btag_corte_400->Fill(pvisible/ST);
	      hist_Cociente_met_raiz_ht_2btag_corte_400->Fill(metpointer->MET/sqrt(HT));
	      hist_MT2BL_2btag_corte_400->Fill(value_mt2bl);
	      if(musolo)hist_Delta_Phi_leading_jet_muon_2btag_corte_400->Fill(Delta_Phi_leading_jet_muon);
	      hist_Delta_phi_hadronic_top_leading_jet_2btag_corte_400->Fill(fabs(phi_th-(leading_jet->Phi+pi)));
	      hist_Delta_phi_leptonic_top_leading_jet_2btag_corte_400->Fill(fabs(phi_tl-(leading_jet->Phi+pi)));
	      if(musolo) hist_Delta_Phi_lepton_bl_2btag_corte_400->Fill(fabs((Muonpointer->Phi+pi)-phi_lb));
	      hist_Delta_phi_w_jets_bjet_hadronic_2btag_corte_400->Fill(fabs(phi_wjets-phi_bjethadronic));
	      hist_Cociente_pt_event_pt_leading_2btag_corte_400->Fill(pt_event/leading_jet->PT);
	      hist_MT_lb_2btag_corte_400->Fill(MT_lb);
	      hist_MT_top_2btag_corte_400->Fill(mt_top);
	      hist_delta_phi_th_leading_jet_2btag_corte_400->Fill(delta_phi_th_leading_jet);
	      hist_delta_eta_th_leading_jet_2btag_corte_400->Fill(delta_eta_th_leading_jet);
	      hist_ST_2btag_corte_400->Fill(rST);
	      hist_Cociente_p_visible_leading_jet_pt_2btag_corte_400->Fill(pvisible/leading_jet->PT);
	      hist_Delta_phi_p_visible_leading_jet_2btag_corte_400->Fill(delta_phi_visible_leading_jet);
	      hist_Delta_eta_p_visible_leading_jet_2btag_corte_400->Fill(delta_eta_visible_leading_jet);
	      hist_rt_2btag_corte_400->Fill(rt);
	      hist_MissingET_2btag_corte_400->Fill(metpointer->MET);
	      hist_Delta_Phi_leading_jet_met_2btag_corte_400->Fill(Delta_Phi_leading_jet_met);
	      hist_Delta_Eta_leading_jet_met_2btag_corte_400->Fill(Delta_Eta_leading_jet_met);
	      hist_rm_2btag_corte_400->Fill(rm);
	      hist_MT2_2btag_corte_400->Fill(value_mt2);
	      hist_MT2BL_2btag_corte_400->Fill(value_mt2bl);
	      hist_MT2W_2btag_corte_400->Fill(value_mt2w);
	      hist_Eta_average_2btag_corte_400->Fill(sum);
	      h_RM_Delta_Phi_leading_jet_met_2btag_corte_400->Fill(rm,Delta_Phi_leading_jet_met);
	      hist_cociente_HT_nobjet_PT_leading_jet_2btag_corte_400->Fill(HT_nobjet/leading_jet->PT);
	      Cociente_PT_top_PT_leading_corte_400->Fill(cociente_PT_top_PT_leading);
	      hist_Cociente_PT_top_lept_PT_leading_jet_2btag_corte_400->Fill(PT_top_lep/leading_jet->PT);
              hist_Cociente_PT_top_lept_PT_top_had_2btag_corte_400->Fill(PT_top_lep/PT_top);
	      if (HT_bjet>0.1) hist_cociente_HT_bjet_PT_leading_jet_2btag_corte_400->Fill(HT_bjet/leading_jet->PT);
	      
	    }

	      if (leading_jet->PT<300) continue;


	    // } // end cuts 



	    // Scan over interesting variables

	      //RT
	      RT_min=0.0;
              RT_max=4.0;
	      deltav=(RT_max-RT_min)/nscans;
              valuev=RT_min;
	      for(int i =0; i<nscans;i++){
		if(rt>valuev) RT_counter[i]++;
		valuev+=deltav;
	      }

	      //MT2
	      MT2_min=0.0;
              MT2_max=350.0;
	      deltav=(MT2_max-MT2_min)/nscans;
              valuev=MT2_min;
	      for(int i =0; i<nscans;i++){
		if(value_mt2>valuev) MT2_counter[i]++;
		valuev+=deltav;
	      }


	    

      } //end btag(2)
      
     
      //Start btag(1)to find the bjet hadronic branch
      
      if(bjet_cntr==1){
	contador[6]++;	
	for(int ijet=0; ijet<NumJets; ijet++){
	  Jetpointer1 = (Jet*)branchJet->At(ijet);
	  if(Jetpointer1->BTag==1){
	    bj = ijet;
	    continue;
	  }
	  
	}
	Jetpointer1 = (Jet*)branchJet->At(bj);
	if (Jetpointer1->PT<30) continue;  // Demand  a threshold on bjet PT
	bool found_wjets=false;
	float chi_w1=9999999999.99;
	for(int ijet=0; ijet<NumJets; ijet++){
	  Jetpointer1 = (Jet*)branchJet->At(ijet);
          if(Jetpointer1->BTag==0) {
	    for(int kjet=0; kjet<NumJets; kjet++){
	      Jetpointer2 = (Jet*)branchJet->At(kjet);
	      if((kjet!=ijet)&& (Jetpointer2->BTag==0)){
		//W Invariant mass
		ej1 = CalculateE(Jetpointer1->Eta, Jetpointer1->PT, Jetpointer1->Mass);
	        ej2 = CalculateE(Jetpointer2->Eta, Jetpointer2->PT, Jetpointer2->Mass);
     		e_j1_j2 = ej1+ej2;
		p2_j1_j2 = CalculateP12_square_2particle(Jetpointer1->Eta, Jetpointer1->Phi, Jetpointer1->PT, Jetpointer2->Eta, Jetpointer2->Phi, Jetpointer2->PT);
		M_12 = TMath::Sqrt(TMath::Power((e_j1_j2),2)-p2_j1_j2);
		chi_dijet = TMath::Power((M_12-mw),2);
		if (chi_dijet<chi_w1) {
		  chi_w1=chi_dijet;
		  jw1=ijet;
		  jw2=kjet;
		  found_wjets=true;
		}
	      }
	    }
	  }
	}
      
	bool bhadronic = false;
	if (found_wjets) {
	  Jetpointerb = (Jet*)branchJet->At(bj);
	  Jetpointer1 = (Jet*)branchJet->At(jw1);
	  Jetpointer2 = (Jet*)branchJet->At(jw2);
	  ej1 = CalculateE(Jetpointer1->Eta, Jetpointer1->PT, Jetpointer1->Mass);
	  ej2 = CalculateE(Jetpointer2->Eta, Jetpointer2->PT, Jetpointer2->Mass);
	  eb = CalculateE(Jetpointerb->Eta, Jetpointerb->PT, Jetpointerb->Mass);
	  ej12b = ej1+ej2+eb;
	  p2_j12b = CalculateP123_square_3particle(Jetpointer1->Eta, Jetpointer1->Phi, Jetpointer1->PT, Jetpointer2->Eta, Jetpointer2->Phi, Jetpointer2->PT, Jetpointerb->Eta, Jetpointerb->Phi, Jetpointerb->PT);
	  M_12b = TMath::Sqrt(TMath::Power((ej12b),2)-p2_j12b);
	  if(fabs(M_12b-mt)<3*sigma_12b){
	    jbh=bj;
            bhadronic = true;
	  }
	
	  if(!bhadronic) {
	    for (int njet=0;njet<NumJets;njet++){
	      if ((njet!=bj)&&(njet!=jw1)&&(njet!=jw2)){
		Jetpointer1 = (Jet*)branchJet->At(jw1);
		Jetpointer2 = (Jet*)branchJet->At(jw2);
		Jetpointer3 = (Jet*)branchJet->At(njet);
		ej1 = CalculateE(Jetpointer1->Eta, Jetpointer1->PT, Jetpointer1->Mass);
		ej2 = CalculateE(Jetpointer2->Eta, Jetpointer2->PT, Jetpointer2->Mass);
		ej3 = CalculateE(Jetpointer3->Eta, Jetpointer3->PT, Jetpointer3->Mass);
		p2_j1_j2_j3 = CalculateP123_square_3particle(Jetpointer1->Eta, Jetpointer1->Phi, Jetpointer1->PT, Jetpointer2->Eta, Jetpointer2->Phi, Jetpointer2->PT, Jetpointer3->Eta, Jetpointer3->Phi, Jetpointer3->PT);
		p2_j1_j2 = CalculateP12_square_2particle(Jetpointer1->Eta, Jetpointer1->Phi, Jetpointer1->PT, Jetpointer2->Eta, Jetpointer2->Phi, Jetpointer2->PT);
		e_j1_j2_j3 = ej1+ej2+ej3;
		e_j1_j2 = ej1+ej2;
		M_12 = TMath::Sqrt(TMath::Power((e_j1_j2),2)-p2_j1_j2);
		M_123 = TMath::Sqrt(TMath::Power((e_j1_j2_j3),2)-p2_j1_j2_j3);
		if(fabs(M_123-mt)<3*sigma_12b){
		  jbh=njet;
		  bhadronic = true;
		}
		
	      }
	    } //end njets
	    
	  } //end no bhadronic
	  //end w jets found
	  
	  
	  if (bhadronic){
	    Jetpointer1 = (Jet*)branchJet->At(jw1);
	    Jetpointer2 = (Jet*)branchJet->At(jw2);
	    Jetpointerb = (Jet*)branchJet->At(jbh);
	    ej1 = CalculateE(Jetpointer1->Eta, Jetpointer1->PT, Jetpointer1->Mass);
	    ej2 = CalculateE(Jetpointer2->Eta, Jetpointer2->PT, Jetpointer2->Mass);
	    e_j1_j2 = ej1+ej2;
	    p2_j1_j2 = CalculateP12_square_2particle(Jetpointer1->Eta, Jetpointer1->Phi, Jetpointer1->PT, Jetpointer2->Eta, Jetpointer2->Phi, Jetpointer2->PT);
	    M_12 = TMath::Sqrt(TMath::Power((e_j1_j2),2)-p2_j1_j2);
	    ejb = CalculateE(Jetpointerb->Eta, Jetpointerb->PT, Jetpointerb->Mass);
	    ej12b = ej1+ej2+ejb;
	    p2_j12b = CalculateP123_square_3particle(Jetpointer1->Eta, Jetpointer1->Phi, Jetpointer1->PT, Jetpointer2->Eta, Jetpointer2->Phi, Jetpointer2->PT, Jetpointerb->Eta, Jetpointerb->Phi, Jetpointerb->PT);
	    M_12b = TMath::Sqrt(TMath::Power((ej12b),2)-p2_j12b);
	    
	    hist_w_recons_1bjet->Fill(M_12);
	    hist_t_recons_1bjet->Fill(M_12b);
	    
	  }
	} // end w jets found

	// Some cuts
	//if((rm>0.5)&& (metpointer->MET>200)&&(Delta_Phi_leading_jet_met<pi)){
	  if (leading_jet->PT>200){
		if (esolo||musolo)  hist_ST_1btag_corte_200->Fill(ST);
                hist_Delta_eta_ave_1btag_corte_200->Fill(delta_eta_ave);
		if(esolo||musolo)hist_Cociente_met_st_1btag_corte_200->Fill(metpointer->MET/ST);
		hist_Cociente_p_visible_leading_jet_pt_1btag_corte_200->Fill(pvisible/leading_jet->PT);
		hist_Delta_phi_p_visible_leading_jet_1btag_corte_200->Fill(delta_phi_visible_leading_jet);
		hist_Delta_eta_p_visible_leading_jet_1btag_corte_200->Fill(delta_eta_visible_leading_jet);
		hist_MissingET_1btag_corte_200->Fill(metpointer->MET);
		hist_Delta_Phi_leading_jet_met_1btag_corte_200->Fill(Delta_Phi_leading_jet_met);
		hist_Delta_Eta_leading_jet_met_1btag_corte_200->Fill(Delta_Eta_leading_jet_met);
		hist_rm_1btag_corte_200->Fill(rm);
		hist_Eta_average_1btag_corte_200->Fill(sum);
		hist_cociente_HT_PT_leading_jet_1btag_corte_200->Fill(HT/leading_jet->PT);
		h_RM_Delta_Phi_leading_jet_met_1btag_corte_200->Fill(rm,Delta_Phi_leading_jet_met);
		hist_cociente_HT_PT_leading_jet_RM_1btag_corte_200->Fill(HT/leading_jet->PT,rm);
		hist_cociente_HT_nobjet_PT_leading_jet_1btag_corte_200->Fill(HT_nobjet/leading_jet->PT);
		if (musolo) hist_MT_Muon_1btag_corte_200->Fill(transverse_mass_Muon);
		if (HT_bjet>0.1) hist_cociente_HT_bjet_PT_leading_jet_1btag_corte_200->Fill(HT_bjet/leading_jet->PT);
		if(musolo||esolo) hist_MT_lepton_1btag_corte_200->Fill(transverse_mass_Lepton);				
		
	    }
	      if (leading_jet->PT>300){
		
		if (esolo||musolo)  hist_ST_1btag_corte_300->Fill(ST);
                hist_Delta_eta_ave_1btag_corte_300->Fill(delta_eta_ave);
		if(esolo||musolo)hist_Cociente_met_st_1btag_corte_300->Fill(metpointer->MET/ST);
		hist_Cociente_p_visible_leading_jet_pt_1btag_corte_300->Fill(pvisible/leading_jet->PT);
		hist_Delta_phi_p_visible_leading_jet_1btag_corte_300->Fill(delta_phi_visible_leading_jet);
		hist_Delta_eta_p_visible_leading_jet_1btag_corte_300->Fill(delta_eta_visible_leading_jet);
		hist_MissingET_1btag_corte_300->Fill(metpointer->MET);
		hist_Delta_Phi_leading_jet_met_1btag_corte_300->Fill(Delta_Phi_leading_jet_met);
		hist_Delta_Eta_leading_jet_met_1btag_corte_300->Fill(Delta_Eta_leading_jet_met);
		hist_rm_1btag_corte_300->Fill(rm);
		hist_Eta_average_1btag_corte_300->Fill(sum);
		hist_cociente_HT_PT_leading_jet_1btag_corte_300->Fill(HT/leading_jet->PT);
		h_RM_Delta_Phi_leading_jet_met_1btag_corte_300->Fill(rm,Delta_Phi_leading_jet_met);
		hist_cociente_HT_PT_leading_jet_RM_1btag_corte_300->Fill(HT/leading_jet->PT,rm);
		hist_cociente_HT_nobjet_PT_leading_jet_1btag_corte_300->Fill(HT_nobjet/leading_jet->PT);
		if (musolo)hist_MT_Muon_1btag_corte_300->Fill(transverse_mass_Muon);
		if (HT_bjet>0.1) hist_cociente_HT_bjet_PT_leading_jet_1btag_corte_300->Fill(HT_bjet/leading_jet->PT);
		if(musolo||esolo) hist_MT_lepton_1btag_corte_300->Fill(transverse_mass_Lepton);				
				
		
	      }
	      
	      if (leading_jet->PT>400){
		if (esolo||musolo)  hist_ST_1btag_corte_400->Fill(ST);
                hist_Delta_eta_ave_1btag_corte_400->Fill(delta_eta_ave);
		if(esolo||musolo)hist_Cociente_met_st_1btag_corte_400->Fill(metpointer->MET/ST);
		hist_Cociente_p_visible_leading_jet_pt_1btag_corte_400->Fill(pvisible/leading_jet->PT);
		hist_Delta_phi_p_visible_leading_jet_1btag_corte_400->Fill(delta_phi_visible_leading_jet);
		hist_Delta_eta_p_visible_leading_jet_1btag_corte_400->Fill(delta_eta_visible_leading_jet);
		hist_MissingET_1btag_corte_400->Fill(metpointer->MET);
		hist_Delta_Phi_leading_jet_met_1btag_corte_400->Fill(Delta_Phi_leading_jet_met);
		hist_Delta_Eta_leading_jet_met_1btag_corte_400->Fill(Delta_Eta_leading_jet_met);
		hist_rm_1btag_corte_400->Fill(rm);
		hist_Eta_average_1btag_corte_400->Fill(sum);
		hist_cociente_HT_PT_leading_jet_1btag_corte_400->Fill(HT/leading_jet->PT);
		h_RM_Delta_Phi_leading_jet_met_1btag_corte_400->Fill(rm,Delta_Phi_leading_jet_met);
		hist_cociente_HT_PT_leading_jet_RM_1btag_corte_400->Fill(HT/leading_jet->PT,rm);
		hist_cociente_HT_nobjet_PT_leading_jet_1btag_corte_400->Fill(HT_nobjet/leading_jet->PT);
		if (musolo)hist_MT_Muon_1btag_corte_400->Fill(transverse_mass_Muon);	  
		if (HT_bjet>0.1) hist_cociente_HT_bjet_PT_leading_jet_1btag_corte_400->Fill(HT_bjet/leading_jet->PT);
		if(musolo||esolo) hist_MT_lepton_1btag_corte_400->Fill(transverse_mass_Lepton);				

				
	      }

	      //  } // end cuts


	      // SCAN over discriminating variables
	      //first decide the leading JET PT

	      if (leading_jet->PT<300) continue;


	      //MET
	      
	      
	      MET_min=0.0;
              MET_max=500.0;
	      deltav=(MET_max-MET_min)/nscans;
              valuev=MET_min;
	      for(int i =0; i<nscans;i++){
		if(metpointer->MET>valuev) MET_counter[i]++;
		valuev+=deltav;
	      }

	      //RM
	      RM_min=0.0;
              RM_max=1.0;
	      deltav=(RM_max-RM_min)/nscans;
              valuev=RM_min;
	      for(int i =0; i<nscans;i++){
		if(rm>valuev) RM_counter[i]++;
		valuev+=deltav;
	      }

	      //Delta_Phi_leading_jet_met
	      Delta_Phi_leading_jet_met_min=0.0;
              Delta_Phi_leading_jet_met_max=5.0;
	      deltav=(Delta_Phi_leading_jet_met_max-Delta_Phi_leading_jet_met_min)/nscans;
              valuev=Delta_Phi_leading_jet_met_min;

	      for(int i =0; i<nscans;i++){
		if(Delta_Phi_leading_jet_met>valuev) Delta_Phi_leading_jet_met_counter[i]++;
		valuev+=deltav;
	      }

	      //Cociente_p_visible_leading_jet_pt
	      Cociente_p_visible_leading_jet_pt_min=0.0;
              Cociente_p_visible_leading_jet_pt_max=1.0;
	      deltav=(Cociente_p_visible_leading_jet_pt_max-Cociente_p_visible_leading_jet_pt_min)/nscans;
              valuev=Cociente_p_visible_leading_jet_pt_min;
	      for(int i =0; i<nscans;i++){
		if(pvisible/leading_jet->PT>valuev) Cociente_p_visible_leading_jet_pt_counter[i]++;
		valuev+=deltav;
	      }

	      //Cociente_HT_bjet_PT_leading_jet
	      Cociente_HT_bjet_PT_leading_jet_min=0.0;
              Cociente_HT_bjet_PT_leading_jet_max=1.0;
	      deltav=(Cociente_HT_bjet_PT_leading_jet_max-Cociente_HT_bjet_PT_leading_jet_min)/nscans;
              valuev=Cociente_HT_bjet_PT_leading_jet_min;
	      for(int i =0; i<nscans;i++){
		if(HT_bjet/leading_jet->PT>valuev) Cociente_HT_bjet_PT_leading_jet_counter[i]++;
		valuev+=deltav;
	      }

	      //PT leading jet

	      PT_leading_jet_min=0.0;
              PT_leading_jet_max=500.0;
	      deltav=(PT_leading_jet_max-PT_leading_jet_min)/nscans;
              valuev=PT_leading_jet_min;
	      for(int i =0; i<nscans;i++){
		if(leading_jet->PT>valuev) PT_leading_jet_counter[i]++;
		valuev+=deltav;
	      }

	      //MT
	      MT_min=0.0;
              MT_max=200.0;
	      deltav=(MT_max-MT_min)/nscans;
              valuev=MT_min;
	      for(int i =0; i<nscans;i++){
		if(transverse_mass_Lepton>valuev) MT_counter[i]++;
		valuev+=deltav;
	      }



	      	 
      } //end Btag(1)
 
       }// end topology analysis
  }//end loop
  ofstream events;
  //events.open("semileptonic_events.txt",std::ios_base::app);
  //events<<nevents<<endl;
  //events.close();

  // events.open("dileptonic_events.txt",std::ios_base::app);
  //events<<nevents<<endl;
  //events.close();

  events.open("stop_400_n1_235_events.txt",std::ios_base::app);
  events<<nevents<<endl;
  events.close();
  
  // events.open("stop_200_n1_20_events.txt",std::ios_base::app);
  //events<<nevents<<endl;
  //events.close();

  cout<<"**********************************************"<<endl;
  cout<<"******** # events processed = "<<nevents<<" *********** "<<endl;
  cout<<"**********************************************"<<endl;
  //cout<<"Number of muons: "<<one_muon<<"Number of electrons: "<<one_electron<<endl;
  //cout<<"More muons: "<<more_muon<<"More electrons: "<<more_electron<<endl;
  //cout<<"Number of muons and electrons: "<<more_muon_electron<<endl;
  //cout<<"Number of bjets: "<<bjet_counter<<endl;
  cout<<"Numero de eventos procesados: "<<contador[0]<<endl;
  cout<<"Topology analysis "<<contador[1]<<endl;
  cout<<"Lepton pt cut:  "<<contador[2]<<endl;
  cout<<"bjet pt cut:   "<<contador[3]<<endl;
  //cout<<"Eventos con lepton pt cut:   "<<contador[4]<<endl;
  //cout<<"Eventos con bjet pt cut   "<<contador[5]<<endl;
  //cout<<"Eventos con 1 bjet:   "<<contador[6]<<endl;
  

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////WRITE SCAN RESULTS//////////////////////////////////////////////////////////////////////////////////
  //MET 
  ofstream scan_out;
  valuev=0;
  deltav=(MET_max-MET_min)/nscans;
  valuev=MET_min;
  scan_out.open("MET_scan_out",std::ios_base::app);
  //scan_out<<" "<<nscans<<" "<<MET_min<<" "<<MET_max<<endl;
  for(int i=0;i<nscans;i++){
    scan_out<<"   "<<valuev<<"  "<<MET_counter[i];
    valuev+=deltav;
  }
  scan_out <<"\n";
  scan_out.close();


  //RM
  deltav=(RM_max-RM_min)/nscans;
  valuev=RM_min;
  scan_out.open("RM_scan_out",std::ios_base::app);
  //scan_out<<" "<<nscans<<" "<<RM_min<<" "<<RM_max<<endl;
  for(int i=0;i<nscans;i++){
    scan_out<<"   "<<valuev<<"   "<<RM_counter[i];
    valuev+=deltav;

  }
  scan_out <<"\n";
  scan_out.close();
 //Delta_Phi_leading_jet_met
  deltav=(Delta_Phi_leading_jet_met_max-Delta_Phi_leading_jet_met_min)/nscans;
  valuev=Delta_Phi_leading_jet_met_min;
  scan_out.open("Delta_Phi_leading_jet_met_scan_out",std::ios_base::app);
  //scan_out<<" "<<nscans<<" "<<Delta_Phi_leading_jet_met_min<<" "<<Delta_Phi_leading_jet_met_max<<endl;
  for(int i=0;i<nscans;i++){
    scan_out<<"   "<<valuev<<"   "<<Delta_Phi_leading_jet_met_counter[i];
    valuev+=deltav;

  }
  scan_out <<"\n";
  scan_out.close();
  //Cociente_p_visible_leading_jet_pt
  deltav=(Cociente_p_visible_leading_jet_pt_max-Cociente_p_visible_leading_jet_pt_min)/nscans;
  valuev=Cociente_p_visible_leading_jet_pt_min;
  scan_out.open("Cociente_p_visible_leading_jet_pt_scan_out",std::ios_base::app);
  //scan_out<<" "<<nscans<<" "<<Cociente_p_visible_leading_jet_pt_min<<" "<<Cociente_p_visible_leading_jet_pt_max<<endl;
  for(int i=0;i<nscans;i++){
    scan_out<<"   "<<valuev<<"   "<<Cociente_p_visible_leading_jet_pt_counter[i];
    valuev+=deltav;

  }
  scan_out <<"\n";
  scan_out.close();
  //Cociente_HT_bjet_PT_leading_jet
  deltav=(Cociente_HT_bjet_PT_leading_jet_max-Cociente_HT_bjet_PT_leading_jet_min)/nscans;
  valuev=Cociente_HT_bjet_PT_leading_jet_min; 
  scan_out.open("Cociente_HT_bjet_PT_leading_jet_scan_out",std::ios_base::app);
  //scan_out<<" "<<nscans<<" "<<Cociente_HT_bjet_PT_leading_jet_min<<" "<<Cociente_HT_bjet_PT_leading_jet_max<<endl;
  for(int i=0;i<nscans;i++){
    scan_out<<"   "<<valuev<<"   "<<Cociente_HT_bjet_PT_leading_jet_counter[i];
    valuev+=deltav;

  }
  scan_out <<"\n";
  scan_out.close();

  //RT
  deltav=(RT_max-RT_min)/nscans;
  valuev=RT_min;
  scan_out.open("RT_scan_out",std::ios_base::app);
  //scan_out<<" "<<nscans<<" "<<RT_min<<" "<<RT_max<<endl;
  for(int i=0;i<nscans;i++){
    scan_out<<"   "<<valuev<<"  "<<RT_counter[i];
    valuev+=deltav;

  }
  scan_out <<"\n";
  scan_out.close();

 //MT2
  deltav=(MT2_max-MT2_min)/nscans;
  valuev=MT2_min;
  scan_out.open("MT2_scan_out",std::ios_base::app);
  //scan_out<<" "<<nscans<<" "<<MT2_min<<" "<<MT2_max<<endl;
  for(int i=0;i<nscans;i++){
    scan_out<<"   "<<valuev<<"  "<<MT2_counter[i];
    valuev+=deltav;

  }
  scan_out <<"\n";
  scan_out.close();

 //MT
  deltav=(MT_max-MT_min)/nscans;
  valuev=MT_min;
  scan_out.open("MT_scan_out",std::ios_base::app);
  //scan_out<<" "<<nscans<<" "<<MT_min<<" "<<MT_max<<endl;
  for(int i=0;i<nscans;i++){
    scan_out<<"   "<<valuev<<"  "<<MT_counter[i];
    valuev+=deltav;

  }
  scan_out <<"\n";
  scan_out.close();
  //PT_leading_jet
  deltav=(PT_leading_jet_max-PT_leading_jet_min)/nscans;
  valuev=PT_leading_jet_min;
  scan_out.open("PT_leading_jet_scan_out",std::ios_base::app);
  //scan_out<<" "<<nscans<<" "<<PT_leading_jet_min<<" "<<PT_leading_jet_max<<endl;
  for(int i=0;i<nscans;i++){
    scan_out<<"   "<<valuev<<"  "<<PT_leading_jet_counter[i];
    valuev+=deltav;

  }
  scan_out <<"\n";
  scan_out.close();

  
  // ROOT Program output
  
  TFile* hfile = new TFile(argv[2], "RECREATE");
  
    histNumJets->Write();
    histNumElectrons->Write();
    histNumMuons->Write();
    histsigma_123->Write();
    histNum_btags->Write();
    hist_Eta_average->Write();
    hist_w_recons_2bjet->Write();
    hist_t_recons_2bjet->Write();
    hist_w_recons_1bjet->Write();
    hist_t_recons_1bjet->Write();
    hist_MT_Muon->Write();
    hist_MT_Electron->Write();
    hist_MT2->Write();
    hist_MT2W->Write();
    hist_MT2BL->Write();
    hist_rt_2btag->Write();
    hist_delta_phi_th_leading_jet_2btag->Write();
    hist_delta_eta_th_leading_jet_2btag->Write();
    hist_cociente_HT_nobjet_PT_leading_jet->Write();
    ////////////////////////////////////////////////////////////
    hist_MT_lepton_1btag_corte_200->Write();
    hist_MT_lepton_1btag_corte_300->Write();
    hist_MT_lepton_1btag_corte_400->Write();


    ///////////////////////////////////////////////////////////

    hist_Cociente_pt_visible_st_2btag_corte_200->Write();
    hist_Cociente_pt_visible_st_2btag_corte_300->Write();
    hist_Cociente_pt_visible_st_2btag_corte_400->Write();
    ////////////////////////////////////////////////////////////
    hist_Cociente_met_st_1btag_corte_200->Write();
    hist_Cociente_met_st_1btag_corte_300->Write();
    hist_Cociente_met_st_1btag_corte_400->Write();

    ////////////////////////////////////////////////////////////
    hist_Delta_eta_ave_1btag_corte_200->Write();
    hist_Delta_eta_ave_1btag_corte_300->Write();
    hist_Delta_eta_ave_1btag_corte_400->Write();

    /////////////////////////////////////////////////////////////
    hist_Cociente_met_raiz_ht_2btag_corte_200->Write();
    hist_Cociente_met_raiz_ht_2btag_corte_300->Write();
    hist_Cociente_met_raiz_ht_2btag_corte_400->Write();

    ////////////////////////////////////////////////////////////

    hist_Cociente_pt_event_pt_leading_2btag_corte_200->Write();
    hist_Cociente_pt_event_pt_leading_2btag_corte_300->Write();
    hist_Cociente_pt_event_pt_leading_2btag_corte_400->Write();

    /////////////////////////////////////////////////////////////
    hist_Delta_phi_w_jets_bjet_hadronic_2btag_corte_200->Write();
    hist_Delta_phi_w_jets_bjet_hadronic_2btag_corte_300->Write();
    hist_Delta_phi_w_jets_bjet_hadronic_2btag_corte_400->Write();

    //////////////////////////////////////////////////////////////
    hist_Delta_Phi_lepton_bl_2btag_corte_200->Write();
    hist_Delta_Phi_lepton_bl_2btag_corte_300->Write();
    hist_Delta_Phi_lepton_bl_2btag_corte_400->Write();
    /////////////////////////////////////////////////////////////
    hist_Delta_phi_hadronic_top_leading_jet_2btag_corte_200->Write();
    hist_Delta_phi_hadronic_top_leading_jet_2btag_corte_300->Write();
    hist_Delta_phi_hadronic_top_leading_jet_2btag_corte_400->Write();

    ////////////////////////////////////////////////////////////////
    hist_Delta_phi_leptonic_top_leading_jet_2btag_corte_200->Write();
    hist_Delta_phi_leptonic_top_leading_jet_2btag_corte_300->Write();
    hist_Delta_phi_leptonic_top_leading_jet_2btag_corte_400->Write();
    /////////////////////////////////////////////////////////////////

    hist_Cociente_PT_top_lept_PT_leading_jet_2btag_corte_200->Write();
    hist_Cociente_PT_top_lept_PT_leading_jet_2btag_corte_300->Write();
    hist_Cociente_PT_top_lept_PT_leading_jet_2btag_corte_400->Write();

    ///////////////////////////////////////////////////////////////
    hist_Cociente_PT_top_lept_PT_top_had_2btag_corte_200->Write();
    hist_Cociente_PT_top_lept_PT_top_had_2btag_corte_300->Write();
    hist_Cociente_PT_top_lept_PT_top_had_2btag_corte_400->Write();
    ///////////////////////////////////////////////////////////
    Cociente_PT_top_PT_leading_corte_200->Write();
    Cociente_PT_top_PT_leading_corte_300->Write();
    Cociente_PT_top_PT_leading_corte_400->Write();

  
  //////////////////////////////////////////////
    hist_Delta_Phi_leading_jet_muon_2btag_corte_200->Write();
    hist_Delta_Phi_leading_jet_muon_2btag_corte_300->Write();
    hist_Delta_Phi_leading_jet_muon_2btag_corte_400->Write();
    ///////////////////////////////////////////
    hist_MT_lb_2btag->Write();
  hist_MT_lb_2btag_corte_200->Write();
  hist_MT_lb_2btag_corte_300->Write();
  hist_MT_lb_2btag_corte_400->Write();
  hist_MT_lb_1btag_corte_200->Write();
  hist_MT_lb_1btag_corte_300->Write();
  hist_MT_lb_1btag_corte_400->Write();

  ///////////////////////////////////////////
  hist_MT_top_2btag->Write();
  hist_MT_top_2btag_corte_200->Write();
  hist_MT_top_2btag_corte_300->Write();
  hist_MT_top_2btag_corte_400->Write();
  ///////////////////////////////////////
  //ST
  histNum_btags->Write();
  hist_ST_2btag_corte_200->Write();
  hist_ST_2btag_corte_300->Write();
  hist_ST_2btag_corte_400->Write();
  hist_ST_1btag_corte_200->Write();
  hist_ST_1btag_corte_300->Write();
  hist_ST_1btag_corte_400->Write();
  //////////////////////////////////////////
  //Cociente pvisible/leading_jet_PT
  hist_Cociente_p_visible_leading_jet_pt->Write();
  hist_Cociente_p_visible_leading_jet_pt_2btag_corte_200->Write();
  hist_Cociente_p_visible_leading_jet_pt_2btag_corte_300->Write();
  hist_Cociente_p_visible_leading_jet_pt_2btag_corte_400->Write();
  hist_Cociente_p_visible_leading_jet_pt_1btag_corte_200->Write();
  hist_Cociente_p_visible_leading_jet_pt_1btag_corte_300->Write();
  hist_Cociente_p_visible_leading_jet_pt_1btag_corte_400->Write();
  /////////////////////////////////////////
  //Delta_Phi_p_visible_leading_jet
  hist_Delta_phi_p_visible_leading_jet->Write();
  hist_Delta_phi_p_visible_leading_jet_2btag_corte_200->Write();
  hist_Delta_phi_p_visible_leading_jet_2btag_corte_300->Write();
  hist_Delta_phi_p_visible_leading_jet_2btag_corte_400->Write();
  hist_Delta_phi_p_visible_leading_jet_1btag_corte_200->Write();
  hist_Delta_phi_p_visible_leading_jet_1btag_corte_300->Write();
  hist_Delta_phi_p_visible_leading_jet_1btag_corte_400->Write();

  //////////////////////////////////////////

  //Delta_Eta_p_visible_leading_jet
  hist_Delta_eta_p_visible_leading_jet->Write();
  hist_Delta_eta_p_visible_leading_jet_2btag_corte_200->Write();
  hist_Delta_eta_p_visible_leading_jet_2btag_corte_300->Write();
  hist_Delta_eta_p_visible_leading_jet_2btag_corte_400->Write();
  hist_Delta_eta_p_visible_leading_jet_1btag_corte_200->Write();
  hist_Delta_eta_p_visible_leading_jet_1btag_corte_300->Write();
  hist_Delta_eta_p_visible_leading_jet_1btag_corte_400->Write();

  /////////////////////////////////////////
  //RT
  hist_rt_2btag_corte_200->Write();
  hist_rt_2btag_corte_300->Write();
  hist_rt_2btag_corte_400->Write();
  /////////////////////////////////
  hist_rm_1btag_corte_200->Write();
  hist_rm_1btag_corte_300->Write();
  hist_rm_1btag_corte_400->Write();
  /////////////////////////////////
  hist_cociente_HT_bjet_PT_leading_jet_1btag_corte_200->Write();
  hist_cociente_HT_bjet_PT_leading_jet_1btag_corte_300->Write();
  hist_cociente_HT_bjet_PT_leading_jet_1btag_corte_400->Write();
  ////////////////////////////////////////////////////////////////
  hist_cociente_HT_bjet_PT_leading_jet_2btag_corte_200->Write();
  hist_cociente_HT_bjet_PT_leading_jet_2btag_corte_300->Write();
  hist_cociente_HT_bjet_PT_leading_jet_2btag_corte_400->Write();

  //////////////////////////////////////////////////////////////
  hist_MT2_2btag_corte_200->Write();
  hist_MT2_2btag_corte_300->Write();
  hist_MT2_2btag_corte_400->Write();

  /////////////////////////////////////////////////////////////

  hist_MT_Muon_1btag_corte_200->Write();
  hist_MT_Muon_1btag_corte_300->Write();
  hist_MT_Muon_1btag_corte_400->Write();

  /////////////////////////////////////////////////////////////

  hist_delta_phi_th_leading_jet_2btag_corte_200->Write();
  hist_delta_phi_th_leading_jet_2btag_corte_300->Write();
  hist_delta_phi_th_leading_jet_2btag_corte_400->Write();
  hist_delta_eta_th_leading_jet_2btag_corte_200->Write();
  hist_delta_eta_th_leading_jet_2btag_corte_300->Write();
  hist_delta_eta_th_leading_jet_2btag_corte_400->Write();

  hist_cociente_HT_bjet_PT_leading_jet->Write();
  hist_cociente_PT_Muon_PT_leading_jet_2btag_corte_200->Write();
  hist_cociente_PT_Muon_PT_leading_jet_2btag_corte_300->Write();
  hist_cociente_PT_Muon_PT_leading_jet_2btag_corte_400->Write();
  hist_cociente_PT_Muon_PT_leading_jet_1btag_corte_200->Write();
  hist_cociente_PT_Muon_PT_leading_jet_1btag_corte_300->Write();
  hist_cociente_PT_Muon_PT_leading_jet_1btag_corte_400->Write();
  hist_cociente_PT_Electron_PT_leading_jet_2btag_corte_200->Write();
  hist_cociente_PT_Electron_PT_leading_jet_2btag_corte_300->Write();
  hist_cociente_PT_Electron_PT_leading_jet_2btag_corte_400->Write();
  hist_cociente_PT_Electron_PT_leading_jet_1btag_corte_200->Write();
  hist_cociente_PT_Electron_PT_leading_jet_1btag_corte_300->Write();
  hist_cociente_PT_Electron_PT_leading_jet_1btag_corte_400->Write();

  hist_cociente_HT_nobjet_PT_leading_jet_2btag_corte_200->Write();
  hist_cociente_HT_nobjet_PT_leading_jet_2btag_corte_300->Write();
  hist_cociente_HT_nobjet_PT_leading_jet_2btag_corte_400->Write();
  hist_cociente_HT_nobjet_PT_leading_jet_1btag_corte_200->Write();
  hist_cociente_HT_nobjet_PT_leading_jet_1btag_corte_300->Write();
  hist_cociente_HT_nobjet_PT_leading_jet_1btag_corte_400->Write();
  hist_cociente_HT_bjet_PT_leading_jet_2btag_corte_200->Write();
  hist_cociente_HT_bjet_PT_leading_jet_2btag_corte_300->Write();
  hist_cociente_HT_bjet_PT_leading_jet_2btag_corte_400->Write();
  hist_cociente_HT_bjet_PT_leading_jet_1btag_corte_200->Write();
  hist_cociente_HT_bjet_PT_leading_jet_1btag_corte_300->Write();
  hist_cociente_HT_bjet_PT_leading_jet_1btag_corte_400->Write();

  h_RM_Delta_Phi_leading_jet_met_2btag_corte_200->Write();
  h_RM_Delta_Phi_leading_jet_met_2btag_corte_300->Write();
  h_RM_Delta_Phi_leading_jet_met_2btag_corte_400->Write();
  h_RM_Delta_Phi_leading_jet_met_1btag_corte_200->Write();
  h_RM_Delta_Phi_leading_jet_met_1btag_corte_300->Write();
  h_RM_Delta_Phi_leading_jet_met_1btag_corte_400->Write();
  hist_cociente_HT_PT_leading_jet_RM_1btag_corte_200->Write();
  hist_cociente_HT_PT_leading_jet_RM_1btag_corte_300->Write();
  hist_cociente_HT_PT_leading_jet_RM_1btag_corte_400->Write();
  hist_cociente_HT_PT_leading_jet_1btag_corte_200->Write();
  hist_cociente_HT_PT_leading_jet_1btag_corte_300->Write();
  hist_cociente_HT_PT_leading_jet_1btag_corte_400->Write();
  hist_MT2_2btag_corte_200->Write();
  hist_MT2BL_2btag_corte_200->Write();
  hist_MT2W_2btag_corte_200->Write();
  hist_Eta_average_2btag_corte_200->Write();
  hist_Eta_average_2btag_corte_300->Write();
  hist_Eta_average_2btag_corte_400->Write();
  hist_MT2_2btag_corte_300->Write();
  hist_MT2BL_2btag_corte_300->Write();
  hist_MT2W_2btag_corte_300->Write();
  hist_MT2_2btag_corte_400->Write();
  hist_MT2BL_2btag_corte_400->Write();
  hist_MT2W_2btag_corte_400->Write();
  hist_Delta_Phi_leading_jet_muon_2btag->Write();
  hist_Delta_Phi_leading_jet_electron_2btag->Write();
  hist_rm_2btag_corte_200->Write();
  hist_Delta_Phi_leading_jet_met_2btag_corte_200->Write();
  hist_Delta_Eta_leading_jet_met_2btag_corte_200->Write();
  hist_MissingET_2btag_corte_200->Write();
  hist_rm_2btag_corte_300->Write();
  hist_Delta_Phi_leading_jet_met_2btag_corte_300->Write();
  hist_Delta_Eta_leading_jet_met_2btag_corte_300->Write();
  hist_MissingET_2btag_corte_300->Write();
  hist_rm_2btag_corte_400->Write();
  hist_Delta_Phi_leading_jet_met_2btag_corte_400->Write();
  hist_Delta_Eta_leading_jet_met_2btag_corte_400->Write();
  hist_MissingET_2btag_corte_400->Write();

  hist_rm_1btag_corte_200->Write();
  hist_Delta_Phi_leading_jet_met_1btag_corte_200->Write();
  hist_Delta_Eta_leading_jet_met_1btag_corte_200->Write();
  hist_MissingET_1btag_corte_200->Write();
  hist_Eta_average_1btag_corte_200->Write();
  hist_Eta_average_1btag_corte_300->Write();
  hist_Eta_average_1btag_corte_400->Write();
  hist_rm_1btag_corte_300->Write();
  hist_Delta_Phi_leading_jet_met_1btag_corte_300->Write();
  hist_Delta_Eta_leading_jet_met_1btag_corte_300->Write();
  hist_MissingET_1btag_corte_300->Write();
  hist_rm_1btag_corte_400->Write();
  hist_Delta_Phi_leading_jet_met_1btag_corte_400->Write();
  hist_Delta_Eta_leading_jet_met_1btag_corte_400->Write();
  hist_MissingET_1btag_corte_400->Write();



  hfile->Close();



}  //end of the program



