////////////////////////////////////////////////////////////////
//                                                            //
// Author: Andrés Flórez, Universidad de los Andes, Colombia  //
//                                                            //
////////////////////////////////////////////////////////////////


#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ROOTFunctions.h"
#include "PhenoAnalyzer.h"
#include "DelphesFunctions.h"


int main(int argc, char *argv[]) {
  
  TChain chain("Delphes");
  chain.Add(argv[1]);
  TFile * HistoOutputFile = new TFile(argv[2], "RECREATE");
  int nDir = 18;
  TDirectory *theDirectory[nDir];
  theDirectory[0]  = HistoOutputFile->mkdir("No_cuts");
  theDirectory[1]  = HistoOutputFile->mkdir("After_jet_pt");
  theDirectory[2]  = HistoOutputFile->mkdir("After_bjet_pt");
  theDirectory[3]  = HistoOutputFile->mkdir("After_bjet_eta");
  theDirectory[4]  = HistoOutputFile->mkdir("After_electron");
  theDirectory[5]  = HistoOutputFile->mkdir("After_muon");
  theDirectory[6]  = HistoOutputFile->mkdir("After_DeltaPhi_el");
  theDirectory[7]  = HistoOutputFile->mkdir("After_DeltaPhi_mu");
  theDirectory[8]  = HistoOutputFile->mkdir("After_met_only_el");
  theDirectory[9]  = HistoOutputFile->mkdir("After_met_only_mu");
  theDirectory[10]  = HistoOutputFile->mkdir("After_mt_only_ele");
  theDirectory[11]  = HistoOutputFile->mkdir("After_mt_only_mu");
  theDirectory[12]  = HistoOutputFile->mkdir("After_mt_top_only_ele");
  theDirectory[13]  = HistoOutputFile->mkdir("After_mt_top_only_mu");
  theDirectory[14]  = HistoOutputFile->mkdir("Hadronic_After_onebjet");
  theDirectory[15]  = HistoOutputFile->mkdir("Hadronic_After_lepton_veto"); 
  theDirectory[16]  = HistoOutputFile->mkdir("Hadronic_After_met");
  theDirectory[17]  = HistoOutputFile->mkdir("Hadronic_After_masstop");
  
  PhenoAnalysis BSM_analysis(chain, HistoOutputFile, theDirectory, nDir);
  
}

using namespace std;
PhenoAnalysis::PhenoAnalysis(TChain& chain, TFile* theFile, TDirectory *cdDir[], int nDir)
{
  ifstream inFile;
  inFile.open ("config.in", ios::in);
  
  if (!inFile)
    {
      cerr << "ERROR: Can't open input file: " << endl;
      exit (1);
    }
  
  string inputType = "";
  
  //This set of lines are used to open and read the "config.in" file. 
  /////////////////////////////////////////////////////////////////////// 
  TEnv *params = new TEnv ("config_file");
  params->ReadFile ("config.in", kEnvChange);
  
  double DR_jet_elec_max  = params->GetValue ("DR_jet_elec_max", 0.);
  double DR_jet_muon_max  = params->GetValue ("DR_jet_muon_max", 0.); 
  double DR_elec_muon_max = params->GetValue ("DR_elec_muon_max", 0.);
  double jet_pt_min       = params->GetValue ("jet_pt_min", 0.);
  double jet_eta_max      = params->GetValue ("jet_eta_max", 0.);
  double b_pt_min         = params->GetValue ("b_pt_min", 0.);
  double b_eta_max        = params->GetValue ("b_eta_max", 0.);
  double elec_pt_min      = params->GetValue ("elec_pt_min", 0.);
  double muon_pt_min      = params->GetValue ("muon_pt_min", 0.);
  double met_min          = params->GetValue ("met_min", 0.);
  double mt_min           = params->GetValue ("mt_min", 0.);
  double mt_top_min       = params->GetValue ("mt_top_min", 0.);
  double blpt_min         = params->GetValue ("blpt_min", 0.);
  
  crateHistoMasps(nDir);
  
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET"); 
 // TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
 // TClonesArray *branchTrack = treeReader->UseBranch("Track");
//  TClonesArray *branchTower = treeReader->UseBranch("Tower");

  MissingET *METpointer; 
  Double_t electron_mass = 0.000510998902;
  Double_t muon_mass = 0.105658369;  
  
  //Loop over entries 
//for(Int_t entry = 0; entry < 4; ++entry)
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
    {
      treeReader->ReadEntry(entry);
      //Missing ET
      METpointer = (MissingET*) branchMissingET->At(0);
      double MET = METpointer->MET;
      double MET_phi = METpointer->Phi;
      
      TLorentzVector jet_i(0., 0., 0., 0.);
      TLorentzVector elec_1(0., 0., 0., 0.);
      TLorentzVector muon_1(0., 0., 0., 0.);
      TLorentzVector bjet_i(0.,0.,0.,0.);
      TLorentzVector Tau_1(0.,0.,0.,0.);
      TLorentzVector Tau_2(0.,0.,0.,0.);
      TLorentzVector jet_1(0., 0., 0., 0.);
      TLorentzVector jet_2(0., 0., 0., 0.);
      TLorentzVector jet_3(0., 0., 0., 0.);
      TLorentzVector jet_4(0., 0., 0., 0.);
      TLorentzVector jet_5(0., 0., 0., 0.);
      TLorentzVector jet_sum(0., 0., 0., 0.);      
      
      TLorentzVector top_e_i(0.,0.,0.,0.);
      TLorentzVector top_m_i(0.,0.,0.,0.);
     
      TLorentzVector top_gen(0.,0.,0.,0.);      
      
      int pass_cuts[nDir] = {0};
      double p_visible_w[2]={0., 0.};
      double p_visible[2]={0., 0.}; 
      
      vector<int> tau_index;
      vector<TLorentzVector> vec_jet(0);
      
      //Particle counter   0=njets, 1=nelectrons, 2=nmuons, 3=bjets, 4=taus
      int particle_counter[5] = {0};
      //Particle identificators
      bool is_jet_elec_overlap = false;
      bool is_jet_muon_overlap = false;
      bool is_b_jet = false;
      bool filled_tau1 = false;
      bool filled_tau2 = false;
      double jets_energy[5] = {0.};
      double jets_mass[5] = {0.};
      double bjet_energy = 0;
      double bjet_mass = 0;
      bool filled_first_jet = false;
      bool filled_second_jet = false;
      bool filled_third_jet = false;
      bool filled_fourth_jet = false;
      bool filled_fifth_jet = false;
      bool passed_mass_w = false;
      bool passed_mass_t = false;
      bool itisbetter = false;
      //Soft initial cuts
      double jet_ipt = 10.0;
      double jet_ieta = 5.0;
      double elec_1pt = 1.0;
      double elec_1eta = 2.1;
      double muon_1pt = 1.0;
      double muon_1eta = 2.1;
      double bjet_pt = 10.0;
      double bjet_eta = 5.0;
      double tau_pt = 10.0;
      double tau_eta = 5.0;
      double jet1_pt = 10.0;
      double jet1_eta = 5.0;
      double jet2_pt = 10.0;
      double jet2_eta = 5.0;
      // Transverse mass
      double MT_el = 0.0;
      double MT_mu = 0.0;
      // Transverse top mass
      double MT_top_el = 0.0;
      double MT_top_mu = 0.0;
      //Invariant mass
      double delta_eta_jj=0.;
      double delta_eta_bj=0.;
      double fmin_deltam = 0.;
      double mass_w = 0.;
      double mass_wfinal = 0.;
      double mass_t = 0.;
      double energy_t = 0.;
      double sum_energyjjb = 0;
      double sum_energy_jjb = 0;
      
      double hbratio = 0.;
      double mass_j1 = 0.;
      double mass_j2 = 0.;
      double energy_j1 = 0.;
      double energy_j2 = 0.;
      double pt_j1 = 0.;
      double pt_j2 = 0.;
      
      int Nbjets = 0;
      int Njets  = 0;
      int NCjets = 0;
      int jet_index = 0;
     
      GenParticle *particle;
      TObject *object;
      Track *track;
      Tower *tower;
      int Numberjets = 0;
      //Loop over jets    
      for (int j = 0; j < branchJet->GetEntriesFast(); j++)
	{ 
      passed_mass_t = false;
      is_jet_elec_overlap = false;
	  is_jet_muon_overlap = false;
 	  is_b_jet = false;
	 
      Jet *jet = (Jet*) branchJet->At(j);

      Numberjets++;
      _hmap_jet_consti[0]->Fill(jet->Constituents.GetEntriesFast());  
      _hmap_jet_consti2[0]->Fill(Numberjets, jet->Constituents.GetEntriesFast());
      


       /* // Loop over all jet's constituents          
           for (int k = 0; k < jet->Constituents.GetEntriesFast(); ++k){
            object = jet->Constituents.At(k);
            // Check if the constituent is accessible
            if(object == 0) continue;
                if(object->IsA() == GenParticle::Class()) {
                 particle = (GenParticle*) object;
                 double PID = TMath::Abs(particle->PID);
                 _hmap_jet_consti[0]->Fill(PID);
                 }
                else if(object->IsA() == Track::Class())
                        {
                        track = (Track*) object;
                        double PID_1 = TMath::Abs(track->PID);
                        _hmap_jet_consti[1]->Fill(PID_1);
                        }
                else if(object->IsA() == Tower::Class())
                        {
                        tower = (Tower*) object;
                        _hmap_jet_consti[2]->Fill(tower->Eta);
                        }
:
            }             
            
	   
	   for(int i = 0; i < branchGenParticle->GetEntriesFast(); i++){
	     GenParticle *t = (GenParticle*) branchGenParticle->At(i);
	        if(TMath::Abs(t->PID) == 6){
	        double top_energy = calculateE(t->Eta, t->PT, t->Mass);
	        top_gen.SetPtEtaPhiE(t->PT, t->Eta, t->Phi, top_energy);
            _hmap_jet_consti[3]->Fill(top_gen.Pt());   
           
           }
	   }    */                                                                                            

      double jet_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
      if((jet->PT > jet_ipt) && (TMath::Abs(jet->Eta) < jet_ieta)){
	particle_counter[0]++;
	jet_i.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
        
	// Find if there are electrons overlapping with the jets
	for (int e = 0; e < branchElectron->GetEntriesFast(); e++){
	  Electron *electron = (Electron*) branchElectron->At(e);
	  double electron_energy = calculateE(electron->Eta, electron->PT, electron_mass);
	  elec_1.SetPtEtaPhiE(electron->PT, electron->Eta, electron->Phi, electron_energy);
	  double DR_jet_electron = jet_i.DeltaR(elec_1);
	  if ((electron->PT > elec_1pt) && (TMath::Abs(electron->Eta) < elec_1eta)){
	    if ( DR_jet_electron < DR_jet_elec_max ){
	      is_jet_elec_overlap = true;
	      elec_1.SetPtEtaPhiE(0.0, 0.0, 0.0, 0.0);
	      break;} 
	    else{
	      elec_1pt = electron->PT; 
	      particle_counter[1]++;
	      p_visible_w[0]=(electron->PT)*(cos(electron->Phi));
	      p_visible_w[1]=(electron->PT)*(sin(electron->Phi));
	      p_visible[0]=(electron->PT)*(cos(electron->Phi));
	      p_visible[1]=(electron->PT)*(sin(electron->Phi));
	      elec_1.SetPtEtaPhiE(electron->PT, electron->Eta, electron->Phi, electron_energy);
	    }
	  }
	}
	
	// Find if there are muons overlapping with the jets
	for (int m = 0; m < branchMuon->GetEntriesFast(); m++){
	  Muon *muon = (Muon*) branchMuon->At(m);
	  double muon_energy = calculateE(muon->Eta, muon->PT, muon_mass);
	  muon_1.SetPtEtaPhiE(muon->PT, muon->Eta, muon->Phi, muon_energy);
	  double DR_jet_muon = jet_i.DeltaR(muon_1);
	  if ((muon->PT > muon_1pt) && (TMath::Abs(muon->Eta) < muon_1eta)){ 
	    if (DR_jet_muon < DR_jet_muon_max){ 
	      is_jet_muon_overlap = true;
	      muon_1.SetPtEtaPhiE(0.0, 0.0, 0.0, 0.0);
	      break;} 
	    else{
	      particle_counter[2]++;
	      muon_1pt = muon->PT;
	      p_visible_w[0]=(muon->PT)*(cos(muon->Phi));
	      p_visible_w[1]=(muon->PT)*(sin(muon->Phi));
	      p_visible[0]=(muon->PT)*(cos(muon->Phi));
	      p_visible[1]=(muon->PT)*(sin(muon->Phi));
	      muon_1.SetPtEtaPhiE(muon->PT, muon->Eta, muon->Phi, muon_energy);
	    }
	  }
	}
      } //soft jets cuts
      
      // If there no jets overlapping with electrons or muons, then fill the info for the b-jets....
      if((!is_jet_elec_overlap) && (!is_jet_muon_overlap)){
	//--------------------------------------------------------------------
	if ((jet->PT > bjet_pt) && (TMath::Abs(jet->Eta) < bjet_eta) && (jet->BTag == 1)){
          bjet_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
          bjet_mass = jet->Mass;
          bjet_i.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, bjet_energy);
          Nbjets++;
	  particle_counter[3]++;
          is_b_jet = true;
	}
	//--------------------------------------------------------------------
	if ((jet->TauTag == 1) && (jet->BTag == 0)){
	  if ((jet->PT > jet_pt_min) && (abs(jet->Eta) < jet_eta_max)){Njets++;}
	  particle_counter[4]++;
	  if ((jet->PT > tau_pt) && (TMath::Abs(jet->Eta) < tau_eta))
	    if (!filled_tau1){
	      Tau_1.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
	      filled_tau1 = true;
	      continue;                  
	    }
	  if (filled_tau1 && !filled_tau2){
	    Tau_2.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
	    filled_tau2 = true;
	  }
	}
        
	//--------------------------------------------------------------------
	if ((jet->TauTag == 0) && ((jet->BTag == 0))){
	  if ((jet->PT > jet_pt_min) && (abs(jet->Eta) < jet_eta_max)){Njets++;}
	  if ((jet->PT > jet1_pt) && (TMath::Abs(jet->Eta) < jet1_eta)){
	    if((filled_first_jet == false) && (filled_second_jet == false) && (filled_third_jet == false) && 
	       (filled_fourth_jet == false) && (filled_fifth_jet == false)){
	      jets_energy[0]=calculateE(jet->Eta, jet->PT, jet->Mass);
	      jets_mass[0]=jet->Mass;
	      //cout << jet->Mass << endl;
	      jet_1.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jets_energy[0]);
	      vec_jet.push_back(jet_1);
	      filled_first_jet = true;
	      continue;
	    }
	    if((filled_first_jet == true) && (filled_second_jet == false) && (filled_third_jet == false) && 
	       (filled_fourth_jet == false) && (filled_fifth_jet == false)){
	      jets_energy[1] = calculateE(jet->Eta, jet->PT, jet->Mass);
	      jets_mass[1]=jet->Mass;
	      jet_2.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jets_energy[1]);
	      vec_jet.push_back(jet_2);
	      filled_second_jet = true;
	      continue;
	    }
	    if((filled_first_jet == true) && (filled_second_jet == true) && (filled_third_jet == false) && 
	       (filled_fourth_jet == false) && (filled_fifth_jet == false)){
	      jets_energy[2] = calculateE(jet->Eta, jet->PT, jet->Mass);
	      jets_mass[2]=jet->Mass;
	      jet_3.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jets_energy[2]);
	      vec_jet.push_back(jet_3);
	      filled_third_jet = true;
	      continue;
	    }
	    if((filled_first_jet == true) && (filled_second_jet == true) && (filled_third_jet == true) && 
	       (filled_fourth_jet == false) && (filled_fifth_jet == false)){
	      jets_energy[3] = calculateE(jet->Eta, jet->PT, jet->Mass);
	      jets_mass[3]=jet->Mass;
	      jet_4.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jets_energy[3]);
	      vec_jet.push_back(jet_4);
	      filled_fourth_jet = true;
	      continue;
	    }
	    if((filled_first_jet == true) && (filled_second_jet == true) && (filled_third_jet == true) && 
	       (filled_fourth_jet == true) && (filled_fifth_jet == false)){
	      jets_energy[4] = calculateE(jet->Eta, jet->PT, jet->Mass);
	      jets_mass[4]=jet->Mass;
	      jet_5.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jets_energy[4]);
	      vec_jet.push_back(jet_5);
	      filled_fifth_jet = true;
	      continue;
	    }
	  }
	}// it's not b or tau 
	
        // At least two jets in the event
	if ( vec_jet.size() >= 2 ){  
	  // Combinatorial calculation
	  double M_jj[5][5] = {0., 0.};
	  int final_pair[1][2] = {0, 0};
	  for (int u = 0; u < vec_jet.size(); ++u){
	    for (int v = u+1; v < vec_jet.size(); ++v){
	      double DR_jet = vec_jet[u].DeltaR(vec_jet[v]); 
	      if(DR_jet > 0.3){
		    double min_deltam = 60.;
		    double sum_energyjj = jets_energy[u]+jets_energy[v];          
		    double p12 = CalculateP12_square_2particle(vec_jet[u].Eta(), vec_jet[u].Pt(), vec_jet[v].Eta(), vec_jet[v].Pt());           
		    M_jj[u][v]= TMath::Sqrt(TMath::Power((sum_energyjj),2)-p12);
            //M_jj[u][v]= (vec_jet[u] + vec_jet[v]).M;
            double deltam = TMath::Abs(M_jj[u][v]-80.4);
		    if (deltam < min_deltam){
		      min_deltam = deltam;
		      mass_wfinal = M_jj[u][v];
		      final_pair[0][0]=u;
		      final_pair[0][1]=v;
		    }
	      }
	    }
	  }
	  
	  
	  mass_w = mass_wfinal; 
	  if(is_b_jet){
	    int p = final_pair[0][0];
	    int q = final_pair[0][1];
	    energy_j1 = jets_energy[p];
	    energy_j2 = jets_energy[q];
	    pt_j1 = vec_jet[p].Pt();
	    pt_j2 = vec_jet[q].Pt();
	    double p_j1 = CalculateP(vec_jet[p].Eta(), vec_jet[p].Pt());
	    mass_j1 = TMath::Sqrt(TMath::Power(jets_energy[p],2)-p_j1);
	    double p_j2 = CalculateP(vec_jet[q].Eta(), vec_jet[q].Pt());
	    mass_j2 = TMath::Sqrt(TMath::Power(jets_energy[q],2)-p_j2);
	    
	    sum_energyjjb = jets_energy[p]+jets_energy[q]+bjet_energy;
	    double p123 = CalculateP123_square_3particle(vec_jet[p].Eta(), vec_jet[p].Phi(), vec_jet[p].Pt(), vec_jet[q].Eta(), vec_jet[q].Phi(), vec_jet[q].Pt(), bjet_i.Eta(), bjet_i.Phi(), bjet_i.Pt());
	    double M_jjb= TMath::Sqrt(TMath::Power((sum_energyjjb),2)-p123);
	    mass_t = M_jjb;
	    hbratio = bjet_energy/sum_energyjjb;
	    energy_t = sum_energyjjb;
	    if (mass_t <= 450. ){passed_mass_t = true;}
	  }
	  
	  
	  
	}// finish vector size >=2
	
	
      } //no muon or electron overlap
      
      
	}    //fisnish loop over branch jets 
      //Other parameters
      double bpt=0.;
      double bpluslep=0.;
      double jet1_met_dphi=0.;
      double jet2_met_dphi=0.;
      jet1_met_dphi=normalizedDphi(jet_1.Phi()-MET_phi);
      jet2_met_dphi=normalizedDphi(jet_2.Phi()-MET_phi);
      //Calculating transverse mass
      double dphi_met_ele=normalizedDphi(elec_1.Phi()-MET_phi);
      double dphi_met_mu=normalizedDphi(muon_1.Phi()-MET_phi);
      MT_el = transverseMass(elec_1.Pt(), MET, dphi_met_ele);    
      MT_mu = transverseMass(muon_1.Pt(), MET, dphi_met_mu);
      //Calculation transverse top mass 
      top_e_i = bjet_i+elec_1;
      top_m_i = bjet_i+muon_1;
      double dphi_met_top_e = normalizedDphi(top_e_i.Phi()-MET_phi);
      double dphi_met_top_m = normalizedDphi(top_m_i.Phi()-MET_phi);   
      MT_top_el = transverseMass(top_e_i.Pt(), MET, dphi_met_top_e);
      MT_top_mu = transverseMass(top_m_i.Pt(), MET, dphi_met_top_m);
      
      
      p_visible_w[0]+=(MET)*(cos(MET_phi));
      p_visible_w[1]+=(MET)*(sin(MET_phi));
      p_visible[0]+=(bjet_i.Pt())*(cos(bjet_i.Phi()));
      p_visible[1]+=(bjet_i.Pt())*(sin(bjet_i.Phi()));
      
      double top_x = p_visible_w[0]+(bjet_i.Pt())*(cos(bjet_i.Phi()));
      double top_y = p_visible_w[1]+(bjet_i.Pt())*(sin(bjet_i.Phi()));
      double pT_top = sqrt(top_x*top_x + top_y*top_y);
      
      double pT_W = sqrt(p_visible_w[0]*p_visible_w[0]+p_visible_w[1]*p_visible_w[1]);
      double pT_bl = sqrt(p_visible[0]*p_visible[0]+p_visible[1]*p_visible[1]);
      
      double top_e = sqrt(pT_top*pT_top + 173.34*173.34);
      //cout <<"top_e "<<bjet_i.E()/top_e<<endl; 
      //cut
      double dphi_bjet_ele=normalizedDphi(elec_1.Phi()-bjet_i.Phi());
      double dphi_bjet_muon=normalizedDphi(muon_1.Phi()-bjet_i.Phi()); 
      
      
      pass_cuts[0]=1;
      
      if(Njets == 0){pass_cuts[1] = 1;}
      if((bjet_i.Pt() > b_pt_min) && (Nbjets ==1) && (pass_cuts[1] == 1)){pass_cuts[2] = 1;}
      if((TMath::Abs(bjet_i.Eta()) < jet_eta_max) && (pass_cuts[2] == 1)){pass_cuts[3] = 1;}
      if((elec_1.Pt() > elec_pt_min) && (abs(elec_1.Eta()) < elec_1eta) && (muon_1.Pt()==0.) && (particle_counter[1]==1) 
	 && (pass_cuts[3] == 1) && (pT_W > 50.) && (!filled_tau1) && (!filled_tau2)){
	
	bpluslep=bjet_i.Pt()+elec_1.Pt();
	bpt = bjet_i.Pt()/bpluslep;
	pass_cuts[4] = 1;
      }  
      if((muon_1.Pt() > muon_pt_min) && (abs(muon_1.Eta()) < muon_1eta) && (elec_1.Pt()==0.) && (particle_counter[2]==1) 
	 && (pass_cuts[3] == 1) && (pT_W > 50.) && (!filled_tau1) && (!filled_tau2)){
	
	bpluslep=bjet_i.Pt()+muon_1.Pt();
	bpt = bjet_i.Pt()/bpluslep;
	pass_cuts[5] = 1;
      } 
      if ( (pass_cuts[4] == 1) && (abs(dphi_bjet_ele) < 1.7) ){pass_cuts[6] =1;} 
      if ( (pass_cuts[5] == 1) && (abs(dphi_bjet_muon) < 1.7) ){pass_cuts[7] =1;} 
      if ((MET > met_min) && (pass_cuts[6] == 1)){pass_cuts[8] = 1;}
      if ((MET > met_min) && (pass_cuts[7] == 1)){pass_cuts[9] = 1;}
      if ((MT_el > mt_min) &&  (pass_cuts[8] == 1) ){pass_cuts[10] = 1;}
      if ((MT_mu > mt_min) &&  (pass_cuts[9] == 1) ){pass_cuts[11] = 1;}
      if ((MT_top_el > mt_top_min) &&  (pass_cuts[8] == 1) ){pass_cuts[12] = 1;}
      if ((MT_top_mu > mt_top_min) &&  (pass_cuts[9] == 1) ){pass_cuts[13] = 1;}
      //Hadronic Monotop
      if((bjet_i.Pt() > b_pt_min) && (TMath::Abs(bjet_i.Eta()) < jet_eta_max) && (Nbjets == 1)){pass_cuts[14] = 1;}      
      if((muon_1.Pt() == 0) && (elec_1.Pt() == 0) && (!filled_tau1) && (!filled_tau2) && (pass_cuts[14] == 1)){pass_cuts[15] = 1;} 
      if((MET > met_min) && (pass_cuts[15] == 1)){pass_cuts[16] = 1;}
      if((passed_mass_t)  && (pass_cuts[16] == 1)){pass_cuts[17] = 1;}  
      
      
      // save some important histograms
      for (int i = 0; i < nDir; i++){
        _hmap_Nevents[i]->Fill(0.0);
        if (pass_cuts[i] == 1){
          _hmap_Nevents[i]->Fill(1.0);
	      _hmap_N_jets[i]->Fill(particle_counter[0]);
          _hmap_N_elec[i]->Fill(particle_counter[1]);
          _hmap_N_muon[i]->Fill(particle_counter[2]);
	      _hmap_N_bjet[i]->Fill(particle_counter[3]);
          _hmap_N_tau[i]->Fill(particle_counter[4]);
          _hmap_pt_w[i]->Fill(pT_W);    
          _hmap_pt_t[i]->Fill(pT_top);
          if (top_e > 0)_hmap_e_tb[i]->Fill(bjet_i.E()/top_e);
          _hmap_e_t[i]->Fill(top_e);
          _hmap_pt_bl[i]->Fill(pT_bl);
	  
	  if(mass_w > 0.){
	    _hmap_invariant_mjj[i]->Fill(mass_w);      
	  }
	  if(mass_t > 0.){
	    _hmap_invariant_mjjb[i]->Fill(mass_t);
	    _hmap_top_energy[i]->Fill(energy_t);
	    _hmap_Hbjet_ratio_pT[i]->Fill(hbratio);
	    if (energy_j1 > 0.){
	      _hmap_invariant_mj1[i]->Fill(mass_j1);
	      _hmap_invariant_ej1[i]->Fill(energy_j1);
	      _hmap_invariant_pj1[i]->Fill(pt_j1);
	    }
	    if (energy_j2 > 0.){
	      _hmap_invariant_mj2[i]->Fill(mass_j2);
	      _hmap_invariant_ej2[i]->Fill(energy_j2);
	      _hmap_invariant_pj2[i]->Fill(pt_j2);
	    }
	    if((energy_j1 > 0.) && (energy_j2 > 0.)){
	      _hmap_invariant_ej1j2[i]->Fill(energy_j1+energy_j2);
	    }
	  }
	  
	  
	  if(jet_1.Pt()>0.){ 
	    _hmap_jet1_pT[i]->Fill(jet_1.Pt());
	    _hmap_jet1_eta[i]->Fill(jet_1.Eta());
	    _hmap_jet1_phi[i]->Fill(jet_1.Phi());
            _hmap_jet_1_met_dphi[i]->Fill(jet1_met_dphi);
            _hmap_met_jet1_met_dphi[i]->Fill(jet1_met_dphi,MET);
	  }
	  
          if(jet_2.Pt()>0.){
            _hmap_jet2_pT[i]->Fill(jet_2.Pt());
            _hmap_jet2_eta[i]->Fill(jet_2.Eta());
            _hmap_jet2_phi[i]->Fill(jet_2.Phi());
            _hmap_jet_2_met_dphi[i]->Fill(jet2_met_dphi);
          }
	  
	  if(Tau_1.Pt()>0.){
            _hmap_tau1_pT[i]->Fill(Tau_1.Pt());
            _hmap_tau1_eta[i]->Fill(Tau_1.Eta());
            _hmap_tau1_phi[i]->Fill(Tau_1.Phi());
	  }
	  if(Tau_2.Pt()>0.){
            _hmap_tau2_pT[i]->Fill(Tau_2.Pt());
            _hmap_tau2_eta[i]->Fill(Tau_2.Eta());
            _hmap_tau2_phi[i]->Fill(Tau_2.Phi());
          }
	  if(bjet_i.Pt()>0.){
	    _hmap_bjet_ratio_pT[i]->Fill(bpt);
            _hmap_bjet_plus_lep[i]->Fill(bpluslep);	
	    _hmap_bjet_pT[i]->Fill(bjet_i.Pt());
	    _hmap_bjet_eta[i]->Fill(bjet_i.Eta());
	    _hmap_bjet_phi[i]->Fill(bjet_i.Phi());
	    _hmap_bjet_energy[i]->Fill(bjet_energy);
	  }
	  if (elec_1.Pt() > 0.0) { 
            _hmap_mt[i]->Fill(MT_el);
            _hmap_mt_top[i]->Fill(MT_top_el);
	    _hmap_elec_pT[i]->Fill(elec_1.Pt());
	    _hmap_elec_eta[i]->Fill(elec_1.Eta());
	    _hmap_elec_phi[i]->Fill(elec_1.Phi());
	  }
	  if (muon_1.Pt() > 0.0) {
            _hmap_mt[i]->Fill(MT_mu);
            _hmap_mt_top[i]->Fill(MT_top_mu);
            _hmap_mt_mt_top[i]->Fill(MT_mu, MT_top_mu); 
	    _hmap_muon_pT[i]->Fill(muon_1.Pt());
	    _hmap_muon_eta[i]->Fill(muon_1.Eta());
	    _hmap_muon_phi[i]->Fill(muon_1.Phi());
	  }
	  _hmap_met[i]->Fill(MET);
	} 
      }
      
    }//Finish loop over entries
  
  
  theFile->cd();
  for (int d = 0; d < nDir; d++)
    {
      cdDir[d]->cd();
      _hmap_Nevents[d]->Write();
      _hmap_N_jets[d]->Write();
      _hmap_jet1_pT[d]->Write();
      _hmap_jet1_eta[d]->Write();
      _hmap_jet1_phi[d]->Write();
      _hmap_jet2_pT[d]->Write();
      _hmap_jet2_eta[d]->Write();
      _hmap_jet2_phi[d]->Write();
      _hmap_N_tau[d]->Write();
      _hmap_tau1_pT[d]->Write();
      _hmap_tau1_eta[d]->Write();
      _hmap_tau1_phi[d]->Write();
      _hmap_tau2_pT[d]->Write();
      _hmap_tau2_eta[d]->Write();
      _hmap_tau2_phi[d]->Write();
      _hmap_N_bjet[d]->Write();
      _hmap_bjet_plus_lep[d]->Write();
      _hmap_bjet_ratio_pT[d]->Write();
      _hmap_Hbjet_ratio_pT[d]->Write();
      _hmap_bjet_pT[d]->Write();
      _hmap_bjet_eta[d]->Write();
      _hmap_bjet_phi[d]->Write();
      _hmap_bjet_energy[d]->Write();
      _hmap_top_energy[d]->Write();
      _hmap_N_elec[d]->Write();
      _hmap_elec_pT[d]->Write();
      _hmap_elec_eta[d]->Write();
      _hmap_elec_phi[d]->Write();
      _hmap_N_muon[d]->Write();
      _hmap_muon_pT[d]->Write();
      _hmap_muon_eta[d]->Write();
      _hmap_muon_phi[d]->Write();
      _hmap_met[d]->Write();
      _hmap_mt[d]->Write();
      _hmap_mt_top[d]->Write();
      _hmap_pt_w[d]->Write();
      _hmap_pt_t[d]->Write();
      _hmap_e_tb[d]->Write();
      _hmap_e_t[d]->Write();
      _hmap_pt_bl[d]->Write();
      _hmap_jet_1_met_dphi[d]->Write();
      _hmap_jet_2_met_dphi[d]->Write();
      _hmap_invariant_mjjb[d]->Write();
      _hmap_invariant_mjj[d]->Write();
      _hmap_invariant_mj1[d]->Write();
      _hmap_invariant_mj2[d]->Write();
      _hmap_invariant_ej1[d]->Write();
      _hmap_invariant_ej2[d]->Write();
      _hmap_invariant_ej1j2[d]->Write();
      _hmap_invariant_pj1[d]->Write();
      _hmap_invariant_pj2[d]->Write();
      _hmap_met_jet1_met_dphi[d]->Write();   
      _hmap_mt_mt_top[d]->Write(); 
      _hmap_jet_consti[d]->Write();
      _hmap_jet_consti2[d]->Write();
      _hmap_top_jet[d]->Write();
    }
  
  theFile->Close();
  
} //finish class

PhenoAnalysis::~PhenoAnalysis()
{
	      
// do anything here that needs to be done at desctruction time
}

double PhenoAnalysis::calculateE(double eta, double pt, double mass){

  double theta = TMath::ATan(TMath::Exp(-eta)); 
  double sin_theta = TMath::Sin(2*theta);
  double p= pt/sin_theta;
  double e = sqrt(pow(p, 2) + pow(mass, 2));
  
  return e;
  
}


double PhenoAnalysis::normalizedDphi(double phi){
  double twoPI = 2.0*TMath::Pi();
  if ( phi < -TMath::Pi() ){phi += twoPI;}
  if ( phi > TMath::Pi() ){phi = twoPI-phi;}
  else phi = TMath::Abs(phi); 
  return phi;
}

double PhenoAnalysis::transverseMass(double ptl, double met, double deltaphi){
  
  double Mt = TMath::Sqrt(2*ptl*met*(1-cos(deltaphi))); 
  return Mt;
}


// approximate Invariant mass of two jets

double PhenoAnalysis::calculateM(double deta, double pt1, double pt2){
  double cosh_deta = TMath::CosH(deta);
  double M = sqrt(2*pt1*pt2*cosh_deta);
  return M;
}


double PhenoAnalysis::CalculateP(double Eta1, double PT1){
double theta1 = 2*TMath::ATan(TMath::Exp(-Eta1));
double pz1 = PT1*TMath::Cos(theta1)/TMath::Sin(theta1);
double p2  = TMath::Power((PT1),2)+TMath::Power((pz1),2);
return p2;
}



double PhenoAnalysis::CalculateP12_square_2particle(double Eta1, double PT1, double Eta2, double PT2){
  double theta1 = 2*TMath::ATan(TMath::Exp(-Eta1));
  double theta2 = 2*TMath::ATan(TMath::Exp(-Eta2));
  double pz1 = PT1*TMath::Cos(theta1)/TMath::Sin(theta1);
  double pz2 = PT2*TMath::Cos(theta2)/TMath::Sin(theta2);
  double p2  = TMath::Power((PT1+PT2),2)+TMath::Power((pz1+pz2),2);
  return p2;
}

double PhenoAnalysis::CalculateP123_square_3particle(double Eta1, double Phi1, double PT1, double Eta2, double Phi2, double PT2, double Eta3, double Phi3, double PT3){
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
  


void PhenoAnalysis::crateHistoMasps (int directories)
{
  for (int i = 0; i < directories; i++)
    {
      _hmap_Nevents[i]		= new TH1F("Nevents", "Nevents", 5,-0.5,4.5); 
      
      _hmap_N_jets[i]		= new TH1F("N_jets",  "N_jets", 11, -0.5, 10.5);    
      
      _hmap_jet1_pT[i]		= new TH1F("jet1_pT",  "P_{T} jet1 ", 100, 0., 1000.);
      _hmap_jet1_eta[i]		= new TH1F("jet1_eta", "#eta jet1 ", 50, -5.0, 5.0);
      _hmap_jet1_phi[i]		= new TH1F("jet1_phi", "#phi jet1 ", 70, -3.6, 3.6);
      
      _hmap_jet2_pT[i]          = new TH1F("jet2_pT",  "P_{T} jet2 ", 100, 0., 1000.);
      _hmap_jet2_eta[i]         = new TH1F("jet2_eta", "#eta jet2 ", 50, -5.0, 5.0);
      _hmap_jet2_phi[i]         = new TH1F("jet2_phi", "#phi jet2 ", 70, -3.6, 3.6);
      
      _hmap_N_tau[i]		= new TH1F("N_tau",  "N_tau", 5, -0.5, 4.5);
      
      _hmap_tau1_pT[i]		= new TH1F("tau1_pT",  "p_{T}(#tau_{1})", 100, 0., 1000.);
      _hmap_tau1_eta[i]		= new TH1F("tau1_eta", "#eta(#tau_{1})", 50, -5.0, 5.0);
      _hmap_tau1_phi[i]		= new TH1F("tau1_phi", "#phi(#tau_{1})", 72, -3.6, 3.6);
      
      _hmap_tau2_pT[i]          = new TH1F("tau2_pT",  "p_{T}(#tau_{2})", 100, 0., 1000.);
      _hmap_tau2_eta[i]         = new TH1F("tau2_eta", "#eta(#tau_{2})", 50, -5.0, 5.0);
      _hmap_tau2_phi[i]         = new TH1F("tau2_phi", "#phi(#tau_{2})", 72, -3.6, 3.6);
      
      _hmap_bjet_plus_lep[i]    = new TH1F("bjet_plus_lep", "p_{T} bjet+lep", 50, 0., 500);
      
      _hmap_N_bjet[i]		= new TH1F("N_bjet",    "N_bjet", 5, -0.5, 4.5);
      _hmap_bjet_ratio_pT[i]	= new TH1F("bjet_ratio_pT", "p_{T} bjet", 50, 0, 1);
      _hmap_Hbjet_ratio_pT[i]    = new TH1F("Hbjet_ratio_pT", "p_{T} bjet ratio", 50, 0, 1);
       _hmap_bjet_pT[i]		= new TH1F("bjet_pT",  "P_{T} bjet ", 100, 0., 1000.);
      _hmap_bjet_eta[i]		= new TH1F("bjet_eta", "#eta bjet ", 50, -5.0, 5.0);
      _hmap_bjet_phi[i]		= new TH1F("bjet_phi", "#phi bjet ", 70, -3.6, 3.6);
      _hmap_bjet_energy[i]	= new TH1F("bjet_energy", "Energy b-quark", 150, 0., 1500.);
      _hmap_top_energy[i]   = new TH1F("top_energy", "Energy top-quark", 150, 0., 1500.);
      
      _hmap_N_elec[i]           = new TH1F("N_elec",    "N_elec", 11, -0.5, 10.5);
      _hmap_elec_pT[i]          = new TH1F("elec_pT", "p_{T} e ", 100, 0., 1000.);
      _hmap_elec_eta[i]         = new TH1F("elec_eta","#eta e ", 50, -5.0, 5.0);
      _hmap_elec_phi[i]         = new TH1F("elec_phi", "#phi e", 70, -3.6, 3.6);      
      
      _hmap_N_muon[i]           = new TH1F("N_muon",    "N_muon", 11, -0.5, 10.5);
      _hmap_muon_pT[i]          = new TH1F("muon_pT", "p_{T} #mu ", 100, 0., 1000.);
      _hmap_muon_eta[i]         = new TH1F("muon_eta","#eta #mu ", 50, -5.0, 5.0);
      _hmap_muon_phi[i]         = new TH1F("muon_phi", "#phi #mu", 70, -3.6, 3.6); 
      _hmap_met[i]		        = new TH1F("met", "MET", 200, 0. , 2000.);
      _hmap_mt[i]               = new TH1F("mt", "Transverse mass", 200, 0. , 2000.);
      _hmap_mt_top[i]           = new TH1F("mt_top", "Transverse top mass", 200, 0. , 2000.);	
      _hmap_pt_w[i]             = new TH1F("pT_w", "pT W", 200, 0. , 2000.);
      _hmap_e_tb[i]             = new TH1F("e_t_over_b", "E(b/t)", 30, 0. , 3.);
      _hmap_e_t[i]              = new TH1F("e_top", "E(t)", 200, 0. , 2000.);
      _hmap_pt_t[i]             = new TH1F("pT_t", "pT top", 200, 0. , 2000.);
      _hmap_pt_bl[i]            = new TH1F("pT_bl", "pT b+l", 200, 0. , 2000.); 
      _hmap_jet_1_met_dphi[i]   = new TH1F("jet_1_met_dphi", "#Delta #phi(jet1, MET)", 32, 0, 3.2); 
      _hmap_jet_2_met_dphi[i]   = new TH1F("jet_2_met_dphi", "#Delta #phi(jet2, MET)", 32, 0, 3.2);     
      _hmap_invariant_mjjb[i]   = new TH1F("invariant_mjjb", "invariant_mjjb", 150, 0, 1500);
      _hmap_invariant_mjj[i]    = new TH1F("invariant_mjj", "invariant_mjj", 150, 0, 1500);
      _hmap_invariant_mj1[i]    = new TH1F("invariant_mj1", "invariant_mj1", 150, 0, 1500);
      _hmap_invariant_pj1[i]    = new TH1F("invariant_pj1", "invariant_pj1", 150, 0, 1500);
      _hmap_invariant_ej1[i]    = new TH1F("invariant_ej1", "invariant_ej1", 150, 0, 1500);
      _hmap_invariant_mj2[i]    = new TH1F("invariant_mj2", "invariant_mj2", 150, 0, 1500);
      _hmap_invariant_pj2[i]    = new TH1F("invariant_pj2", "invariant_pj2", 150, 0, 1500);
      _hmap_invariant_ej2[i]    = new TH1F("invariant_ej2", "invariant_ej2", 150, 0, 1500);
      _hmap_invariant_ej1j2[i]  = new TH1F("invariant_ej1j2", "invariant_ej1j2", 150, 0, 1500);
      _hmap_met_jet1_met_dphi[i]= new TH2F("met_jet1_met_dphi", "#Delta #phi(jet1,MET), MET", 32, 0., TMath::Pi(), 50, 0., 500.);
      _hmap_mt_mt_top[i]        = new TH2F("mt_mt_top", "M_{T}, M_{T}top", 100, 0., 1000., 100, 0., 1000.);	
      _hmap_jet_consti[i]       = new TH1F("jet_consti", "jets constituens", 50, -0.5, 49.5); 
      _hmap_jet_consti2[i]      = new TH2F("jet_consti2", "Jets-constituens", 20, -0.5, 19.5, 50, -0.5, 49.5);
      _hmap_top_jet[i]          = new TH2F("top_jet", "top test", 100, 0., 1000., 100, 0., 1000.);
}
}

