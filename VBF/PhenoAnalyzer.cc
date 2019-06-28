////////////////////////////////////////////////////////////////
//                                                            //
// Author: Andrés Flórez, Universidad de los Andes, Colombia  //
//                                                            //
////////////////////////////////////////////////////////////////


#include <iostream>
#include "ROOTFunctions.h"
#include "PhenoAnalyzer.h"
//#include "DelphesFunctions.h"

//FILE *docDelphes;

int main(int argc, char *argv[]) {

//docDelphes = fopen("pruebaDelphes.txt","w");

  //TApplication app("App",&argc, argv);
  TChain chain("Delphes");
  chain.Add(argv[1]);
  TFile * HistoOutputFile = new TFile(argv[2], "RECREATE");
  int nDir = 28;
  TDirectory *theDirectory[nDir];
  theDirectory[0]  = HistoOutputFile->mkdir("No_cuts");
  theDirectory[1]  = HistoOutputFile->mkdir("After_tau_pt");
  theDirectory[2]  = HistoOutputFile->mkdir("After_b_jet_veto");
  theDirectory[3]  = HistoOutputFile->mkdir("After_MET");
  theDirectory[4]  = HistoOutputFile->mkdir("After_Ntau_1");
  theDirectory[5]  = HistoOutputFile->mkdir("After_Ntau_2");
  theDirectory[6]  = HistoOutputFile->mkdir("After_Ntau_3");
  theDirectory[7]  = HistoOutputFile->mkdir("After_OnlyJets_forward");
  theDirectory[8]  = HistoOutputFile->mkdir("After_OnlyJets_pt");
  theDirectory[9]  = HistoOutputFile->mkdir("After_OnlyJets_pt_tau1");
  theDirectory[10]  = HistoOutputFile->mkdir("After_OnlyJets_pt_tau2");
  theDirectory[11] = HistoOutputFile->mkdir("After_OnlyJets_pt_tau3");
  theDirectory[12]  = HistoOutputFile->mkdir("After_OnlyJets_eta");
  theDirectory[13]  = HistoOutputFile->mkdir("After_OnlyJets_eta_tau1");
  theDirectory[14]  = HistoOutputFile->mkdir("After_OnlyJets_eta_tau2");
  theDirectory[15]  = HistoOutputFile->mkdir("After_OnlyJets_eta_tau3");
  theDirectory[16]  = HistoOutputFile->mkdir("After_OnlyJets_backToback");
  theDirectory[17]  = HistoOutputFile->mkdir("After_OnlyJets_backToback_tau1");
  theDirectory[18]  = HistoOutputFile->mkdir("After_OnlyJets_backToback_tau2");
  theDirectory[19]  = HistoOutputFile->mkdir("After_OnlyJets_backToback_tau3");
  theDirectory[20]  = HistoOutputFile->mkdir("After_OnlyJets_deltaEta");
  theDirectory[21]  = HistoOutputFile->mkdir("After_OnlyJets_deltaEta_tau1");
  theDirectory[22]  = HistoOutputFile->mkdir("After_OnlyJets_deltaEta_tau2"); 
  theDirectory[23]  = HistoOutputFile->mkdir("After_OnlyJets_deltaEta_tau3"); 
  theDirectory[24] = HistoOutputFile->mkdir("After_OnlyJets_mjj");
  theDirectory[25] = HistoOutputFile->mkdir("After_Ntau_1_pass_VBF");
  theDirectory[26] = HistoOutputFile->mkdir("After_Ntau_2_pass_VBF");
  theDirectory[27] = HistoOutputFile->mkdir("After_Ntau_3_pass_VBF");
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

   // This set of lines are used to open and read the "config.in" file. 
   /////////////////////////////////////////////////////////////////////// 
   TEnv *params = new TEnv ("config_file");
   params->ReadFile ("config.in", kEnvChange);
   
   double lead_jet_pt      = params->GetValue ("lead_jet_pt", 100.);
   double slead_jet_pt     = params->GetValue ("slead_jet_pt",100.);
   double jet_eta_max      = params->GetValue ("jet_eta_max", 5.0);
   double jj_delta_eta_min = params->GetValue ("jj_delta_eta_min", 3.5);  
   double jj_mass_min      = params->GetValue ("jj_mass_min", 100.);  
   double tau1_pt_min      = params->GetValue ("tau1_pt_min", 10.);
   double tau2_pt_min      = params->GetValue ("tau2_pt_min", 10.);
   double tau3_pt_min      = params->GetValue ("tau3_pt_min", 10.);
   double tau_eta_max      = params->GetValue ("tau_eta_max", 2.5);
   double b_jet_pt_min     = params->GetValue ("b_jet_pt_min", 20.0);
   double DR_jet_elec_max  = params->GetValue ("DR_jet_elec_max", 0.3);
   double met_min          = params->GetValue ("met_min", 10.);
   double jet_eta_min      = params->GetValue ("jet_eta_min", 2.0);
    
   crateHistoMasps(nDir);

   TH1 *hmap_numero_taus = new TH1F("Numero_de_taus", "Ntaus", 5 ,0, 5);
   TH1 *hmap_numero_taus_gen = new TH1F("Numero_de_taus_GenParticle", "N taus GenParticle", 5 ,0, 5);
   //Prueba masa W:
   //TH1 *hmasa_w = new TH1F("masa w", "masa w", 32 ,0, 400 );
   TH1 *hmap_numero_elec= new TH1F("Numero_de_electrones", "Nelec",5,0,5);
  TH1 *hmap_numero_elec_gen= new TH1F("Numero_de_electrones GenParticle", "N elec GenParticle",5,0,5); 

  //TH1 *hmasam_w = new TH1F("masa de w", "masa de w", 32 ,0, 400 );
   TH1 *hmap_numero_muones= new TH1F("Numero_de_Muones", "Nmuones",5,0,5);
   TH1 *hmap_numero_muones_gen= new TH1F("Numero_de_Muones_GenParticle", "N muones GenParticle",5,0,5);

   ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
   Long64_t numberOfEntries = treeReader->GetEntries();
   
   TClonesArray *branchJet = treeReader->UseBranch("Jet");
   TClonesArray *branchElectron = treeReader->UseBranch("Electron");
   TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");  
   TClonesArray *branchMuon = treeReader->UseBranch("Muon");
   TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
   MissingET *METpointer; 
   
   //numberOfEntries = 10000;
   for(Int_t entry = 0; entry < numberOfEntries; entry++)
     {
       int pass_cuts[nDir] = {0};
       pass_cuts[0] = 1;
       TLorentzVector Jet_leading_vec(0.,0., 0., 0.);
       TLorentzVector Jet_subleading_vec(0., 0., 0., 0.);
       TLorentzVector Tau1_vec (0., 0., 0., 0.);
       TLorentzVector Tau2_vec (0., 0., 0., 0.);
       TLorentzVector Tau3_vec (0., 0., 0., 0.);
       vector<int> tau_index;
       bool is_b_jet = false;
       treeReader->ReadEntry(entry);
       METpointer = (MissingET*) branchMissingET->At(0);
       Double_t MET = METpointer->MET;
       int tau_counter=0;
       int elec_counter=0;
       int muon_counter=0;
       int tau_counter_gen=0;
       int elec_counter_gen=0;
       int muon_counter_gen=0;
       int jet_counter=0;
       double Mt; 
       double Mtm;
       
//CONTADORES DE MUONES,ELECTRONES Y TAUS DESDE GEN_PARTICLE
	for(int l=0; l < branchGenParticle->GetEntriesFast(); l++)
	{
	GenParticle *particle = (GenParticle*) branchGenParticle->At(l);
	GenParticle *particle_mother =(GenParticle*) branchGenParticle->At(l);
	GenParticle *particle_daughter = (GenParticle*) branchGenParticle->At(l);
	if((particle->PID == abs(15)))
		{
		tau_counter_gen++;
		}
	if((particle->PID == abs(11)))
		{
		elec_counter_gen++;
		}
     
        if((particle->PID == abs(13)))
                {
                muon_counter_gen++;
                }
	}	

       if(branchJet->GetEntries() > 0)
	 {
           // For Jets
           double jet_highest_pt = 0.;
           double jet_Secondhighest_pt = 0.;
           int lead_ljet_index = 10000;
	   
	   // if (branchJet->GetEntriesFast() < 3) continue;
	   bool is_jet_elec_overlap = false;
           TLorentzVector jet_i(0., 0. , 0., 0.);
           TLorentzVector elec_i(0., 0., 0., 0.);
	   
           for (int k= 0; k < branchMuon->GetEntriesFast(); k++){
	     Muon *muon = (Muon*) branchMuon->At(k);
	     //if ((muon->PT > 10.) && (abs(muon->Eta) < 2.5)){muon_counter++;}
	     muon_counter++;
	   }
	   
           for (int j = 0; j < branchJet->GetEntriesFast(); j++)
             {
               jet_counter++;
	       //cout << "numero de jets "<< jet_counter <<" en el evento " << entry <<endl;
               Jet *jet = (Jet*) branchJet->At(j);
	       double jet_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
               jet_i.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
               for (int el = 0; el < branchElectron->GetEntriesFast(); el++){
		 Electron *elec = (Electron*) branchElectron->At(el);
		 double elec_energy = calculateE(elec->Eta, elec->PT, 0.000510998902);
		 elec_i.SetPtEtaPhiE(elec->PT, elec->Eta, elec->Phi, elec_energy);
		 //if ((elec->PT > 10.) && (abs(elec->Eta)) < 2.5){elec_counter++;} 
		 elec_counter++;
		 double DR_jet_elec = jet_i.DeltaR(elec_i);
		 if (DR_jet_elec < DR_jet_elec_max){
		   is_jet_elec_overlap = true;
		   break;
		 }
		// if ((elec->PT > 10.) && (abs(elec->Eta)) < 2.5){elec_counter++;}
               }
	      if (!is_jet_elec_overlap){	       
               //if ((!is_jet_elec_overlap) && (muon_counter == 0) && (elec_counter == 0) ){
		 if ((jet->PT > b_jet_pt_min) && (jet->BTag == 1)){is_b_jet = true;}
                 cout << " TauTag " << jet->TauTag << endl;
		 if (jet->TauTag == 1){
	           tau_counter++;
		   if ((abs(jet->Eta) < tau_eta_max)){
		     tau_index.push_back(j);
			}
		 }
		 
		 if (jet->TauTag == 0){
		   if (jet->PT > jet_highest_pt){
                     Jet_leading_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy); 
		     jet_highest_pt = jet->PT; 
		     lead_ljet_index = j;
		   }
		 }
	       }
	     }
	   
           is_jet_elec_overlap = false;
           for (int j = 0; j < branchJet->GetEntriesFast(); j++)
             {
               Jet *jet = (Jet*) branchJet->At(j);
	       double jet_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
	       jet_i.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
               for (int el = 0; el < branchElectron->GetEntriesFast(); el++){
                 Electron *elec = (Electron*) branchElectron->At(el);
                 double elec_energy = calculateE(elec->Eta, elec->PT, 0.000510998902);
                 elec_i.SetPtEtaPhiE(elec->PT, elec->Eta, elec->Phi, elec_energy);
                 double DR_jet_elec = jet_i.DeltaR(elec_i);
                 if (DR_jet_elec < DR_jet_elec_max){
                   is_jet_elec_overlap = true;
                   break;
		 }
		 //if ((elec->PT > 10.) && (abs(elec->Eta)) < 2.5){elec_counter2++;}	      		
	       }
               if (!is_jet_elec_overlap){
	       //if ((!is_jet_elec_overlap) && (muon_counter == 0) && (elec_counter2 == 0) ){
		 if ((j != lead_ljet_index) && (jet->TauTag != 1)){
		   if (jet->PT > jet_Secondhighest_pt){
		     jet_Secondhighest_pt = jet->PT; 
		     Jet_subleading_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
		   }
                 }
	       }
	     }	
	   
	   if (tau_index.size() > 0)
	     {
	       for (int i = 0; i < tau_index.size(); i++)
		 {
		   int tau_i = tau_index.at(i);
		   Jet *jet = (Jet*) branchJet->At(tau_i); 
		   double jet_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
		   if (tau_index.size() == 1){
		     Tau1_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
		   }
		   if (tau_index.size() == 2){
		     if(i == 0) Tau1_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
		     if(i == 1) Tau2_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
		   }
		   if (tau_index.size() == 3){
		     if(i == 0) Tau1_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
		     if(i == 1) Tau2_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
		     if(i == 2) Tau3_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
		   }			 
		   //tau_counter++;	 
		 }
	     }
	   
	   // Apply central event selection criteria
	   bool pass_lead_jet_cuts = false;
	   bool pass_sublead_jet_cuts = false;
	   
           if ((Tau1_vec.Pt() > tau1_pt_min) || (Tau2_vec.Pt() > tau2_pt_min) || (Tau3_vec.Pt() > tau3_pt_min)){
	     pass_cuts[1] = 1; 
	   }
	   if ((pass_cuts[1] == 1) && (!is_b_jet)){pass_cuts[2] = 1;}
	   
           if ((pass_cuts[2] == 1) && (MET > met_min)){pass_cuts[3] = 1;}
	   if ((pass_cuts[3] == 1) && (tau_index.size() == 1)){pass_cuts[4] = 1;}
	   if ((pass_cuts[3] == 1) && (tau_index.size() == 2)){pass_cuts[5] = 1;}
	   if ((pass_cuts[3] == 1) && (tau_index.size() > 2)){pass_cuts[6] = 1;}
	   
	   // Apply VBF selections
	   if ((abs(Jet_leading_vec.Eta()) > jet_eta_min) && (abs(Jet_subleading_vec.Eta()) > jet_eta_min)){pass_cuts[7] = 1;}	
	   
	   if ((pass_cuts[7] == 1) && (Jet_leading_vec.Pt() > lead_jet_pt) && (Jet_subleading_vec.Pt() > slead_jet_pt)){pass_cuts[8] = 1;} 
	   if ((pass_cuts[4] == 1) && (pass_cuts[8]==1)){pass_cuts[9] = 1;}
	   if ((pass_cuts[5] == 1) && (pass_cuts[8]==1)){pass_cuts[10] = 1;}
	   if ((pass_cuts[6] == 1) && (pass_cuts[8]==1)){pass_cuts[11] = 1;}
           if ((pass_cuts[8] == 1) && (abs(Jet_leading_vec.Eta()) < jet_eta_max) && (abs(Jet_subleading_vec.Eta()) < jet_eta_max)){pass_cuts[12] = 1;}
	   if ((pass_cuts[4] == 1) && (pass_cuts[12] == 1)){pass_cuts[13] = 1;}
	   if ((pass_cuts[5] == 1) && (pass_cuts[12] == 1)){pass_cuts[14] = 1;}
	   if ((pass_cuts[6] == 1) && (pass_cuts[12] == 1)){pass_cuts[15] = 1;}
	   if (pass_cuts[12] == 1){
	     double eta_jj_product = Jet_leading_vec.Eta()*Jet_subleading_vec.Eta();
	     if (eta_jj_product < 0.){pass_cuts[16] = 1;}
	   }
	   if ((pass_cuts[4] == 1) && (pass_cuts[16] == 1)){pass_cuts[17] = 1;}
	   if ((pass_cuts[5] == 1) && (pass_cuts[16] == 1)){pass_cuts[18] = 1;}
	   if ((pass_cuts[6] == 1) && (pass_cuts[16] == 1)){pass_cuts[19] = 1;}
	   
	   if (pass_cuts[16] == 1){
	     double Delta_eta_jj = abs(Jet_leading_vec.Eta() - Jet_subleading_vec.Eta());
	     if (Delta_eta_jj > jj_delta_eta_min){pass_cuts[20] = 1;}
	   }
	   if ((pass_cuts[4] == 1) && (pass_cuts[20] == 1)){pass_cuts[21] = 1;}
	   if ((pass_cuts[5] == 1) && (pass_cuts[20] == 1)){pass_cuts[22] = 1;}
	   if ((pass_cuts[6] == 1) && (pass_cuts[20] == 1)){pass_cuts[23] = 1;}
	   
	   if (pass_cuts[20] == 1){
	     //  TLorentzVector jj_LV_sum = Jet_leading_vec + Jet_subleading_vec;
	     // double mjj = jj_LV_sum.M();
	     double Delta_eta_jj = abs(Jet_leading_vec.Eta() - Jet_subleading_vec.Eta());
	     double mjj = calculateM(Delta_eta_jj, Jet_leading_vec.Pt(), Jet_subleading_vec.Pt());
	     if (mjj > jj_mass_min){pass_cuts[24] = 1;}
	   }
	   if ((pass_cuts[4] == 1) && (pass_cuts[24] == 1)){pass_cuts[25] = 1;}
	   if ((pass_cuts[5] == 1) && (pass_cuts[24] == 1)){pass_cuts[26] = 1;}
	   if ((pass_cuts[6] == 1) && (pass_cuts[24] == 1)){pass_cuts[27] = 1;}

	 //LLENA HISTOGRAMAS DE TAUS, MUONES, ELECTRONES DESDE DETECTOR
	   hmap_numero_taus->Fill(tau_counter);
           hmap_numero_elec->Fill(elec_counter);
           hmap_numero_muones->Fill(muon_counter);
         //LLENTA HISTOGRAMAS DE TAUS, MUONES ELECTRONES DESDE GEN_PARTICLE
	   hmap_numero_taus_gen->Fill(tau_counter_gen);
           hmap_numero_elec_gen->Fill(elec_counter_gen);
           hmap_numero_muones_gen->Fill(muon_counter_gen);	
	 }
       
       
	for (int i = 0; i < nDir; i++){
            _hmap_Nevents[i]->Fill(0.0);
	   if ( pass_cuts[i] == 1){
	   if(Jet_leading_vec.Pt()>2){        // Se exige jets con momento mayor o igual a 2 GeV
	     _hmap_lead_jet_pT[i]->Fill(Jet_leading_vec.Pt());
	     _hmap_lead_jet_eta[i]->Fill(Jet_leading_vec.Eta());
	     _hmap_lead_jet_phi[i]->Fill(Jet_leading_vec.Phi());		   
	     _hmap_Nevents[i]->Fill(1.0);
	   }
	   if(Jet_subleading_vec.Pt()>2){     // Se exige jets con momento mayor o igual a 2 GeV
	     _hmap_sublead_jet_pT[i]->Fill(Jet_subleading_vec.Pt());
	     _hmap_sublead_jet_eta[i]->Fill(Jet_subleading_vec.Eta());
	     _hmap_sublead_jet_phi[i]->Fill(Jet_subleading_vec.Phi());
	   }
	   
	   if(Jet_leading_vec.Pt()>2 && Jet_subleading_vec.Pt()>2){
	     double delta_phi= abs(Jet_leading_vec.Phi()-Jet_subleading_vec.Phi());
             _hmap_delta_phi_jj[i]->Fill(delta_phi);	
       	   }	
           _hmap_MET[i]->Fill(MET);
           if(Jet_leading_vec.Pt()>2 && Jet_subleading_vec.Pt()>2){
	     double delta_eta_jj = abs(Jet_leading_vec.Eta() - Jet_subleading_vec.Eta());
	     _hmap_delta_eta_jj[i]->Fill(delta_eta_jj);
	     double mjj = (Jet_leading_vec + Jet_subleading_vec).M();
	     _hmap_jet_mjj[i]->Fill(mjj);
	     double mjj_formula = calculateM(delta_eta_jj, Jet_leading_vec.Pt(), Jet_subleading_vec.Pt());
	     _hmap_jet_mjj_formula[i]->Fill(mjj_formula);
	   }
	   if(Tau1_vec.Pt() > 2){
	     _hmap_tau1_pT[i]->Fill(Tau1_vec.Pt());
	     _hmap_tau1_eta[i]->Fill(Tau1_vec.Eta());
	     _hmap_tau1_phi[i]->Fill(Tau1_vec.Phi());
  	     double Delta_phi_tau_met=abs(Tau1_vec.Phi()-METpointer->Phi);
             _hmap_delta_phi_tau_met[i]->Fill(Delta_phi_tau_met);	
       	     double razon_Delta_phi_Ptau=Delta_phi_tau_met/Tau1_vec.Pt();
             _hmap_razon_delta_phi_Pt[i]->Fill(razon_Delta_phi_Ptau);
           }	
	   if(Tau2_vec.Pt() > 2){
	     _hmap_tau2_pT[i]->Fill(Tau2_vec.Pt());
	     _hmap_tau2_eta[i]->Fill(Tau2_vec.Eta());
	     _hmap_tau2_phi[i]->Fill(Tau2_vec.Phi());
	   }
	   if(Tau3_vec.Pt() > 2){
	     _hmap_tau3_pT[i]->Fill(Tau3_vec.Pt());
	     _hmap_tau3_eta[i]->Fill(Tau3_vec.Eta());
	     _hmap_tau3_phi[i]->Fill(Tau3_vec.Phi());
	   }
	   
         }
       }  
     }
   theFile->cd();
   
//DELPHES
           hmap_numero_taus->Write();
           hmap_numero_elec->Write();
           hmap_numero_muones->Write();
//GEN_PARTICLE
	   hmap_numero_taus_gen->Write();
           hmap_numero_elec_gen->Write();
           hmap_numero_muones_gen->Write();

	  for (int d = 0; d < nDir; d++)
	     {
       cdDir[d]->cd();
       _hmap_Nevents[d]->Write();
       _hmap_lead_jet_pT[d]->Write();
       _hmap_lead_jet_eta[d]->Write();
       _hmap_lead_jet_phi[d]->Write();
       _hmap_sublead_jet_pT[d]->Write();
       _hmap_sublead_jet_eta[d]->Write();
       _hmap_sublead_jet_phi[d]->Write();
       _hmap_delta_phi_jj[d]->Write();
       _hmap_MET[d]->Write();
       _hmap_delta_eta_jj[d]->Write();
       _hmap_jet_mjj[d]->Write();
       _hmap_jet_mjj_formula[d]->Write();
       _hmap_tau1_pT[d]->Write();
       _hmap_tau1_eta[d]->Write();
       _hmap_tau1_phi[d]->Write();
       _hmap_delta_phi_tau_met[d]->Write();
       _hmap_razon_delta_phi_Pt[d]->Write();
       _hmap_tau2_pT[d]->Write();
       _hmap_tau2_eta[d]->Write();
       _hmap_tau2_phi[d]->Write();
       _hmap_tau3_pT[d]->Write();
       _hmap_tau3_eta[d]->Write();
       _hmap_tau3_phi[d]->Write();
	     }
   theFile->Close();
   
}

PhenoAnalysis::~PhenoAnalysis()
{
  // do anything here that needs to be done at desctruction time
}

double PhenoAnalysis::calculateE(double eta, double pt, double mass){
  
  double theta = TMath::ATan(TMath::Exp(-eta));
  double cos_theta = TMath::Cos(2*theta);
  double p= pt/cos_theta;
  double e = sqrt(pow(p, 2) + pow(mass, 2));
  
  return e;
}

//Calcular masa reconstruida de los dos jets de mas alto momento transverso.

double PhenoAnalysis::calculateM(double deta, double pt1, double pt2){
  double cosh_deta = TMath::CosH(deta);
  double M = sqrt(2*pt1*pt2*cosh_deta);
  return M;
}

//Calcular la masa transversa

double PhenoAnalysis::calculateMt(double phie, double phimet, double met, double pt){
  double Delta =TMath::Cos(phie-phimet);
  double mt = sqrt(2*pt*met*(1-Delta));
  return mt;
}

double PhenoAnalysis::normalizedDphi(double phi){
  const double PI  = 3.141592653589793238463;
  double twoPI = 2.0*PI;
  if ( phi < -PI ){phi += twoPI;}
  if ( phi > PI ){phi -= twoPI;}
  return phi;
}

   
void PhenoAnalysis::crateHistoMasps (int directories)
{
   for (int i = 0; i < directories; i++)
    {
      _hmap_Nevents[i]         = new TH1F("Nevents",          "Nevents", 3, 0., 3);
      _hmap_lead_jet_pT[i]     = new TH1F("jet_lead_pT",     " PT leading Jet", 100, 0., 1000.);
      _hmap_lead_jet_eta[i]    = new TH1F("jet_lead_eta",     "#eta leading Jet", 160, -6.0, 6.0);
      _hmap_lead_jet_phi[i]    = new TH1F("jet_lead_phi",     "#phi leading Jet", 64, -3.6, 3.6); 
      _hmap_sublead_jet_pT[i]  = new TH1F("jet_sublead_pT",   "PT subleading Jet", 100, 0., 1000.);
      _hmap_sublead_jet_eta[i] = new TH1F("jet_sublead_eta",  "#eta subleading Jet", 160, -6.0, 6.0);
      _hmap_sublead_jet_phi[i] = new TH1F("jet_sub_lead_phi", "#phi subleading Jet", 64, -3.6, 3.6);
      _hmap_MET[i]             = new TH1F("MET",   "E_{T}^{miss}", 1000, 0., 10000);
      _hmap_delta_phi_jj[i]     = new TH1F("Delta_phi_jj", "#Delta#phi_{jj}",64,-0., 7.6);
      _hmap_delta_eta_jj[i]     = new TH1F("delta_eta_jj",     "#Delta #eta_{jj}", 100, 0., 10.0);
      _hmap_jet_mjj[i]         = new TH1F("dijet_mjj",        "Masa invariante M_{jj}", 1000, 0., 10000.0);
      _hmap_jet_mjj_formula[i] = new TH1F("dijet_mjj_formula", "Masa invariante M_{jj}", 1000, 0., 10000.0);
      _hmap_tau1_pT[i]         = new TH1F("tau1_pT",          "PT (#tau_{1})", 100, 0., 500.);
      _hmap_tau1_eta[i]        = new TH1F("tau1_eta",         "#eta (#tau_{1})", 150, -6, 6);
      _hmap_tau1_phi[i]        = new TH1F("tau1_phi",         "#phi (#tau_{1})", 64, -3.6, 3.6);
      _hmap_delta_phi_tau_met[i] = new TH1F("Delta_phi_tau_met", "#Delta#phi_{#tau E_{t}^{miss}}",64,0., 7.6);
      _hmap_razon_delta_phi_Pt[i] =new TH1F("Razon_Delta_phi_Pt", "#Delta#phi_{#tau E_{t}^{miss}}/PT_{#tau_{1}}",64,0., 1);
      _hmap_tau2_pT[i]         = new TH1F("tau2_pT",          "PT (#tau_{2})", 100, 0., 500.);
      _hmap_tau2_eta[i]        = new TH1F("tau2_eta",         "#eta (#tau_{2})", 160, -6., 6.);
      _hmap_tau2_phi[i]        = new TH1F("tau2_phi",         "#phi (#tau_{2})", 64, -3.6, 3.6);
      _hmap_tau3_pT[i]         = new TH1F("tau3_pT",          "PT (#tau_{3})", 100, 0., 350.);
      _hmap_tau3_eta[i]        = new TH1F("tau3_eta",         "#eta (#tau_{3})", 160, -3.6, 3.6);
      _hmap_tau3_phi[i]        = new TH1F("tau3_phi",         "#phi (#tau_{3})", 64, -3.6, 3.6);       
 }
}
