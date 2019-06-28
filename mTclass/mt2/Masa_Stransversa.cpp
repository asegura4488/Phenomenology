#include <iostream>
#include "ROOTFunctions.h"
#include "DelphesFunctions.h"
#include "mt2_bisect.h"


using namespace std;

int main(){


  // Load shared library
  
  gSystem->Load("lib/libExRootAnalysis.so");
  gSystem->Load("libPhysics");
  
  
  // Create chains of root trees
  
  //Susy1: Decaimiento de pares de stops en pares top-antitop
  
  TChain chainsusy1_MadGraph("LHEF");
  chainsusy1_MadGraph.Add("/home/rrodriguez/Simulación/MG_pythia8_delphes_parallel/Runs/Susy_stops1_13Tev_con_1_ISR_001/Events/run_01/unweighted_events.root");
  
  TChain chainsusy1_Pythia("STDHEP");
  chainsusy1_Pythia.Add("/home/rrodriguez/Simulación/MG_pythia8_delphes_parallel/Runs/Susy_stops1_13Tev_con_1_ISR_001/Events/run_01/output_pythia8.root");
  
  TChain chainsusy1_Delphes("Delphes");
  chainsusy1_Delphes.Add("/home/rrodriguez/Simulación/MG_pythia8_delphes_parallel/Runs/Susy_stops1_13Tev_con_1_ISR_001/Events/run_01/output_delphes.root");

 // Create object of class ExRootTreeReader

  ExRootTreeReader *treeReader_susy1_MadGraph = new ExRootTreeReader(&chainsusy1_MadGraph);
  ExRootTreeReader *treeReader_susy1_Pythia = new ExRootTreeReader(&chainsusy1_Pythia);
  ExRootTreeReader *treeReader_susy1_Delphes = new ExRootTreeReader(&chainsusy1_Delphes);

  Long64_t numberOfEntries_susy1 = treeReader_susy1_MadGraph->GetEntries();

  // Get pointers to branches used in this analysis
    
  TClonesArray *branchJet_susy1_Delphes = treeReader_susy1_Delphes->UseBranch("Jet");
  TClonesArray *branchElectron_susy1_Delphes = treeReader_susy1_Delphes->UseBranch("Electron");
  TClonesArray *branchMuon_susy1_Delphes = treeReader_susy1_Delphes->UseBranch("Muon");  
  TClonesArray *branchPhoton_susy1_Delphes = treeReader_susy1_Delphes->UseBranch("Photon"); 
  TClonesArray *branchMissingET_susy1_Delphes = treeReader_susy1_Delphes->UseBranch("MissingET");


  // Book histograms
  
  TH1 *histMT2_susy1 = new TH1F("MT2_susy1","mt2_susy1",100, 100.0, 490.0);
  TH1 *histNumJets_susy1  = new TH1F("num_jets","jets",100, 100.0, 10.0);
  double pa[3]={0};
  double pb[3]={0};
  double pmiss[3]={0};
  double mn    = 50.;


  // Loop over all events (susy)


  int nbjet;

  for ( Int_t entry = 0; entry < numberOfEntries_susy1; ++entry )

    {

      // Load selected branches with data from specified event


      treeReader_susy1_MadGraph->ReadEntry(entry);
      treeReader_susy1_Pythia->ReadEntry(entry);
      treeReader_susy1_Delphes->ReadEntry(entry);

      // TLorentz vector to calculate the di-jet invariant mass and tri-jet invariant mass
      
      TLorentzVector jetvec1, jetvec2, jetvec3;
      TLorentzVector jetvecb1, jetvecb2, jetvecb3;
      TLorentzVector vec1, vec2;
      TLorentzVector vec3, vec4;
      bool  esolo=false;
      bool musolo=false;
      MissingET *METpointer;
      Muon *Muonpointer;
      Jet *Jetpointer;
            
      int NumJets = branchJet_susy1_Delphes->GetEntries();
      histNumJets_susy1->Fill(NumJets);
      
      // Initializing cuadrimomentum
      for (int i=0; i<3; i++) {
	pa[i]=0;
	pb[i]=0;
	if (i<3) pmiss[i]=0;
      }	 


      
      if(branchJet_susy1_Delphes->GetEntries()>0 )
      {
	nbjet=0;
	  for(int i=0; i<branchJet_susy1_Delphes->GetEntries();i++)
	    {

	  Jetpointer = (Jet*) branchJet_susy1_Delphes->At(i);
	  if (Jetpointer->BTag>0) {
	    nbjet++;
	    if (nbjet==1){
	      double theta_jetb1 = TMath::ATan(TMath::Exp(-Jetpointer->Eta));
	      double cos_theta_jetb1 = TMath::Cos(2*theta_jetb1);
	      double P_jetb1 = Jetpointer->PT/cos_theta_jetb1;
	      double E_jetb1 = sqrt(pow(P_jetb1, 2) + pow(105.65, 2));
	      pa[0]= 0.106;
	      pa[1]= Jetpointer->PT*sin(theta_jetb1)*cos(Jetpointer->Phi);
	      pa[2]= Jetpointer->PT*sin(theta_jetb1)*sin(Jetpointer->Phi);
	    }
	    if (nbjet==2) {
	      double theta_jetb2 = TMath::ATan(TMath::Exp(-Jetpointer->Eta));
	      double cos_theta_jetb2 = TMath::Cos(2*theta_jetb2);
	      double P_jetb2 = Jetpointer->PT/cos_theta_jetb2;
	      double E_jetb2 = sqrt(pow(P_jetb2, 2) + pow(105.65, 2));
	      pb[0]= 0.106;
	      pb[1]= Jetpointer->PT*sin(theta_jetb2)*cos(Jetpointer->Phi);
	      pb[2]= Jetpointer->PT*sin(theta_jetb2)*sin(Jetpointer->Phi);
	      
	    }
	  }
	    }
      
      }
      
      bool wehave_etmiss=false;
    
	  if (branchMissingET_susy1_Delphes->GetEntries()>0)
	 {
	      for(Int_t i=0; i<branchMissingET_susy1_Delphes->GetEntries();i++)
		{
	      METpointer = (MissingET*) branchMissingET_susy1_Delphes->At(i);
	      double missET_eta = METpointer->Eta;
	      double theta_missET = 2*atan(exp(-missET_eta));
	      pmiss[0]=0;
	      pmiss[1]= METpointer->MET*sin(theta_missET)*cos(METpointer->Phi);
	      pmiss[2]= METpointer->MET*sin(theta_missET)*sin(METpointer->Phi);
	      wehave_etmiss=true;
		}	 
	      } 

	  if (wehave_etmiss&&nbjet>1){
	      mt2_bisect::mt2 mt2_event;
	      mt2_event.set_momenta(pa,pb,pmiss);
	      mt2_event.set_mn(mn);
	      double value_mt2=mt2_event.get_mt2();
	      histMT2_susy1->Fill(value_mt2);
	            
	      cout << " P_bjet1 ="<<pa[0]<<" , "<<pa[1]<<" , "<<pa[2]<<endl;
	      cout << " P_bjet2 ="<<pb[0]<<" , "<<pb[1]<<" , "<<pb[2]<<endl;
	      cout << " P- Miss Et =  0 , "<<pmiss[1]<<" , "<<pmiss[2]<<endl;


	      
	      cout << " mt2 = " << value_mt2 << endl;
	      
	  }

    }

	  
  
	  
	  TFile* hfile = new TFile("mt2_prueba.root", "RECREATE");
	  
	  histMT2_susy1->Write();
	  histNumJets_susy1->Write();
	  
}
