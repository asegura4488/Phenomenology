#include <iostream>
#include "ROOTFunctions.h"
#include "DelphesFunctions.h"
#include "mt2w_bisect.h"


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
  
  TH1 *histMT2W_susy1 = new TH1F("MT2W_susy1","mt2w_susy1",100, 100.0, 490.0);
  double pl[4]={0}; // El, plx, ply, plz,     (visible lepton)
  double pb1[4]={0};
  double pb2[4]={0};
  double pmiss[4]={0};

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
      TRootMuon *muon1, *muon2;
      TLorentzVector vec1, vec2;
      TRootElectron *elec1, *elec2;
      TLorentzVector vec3, vec4;
      TRootElectron *electron;
      TRootMuon *muon;
      TRootMissingET *met;
      bool  esolo=false;
      bool musolo=false;
      MissingET *METpointer;
      Muon *Muonpointer;
      Jet *Jetpointer;
         
      
      // Initializing cuadrimomentum
      for (int i=0; i<4; i++) {
	pl[i]=0;
	pb1[i]=0;
	pb2[i]=0;
	if (i<3) pmiss[i]=0;
      }	 

      bool thereisamuon=false;
	  
      if(branchMuon_susy1_Delphes->GetEntries()>0)
      {
	  for(int i=0; i<branchMuon_susy1_Delphes->GetEntries();i++)
	    {

	  Muonpointer = (Muon*) branchMuon_susy1_Delphes->At(i);
	  double theta_Muon = TMath::ATan(TMath::Exp(-Muonpointer->Eta));
	  theta_Muon=theta_Muon*2;
	  double cos_theta_Muon = TMath::Cos(2*theta_Muon);
	  double P_Muon = Muonpointer->PT/cos_theta_Muon;
	  double E_Muon = sqrt(pow(P_Muon, 2) + pow(105.65, 2));
	  pl[0]= E_Muon;
	  pl[1]= Muonpointer->PT*sin(theta_Muon)*cos(Muonpointer->Phi);
	  pl[2]= Muonpointer->PT*sin(theta_Muon)*sin(Muonpointer->Phi);
	  pl[3]= Muonpointer->PT*cos(theta_Muon);
	  thereisamuon=true;
	  }

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
	      pb1[0]= E_jetb1;
	      pb1[1]= Jetpointer->PT*sin(theta_jetb1)*cos(Jetpointer->Phi);
	      pb1[2]= Jetpointer->PT*sin(theta_jetb1)*sin(Jetpointer->Phi);
	      pb1[3]= Jetpointer->PT*cos(theta_jetb1);
	    }
	    if (nbjet==2) {
	      double theta_jetb2 = TMath::ATan(TMath::Exp(-Jetpointer->Eta));
	      double cos_theta_jetb2 = TMath::Cos(2*theta_jetb2);
	      double P_jetb2 = Jetpointer->PT/cos_theta_jetb2;
	      double E_jetb2 = sqrt(pow(P_jetb2, 2) + pow(105.65, 2));
	      pb2[0]= E_jetb2;
	      pb2[1]= Jetpointer->PT*sin(theta_jetb2)*cos(Jetpointer->Phi);
	      pb2[2]= Jetpointer->PT*sin(theta_jetb2)*sin(Jetpointer->Phi);
	      pb2[3]= Jetpointer->PT*cos(theta_jetb2);
	      
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

	  if (wehave_etmiss&&thereisamuon&&nbjet>1){
	      mt2w_bisect::mt2w mt2w_event;
	      
	      mt2w_event.set_momenta(pl,pb1,pb2,pmiss);
	      double value_mt2w=mt2w_event.get_mt2w();
	      histMT2W_susy1->Fill(value_mt2w);
	      /*      
	      cout << " Plepton ="<<pl[0]<<" , "<<pl[1]<<" , "<<pl[2]<<" , "<<pl[3]<<endl;
	      cout << " P_bjet1 ="<<pb1[0]<<" , "<<pb1[1]<<" , "<<pb1[2]<<" , "<<pb1[3]<<endl;
	      cout << " P_bjet2 ="<<pb2[0]<<" , "<<pb2[1]<<" , "<<pb2[2]<<" , "<<pb2[3]<<endl;
	      cout << " P- Miss Et =  0 , "<<pmiss[1]<<" , "<<pmiss[2]<<endl;


	      
	      cout << " mt2w = " << value_mt2w << endl;
	      */
	  }

    }

	  
  
	  
	  TFile* hfile = new TFile("mt2w_prueba.root", "RECREATE");
	  
	  histMT2W_susy1->Write();
	  
	  
}
