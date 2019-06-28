#include <iostream>
#include "ROOTFunctions.h"
#include "DelphesFunctions.h"
using namespace std;

// Started
int main(){

 // Root ntuplas semileptonic background (prueba)

  TChain chain_Pythia("STDHEP");
  chain_Pythia.Add("/Disco1/Pheno/BackgroudSamples/ttbar_semileptonico/ttbar_semileptonico_prueba_001/Events/run_01/output_pythia8.root");

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader_Pythia = new ExRootTreeReader(&chain_Pythia);

 // Number of entries access
  long int numberOfEntries = treeReader_Pythia->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchGenParticle = treeReader_Pythia->UseBranch("GenParticle");

  // Book histograms

 TH1 *histNumElectrons_pythia = new TH1F("num_electrons_pythia","numelectrons_pythia",21,0.0,20.0);
 TH1 *histNumMuons_pythia = new TH1F("num_muons_pythia","nummuons_pythia",21,0.0,20.0);
 TH1 *histNumLeptons_pythia = new TH1F("num_leptons_pythia","num_Leptons_pythia",21,0.0,20.0);

 int imuon_pythia=0;
 int ielectron_pythia=0;
 int ielectronw=0;
 int counter_mother=0;
 int imuon_electron_pythia=0;

 // Loop over all events
  cout<<"numberofEntries="<<numberOfEntries<<endl;

  for ( Int_t entry = 0; entry < numberOfEntries; ++entry ){

     imuon_pythia=0;
     ielectron_pythia=0;
     ielectronw=0;
     counter_mother=0;
     imuon_electron_pythia=0;

    // Load selected branches with data from specified event
    treeReader_Pythia->ReadEntry(entry);
   //Definition of particles
    GenParticle *Particlepointer, *Motherpointer, *Daughterpointer;  

    //Topology Analysis in Pythia 8
    int motherPID, daughterPID;
    //cout<<"entries in branchGenParticle="<<branchGenParticle->GetEntries()<<endl;
    for (Int_t k=0; k<branchGenParticle->GetEntries(); k++){

      Particlepointer = (GenParticle*)branchGenParticle->At(k);
      cout<<"particle ="<<k<<" ID ="<<Particlepointer->PID<<": Status ="<<Particlepointer->Status<<"; Mother ID1="<<Particlepointer->M1<<"; Mother ID2="<<Particlepointer->M2<<"; daughter ID1="<<Particlepointer->D1<<"; daughter ID2="<<Particlepointer->D2<<endl;
 
      //cout<<"mother ="<<Particlepointer->M1<<endl;
      //Motherpointer   = (GenParticle*)branchGenParticle->At(Particlepointer->M1);
      Motherpointer   = (GenParticle*)branchGenParticle->At(k);
      Daughterpointer = (GenParticle*)branchGenParticle->At(k);
      
      motherPID = Motherpointer->PID;
      daughterPID = Daughterpointer->PID;

      if((Particlepointer->Status == 23) || (Particlepointer->Status == 24) ){ //Particle of the hardest process (outgoing)

	if((Particlepointer->PID == 13) ||(Particlepointer->PID == -13)){ //Checking muons
	  imuon_pythia++;
	  imuon_electron_pythia++;

	  counter_mother = Motherpointer->M1;
	  //cout<<"muon: Status ="<<Particlepointer->Status<<"; Mother ID="<<Particlepointer->M1<<"; daughter="<<Particlepointer->D1<<endl;
   // cout<<"Muons_mother"<<counter_mother<<endl;
	}

	if((Particlepointer->PID == 11)||(Particlepointer->PID == -11)){ //Checking electrons
	  ielectron_pythia++;
	  imuon_electron_pythia++;
	  counter_mother = Motherpointer->M1;
	  //  cout<<"Electrons_mother"<<counter_mother<<endl;

	}

	
	 }//hardest process
      /*
	 if((Particlepointer->PID == 13) ||(Particlepointer->PID == -13)){ 
	   cout<<"muon: Status ="<<Particlepointer->Status<<"; Mother= "<<Particlepointer->M1<<" ; daughter = "<<Particlepointer->D1<<endl;
	 }
	 if((Particlepointer->PID == 24) ||(Particlepointer->PID == -24)){ 
	   cout<<"Wboson: Status ="<<Particlepointer->Status<<"; Mother= "<<Particlepointer->M1<<" ; daughter = "<<Particlepointer->D1<<endl;
	 }
      */
    } //end branch particle
 
    histNumMuons_pythia->Fill(imuon_pythia++);
    histNumElectrons_pythia->Fill(ielectron_pythia++);
    histNumLeptons_pythia->Fill(imuon_electron_pythia++);

  } //end cicle of the program

  cout<<"events with muons pythia="<<imuon_pythia<<endl;
  cout<<"events with electrons pythia="<<ielectron_pythia<<endl;

 

  // ROOT Program output

  TFile* hfile = new TFile("prueba_semileptonic_pythia_leading_jet.root", "RECREATE");
  histNumElectrons_pythia->Write();
  histNumMuons_pythia->Write();
  histNumLeptons_pythia->Write();
  hfile->Close();


}//end program
