#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TROOT.h>
#include <TF1.h>
#include <TObject.h>
#include <TStorage.h>
#include <TRint.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include "TClonesArray.h"


#include <classes/DelphesClasses_1.h>
#include "ExRootTreeReader.h"

using namespace std;

int main(){

  TChain chain("Delphes");
  chain.Add("/Disco2/Pheno/Backgrounds/W+jets/W+jets_1/Events/run_01/output_delphes.root");
  TFile * thefile = new TFile("salida.root", "RECREATE");



  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");

  TH1 *histJetPT = new TH1F("jet_pt", "jet P_{T}", 100, 0.0, 500.0);
  TH1 *histNtaus = new TH1F("Ntaus", "N_{taus}", 5, 0, 5);

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    treeReader->ReadEntry(entry);

    if(branchJet->GetEntries() > 0)
    {
      int Jet_size = branchJet->GetEntries();
      int Ntaus = 0;

      Jet *jet = (Jet*) branchJet->At(0);

      // Plot jet transverse momentum
  if(jet->PT > 70.0){
      histJetPT->Fill(jet->PT);
	}
      for(Int_t jeti = 0; jeti < Jet_size; ++jeti)
	{
	  Jet *jet1 = (Jet*) branchJet->At(jeti);
	  if (jet1->TauTag==1) Ntaus++;
	}
      histNtaus->Fill(Ntaus);
    }

  }

   thefile->cd();
   histJetPT->Write();
   histNtaus->Write();  
   thefile->Close();
  
}
