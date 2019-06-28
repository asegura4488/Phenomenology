/*
Simple macro showing how to access branches from the delphes output root file,
loop over events, and plot simple quantities such as the jet pt and the di-electron invariant
mass.

root -l examples/Example1.C'("delphes_output.root")'
*/


//#ifdef __CLING__
//R__LOAD_LIBRARY(libDelphes);
#include "DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
//#endif

//------------------------------------------------------------------------------

void TausExample(const char *inputFile)
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  //TClonesArray *branchElectron = treeReader->UseBranch("Electron");

  // Book histograms
  TH1 *histJetPT = new TH1F("jet_pt", "jet P_{T}", 100, 0.0, 100.0);
  TH1 *histNtaus = new TH1F("Ntaus", "N_{taus}", 5, 0, 5);

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // If event contains at least 1 jet
    if(branchJet->GetEntries() > 0)
    {
      int Jet_size = branchJet->GetEntries();

      //Counting tau tagged jets per event
      int Ntaus = 0;

      // Take first jet
      Jet *jet = (Jet*) branchJet->At(0);

      // Plot jet transverse momentum
      histJetPT->Fill(jet->PT);

      // Print jet transverse momentum
//      cout << "Jet pt: "<<jet->PT << endl;

      for(Int_t jeti = 0; jeti < Jet_size; ++jeti)
	{
	  Jet *jet1 = (Jet*) branchJet->At(jeti);
	  if (jet1->TauTag==1) Ntaus++;
	}
      histNtaus->Fill(Ntaus);

    }

    //Electron *elec1, *elec2;

    // If event contains at least 2 electrons
    /*if(branchElectron->GetEntries() > 1)
    {
      // Take first two electrons
      elec1 = (Electron *) branchElectron->At(0);
      elec2 = (Electron *) branchElectron->At(1);

      // Plot their invariant mass
      histMass->Fill(((elec1->P4()) + (elec2->P4())).M());
      }*/
  }

  // Show resulting histograms
  histJetPT->Draw();
 // histNtaus->Draw();
}

