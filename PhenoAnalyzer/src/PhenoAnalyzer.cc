/*
  @file PhenoAnalyzer.cc
  @author Andres Florez
  
  Code used to perform phenomenological analysis of Heavy Neutrinos in the tau channel
*/

#include "PhenoAnalyzer.h"
#include "ElectronA.h"

int main(int argc, char *argv[]) {
  
  //TApplication app("App",&argc, argv);
  // gSystem->Load("libDelphes.so");
  TChain chain("Delphes");
  chain.Add(argv[1]);
  TFile * HistoOutputFile = new TFile(argv[2], "RECREATE");
  int nDir = 7;
  TDirectory *theDirectory[nDir];
  theDirectory[0]  = HistoOutputFile->mkdir("no_cuts");//We make a directory for each historams.
  theDirectory[1]  = HistoOutputFile->mkdir("jet_kinematics");
  theDirectory[2]  = HistoOutputFile->mkdir("met_selection");
  theDirectory[3]  = HistoOutputFile->mkdir("muon_cuts");
  theDirectory[4]  = HistoOutputFile->mkdir("electron_cuts");
  theDirectory[5]  = HistoOutputFile->mkdir("tau_cuts");
  theDirectory[6]  = HistoOutputFile->mkdir("bjet_cuts");
  
  
  PhenoAnalysis BSM_analysis(chain, HistoOutputFile, theDirectory, nDir);
  
}

using namespace std;
PhenoAnalysis::PhenoAnalysis(TChain& chain, TFile* theFile, TDirectory *cdDir[], int nDir)
{
  ifstream inFile;
  inFile.open ("config.in", ios::in);//open the file
  
  if (!inFile)
    {
      cerr << "ERROR: Can't open input file: " << endl;
      exit (1);//check if we opened the file
    }
  
  string inputType = "";
  
  //This set of lines are used to open and read the "config.in" file.
  ///////////////////////////////////////////////////////////////////////
  TEnv *params = new TEnv ("config_file");
  params->ReadFile ("config.in", kEnvChange);
  
  double b_jet_pt_min              = params->GetValue ("b_jet_pt_min",      30.0);
  double lead_jet_ISR_pt           = params->GetValue ("lead_jet_ISR_pt",   100.0);
  double lead_jet_eta 		   = params->GetValue ("lead_jet_eta", 2.4);
  double lead_jet_pt_min           = params->GetValue ("lead_jet_pt_min", 30.);
  double muon_pt_min               = params->GetValue ("muon_pt_min",      20.0);
  double muon_pt_max               = params->GetValue ("muon_pt_max",   40.0);
  double muon_eta_min              = params->GetValue ("muon_eta_min",      2.1);
  double lepton_pt_min             = params->GetValue ("lepton_pt_min", 10.);
  double lepton_pt_max             = params->GetValue ("lepton_pt_max",    40.);
  double lepton_eta_max            = params->GetValue ("lepton_eta_max",    2.4);
  double et_miss                   = params->GetValue ("et_miss", 200);
  double b_jet_eta                 = params->GetValue ("b_jet_eta", 0.);
  double delta_R_jet_electron      = params->GetValue ("delta_R_jet_electron", 0.4);
  double delta_R_jet_muon          = params->GetValue ("delta_R_jet_muon", 0.4);
  double delta_R_jet_tau           = params->GetValue ("delta_R_jet_tau", 0.4);
  crateHistoMasps(nDir);
  
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
  
  MissingET *METpointer;//Energia transversal perdida
  auto el_pt_min = 10000.;
  int n_elec_min = 0;

  for(Int_t entry = 0; entry < numberOfEntries; ++entry){
    
    
    treeReader->ReadEntry(entry);
    int pass_cuts[nDir];
    TLorentzVector jetLeadingVec (0., 0., 0., 0.); //tc stands for tau channel
    TLorentzVector jetLeadingVec_cuts (0.,0.,0.,0.);
    TLorentzVector muonVec (0.,0.,0.,0.);
    TLorentzVector muonVec_cuts (0.,0.,0.,0.);
    TLorentzVector electronVec (0.,0.,0.,0.);
    TLorentzVector electronVec_min (0.,0.,0.,0.);  
    TLorentzVector tauVec (0.,0.,0.,0.); 
    TLorentzVector bjetVec(0.,0.,0.,0.);   
    double jet_pt_max = 0.0; 
    double met = 0.0;
    double met_phi = 0.0;
    double met_mass = 0.0;
    double met_pt = 0.0;
    double norm_hist = 1.;
    //Search for taus and bjets
    int n_taus = 0;
    int n_bjets = 0;
    int n_muons = 0; 
    int n_electrons = 0;
   
    ElectronA *Electron_ = new ElectronA(branchElectron, 0.00051099802);
    TLorentzVector PTmin = Electron_->GetMinVectorPt();
    
   if(PTmin.Pt() != 0.){
    TLorentzVector Final=Cuts(Electron_);
    std::cout << Final.Pt() << endl;

    std::cout << PTmin.Pt() << std::endl;
 
   }
    for (int j = 0; j < branchJet->GetEntriesFast(); j++){
      Jet *jet = (Jet*) branchJet->At(j);
      
      double jet_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
      // choose leading jet
      if ((jet->PT > lead_jet_pt_min) && (abs(jet->Eta) < lead_jet_eta)){
        if((jet->TauTag!=1) && (jet->BTag!=1) && (jet->PT > jet_pt_max)){
          jet_pt_max = jet->PT;
          jetLeadingVec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);	  
        }
        if (jet->TauTag == 1){
          n_taus++;
          tauVec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy); 
        }
        if (jet->BTag == 1){ 
          bjetVec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
          n_bjets++;         
        }
      } 
    }
    //Check if taus overlap with muons
    //we check the number of muons
    for (int muo = 0; muo < branchMuon->GetEntriesFast(); muo++){
      Muon *muon = (Muon*) branchMuon->At(muo);
      
      double muon_energy = calculateE(muon->Eta, muon->PT, 0.1056583745);
      if((muon->PT > muon_pt_min) && (muon->PT < muon_pt_max) && (abs(muon->Eta) > muon_eta_min)){
        n_muons++; 
	muonVec_cuts.SetPtEtaPhiE(muon->PT,muon->Eta, muon->Phi, muon_energy);
      }
    } // for (int muo = 0; muo < branchMuon->GetEntriesFast(); muo++)
    
    METpointer = (MissingET*) branchMissingET->At(0);//we fill the MET
    met = METpointer->MET;
    met_phi = METpointer->Phi;
    //    met_mass = METpointer->Mass; 
    //    met_pt = METpointer->Pt;
   // int n_elec_min = 0;
    int elec_pt_min =0;
//    auto el_pt_min = 10000.;
    vector<double> a;
    //loop for electrons
   // cout << entry-> << endl;
    for (int el = 0; el < branchElectron->GetEntriesFast(); el++){
      Electron* electron = (Electron*) branchElectron->At(el);   
      double	electron_energy = calculateE(electron->PT, electron->Phi, 0.00051099802);
      if((electron->PT > lepton_pt_min) && (electron->PT < lepton_pt_max) && (abs(electron->Eta) < lepton_eta_max) ){
	electronVec.SetPtEtaPhiE(electron->PT,electron->Eta, electron->Phi, electron_energy);
     
     
	if(electronVec.Pt() < el_pt_min){
	 el_pt_min = electronVec.Pt();
         n_elec_min++;
        
  }
}
// method for finding the least pt lepton
// yet to find a more general method
//     if((lepton_pt_min < electronVec.Pt()) && (electronVec.Pt() <= lepton_pt_min + 0.2) ){
//   electronVec_min.SetPtEtaPhiE(electron->PT,electron->Eta, electron->Phi, electron_energy); 
// cout << electronVec.Pt() << endl;  
    n_electrons++;
    }
   // cout << el_pt_min << endl;
      /*
       double first = a[0];
       for(int p = 0; p < size ; p++){
        if(a[p] < first){
           first = a[p];
         }
        // cout << first << endl;
       }
  */   
    
    //Check for overlapping of leptons and jets
    electronVec.DeltaR(jetLeadingVec) > delta_R_jet_electron; 
    muonVec.DeltaR(jetLeadingVec) >  delta_R_jet_muon;
    tauVec.DeltaR(jetLeadingVec) >  delta_R_jet_tau;
    
    double tr_mass = sqrt((electronVec.M()*electronVec.M()) +  + 2*(electronVec.Et()*met) );
    ///cout << tr_mass << endl;
    //////// Apply cuts /////////
    // Events with no cuts
    pass_cuts[0] = 1;
    
    // Events with cuts
    if (jetLeadingVec.Pt() > lead_jet_ISR_pt){
      pass_cuts[1] = 1;
    }
    if ((pass_cuts[1] == 1) && (met > et_miss )){
      pass_cuts[2] = 1;
    }
    if ((pass_cuts[2] == 1) && (n_muons== 1) && (n_taus == 0) && (n_electrons == 0) ) {
      pass_cuts[3] = 1;
    } 
    if ((pass_cuts[2] == 1) && (n_electrons==1) && (n_taus == 0) && (n_muons == 0)) {
      pass_cuts[4] = 1;
    }
    if ((pass_cuts[2] == 1) && (n_taus==1) && (n_muons == 0) && (n_electrons == 0)) {
      pass_cuts[5] = 1;
    }
    if((pass_cuts[2]== 1) && (bjetVec.Eta() == b_jet_eta)){
      pass_cuts[6] = 1;
    }  
    //Fill histograms (MM: No sé como llenar los histogramas de masa...)
    for (Int_t i = 0; i < nDir; i++){
      if (pass_cuts[i] == 1){
	_hmap_Nevents[i]->Fill(0.0);
//Binning histograms
    for (int bin = 1 ; bin < nDir; bin++){
        int  bin_no_cuts = _hmap_Nevents[0]->GetBinContent(1);
        _hmap_Nev_bins[bin]->SetBinContent(1,bin_no_cuts);
        int  bin_cuts  = _hmap_Nevents[bin]->GetBinContent(1);
        _hmap_Nev_bins[bin]->SetBinContent(2,bin_cuts);
    }
        _hmap_jet_pt[i]->Fill(jetLeadingVec.Pt());
        _hmap_jet_pt_cuts[i]->Fill(jetLeadingVec_cuts.Pt());
        _hmap_jet_eta_cuts[i]->Fill(jetLeadingVec_cuts.Eta());
        _hmap_jet_phi_cuts[i]->Fill(jetLeadingVec_cuts.Phi());
        _hmap_MET[i]->Fill(met);
        _hmap_muon_pt[i]->Fill(muonVec.Pt());
        _hmap_muon_pt_cuts[i]->Fill(muonVec_cuts.Pt());
        _hmap_muon_eta_cuts[i]->Fill(muonVec_cuts.Eta());
        _hmap_muon_phi_cuts[i]->Fill(muonVec_cuts.Phi());
	_hmap_electron_pt[i]->Fill(electronVec.Pt());
	_hmap_electron_eta[i]->Fill(electronVec.Eta());
	_hmap_electron_phi[i]->Fill(electronVec.Phi()); 
	_hmap_electron_met_phi[i]->Fill(abs(electronVec.Phi() - met_phi));
	_hmap_muon_met_phi[i]->Fill(abs(muonVec.Phi() - met_phi)) ;
	_hmap_tau_met_phi[i]->Fill(abs(tauVec.Phi() - met_phi));
       	_hmap_electron_jet_phi[i]->Fill(abs(electronVec.Phi() - jetLeadingVec.Phi()));
	_hmap_muon_jet_phi[i]->Fill(abs(muonVec.Phi() -  jetLeadingVec.Phi())) ;
	_hmap_tau_jet_phi[i]->Fill(abs(tauVec.Phi() -  jetLeadingVec.Phi()));
        _hmap_jet_met_phi[i]->Fill(abs(jetLeadingVec.Phi() - met_phi));
        _hmap_pt_ratio_electron[i]->Fill(electronVec.Pt()/(met + jetLeadingVec.Phi()));
        _hmap_pt_ratio_muon[i]->Fill(muonVec.Pt()/(met + jetLeadingVec.Phi()));
        _hmap_pt_ratio_tau[i]->Fill(tauVec.Pt()/(met + jetLeadingVec.Phi()));
//Transverse mass, very small values.        
        _hmap_mt[i]->Fill(tr_mass);
//2D of jetPT vs MET
        _hmap_ptlead_vs_met[i]->Fill(jetLeadingVec.Pt(),met);	
//Scaling histogram.
        _hmap_Nevents[i]->Scale(norm_hist);
        _hmap_jet_pt[i]->Scale(norm_hist); 	
      
      }
    }
  }
//---------------------- 

//  cout << el_pt_min << endl;
//  cout << n_elec_min << endl;
// end entry loop for tau channel
  theFile->cd();
  for (int d = 0; d < nDir; d++)
    {
      cdDir[d]->cd();
      _hmap_Nevents[d]->Write();
      _hmap_Nev_bins[d]->Write();
      _hmap_jet_pt[d]->Write();
      _hmap_jet_pt_cuts[d]->Write();
      _hmap_jet_eta_cuts[d]->Write();
      _hmap_jet_phi_cuts[d]->Write();
      _hmap_MET[d]->Write();
      _hmap_muon_pt[d]->Write();
      _hmap_muon_pt_cuts[d]->Write();
      _hmap_muon_eta_cuts[d]->Write();
      _hmap_muon_phi_cuts[d]->Write();
      _hmap_electron_pt[d]->Write();
      _hmap_electron_eta[d]->Write();
      _hmap_electron_phi[d]->Write();
      _hmap_electron_met_phi[d]->Write();
      _hmap_muon_met_phi[d]->Write() ;
      _hmap_tau_met_phi[d]->Write();
      _hmap_electron_jet_phi[d]->Write();
      _hmap_muon_jet_phi[d]->Write() ;
      _hmap_tau_jet_phi[d]->Write();	
      _hmap_jet_met_phi[d]->Write();
      _hmap_pt_ratio_electron[d]->Write();
      _hmap_pt_ratio_muon[d]->Write();
      _hmap_pt_ratio_tau[d]->Write();
      _hmap_mt[d]->Write();
      _hmap_ptlead_vs_met[d]->Write();
  
         }
  theFile->Close();
}

PhenoAnalysis::~PhenoAnalysis()
{
  // do anything here that needs to be done at desctruction time
}

TLorentzVector PhenoAnalysis::Cuts(ElectronA *Electron_){

  TLorentzVector tmp;
  bool passCut = true;
  
  passCut = passCut && (Electron_->p4(0).Pt() > 10. && Electron_->p4(0).Pt() < 40.);
  passCut = passCut && (Electron_->p4(0).Eta() > 0. && Electron_->p4(0).Eta() < 2.4);
  
  if(passCut) tmp = Electron_->p4(0);

  return tmp;
}

double PhenoAnalysis::calculateE(double eta, double pt, double mass){
  
  double theta = TMath::ATan(TMath::Exp(-eta));
  double sin_theta = TMath::Sin(2*theta);
  double p= pt/sin_theta;
  double e = sqrt(pow(p, 2) + pow(mass, 2));
  
  return e;
}
double PhenoAnalysis::calculate_deltaR(TLorentzVector vector,  Track* track){
  
  double eta1 = vector.Eta();
  double phi1 = vector.Phi();
  double eta2 = track->Eta;
  double phi2 = track->Phi;
  double deltaR = sqrt(pow(eta1-eta2,2) + pow(phi1-phi2,2));
  return deltaR;
}

double PhenoAnalysis::normalizedDphi(double phi){
  const double PI  = 3.141592653589793238463;
  double twoPI = 2.0*PI;
  if ( phi < -PI ){phi += twoPI;}
  if ( phi > PI ){phi = twoPI-phi;}
  else phi = TMath::Abs(phi);
  return phi;
}
void PhenoAnalysis::crateHistoMasps (int directories)
{
 for (Int_t i = 0; i < directories; i++)
    {
      _hmap_Nevents[i]       = new TH1F("Nevents", "Nevents", 3,0,3);
      _hmap_Nev_bins[i]      = new TH1F("Nevents_cuts" , "Nevents_cuts" , 3 , 0 , 3 );
      _hmap_jet_pt[i]        = new TH1F("jet_pt_max", "jet pT max", 500, 0., 500.);
      _hmap_jet_pt_cuts[i]   = new TH1F("jet_pt_max_cuts", "jet pT max cuts", 500, 100., 500.);
      _hmap_jet_eta_cuts[i]  = new TH1F("jet_pt_max_cuts_eta", "jet pT max cuts_eta", 500, 0.1, 5.);
      _hmap_jet_phi_cuts[i]  = new TH1F("jet_pt_max_cuts_phi", "jet pT max cuts_phi", 500,0., 180.);
      _hmap_muon_pt[i]       = new TH1F("muon_pt", "muon pT", 500,0.1,500.);
      _hmap_muon_pt_cuts[i]  = new TH1F("muon_pt_cuts", "muon pT w/ cuts", 500,20,40);
      _hmap_muon_eta_cuts[i] = new TH1F("muon_eta_cuts"," muon eta w/ cuts", 500,0.,5);
      _hmap_muon_phi_cuts[i] = new TH1F("muon_phi_cuts"," muon phi w/ cuts",1000,0.001,180);
      _hmap_MET[i]           = new TH1F("MET", "Missing Transverse Energy", 500, 0., 500.);
      _hmap_electron_pt[i]   = new TH1F("electron_pt","electron Pt", 500, 0.,100.);
      _hmap_electron_eta[i]  = new TH1F("electron_eta","electron Eta", 500, 0.1,10.);
      _hmap_electron_phi[i]  = new TH1F("electron_phi","electron Phi", 500, 0.01,5.);
      _hmap_electron_met_phi[i]   = new TH1F("electron_met_phi","electron_met_phi", 500, 0.,200.);
      _hmap_muon_met_phi[i]       = new TH1F("muon_met_phi","muon_met_phi", 500, 0.,200.);
      _hmap_tau_met_phi[i]        = new TH1F("tau_met_phi","tau_met_phi", 500, 0.,200.);
      _hmap_electron_jet_phi[i]   = new TH1F("electron_jet_phi","electron_jet_phi", 500, 0.,200.);
      _hmap_muon_jet_phi[i]       = new TH1F("muon_jet_phi","muon_jet_phi", 500, 0.,200.);
      _hmap_tau_jet_phi[i]        = new TH1F("tau_jet_phi","tau_jet_phi", 500, 0.,200.);
      _hmap_jet_met_phi[i]        = new TH1F("jet_met_phi","jet_met_phi", 500, 0.,200.); 
      _hmap_pt_ratio_electron[i]  = new TH1F("pt_ratio_electron", "pt ratio electron", 500, 0., 500.);
      _hmap_pt_ratio_muon[i]      = new TH1F("pt_ratio_muon", "pt ratio muon", 500, 0., 500.);
      _hmap_pt_ratio_tau[i]       = new TH1F("pt_ratio_tau", "pt ratio tau", 500, 0., 500.);
      _hmap_mt[i]                 = new TH1F("transverse_mass", "trasverse mass", 500, 0.000001, 5.);
      _hmap_ptlead_vs_met[i]      = new TH2F("Pt_lead_jet_vs_MET", "Pt_lead_jet vs MET", 500, 0., 500., 500, 0.,500.);  
      
  }
}
