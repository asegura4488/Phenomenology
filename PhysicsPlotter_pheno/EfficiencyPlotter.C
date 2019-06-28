{
  TFile* f2 = (TFile*) TFile::Open ("Normalized_DYLL/AfterTau_IsoLoosDeltaBetaCorr3Hits/normalizedHistos_DYToLL_0j.root");
  //f1->cd("AfterMuonChargeProduct");
  TH1F* den_pt = muon_pT_lumi;
  TH1F* den_eta = muon_eta_lumi;
  TH1F* den_m   = diMuonMass_lumi;
  TFile* f1 = (TFile*) TFile::Open ("Normalized_DYLL/AfterTau_AgainstMuonLoose2/normalizedHistos_DYToLL_0j.root");
  //f->cd("AfterTau_AgainstMuonTight2");
  TH1F* num_pt = muon_pT_lumi;
  TH1F* num_eta = muon_eta_lumi;
  TH1F* num_m   = diMuonMass_lumi;

  den_pt->Rebin(5);
  num_pt->Rebin(5);

  den_m->Rebin(10);
  num_m->Rebin(10);
  
  TCanvas* c1 = new TCanvas("c1", "c1");   
  TGraphAsymmErrors* gr1 = new TGraphAsymmErrors( num_pt, den_pt, "b(1,1) mode" );
  gr1->SetMarkerStyle(20);
  gr1->SetLineColor(kBlack);
  gr1->SetFillStyle(0);
  gr1->SetMaximum(1.05);
  gr1->GetXaxis()->SetTitle("p_{T}(#mu)");
  gr1->GetYaxis()->SetTitle("Efficiency");

  TGraphAsymmErrors* gr2 = new TGraphAsymmErrors( num_eta, den_eta, "b(1,1) mode" );
  gr2->SetMarkerStyle(20);
  gr2->SetLineColor(kBlack);
  gr2->SetFillStyle(0);
  gr2->SetMaximum(1.05);
  gr2->GetXaxis()->SetTitle("#eta(#mu)");
  gr2->GetYaxis()->SetTitle("Efficiency");
  
  TGraphAsymmErrors* gr3 = new TGraphAsymmErrors( num_m, den_m, "b(1,1) mode" );
  gr3->SetMarkerStyle(20);
  gr3->SetLineColor(kBlack);
  gr3->SetFillStyle(0);
  gr3->SetMaximum(1.1);
  gr3->GetXaxis()->SetTitle("m(#mu, #mu) [GeV]");
  gr3->GetYaxis()->SetTitle("Efficiency");

  TLegend* legend = new TLegend(0.5287356,0.8224101,0.8505747,0.8858351,NULL,"brNDC");
  legend->SetTextFont(42);
  legend->SetLineColor(1);
  legend->SetLineStyle(1);
  legend->SetLineWidth(1);
  legend->SetFillColor(0);
  legend->SetFillStyle(1001);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(gr1, "Z #rightarrow #mu#mu + jets");
   
  gr1->Draw("AP");  
  legend->Draw();

  TCanvas* c2 = new TCanvas("c2", "c2");
  c2->cd();
  gr2->Draw("AP");  
  legend->Draw();

  TCanvas* c3 = new TCanvas("c3", "c3");
  c3->cd();
  gr3->Draw("AP");
  legend->Draw();

}