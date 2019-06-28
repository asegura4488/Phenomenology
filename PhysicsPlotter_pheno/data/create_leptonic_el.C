void create_leptonic_el()
{
TFile *Lepton = new TFile("Leptonic_histograms_el.root","new");


TDirectory *theDirectory[8];
theDirectory[0]  = Lepton->mkdir("monotop_left");
theDirectory[1]  = Lepton->mkdir("monotop_right");
theDirectory[2]  = Lepton->mkdir("monotop_full");
theDirectory[3]  = Lepton->mkdir("Drell-yan");
theDirectory[4]  = Lepton->mkdir("W+jets");
theDirectory[5]  = Lepton->mkdir("single-top");
theDirectory[6]  = Lepton->mkdir("ttbar-semileptonic");
theDirectory[7]  = Lepton->mkdir("ttbar-dileptonic");


TFile* f1 = new TFile("normalizedHistos_signal_1jet_left.root");
TH1F* H1 = (TH1F*)f1->Get("bjet_ratio_pT_lumi");
TH1F* H2 = (TH1F*)f1->Get("mt_top_lumi");

TFile* f2 = new TFile("normalizedHistos_signal_1jet_right.root");
TH1F* H3 = (TH1F*)f2->Get("bjet_ratio_pT_lumi");
TH1F* H4 = (TH1F*)f2->Get("mt_top_lumi");

TFile* f3 = new TFile("normalizedHistos_DYToLL.root");
TH1F* H5 = (TH1F*)f3->Get("bjet_ratio_pT_lumi");
TH1F* H6 = (TH1F*)f3->Get("mt_top_lumi");

TFile* f4 = new TFile("normalizedHistos_wjets.root");
TH1F* H7 = (TH1F*)f4->Get("bjet_ratio_pT_lumi");
TH1F* H8 = (TH1F*)f4->Get("mt_top_lumi");


TFile* f5 = new TFile("normalizedHistos_singletop.root");
TH1F* H9 = (TH1F*)f5->Get("bjet_ratio_pT_lumi");
TH1F* H10 = (TH1F*)f5->Get("mt_top_lumi");

TFile* f6 = new TFile("normalizedHistos_ttbar_semi.root");
TH1F* H11 = (TH1F*)f6->Get("bjet_ratio_pT_lumi");
TH1F* H12 = (TH1F*)f6->Get("mt_top_lumi");

TFile* f7 = new TFile("normalizedHistos_ttbar_di.root");
TH1F* H13 = (TH1F*)f5->Get("bjet_ratio_pT_lumi");
TH1F* H14 = (TH1F*)f5->Get("mt_top_lumi");

TFile* f8 = new TFile("normalizedHistos_signal_1jet_full.root");
TH1F* H15 = (TH1F*)f8->Get("bjet_ratio_pT_lumi");
TH1F* H16 = (TH1F*)f8->Get("mt_top_lumi");


Lepton->cd();
{
theDirectory[0]->cd();
H1->Write();
H2->Write();
theDirectory[1]->cd();
H3->Write();
H4->Write();
theDirectory[2]->cd();
H15->Write();
H16->Write();
theDirectory[3]->cd();
H5->Write();
H6->Write();
theDirectory[4]->cd();
H7->Write();
H8->Write();
theDirectory[5]->cd();
H9->Write();
H10->Write();
theDirectory[6]->cd();
H11->Write();
H12->Write();
theDirectory[7]->cd();
H13->Write();
H14->Write();

}
}
