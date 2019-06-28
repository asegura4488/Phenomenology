void create_hadronic()
{
TFile *Hadro = new TFile("Hadronic_histograms.root","RECREATE");


TDirectory *theDirectory[5];
theDirectory[0]  = Hadro->mkdir("monotop_left");
theDirectory[1]  = Hadro->mkdir("monotop_right");
theDirectory[2]  = Hadro->mkdir("monotop_full");
theDirectory[3]  = Hadro->mkdir("Drell-yan");
theDirectory[4]  = Hadro->mkdir("single-top");


TFile* f1 = new TFile("normalizedHistos_signal_1jet_left.root");
TH1F* H1 = (TH1F*)f1->Get("Hbjet_ratio_pT_lumi");
TH1F* H2 = (TH1F*)f1->Get("invariant_mjj_lumi");
TH1F* H3 = (TH1F*)f1->Get("invariant_mjjb_lumi");

TFile* f2 = new TFile("normalizedHistos_signal_1jet_right.root");
TH1F* H4 = (TH1F*)f2->Get("Hbjet_ratio_pT_lumi");
TH1F* H5 = (TH1F*)f2->Get("invariant_mjj_lumi");
TH1F* H6 = (TH1F*)f2->Get("invariant_mjjb_lumi");


TFile* f3 = new TFile("normalizedHistos_signal_1jet_full.root");
TH1F* H7 = (TH1F*)f3->Get("Hbjet_ratio_pT_lumi");
TH1F* H8 = (TH1F*)f3->Get("invariant_mjj_lumi");
TH1F* H9 = (TH1F*)f3->Get("invariant_mjjb_lumi");


TFile* f4 = new TFile("normalizedHistos_DYToLL.root");
TH1F* H10 = (TH1F*)f4->Get("Hbjet_ratio_pT_lumi");
TH1F* H11 = (TH1F*)f4->Get("invariant_mjj_lumi");
TH1F* H12 = (TH1F*)f4->Get("invariant_mjjb_lumi");

TFile* f5 = new TFile("normalizedHistos_singletop.root");
TH1F* H13 = (TH1F*)f5->Get("Hbjet_ratio_pT_lumi");
TH1F* H14 = (TH1F*)f5->Get("invariant_mjj_lumi");
TH1F* H15 = (TH1F*)f5->Get("invariant_mjjb_lumi");


Hadro->cd();
{
theDirectory[0]->cd();
H1->Write();
H2->Write();
H3->Write();
theDirectory[1]->cd();
H4->Write();
H5->Write();
H6->Write();
theDirectory[2]->cd();
H7->Write();
H8->Write();
H9->Write();
theDirectory[3]->cd();
H10->Write();
H11->Write();
H12->Write();
theDirectory[4]->cd();
H13->Write();
H14->Write();
H15->Write();

}
}
