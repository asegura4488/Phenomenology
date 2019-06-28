//*******************************************************************************
//      Filename:       EfficiencyPlotter.cpp
//      Authors: Andres Florez <andres.florez@cern.ch>
//      Last update:    Mar 28, 2013
//
//      This function is called by the shell script RunDigiHistos.sh. It opens ROOT files listed in
//      the configuration file EfficiencyPlotter.in and extracts the 1D histograms.
//      Then, normalizes the histograms to unity and overlay data and MC in one canvas.
//*******************************************************************************

//--------------------------------------------------------------------------------------------------------------
// Load supplied header file and source codes
#include "include/EfficiencyPlotter.hpp"
//--------------------------------------------------------------------------------------------------------------

// Main function to normalize histograms
void EfficiencyPlotter (string inRootFilePath, string configFileList, string configFilePath, string saveHistosPath)

{


  // use ATLAS style
  gROOT->LoadMacro("src/cmsStyle.C");
  SetcmsStyle();

  // Store global variables
  readConfigFile (configFileList, configFilePath);
  readConfigPlotStyle (configFilePath); 
  nHistList = 0;
  getHistoInfo ();
  makePlots(configFilePath);
}

void readConfigFile (string configFileList, string configFilePath)
{
  // Input stream for data file
  ifstream inFile;       
  // Open input file
  inFile.open (configFileList.c_str(), ios::in);
  
  // If can't open input file, exit the code
  if (!inFile)
    {
      cerr << "ERROR: Can't open input file: " << configFileList << endl;
      exit (1);
    }

  // ROOT configuration file


  while (inFile >> inputType) 
    {

      string configfile=configFilePath + "/" + inputType;

      TEnv *params = new TEnv ("config_file");
      params->ReadFile (configfile.c_str(), kEnvChange);
      
      // Get Root file name
      string rootFileNameNum = params->GetValue ("NormalizedrootfileNum", "file");
      inputRootFileNum.push_back (rootFileNameNum.c_str() ); 

      string rootFileNameDen = params->GetValue ("NormalizedrootfileDen", "file");
      inputRootFileDen.push_back (rootFileNameDen.c_str() );

      // Get process name: data, ttbar, wjets etc.
      string process = params->GetValue ("process", "theprocess"); 
      inputProcess.push_back(process.c_str());

      // Marker Style
      int themarker = params->GetValue("marketStyle", 20);
      inputMarkerStyle.push_back(themarker); 

      // Line Color
      int thecolor = params->GetValue("color" , kRed);
      inputColor.push_back (thecolor);

      // Store TLegend text
      string thelegend = params->GetValue("legend" , "histo");
      inputLegend.push_back (thelegend.c_str());

      // Use collision data in your plots? --> 99% of the time yes...just in case you do not want to
      
      string isSignal = params->GetValue("isSignal", "False");
      inputIsSignal.push_back ( isSignal.c_str() );

      // Output root file name
      HistosOutputRootFile = params->GetValue("outputRootFileHistos", "Histos.root");

    }
  
  inFile.close();
}

void readConfigPlotStyle (string configFilePath)
{
  // Input stream for data file
  ifstream inFile;

  string configPlotList=configFilePath + "/style.in"; 
 
  // Open input file
  inFile.open (configPlotList.c_str(), ios::in);
  
  // If can't open input file, exit the code
  if (!inFile)
    {
      cerr << "ERROR: Can't open input file: " << configPlotList << endl;
      exit (1);
    }
  
  // ROOT configuration file
  while (inFile >> inputStyle >> inputString)
    {
      
      if (inputStyle == "histoName")
	inputHistoName.push_back (inputString );
      
      else if (inputStyle == "xtitle")
	inputXtitle.push_back(inputString);
      
      else if (inputStyle == "ytitle")
	inputYtitle.push_back(inputString);
      
      else if (inputStyle == "rebin")
	inputRebin.push_back(atoi (inputString.c_str()) );
     
      else if (inputStyle == "differentBinSize")
        inputdifferentBinSize.push_back(atoi (inputString.c_str()) ); 

      else if (inputStyle == "binSize")
        inputBinSize.push_back(atoi (inputString.c_str()) ); 

      else if (inputStyle == "Nbins")
        inputNbins.push_back(atoi (inputString.c_str()) );

      else if (inputStyle == "NbinMax")
        inputNbinMax.push_back(atoi (inputString.c_str()) );
      else if (inputStyle == "yrange_max")
        inputYrangeMax.push_back (atof (inputString.c_str()) );
    }
  inFile.close();
}

TH1* Rebinh (TH1* h_j)
{
  // set style
  for (int o = 0; o < inputHistoName.size(); o++)
    {
      if ( h_j->GetName() == inputHistoName.at(o) )
        {
          h_j->Rebin(inputRebin.at(o));
          break;
        }
    }
  return h_j;
}

void getHistoInfo ()
{
  listHistosNum = new TList(); 
  listHistosDen = new TList(); 
  bool boolNameStrings = true;
  listHistosNum = readHistos(inputRootFileNum, boolNameStrings);
  boolNameStrings = false;
  listHistosDen = readHistos(inputRootFileDen, boolNameStrings); 

} // close getHistoInfo function


TList* readHistos(vector<string> inputRootFiles, bool fillNameStrings){
  TList* thelistHistos = new TList();
  for (int j=0; j< inputRootFiles.size(); j++) 
    {
      int numMax = 0;
      //cout <<"inputProcess: "<<inputProcess.at(j)<<endl;
      // Access input ROOT file (can access over secure shell)
      cout <<"inputRootFiles.at(j) "<<inputRootFiles.at(j)<<endl; 
      TFile *rootTFile = new TFile((inputRootFiles.at(j)).c_str()); // open root file
      TDirectory *current_sourcedir = gDirectory;
      TIter nextkey( current_sourcedir->GetListOfKeys() );
      TKey *key;
      
      // loop over keys(ROOT objects) within the root file
      while ((key = (TKey*)nextkey())) 
        {
	  TObject *obj = key->ReadObj();
	  
	  TH1 *histoObj = (TH1*)obj;
	  TH1* h1 = 0; //points to the new copied hist we will make
          int totalbins = 0;
	  // rename histogram
          if ( obj->IsA()->InheritsFrom( "TH1" ) ) 
            {
              string histname = obj->GetName();
	      size_t found = histname.find("lumi");
              string theProcess = "lumi";

              string new_theProcess = "_" + theProcess;
              histname.erase(found-1, new_theProcess.size());
              string histName = histname;
              cout << histName << endl;
              for (int o = 0; o < inputHistoName.size(); o++)
                {   
                  if ( histName == inputHistoName.at(o) )
		    {
                      int numBins = histoObj->GetXaxis()->GetNbins();
		      double lowerEdge = histoObj->GetXaxis()->GetXmin();
		      double upperEdge = histoObj->GetXaxis()->GetXmax();
		      
		      h1 = new TH1F (histName.c_str(), histName.c_str(),
			             numBins, lowerEdge, upperEdge);
                      
                      TH1F* tmp2 = static_cast<TH1F*>(h1);
                      totalbins = tmp2->GetSize();

                      for (int k=0; k <= totalbins; k++)
                        {
                          h1->SetBinContent(k,(histoObj->GetBinContent(k)));
                          h1->SetBinError(k,(histoObj->GetBinError(k)));
                        }

                      h1 = Rebinh(h1);
                      //cout<<"h1 "<<h1->GetName()<<"   Integral  "<<h1->Integral()<<endl; 
		      thelistHistos->Add(h1);
		      if(fillNameStrings) nHistList++;
		      if(fillNameStrings) theHistNameStrings.push_back(histName);
		      numMax++;
		    } // close if loop on selected histograms
		} // close for loop on selected histograms
	    } // close if loop InheritsFrom("TH1")
	} // close while loop
      // If first input ROOT file (Data), store total number of histograms
      if (j == 0)
          numHistos = numMax;
      // Else exit the code if a ROOT file has a different number of histograms
      else
        {
          if (numMax != numHistos)
            {
              cerr << "ERROR: Input Root Files DO NOT have the same number"
                   << " of histograms." << endl;
              exit (1);
            }
        }
    } // close for loop over root files

    return thelistHistos;
}

void makePlots (string configFilePath){
  if ((  inputRootFileNum.size() > 0 )) 
    {
      for (int i=0; i<numHistos; i++) 
        {
          double y_max = 0.0;	 
          TCanvas *c = new TCanvas(theHistNameStrings.at(i).c_str(), "", 81,58,500,602);
 
          TH1* h_num;        // histograms for Numerator
          TH1* h_den;        // histohrams for Denominator
          int whichHisto = 0;

	  TMultiGraph *mg = new TMultiGraph();
          for (int j=0;j<inputRootFileNum.size();j++)
	    {
             whichHisto=j*numHistos; 
	      
	      int a = 0;
	      bool foundHisto = false;
	      if (theHistNameStrings.at(i) == theHistNameStrings.at(i+whichHisto)) 
		{
		  foundHisto=true; 
		  a=i;
		}
	      else 
		{
		  for (a = 0; a < numHistos; a++)
		    {
		      if (theHistNameStrings.at(i) == theHistNameStrings.at(a+whichHisto))
			{
			  foundHisto = true;
			  break;
			}
		    }
		}
	      //-----------------------------------------------------------
	      string hist2name = listHistosNum->At(a+whichHisto)->GetName();
	      if (foundHisto == true)
		{
                  h_num = (TH1*)listHistosNum->At(i+whichHisto);
                  h_den = (TH1*)listHistosDen->At(a+whichHisto);

                  TH1F* H_bins_num;
                  TH1F* H_bins_den;
                  int id = 0;
                  for (int k = 0; k < inputXtitle.size(); k++)
                   {
                     if (hist2name == inputHistoName.at(k)){
                      id=k;
                    }
                   }
                  if (inputdifferentBinSize.at(id)==1){
                     int Nbins = h_num->GetXaxis()->GetNbins();
                     int arrayD = (int*)inputNbins.at(id);
                     float *BINS = new float [arrayD];
                     for (int b = 0; b <= arrayD; b++ )
                       {
                          if (b < (arrayD-1))
                           {
                             BINS[b] = (float)b*(h_num->GetXaxis()->GetBinWidth(b));
                           }
                          else {
                             BINS[b] = (float)inputBinSize.at(id);
                          }
                       }
                     H_bins_num = new TH1F("h_NUM", "h_NUM", inputNbinMax.at(id), BINS);
                     H_bins_den = new TH1F("h_DEN", "h_DEN", inputNbinMax.at(id), BINS);
                     for (int b = 0; b <= (Nbins+1); b++ )
                       {
                         H_bins_num->SetBinContent(b, h_num->GetBinContent(b));
                         H_bins_den->SetBinContent(b, h_den->GetBinContent(b));
                       }
                  }

                  if (inputdifferentBinSize.at(id)==1){
                     TGraphAsymmErrors* gr1 = new TGraphAsymmErrors( H_bins_num, H_bins_den, "b(1,1) mode" );
                  } else{
                    cout << "h_num "<<h_num->GetXaxis()->GetNbins() << "  h_den "<<h_den->GetXaxis()->GetNbins()<<endl;
                    cout << "h_num "<<h_num->Integral() << "  h_den "<<h_den->Integral()<<endl;
                    TGraphAsymmErrors* gr1 = new TGraphAsymmErrors( h_num, h_den, "b(1,1) mode" );
                  }
		  int setcolor = inputColor.at(j);
		  int marker   = inputMarkerStyle.at(j);
                  string processlegend = (inputLegend.at(j)).c_str();
                  gr1->SetMarkerStyle(marker);
                  gr1->SetMarkerColor(setcolor);
                  gr1->SetTitle(processlegend.c_str());
                  gr1->SetFillStyle(0);
                  mg->SetTitle(hist2name.c_str());
                  mg->Add(gr1);
		}// close if (foundHisto)
 
	      c->cd();
 	      c->SetGrid(); 
	    } // close for loop inputRootFile

          mg->Draw("APsame");
          for (int k = 0; k < inputXtitle.size(); k++)
           {
              if (hist2name == inputHistoName.at(k))
               {
                 mg->GetXaxis()->SetTitle((inputXtitle.at(k)).c_str());
                 mg->GetYaxis()->SetTitle((inputYtitle.at(k)).c_str());
                 mg->GetXaxis()->SetTitleSize(0.05);
                 mg->GetYaxis()->SetTitleSize(0.05);
                 mg->SetMaximum((float)inputYrangeMax.at(k));
               }
           }
          c->Update();
          TLegend* legend = c->BuildLegend(0.55,0.8,0.90,0.94);
          //TLegend *legend = new TLegend(0.2, 0.85-.035*inputLegend.size(), 0.5, 0.90,NULL,"brNDC");
          legend->SetTextFont(42);
          legend->SetLineColor(1);
          legend->SetLineStyle(1);
          legend->SetLineWidth(1);
          legend->SetFillColor(0);
          legend->SetFillStyle(1001);
          legend->SetBorderSize(0);
          legend->SetFillColor(kWhite);

	  TFile *hfile = (TFile*)gROOT->FindObject(HistosOutputRootFile.c_str());
	  if (hfile) {hfile->Close();}
	  hfile = new TFile(HistosOutputRootFile.c_str(),"UPDATE");

          for (int o = 0; o < inputHistoName.size(); o++)
             {
                if ( c->GetName() == inputHistoName.at(o) )
                  {
                     string save = configFilePath+"/"+inputHistoName.at(o)+".pdf";  
                     c->SaveAs(save.c_str());
                     c->Write();    
                     break;
                  }
             }

	  hfile->Close();
          c->Close();
	  
	}  // close foor loop numHistos     
      
    } // close if numHistos % inputRootFile.size() == 0
} // close  makePlots function

