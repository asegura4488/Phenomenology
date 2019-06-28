//*******************************************************************************
//      Filename:       StackPlotter.cpp
//      Authors:        Andres Florez (Universidad de los Andes) 
//                      Alfredo Gurrola (Vanderbilt University)
//      Last update:    Jan. 31, 2014
//
//      This function is called by the shell script runPlotter.sh. It opens ROOT files listed in
//      the configuration file NormalizeHistos_${Directory}_SAFE.in and extracts the 1D
//      histograms. The probability and normalized histograms are output to ROOT files
//      and a log file is generated.
//*******************************************************************************

//--------------------------------------------------------------------------------------------------------------
// Load supplied header file and source codes
#include "include/StackPlotter.hpp"
#include "src/LogFile.cpp"
//--------------------------------------------------------------------------------------------------------------

// Main function to normalize histograms
void StackPlotter (string inRootFilePath, string configFileList, string configFilePath)

{

  // use CMS Style 
  gROOT->LoadMacro("src/cmsStyle.C");
  SetcmsStyle();

  // Store global variables
  inputRootFilePath = inRootFilePath; 
  readConfigFile (configFileList, configFilePath);
  readConfigPlotStyle (configFilePath); 
            
  // Create log file
  string logFilePath = inputRootFilePath + "/" + logFilename;
  LogFile* logTextFile = new LogFile (logFilePath);

  getHistoInfo (logTextFile);
  makeStackPlots(configFilePath);
  delete logTextFile;           // Finished writing log file 
   
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
      string rootFileName = params->GetValue ("Normalizedrootfile", "file");
      inputRootFile.push_back (rootFileName.c_str() ); 

      // Get process name: data, ttbar, wjets etc.
      string process = params->GetValue ("process", "theprocess"); 
      inputProcess.push_back(process.c_str());

      // Apply scale factor?
      string applyScaleFactor = params->GetValue("applyScaleFactor", "False");
      inputApplyScaleFactor.push_back (applyScaleFactor.c_str());

      // Scale Factor
      float scaleFactor = params->GetValue("scaleFactor", 1.0);
      inputScaleFactor.push_back (scaleFactor);
    
      // Scale factor error
      float scaleFactorError = params->GetValue("scaleFactorError", 0.0);   
      inputScaleFactorError.push_back (scaleFactorError);

      // Rebin histos? 
      string rebinHistoName = params->GetValue("rebinHistoName", "h_name");
      inputRebinHisto.push_back (rebinHistoName.c_str());

      // rebin Factor
      int rebin = params->GetValue("rebinFactor", 1);
      inputRebinFactor.push_back ( rebin ); 

      // Line style 
      int lineStyle = params->GetValue("lineStyle", 1);
      inputSetLineStyle.push_back ( lineStyle );

      // Line width
      int lineWidth = params->GetValue("lineWidth", 1);
      inputSetLineWidth.push_back ( lineWidth );

      // Line Color
      int thecolor = params->GetValue("color" , kRed); 
      inputColor.push_back (thecolor);

      // Store TLegend text
      string thelegend = params->GetValue("legend" , "histo");
      inputLegend.push_back (thelegend.c_str());

      // Use collision data in your plots? --> 99% of the time yes...just in case you do not want to
      
      use_data = params->GetValue("usedata", "True"); 

      use_log  = params->GetValue("setLogScale", "False");

      string isSignal = params->GetValue("isSignal", "False");
      inputIsSignal.push_back ( isSignal.c_str() );

      // Output root file name
      stackedHistosOutputRootFile = params->GetValue("outputRootFileStackedHistos", "stackedHistos.root");
      
      // log file name
      logFilename  = params->GetValue("logFilename", "outputLog_stackPlots.log");

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
      
      else if (inputStyle == "xrange_min")
	inputXrangeMin.push_back (atof (inputString.c_str()) );
      
      else if (inputStyle == "xrange_max")
	inputXrangeMax.push_back (atof (inputString.c_str()) );
      
      else if (inputStyle == "yrange_min")
	inputYrangeMin.push_back (atof (inputString.c_str()) );
      
      else if (inputStyle == "yrange_max")
	inputYrangeMax.push_back (atof (inputString.c_str()) );
      
      else if (inputStyle == "xtitle")
	inputXtitle.push_back(inputString);
      
      else if (inputStyle == "ytitle")
	inputYtitle.push_back(inputString);
      
      else if (inputStyle == "rebin")
	inputRebin.push_back(atoi (inputString.c_str()) );
      
    }
  inFile.close();
}

TH1F* SetStyleh(TH1F* h_j)
{
  // set style
  for (int o = 0; o < inputHistoName.size(); o++)
    {
      if ( h_j->GetName() == inputHistoName.at(o) )
        {
          h_j->GetXaxis()->SetRangeUser(inputXrangeMin.at(o), inputXrangeMax.at(o));
          h_j->SetMinimum(inputYrangeMin.at(o));
          h_j->SetMaximum(inputYrangeMax.at(o));
          h_j->GetXaxis()->SetLabelFont(42);
          h_j->GetXaxis()->SetLabelSize(0.04);
          h_j->GetXaxis()->SetTitleSize(0.09);
          h_j->GetXaxis()->SetTitleOffset(0.70);
          h_j->SetYTitle(inputYtitle.at(o).c_str());
          h_j->SetXTitle(inputXtitle.at(o).c_str());
          break;
        }
    }
  return h_j;
}

TH1F* Rebinh(TH1F* h_j)
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


void getHistoInfo (LogFile* logFile)
{
 
  listHistos = new TList();    // Histograms normalized to cross-section & luminosity
  nHistList = 0;

  // loop over root files
  for (int j=0; j< inputRootFile.size(); j++) 
    {

      int numMax = 0;

      // Access input ROOT file (can access over secure shell)
      TFile *rootTFile = new TFile((inputRootFile.at(j)).c_str()); // open root file
      
      // Change directory accessed by ROOT to cutflow path
      TDirectory *current_sourcedir = gDirectory;
      TIter nextkey( current_sourcedir->GetListOfKeys() );
      TKey *key;
      
      // loop over keys(ROOT objects) within the root file
      while ((key = (TKey*)nextkey())) 
        {
	  TObject *obj = key->ReadObj();
	  
	  // Get only 1D histograms from ROOT file and store in a TList 
	  if ( obj->IsA()->InheritsFrom( "TH1" ) ) 
	    {
	      string histname = obj->GetName();       
	      size_t found = histname.find("lumi");
	      string theProcess = "lumi";
	      
	      if (found!=string::npos)
		{
		  
		  TH1 *histoObj = (TH1*)obj;
		  // rename histogram
		  string new_theProcess = "_" + theProcess;
		  histname.erase(found-1, new_theProcess.size());
		  string histoName = histname;
		 
                  // Write info to log file
                  if (j == 0)
                    logFile->writeHistoName (histoName);
 
                  // Get bin & edge values of the 1D histogram
                  int numBins = histoObj->GetXaxis()->GetNbins();
                  double lowerEdge = histoObj->GetXaxis()->GetXmin();
                  double upperEdge = histoObj->GetXaxis()->GetXmax();
	
                  TH1F* h1 = new TH1F (histoName.c_str(), histoName.c_str(),
				       numBins, lowerEdge, upperEdge);
                  // Apply scale factor?
		  if (inputApplyScaleFactor.at(j) == "True")
		    {
		      for (int k=0; k <= (h1->GetXaxis()->GetNbins() + 1); k++)
			{
			  double scf = inputScaleFactor.at(j);
			  h1->SetBinContent(k,(histoObj->GetBinContent(k)*scf));
			}   
		    } else {

                       for (int k=0; k <= (h1->GetXaxis()->GetNbins() + 1); k++)
                        {
                          h1->SetBinContent(k,(histoObj->GetBinContent(k)));
                        }
                    }
		  
                  // Rebin historgams ?
                  for (int l=0; l < inputRebinHisto.size(); l++)
                    { 
		      if (histname == (inputRebinHisto.at(l)).c_str())  
		        {
		          h1->Rebin(inputRebinFactor.at(l));
                          break;
		        } 
                    }
                  h1 = Rebinh(h1);
                  listHistos->Add(h1);
		  nHistList++;
		  theHistNameStrings.push_back(histname);
		  numMax++;
		  
		}  
	    } // close if loop InheritsFrom("TH1")
	} // close while loop 
      cout<<"Number histos "<<numMax<<"  File "<<inputRootFile.at(j).c_str()<<endl; 
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
} // close getHistoInfo function

void makeStackPlots (string configFilePath){

  if (( nHistList % inputRootFile.size() == 0 )) 
    {
      for (int i=0; i<numHistos; i++) 
        {
	 
          TCanvas *c = new TCanvas(theHistNameStrings.at(i).c_str(), "", 81,58,500,602);
          if (use_log == "True" && use_data == "False" ) gPad->SetLogy();
          if (use_data == "True")
           {
             TPad *pad1 = new TPad("pad1", "pad1",0,0.2176166,1,1);
             pad1->Draw();
             pad1->cd();

             if (use_log == "True" ) gPad->SetLogy(); 
           }
          gStyle->SetOptTitle(0); 
          TH1F* h_clone;
          TH1F* h_mc;      // histograms for BGs
          TH1F* h_data;    // histogram to store data points 
          TH1F* h_signal;   // histogram for signal 
          TH1F* h_allBG;   // histogram with the sum of all BGs
	  TH1F* h_data_clone;
          TH1F* h_allBG_clone;  
          THStack *hs = new THStack(theHistNameStrings.at(i).c_str(),theHistNameStrings.at(i).c_str());
          TLegend *legend = new TLegend(0.60, 0.85-.035*inputLegend.size(), 0.9, 0.9,NULL,"brNDC");

          h_data = (TH1F*)listHistos->At(i);
          h_data = SetStyleh(h_data);
          h_clone = (TH1F*)h_data->Clone();
          if (use_data == "True")
            {
              // Get data histogram from the list
              
              h_data->SetMarkerStyle(20);
              h_data->Sumw2();

              h_data->Draw("ep1");
              h_data->SetLineColor(1);
              legend->AddEntry(h_data,inputLegend.at(0).c_str(), "lep");
              
            }   
          
          int whichHisto = 0;
          int k = 0;
          if(use_data == "True") k=1;    
	  
          int Nbins = h_clone->GetXaxis()->GetNbins();
	  
          h_allBG = new TH1F("h_allBG","allBG", Nbins, h_clone->GetXaxis()->GetXmin(),h_clone->GetXaxis()->GetXmax());
          for (int bg=0; bg <= (Nbins + 1); bg++)
            {  
              h_allBG->SetBinContent(bg,0.0);
            }
	 
          h_allBG->UseCurrentStyle(); 
          h_allBG->SetFillStyle(3004);
          h_allBG->SetFillColor(kRed+2);

          int rootFileCounter = 0;	 
          int rootFileSize = inputRootFile.size();

          for (int j=k;j<inputRootFile.size();j++)
	    {
	      
	      whichHisto=j*numHistos; 
              
	      //---------------- find histo -------------------------------
	      // Make sure that the histograms to be stacked (e.g Mass) 
	      // is the same for every process (e.g z-> mumu, ttbar, Wjets) 
	      // ----------------------------------------------------------
	      
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
              
	      string hist2name = listHistos->At(a+whichHisto)->GetName();
              // Style
              int thelineStyle = inputSetLineStyle.at(j);
              int thelineWidth = inputSetLineWidth.at(j);
              string processlegend = (inputLegend.at(j)).c_str();
              int setcolor = inputColor.at(j);

	      if (foundHisto == true)
		{
		  
		  h_mc = (TH1F*)listHistos->At(a+whichHisto);
                  h_signal = (TH1F*)listHistos->At(a+whichHisto);
                  h_mc = SetStyleh(h_mc);
		  h_mc->SetName(hist2name.c_str());
		  
		  h_mc->SetLineStyle(thelineStyle);
		  h_mc->SetLineWidth(thelineWidth);
		  h_mc->SetLineColor(setcolor);
                  h_mc->SetFillColor(setcolor);                 
 
                  string processlegend = (inputLegend.at(j)).c_str(); 
                  h_mc->SetLineStyle(1);
                  h_mc->SetLineWidth(2);

                  if (inputIsSignal.at(j) == "False")
		    {
                      h_allBG->Add(h_mc);
                      hs->SetHistogram(h_mc);
		      hs->Add(h_mc);
		      legend->AddEntry(h_mc, processlegend.c_str(), "f");
		    } // close if IsSignal
		  else
		    {
		      h_signal->SetLineStyle(thelineStyle);
                      h_signal->SetLineWidth(thelineWidth);
                      h_signal->SetLineColor(setcolor);   
                      legend->AddEntry(h_signal, processlegend.c_str(), "l");
		      //h_signal->Draw("lsame");
		    }
		  
		}// close if (foundHisto)
	      if ((use_data == "True") && (j-1) > 0 && (inputIsSignal.at(j) == "False") )
		{
                   hs->Draw("HISTsame");
                   Double_t hs_ymax = hs->GetMaximum();
                   Double_t ymax = hs_ymax + 0.6*hs_ymax;
                   h_data->SetMaximum(ymax);
                   h_data->Draw("e1psame");
		} 
	      if ( (use_data == "False") && j==0 && (inputIsSignal.at(j) == "False") )
                {
                  hs->SetMaximum(1000.);
                  hs->Draw("HIST");
                  h_allBG->Draw("e2same");
                  gPad->RedrawAxis();
                }
              rootFileCounter++;
	      //if ((inputIsSignal.at(j) == "True") && ((rootFileCounter) == rootFileSize))
              if ((inputIsSignal.at(j) == "True"))
		{
                  h_signal->SetFillColor(0);
                  h_signal->Draw("same");
		}
              if (inputIsSignal.at(j) == "False" && ((rootFileCounter + 1) == rootFileSize))
                {
                  h_allBG->Draw("e2same");
                }
 
	    } // close for loop inputRootFile
          legend->SetTextFont(42);
          legend->SetLineColor(1);
          legend->SetLineStyle(1);
          legend->SetLineWidth(1);
          legend->SetFillColor(0);
          legend->SetFillStyle(1001);
          legend->SetBorderSize(0);
          legend->SetFillColor(kWhite);
          legend->AddEntry(h_allBG, "All BGs (stat uncertainties)", "f");
          legend->Draw();

          c->cd();

          if (use_data == "True")
           {
	     TPad* pad2 = new TPad("pad2", "pad2",0,0.04181185,0.9959677,0.3013937);
	     pad2->Draw();
	     pad2->cd();
	     pad2->Range(-82.03846,-1.26,430.9615,3.21);
	     pad2->SetFillColor(0);
	     pad2->SetBorderMode(0);
	     pad2->SetTickx(1);
	     pad2->SetTicky(1);
	     pad2->SetLeftMargin(0.159919);
	     pad2->SetRightMargin(0.05060729);
	     pad2->SetTopMargin(0.04697987);
	     pad2->SetBottomMargin(0.2818792);
	     pad2->SetFrameBorderMode(0);
	     pad2->SetFrameBorderMode(0); 
	     
	     h_data_clone = (TH1F*)h_data->Clone(); 
	     h_allBG_clone = (TH1F*)h_allBG->Clone();
	     
	     h_data_clone->Divide(h_allBG_clone);
	     h_data_clone->GetYaxis()->SetRangeUser(0, 3);
	     h_data_clone->GetYaxis()->SetTitle("ratio (Data/MC)");
	     h_data_clone->GetYaxis()->SetTitleOffset(0.35);
	     h_data_clone->GetXaxis()->SetLabelFont(42);
	     h_data_clone->GetXaxis()->SetLabelSize(0.14);
	     h_data_clone->GetXaxis()->SetTitleSize(0.14);
	     h_data_clone->GetXaxis()->SetTitleOffset(0.90);
	     h_data_clone->GetXaxis()->SetTitleFont(42);
	     h_data_clone->GetYaxis()->SetLabelFont(42);
	     h_data_clone->GetYaxis()->SetLabelSize(0.13);
	     h_data_clone->GetYaxis()->SetTitleSize(0.14);
	     h_data_clone->GetYaxis()->SetTitleFont(42);
	     
	     h_data_clone->Draw("ep");
	     
	     float startxbin = h_data_clone->GetXaxis()->GetBinWidth(h_data_clone->GetXaxis()->GetFirst()) * 
	       (h_data_clone->GetXaxis()->GetFirst() - 1);
	     float lastxbin  = h_data_clone->GetXaxis()->GetBinWidth(h_data_clone->GetXaxis()->GetFirst()) *
	       (h_data_clone->GetXaxis()->GetLast() - 1);
	     
	     double x1 = 0;
	     double x2 = lastxbin + h_data_clone->GetXaxis()->GetBinWidth(h_data_clone->GetXaxis()->GetNbins())/2.;
	     
	     if (h_data_clone->GetXaxis()->GetBinCenter(0) <= 0)
               x1 = h_data_clone->GetXaxis()->GetBinCenter(0) + h_data_clone->GetXaxis()->GetBinWidth(0)/2. + startxbin;
	     else 
               x1 = h_data_clone->GetXaxis()->GetBinCenter(0) - h_data_clone->GetXaxis()->GetBinWidth(0)/2. - startxbin;
	     
	     TLine *line = new TLine(x1,1.0,x2,1.0);
	     line->SetLineStyle(2);
	     line->SetLineWidth(3);
	     line->SetLineColor(kRed);
	     line->Draw();
	     
	     pad1->cd();
	   }
	  
	  TFile *hfile = (TFile*)gROOT->FindObject(stackedHistosOutputRootFile.c_str());
	  if (hfile) {hfile->Close();}
	  hfile = new TFile(stackedHistosOutputRootFile.c_str(),"UPDATE");
	  
          for (int o = 0; o < inputHistoName.size(); o++)
	    {
	      if ( c->GetName() == inputHistoName.at(o) )
		{
		  string save = configFilePath+"/"+inputHistoName.at(o)+".pdf"; 
		  c->SaveAs(save.c_str());
		  break;
		}
	    }
	  c->Update();
          //pad1->GetUymax();
           
       	  c->Write();
	  hfile->Close();
          c->Close();
	  
	}  // close foor loop numHistos     
      
    } // close if numHistos % inputRootFile.size() == 0
} // close  makeStackPlots function


