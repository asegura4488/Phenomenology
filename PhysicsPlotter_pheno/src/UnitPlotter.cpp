//*******************************************************************************
//      Filename:       UnitPlotter.cpp
//      Authors: Andres Florez <andres.florez@cern.ch>
//      Last update:    Mar 28, 2013
//
//      This function is called by the shell script RunDigiHistos.sh. It opens ROOT files listed in
//      the configuration file UnitPlotter.in and extracts the 1D histograms.
//      Then, normalizes the histograms to unity and overlay data and MC in one canvas.
//*******************************************************************************

//--------------------------------------------------------------------------------------------------------------
// Load supplied header file and source codes
#include "include/UnitPlotter.hpp"
//--------------------------------------------------------------------------------------------------------------

// Main function to normalize histograms
void UnitPlotter (string inRootFilePath, string configFileList, string configFilePath, string saveHistosPath)

{

  // use ATLAS style
  gROOT->LoadMacro("src/cmsStyle.C");
  SetcmsStyle();

  // Store global variables
  inputRootFilePath = inRootFilePath; 
  readConfigFile (configFileList, configFilePath);
  readConfigPlotStyle (configFilePath); 

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
      string rootFileName = params->GetValue ("Normalizedrootfile", "file");
      inputRootFile.push_back (rootFileName.c_str() ); 

      string isdata = params->GetValue ("isdata", "False");
      inputIsData.push_back (isdata.c_str());  

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

      // Line style 
      int lineStyle = params->GetValue("lineStyle", 1);
      inputSetLineStyle.push_back ( lineStyle );

      // Line width
      int lineWidth = params->GetValue("lineWidth", 1);
      inputSetLineWidth.push_back ( lineWidth );

      // Line Color
      int thecolor = params->GetValue("color" , kRed); 
      inputColor.push_back (thecolor);

      // Marker Style
      int themarker = params->GetValue("marketStyle", 20);
      inputMarkerStyle.push_back(themarker); 

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

  string configPlotList=configFilePath + "/style_unit.in"; 
 
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

TH1* SetStyleh(TH1* h_j)
{
   h_j->SetTitle("");
  // set style
  for (int o = 0; o < inputHistoName.size(); o++)
    {
      if ( h_j->GetName() == inputHistoName.at(o) )
        {
	  if(inputRebin.at(o)!=-1) h_j = normalizeToOne (h_j);
          h_j->GetXaxis()->SetRangeUser(inputXrangeMin.at(o), inputXrangeMax.at(o));
          h_j->GetYaxis()->SetRangeUser(inputYrangeMin.at(o), inputYrangeMax.at(o));
          h_j->SetXTitle(inputXtitle.at(o).c_str());
          h_j->SetYTitle(inputYtitle.at(o).c_str());
          h_j->GetXaxis()->SetTitleSize(0.05);
          h_j->GetYaxis()->SetTitleSize(0.05);
          break;
        }
    }
  return h_j;
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

//=======================================================
// Normalize histograms to one
//=======================================================

TH1* normalizeToOne ( TH1* hp)
{

   hp->Sumw2();

   if (hp->Integral() > 0) //by default Integral() does not consider
			   //overflow
     {
       hp->Scale (1.0 / hp->Integral ());
     }

   return hp;

}

void getHistoInfo ()
{
 
  listHistos = new TList();    // Histograms normalized to one
  nHistList = 0;
  
  // loop over root files
  for (int j=0; j< inputRootFile.size(); j++) 
    {
      int numMax = 0;

      // Access input ROOT file (can access over secure shell)
      
      TFile *rootTFile = new TFile((inputRootFile.at(j)).c_str()); // open root file
      std::cout <<"inputRootFile.at(j)).c_str()  "<<inputRootFile.at(j).c_str()<<std::endl;
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
 
	      int numBins = histoObj->GetXaxis()->GetNbins();
	      double lowerEdge = histoObj->GetXaxis()->GetXmin();
	      double upperEdge = histoObj->GetXaxis()->GetXmax();

	      if ( obj->IsA()->InheritsFrom( "TH2" ) ) 
		{
		  int numBinsY = histoObj->GetYaxis()->GetNbins();
		  double lowerEdgeY = histoObj->GetYaxis()->GetXmin();
		  double upperEdgeY = histoObj->GetYaxis()->GetXmax();
		  h1 = new TH2F (histName.c_str(), histName.c_str(),
				 numBins, lowerEdge, upperEdge,
				 numBinsY, lowerEdgeY, upperEdgeY);
	      
		  //here (and below) GetSize() will allow us to
		  //iterate over 1 or 2d hists.  We pick up overflow
		  //bins too, but that's fine
		  TH2F* tmp1 = static_cast<TH2F*>(h1);
		  totalbins = tmp1->GetSize();
		} else {
	      h1 = new TH1F (histName.c_str(), histName.c_str(),
				   numBins, lowerEdge, upperEdge);
	      TH1F* tmp2 = static_cast<TH1F*>(h1);
	      totalbins = tmp2->GetSize();
	      }

	      // Apply scale factor?
	      if (inputApplyScaleFactor.at(j) == "True")
		{
		  for (int k=0; k <= totalbins; k++)
		    {
		      double scf = inputScaleFactor.at(j);
		      h1->SetBinContent(k,(histoObj->GetBinContent(k)*scf));
		      h1->SetBinError(k,(histoObj->GetBinError(k))*scf);
		    }
		} else
		{
		  for (int k=0; k <= totalbins; k++)
		    {
		      h1->SetBinContent(k,(histoObj->GetBinContent(k)));
		      h1->SetBinError(k,(histoObj->GetBinError(k)));
		    }
		}
              h1 = Rebinh(h1);
	      listHistos->Add(h1);
	      nHistList++;
	      theHistNameStrings.push_back(histName);
	      numMax++;
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
} // close getHistoInfo function

void makePlots (string configFilePath){

  if (( nHistList % inputRootFile.size() == 0 )) 
    {
      for (int i=0; i<numHistos; i++) 
        {
	 
          TCanvas *c = new TCanvas(theHistNameStrings.at(i).c_str(), "", 81,58,500,602);
 
          TH1* h_mc;        // histograms for MC
          TH1* h_mc_clone;  // histograms for MC 
          TH1* h_data;      // histogram to store data points 
	  TH1* h_data_clone;
          TLegend *legend = new TLegend(0.2, 0.85-.035*inputLegend.size(), 0.5, 0.90,NULL,"brNDC");
          int whichHisto = 0;
	  
          for (int j=1;j<inputRootFile.size();j++)
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
              
	      string hist2name = listHistos->At(a+whichHisto)->GetName();
	      if (foundHisto == true)
		{
                  h_data = (TH1*)listHistos->At(i);
                  h_data->SetMarkerStyle(inputMarkerStyle.at(0));
                  h_data = SetStyleh(h_data);
                  h_data->Sumw2();
                  int thelineStyle = inputSetLineStyle.at(j);
                  int thelineWidth = inputSetLineWidth.at(j);
                  int setcolor = inputColor.at(j);
                  int marker   = inputMarkerStyle.at(j);
                  h_data->SetMarkerColor(inputColor.at(0));
                  h_data->SetLineColor(inputColor.at(0));
                  if (j == 1) 
                    { 
		       if ( h_data->IsA()->InheritsFrom( "TH2" ) ) {
			 //n.b. that ProfileX will do nothing if the
			 //name you choose is alreayd taken!  so each
			 //profileX needs a unique name
			 TString profileName = TString::Format("%s_data_pfx",h_data->GetName());
			 TProfile *dataProf = ((TH2*)h_data)->ProfileX(profileName.Data());
			 dataProf->SetMinimum(h_data->GetMinimum()); //min/max values were stored in the original hist
			 dataProf->SetMaximum(h_data->GetMaximum());
			 dataProf->Draw("");
			 
			 //profileX cretaes a clone, so no need to
			 //clone again here.  Need to set the hist
			 //name so root doesn't complain
			 TString cloneName = TString::Format("%s_data_pfx_clone",h_data->GetName());
			 TString cloneName2 = TString::Format("%s_data_pfx_clone_px",h_data->GetName());
			 h_data_clone = ((TH2*)h_data)->ProfileX(cloneName.Data())->ProjectionX(cloneName2.Data());
		       }
		       else { //1d case
                         Float_t ymaxc = gPad->GetUymax();
                         //h_data->SetMaximum(ymaxc+0.15*ymaxc);
			 h_data->Draw("ep1");
			 h_data_clone = (TH1*)h_data->Clone();
		       }
                       //h_data->SetLineColor(1);
                       h_data->SetLineStyle(thelineStyle);
                       h_data->SetLineWidth(thelineWidth);
                       float rms_data = h_data->GetRMS();
                       float mean_data = h_data->GetMean();
                       std::stringstream legend_data;
                       legend_data.str("");
                       legend_data << inputLegend.at(0).c_str() << ": mean=" << setprecision(3)
                                   <<mean_data <<", rms=" << setprecision(3) <<rms_data<< endl;
                       legend->AddEntry(h_data,legend_data.str().c_str(), "lep");

                    }
		  h_mc = (TH1*)listHistos->At(a+whichHisto);
                  h_mc = SetStyleh(h_mc);
		  h_mc->SetName(hist2name.c_str());
                  //cout<<"NAME "<<h_mc->GetName()<<endl;
		  // Style
                  
		  string processlegend = (inputLegend.at(j)).c_str();
 
		  h_mc->SetLineStyle(thelineStyle);
		  h_mc->SetLineWidth(thelineWidth);
		  h_mc->SetLineColor(setcolor);
                  h_mc->SetMarkerStyle(marker);
                  h_mc->SetMarkerColor(setcolor);                 
 
                  string processlegend = (inputLegend.at(j)).c_str(); 

		  if ( h_mc->IsA()->InheritsFrom( "TH2" ) ) {
		    TString profileName = TString::Format("%s_mc_pfx", h_mc->GetName());
		    ((TH2*)h_mc)->ProfileX(profileName.Data())->Draw("SAME");
		    //profileX cretaes a clone, so no need to
		    //clone again here.  Need to set the hist
		    //name so root doesn't complain
		    TString cloneName = TString::Format("%s_mc_pfx_clone",h_mc->GetName());
		    TString cloneName2 = TString::Format("%s_mc_pfx_clone_px",h_mc->GetName());
		    h_mc_clone = ((TH2*)h_mc)->ProfileX(cloneName.Data())->ProjectionX(cloneName2.Data());
		  } else { //1d case
                    Float_t ymaxc = gPad->GetUymax();
                    h_mc->SetMaximum(ymaxc+0.15*ymaxc);
		    h_mc->Draw("epsame");
		    h_mc_clone = (TH1*)h_mc->Clone();
		  }
                  h_mc_clone->SetMarkerStyle(marker);
                  h_mc_clone->SetMarkerColor(setcolor);
                  float rms_mc = h_mc->GetRMS();
                  float mean_mc = h_mc->GetMean();
		  std::stringstream legend_mc;
                  legend_mc.str("");
                  legend_mc << processlegend.c_str() << " - "<< "mean: " << setprecision(3)
                              <<mean_mc <<" , rms: " << setprecision(3) <<rms_mc<< endl;
                  legend->AddEntry(h_mc, legend_mc.str().c_str(), "lep"); 
		}// close if (foundHisto)
 
	      legend->SetTextFont(42);
              legend->SetTextSize(.04);
              legend->SetMargin(0.15);
              legend->SetLineColor(1);
              legend->SetLineStyle(1);
              legend->SetLineWidth(1);
              legend->SetFillColor(0);
              legend->SetFillStyle(0);
              legend->SetBorderSize(0);
              legend->SetFillColor(kWhite);
              legend->Draw();
 
	      c->cd();
 	      
              if (j == 1)
                {
	           float startxbin = h_data_clone->GetXaxis()->GetBinWidth(h_data_clone->GetXaxis()->GetFirst()) * 
		      (h_data_clone->GetXaxis()->GetFirst() - 1);
	           float lastxbin  = h_data_clone->GetXaxis()->GetBinCenter(h_data_clone->GetXaxis()->GetNbins());
	      
	            double x1 = 0;
	            double x2 = lastxbin + h_data_clone->GetXaxis()->GetBinWidth(h_data_clone->GetXaxis()->GetNbins())/2.;
	      
	            if (h_data_clone->GetXaxis()->GetBinCenter(0) <= 0)
		        x1 = h_data_clone->GetXaxis()->GetBinCenter(0) + h_data_clone->GetXaxis()->GetBinWidth(0)/2. + startxbin;
	            else 
		        x1 = h_data_clone->GetXaxis()->GetBinCenter(0) - h_data_clone->GetXaxis()->GetBinWidth(0)/2. - startxbin;
	        } 
	  
	    } // close for loop inputRootFile
	  
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


