//*******************************************************************************
//      Filename:       LimitsRootFileMaker.cpp
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
#include "include/LimitsRootFileMaker.hpp"
#include "src/LogFile.cpp"
//--------------------------------------------------------------------------------------------------------------

// Main function to normalize histograms
void LimitsRootFileMaker (string inRootFilePath, string configFileList, string configFilePath)

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
  makePlotsForLimitsTool(configFilePath);
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

      // systematics on the shape?
      string bg_sys_shape = params->GetValue("inputDataDrivenBG", "False");
      inputDataDrivenBG.push_back(bg_sys_shape.c_str());

      // Store TLegend text
      string thelegend = params->GetValue("legend" , "histo");
      inputLegend.push_back (thelegend.c_str());

      // Use collision data in your plots? --> 99% of the time yes...just in case you do not want to
      
      use_data = params->GetValue("usedata", "True"); 

      use_log  = params->GetValue("setLogScale", "False");

      string isSignal = params->GetValue("isSignal", "False");
      inputIsSignal.push_back ( isSignal.c_str() );

      // log file name
      logFilename  = params->GetValue("logFilename", "outputLog_stackPlots.log");

    }
  
  inFile.close();
}

void readConfigPlotStyle (string configFilePath)
{
  // Input stream for data file
  ifstream inFile;

  string configPlotList=configFilePath + "/style_limits.in"; 
 
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
	inputHistoName = inputString;
      
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
  if ( h_j->GetName() == inputHistoName)
    {
      h_j->Rebin(inputRebin.at(0));
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
		  
                  if (histoName == inputHistoName)
                    {  
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

void makePlotsForLimitsTool (string configFilePath){
  

  TH1F* h_clone;
  TH1F* h_mc;      // histograms for BGs
  TH1F* h_data;    // histogram to store data points 
  TH1F* h_signal;   // histogram for signal 
  TH1F* h_allBG;   // histogram with the sum of all BGs

  h_data = (TH1F*)listHistos->At(0);
  cout<<"Bins "<<h_data->GetXaxis()->GetNbins()<<endl;
  h_data = SetStyleh(h_data);

  h_signal = (TH1F*)listHistos->At(inputRootFile.size() - 1);
  h_signal = SetStyleh(h_signal);

  int Nbins = h_data->GetXaxis()->GetNbins();
  
  h_allBG = new TH1F("h_allBG","allBG", Nbins, h_data->GetXaxis()->GetXmin(),h_data->GetXaxis()->GetXmax());

  vector<int> histoId;
  int hID = 0;
  for (int bg=0; bg <= (Nbins); bg++)
    {
      h_allBG->SetBinContent(bg,0.0);
    }
 
  for (int j=1;j< (inputRootFile.size() - 1);j++)
    { 
      if (use_log == "True" && use_data == "False" ) gPad->SetLogy();
      
      h_mc = (TH1F*)listHistos->At(j);
      h_allBG->Add(h_mc); 
      if (inputDataDrivenBG.at(j) == "True")
       {
         histoId.push_back(j);
         
       }      
    }  // 
    
   h_allBG = SetStyleh(h_allBG);
   
   h_allBG->SetName("background");
   
   TH1F *h_background = new TH1F("allbkg","allbkg", Nbins, h_data->GetXaxis()->GetXmin(),h_data->GetXaxis()->GetXmax());
   TH1F* h_signal_new = new TH1F("sig","sig", Nbins, h_data->GetXaxis()->GetXmin(),h_data->GetXaxis()->GetXmax());
   TH1F* h_data_new   = new TH1F("data_obs","data_obs", Nbins, h_data->GetXaxis()->GetXmin(),h_data->GetXaxis()->GetXmax()); 

   hID = histoId.size();
   cout<<"histoId.size() "<<histoId.size()<<endl;
   TFile* f = new TFile ("InputHistosForLimitsVBF.root", "RECREATE");
   
   for (int h = 0; h < hID; h++)
     { 
       int id = histoId.at(h);
       h_mc = (TH1F*)listHistos->At(id);
       cout <<"Nbins "<<Nbins <<endl;
       for (int bg=0; bg <= (Nbins); bg++)
	 {
	   TH1F *h_background_Up = new TH1F("background_Up","background_Up", Nbins, h_data->GetXaxis()->GetXmin(),h_data->GetXaxis()->GetXmax());
	   TH1F *h_background_Down = new TH1F("background_Down","background_Down", Nbins, h_data->GetXaxis()->GetXmin(),h_data->GetXaxis()->GetXmax());
           h_background->SetBinContent(bg, h_allBG->GetBinContent(bg));         
           h_signal_new->SetBinContent(bg, h_signal->GetBinContent(bg)); 
           h_data_new->SetBinContent(bg, h_data->GetBinContent(bg));       
	   for (int bg_up=0; bg_up <= (Nbins); bg_up++)
             {
	       if (bg == bg_up)
		 { 
		   float bin_up   = h_allBG->GetBinContent(bg_up) + h_mc->GetBinError(bg_up);
                   float bin_down = h_allBG->GetBinContent(bg_up) - h_mc->GetBinError(bg_up);
		   h_background_Up->SetBinContent(bg_up, bin_up);
                   h_background_Down->SetBinContent(bg_up, bin_down);  
		 } else {
		   h_background_Up->SetBinContent(bg_up, h_allBG->GetBinContent(bg_up));
	           h_background_Down->SetBinContent(bg_up, h_allBG->GetBinContent(bg_up));
                 }
             }
           std::stringstream newHistoName;
           newHistoName <<"allbkg_bgk"<<h<<"bin"<<bg<<"Up";    
           h_background_Up->SetName(newHistoName.str().c_str());
	   h_background_Up->Write();
           newHistoName.str("");
           newHistoName <<"allbkg_bgk"<<h<<"bin"<<bg<<"Down";
           h_background_Down->SetName(newHistoName.str().c_str());
           h_background_Down->Write();
	 }
         if (h==(hID - 1))
           {
             h_background->Write();
             h_signal_new->Write();
             h_data_new->Write();
           }
     }
   f->Clone(); 
} // close  makePlotsForLimitsTool function


