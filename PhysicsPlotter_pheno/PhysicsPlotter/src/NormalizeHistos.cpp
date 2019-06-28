//*******************************************************************************
//	Filename:	NormalizeHistos.cpp
//	Authors:	Andres Florez 
//	Last update:	Aug. 23, 2012
//
//	This function is called by the shell script runPlotter.sh. It opens ROOT files listed in
//	the configuration file NormalizeHistos_${Directory}_SAFE.in and extracts the 1D
//	histograms. The probability and normalized histograms are output to ROOT files
//	and a log file is generated.
//*******************************************************************************

//--------------------------------------------------------------------------------------------------------------
// Load supplied header file and source codes
#include "include/NormalizeHistos.hpp"
#include "src/calcBayEffError.cpp"
#include "src/LogFile.cpp"
//--------------------------------------------------------------------------------------------------------------

// Main function to normalize histograms
void NormalizeHistos (string inRootFilePath, string configFileList, string configFilePath, string dirName)
{

  // use CMS style
  gROOT->LoadMacro("src/cmsStyle.C");
  SetcmsStyle();
 
  // Store global variables
  inputRootFilePath = inRootFilePath; 
  readConfigFile (configFileList, configFilePath);
  
  // Create log file
  string logFilename_dir = logFilename+"_"+dirName+".log";
  string logFilePath = "LogFiles/" + logFilename_dir;
  
  LogFile* logTextFile = new LogFile (logFilePath);
  
  getHistoInfo (logTextFile);
  delete logTextFile;		// Finished writing log file
  
}

//=======================================================
// Store information from configuration file (e.g. root files, scale factors, cross-sections).
// The config. file lists each piece of information (its type and content) on a new line.
//=======================================================
void readConfigFile (string configFileList, string configFilePath)
{
  ifstream inFile;		// Input stream for data file
  
  // Open input file
  inFile.open (configFileList.c_str(), ios::in);
 
  // If can't open input file, exit the code
  if (!inFile)
    {
      cerr << "ERROR: Can't open input file: " << configFilePath << endl;
      exit (1);
    }
  
  // Get all relevant information from the input file
  while (inFile >> inputType)
    {

      string configfile=configFilePath + "/" + inputType;
      
      TEnv *params = new TEnv ("config_file");
      params->ReadFile (configfile.c_str(), kEnvChange);

      // Get Root file name
      string rootFileName = params->GetValue ("rootfile", "thefile");
      inputRootFile.push_back (rootFileName.c_str() );

      // Get directory where the histograms are stored
      string directory = params->GetValue ("directory", "event");
      inputDirectory.push_back ( directory.c_str() );

      // Get process name: data, ttbar, wjets etc.
      string process = params->GetValue ("process", "theprocess");
      inputProcess.push_back( process.c_str() );

      // Skimming Efficiency 
      float skimmingEff = params->GetValue("skimmingEff", 1.0);
      inputSkimmingEff.push_back (skimmingEff);

      // Skimming Efficiency error
      float skimmingEffError = params->GetValue("skimmingEffError", 0.0);
      inputSkimmingEffError.push_back ( skimmingEffError ); 

      // Skimming Efficiency 
      float scaleFactor = params->GetValue("scaleFactor", 1.0);
      inputScaleFactor.push_back (scaleFactor);

      // Skimming Efficiency error
      float scaleFactorError = params->GetValue("scaleFactorError", 0.0);
      inputScaleFactorError.push_back ( scaleFactorError );

      // Cross section for each process       
      float xsection = params->GetValue("xsection", 1.0); 
      xSection.push_back ( xsection );

      // luminosity 
      luminosity  = params->GetValue("luminosity", 1.0 );
      ostringstream  lumi;
      lumi << luminosity;
      luminosityString = lumi.str();
  
      // Is collision data? no normalization is applied, but ATLAS style is used for the plots
      string is_data = params->GetValue("isdata", "True"); 
      inputIsData.push_back ( is_data.c_str() );

      // Output path
      string outputfilepath = params->GetValue("OutputFilePath" , "output"); 
      outputRootFilePath.push_back( outputfilepath.c_str() ); 
 
      // Output root file name
      string normHistosOutputName = params->GetValue("NormOutputFileName", "norm.root");
      normHistosOutputRootFile.push_back( normHistosOutputName.c_str() ); 
     
      probHistosOutputRootFile = params->GetValue("NormToOneOutputFileName", "normToOne.root" );
      normOneHistosOutputRootFile.push_back( probHistosOutputRootFile.c_str() );
  
      // log file name 
      logFilename  = params->GetValue("logFilename", "outputLog_stackPlots.log"); 
      
    }
  
  inFile.close();
}

//=======================================================
// 
//=======================================================
void getHistoInfo (LogFile* logFile)
{
  // Initialize variables
  
  // Loop over input ROOT files
  for (int j = 0; j < inputRootFile.size(); j++)
    {
      int numMax = 0;			// Total entries of maxHistoContent
      
      // Write info to log file
      logFile->writeIdentifier (j, inputRootFile.at(j), inputDirectory.at(j), inputProcess.at(j));
 
      // Access input ROOT file (can access over secure shell)
      TFile* rootTFile = (TFile*) TFile::Open (inputRootFile.at(j).c_str());
 
      // Change directory accessed by ROOT to cutflow path
      rootTFile->cd(inputDirectory.at(j).c_str());

      TH1F* n_events = Events_0;

      // Skimming efficiency
      logFile->writeSkimEff(inputSkimmingEff.at(j), inputSkimmingEffError.at(j));
      calcCumlEff (j, n_events, logFile);

      // Create iterator for ROOT objects
      TDirectory *currentSourceDir = gDirectory;
      TIter nextKey (currentSourceDir->GetListOfKeys());
      TKey* key;

      // Loop over keys (ROOT objects) from ROOT file
      while ( (key = (TKey*) nextKey()) )			// Need double parentheses
	{
	  TObject* obj = key->ReadObj();
	  TH1* histoObj = (TH1*) obj;
	  // Get only 1D histograms from ROOT file and store in a TList
          if ( (histoObj->GetYaxis()->GetNbins() == 1) )   
	    {
              double maxContent = 0;			// Max histogram content
	      TH1F* h1;     // Histogram normalized to cross section & luminosity
	      TH1F* h2;     // Histogram normalized to 1 (Must declare separately)
	    
              // Store object as a 1D histogram
              string histoName = histoObj->GetName(); 
                
	      // Get name of histogram and create (name+cut) & (name+probability)
              string histoNameLumi = histoName + "_lumi"; 
	      string histoNameProb = histoName + "_pro";
	      
	      // Store default histogram names
	      histoNameDefault.push_back (histoName);
	      
	      // Write info to log file
	      //if (j == 0) 
                //logFile->writeHistoName (histoName);
	      
	      // Create histograms based on the 1D histogram
              h1 = (TH1F*)histoObj;
              h1->SetName( histoNameLumi.c_str() );
              normalizeHistosToLumi( h1, j, logFile );
             
              h2 = (TH1F*)histoObj;
              h2->SetName( histoNameProb.c_str() );
              normalizeToOne ( h2, j );
 
	      // Store histograms and increment total count
	      numMax++;

	    } // Close if()
	} // Close while loop
      
      // Calculate cumulative efficiency
      cout<<"ROOT FILE "<<j <<"  "<< inputRootFile.at(j)<<endl;  
    
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
    } // Close for loop over root files
} // Close function

//=======================================================
// Calculate skimming efficiency using number of events before & after cuts
//=======================================================
void calcCumlEff (int i, TH1F* N_events, LogFile* logFile1)
{
  float cumlEfficiency = 0;			        // Skimming efficiency
  
  // If there were zero events analysed, exit code
  if (N_events->GetBinContent(1) == 0)
    {
      cerr << "ERROR: 'Events' histogram contains zero entries in bin 1:"
	   << " Zero events were analyzed." << endl;
      exit (1);
    }
  // Else extract number of events from histogram
  else
    {
      float numEventsBeforeCuts = (float)N_events->GetBinContent(1);
      float numEventsAfterCuts =  (float)N_events->GetBinContent(2);
      
      // Calculate cuml efficiency and warn if zero
      cumlEfficiency = (numEventsAfterCuts / numEventsBeforeCuts);
      if (cumlEfficiency == 0)
	cerr << "WARNING: Cuml efficiency is zero." << endl;
      
      // Write info to log file
      logFile1->writeCumEff (N_events->GetBinContent(1), N_events->GetBinContent(2),
			      cumlEfficiency);
      
      // Store event numbers in global variables
      numEventsAnalysed.push_back(numEventsBeforeCuts);
      numEventsPassing.push_back(numEventsAfterCuts);
    }
}

//=======================================================
// Normalize historgrams to luminosity
//=======================================================
void normalizeHistosToLumi (TH1* h, int k, LogFile* logFile)
{

   h->Sumw2();

   if ( inputIsData.at(k) == "False" )
     {
       // Calculate scale factor for number of events
       float normFactor = (1.0 / numEventsAnalysed.at(k)) * xSection.at (k) * inputSkimmingEff.at(k) * luminosity;

       //cout<<"normFactor "<< normFactor <<" numEventsAnalysed.at(k) "<<numEventsAnalysed.at(k)<< " xSection.at (k)  " << xSection.at (k) <<" luminosity "<<luminosity <<" inputSkimmingEff.at(k) "<<inputSkimmingEff.at(k)<<endl;

       h->Scale (normFactor);

       TH1F* h_sf = new TH1F ("h_sf", "h_sf", h->GetXaxis()->GetNbins(), 
			      h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
       
       for (int ibin=0; ibin <= (h->GetXaxis()->GetNbins() + 1); ibin++)
	 {
	   h_sf->SetBinContent (ibin, inputScaleFactor.at(k));
	   h_sf->SetBinError (ibin, inputScaleFactorError.at(k));
	 }
       
       h_sf->Sumw2();
       h->Multiply (h, h_sf, 1, 1, "");
       
     } 

   string YtitleString = "N / " + luminosityString + " fb^{-1}";
   h->GetYaxis()->SetTitle (YtitleString.c_str());
   
   string histoname = h->GetName();
   if (histoname == "Met_0_lumi")
    {
       float NormNevents = h->Integral();
       logFile->writeNumberOfEvents(NormNevents);
       logFile->writeFactors(xSection.at(k), luminosity, inputScaleFactor.at(k), inputScaleFactorError.at(k));
    } 
   TFile* hfile;
   string outputFilePath = outputRootFilePath.at(k) + "/" + normHistosOutputRootFile.at(k);
   hfile = (TFile*) gROOT->FindObject (outputFilePath.c_str());
   
   if (hfile)
     {
       hfile->Close();
     }
   
   hfile = new TFile (outputFilePath.c_str(), "UPDATE");
   h->Write();
   hfile->Close();
   
}

//=======================================================
// Normalize histograms to one
//=======================================================

void normalizeToOne ( TH1* hp, int k) 
{

   hp->Sumw2();

   if (hp->Integral (0, (hp->GetXaxis()->GetNbins() + 1)) > 0)
     {
        hp->Scale (1.0 / hp->Integral (0, (hp->GetXaxis()->GetNbins() + 1)));
     }

   TFile* hfile;
   string outputFilePath = outputRootFilePath.at(k) + "/" + normOneHistosOutputRootFile.at(k);
   hfile = (TFile*) gROOT->FindObject (outputFilePath.c_str());

   if (hfile)
     {
       hfile->Close();
     }

   hfile = new TFile (outputFilePath.c_str(), "UPDATE");
   hp->Write();
   hfile->Close();     

} 

