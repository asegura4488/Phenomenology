//*******************************************************************************
//	Filename:	LogFile.cpp
//	Last update:	Aug. 23, 2012
//
//	This class creates and adds information to a log file.
//
//	To debug, type the following:
//
//		root -l -L LogFile.cpp+
//*******************************************************************************

#include <iostream>			// Object cerr, Class ios
#include <fstream>			// Class ofstream
#include <stdlib.h>		        // Function exit
using namespace std;			// Class ofstream, Class string

//=======================================================
// Constructor: Create and open output log file
//=======================================================
LogFile::LogFile (string filePath)
{
  // Create an output log file that will contain the cut flow efficiency table
  outFile.open (filePath.c_str(), ios::out);
  
  // If output log file cannot be opened, exit the code
  if (!outFile)
    {
      cerr << "Can't open output file: " << filePath << endl;
      exit (1);
    }
  else
    {
      outFile << endl;
    }
}

//=======================================================
// Write name of ROOT file, path containing histograms, and name of cut
//=======================================================
void LogFile::writeIdentifier (int num, string file, string directory, string process)
{
  outFile << endl << "ROOTFile " << file << endl;
  outFile << "Directory " << directory << endl;
  outFile << "Process " <<process << endl;
}

//=======================================================
// Write number of events before & after cuts and skimming efficiency
//=======================================================
void LogFile::writeSkimEff (float skimEff, float skimEffError)
{
  outFile << "SkimmingEfficiency " <<skimEff<< endl;
  outFile << "SkimmingEfficiencyError " <<skimEffError<< endl;
}

//=======================================================
// Write name of histogram
//=======================================================
void LogFile::writeHistoName (string histogramName)
{
  outFile << histogramName << endl;
}

//=======================================================
// Write number of histograms and cumulative efficiency
//=======================================================

void LogFile::writeCumEff (float eventsBeforeCuts, float eventsAfterCuts, float cumlEff)
{
  outFile << "NumberOfHistorgams " << numHistos << endl;
  outFile << "EventsAnalysed " << eventsBeforeCuts << endl;
  outFile << "EventsPassing " << eventsAfterCuts << endl;
  outFile << "CumulativeEfficiency " << cumlEff << endl;
}
//=======================================================
// Write number of events 
//=======================================================
void LogFile::writeNumberOfEvents (float events)
{
  outFile << "Events " <<events << endl;
}

//=======================================================
// Write number of histograms and "K"
//=======================================================
void LogFile::writeHistoNameK (string histogramName, int kVal)
{
  outFile << endl << "HistoName :\t" << histogramName << " K :\t" << kVal << endl;
}

//=======================================================
// Write histogram integral, events analysed, and cross-section
//=======================================================
void LogFile::writeFactors (float xSect, float lumi, float scaleFactor, float scaleFactorError)
{
  outFile << "CrossSection " << xSect << endl;
  outFile << "luminosity "   << lumi << endl;
  outFile << "ScaleFactor " << scaleFactor << endl;
  outFile << "ScaleFactorError " << scaleFactorError << endl;
}

//=======================================================
// Write scale factor and histogram integral scaled by it
//=======================================================
void LogFile::writeScaledIntegral (float scalingFactor, float integral)
{
  outFile << "Scale Factor :\t" << scalingFactor << endl;
  outFile << "Integral :\t" << integral << endl;
}

//=======================================================
// Close output log file
//=======================================================
void LogFile::close()
{
  outFile.close();
}

