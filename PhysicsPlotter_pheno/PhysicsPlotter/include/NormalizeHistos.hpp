//*******************************************************************************
//	Filename:	NormalizeHistos.hpp
//	Authors:	Andres Florez 
//	Last update:	Aug. 23, 2012
//
//	This is the header file for the source file NormalizeHistos.cpp.
//*******************************************************************************

// Header guard
#ifndef NORMALIZEHISTOS_H
#define NORMALIZEHISTOS_H

//=======================================================
// Contents of header file
//=======================================================

// Load C++ libraries
#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <fstream>

// Load ROOT library
#include <TError.h>

// Load libraries (automatically included by ROOT)
#include <TFile.h>
#include <TH1.h>
#include <TCollection.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TGraphAsymmErrors.h>
#include <TClass.h>

using namespace std;


// Set default style options using global variables
gROOT->Reset();						// Reset all global variables
//gStyle->SetOptStat(0);					// No histogram statistics box
gStyle->SetOptTitle(0);					// No histogram title

// Stop printing of error messages. Used to suppress the warning for function
// Sumw2() that says "sums of squares of weights structure already created"
gErrorIgnoreLevel = kError;


// Declare global variables
vector<string> inputRootFile;			// ROOT file storing cut information
vector<string> inputDirectory;			// Path containing histograms
vector<string> inputProcess;			// Name of cut
vector<string> histoNameDefault;		// Default name of histogram
vector<string> inputIsData;

vector<float> inputSkimmingEff;			// Skimming Efficiency
vector<float> inputSkimmingEffError;		// Skimming Efficiency error
vector<float> inputScaleFactor;                 // Scale factor
vector<float> inputScaleFactorError;            // SScale factor error
vector<float> xSection;				// Cross section
vector<float> numEventsAnalysed;		// Number of events before cuts
vector<float> numEventsPassing;		        // Number of events after cuts
vector<float> cumulativeEfficiency;		// Cumulative efficiency

string inputRootFilePath;			// Path to ROOT file storing cut information
string inputType;				// Type of info from input file
string inputString;				// Content of info from input file
string luminosityString;			// Luminosity
string use_data;				// Process data? ("True" or "False")
string probHistosOutputRootFile;		// Output ROOT file
string logFilename;				// Name of log file

vector<string> outputRootFilePath;              // Path to output file
vector<string> normHistosOutputRootFile;        // Name of root file for normalized histos
vector<string> normOneHistosOutputRootFile;     // Name of root file for normalized to 1 histos

float luminosity;				// Luminosity (a constant)
float branchingRatio;				// Branching ratio
float branchingRatioError;			// Branching ratio error

int numHistos;					// Number of histograms in input ROOT file
int numHistos1;					// Number of histograms in listHistos1
int numHistos2;					// Number of histograms in listHistos2

// Declare class
class LogFile
{
  ofstream outFile;				// Output stream for log file
  
public:
  LogFile (string filePath);
  void writeIdentifier (int num, string file, string directory, string process);
  void writeSkimEff (float eventsBeforeCuts, float eventsAfterCuts, float skimEff);
  void writeHistoName (string histogramName);
  void writeCumEff (int numHistos, float cumEff);
  void writeHistoNameK (string histogramName, int kVal);
  void writeIntegral (float integral, float eAnalysed, float xSect);
  void writeScaledIntegral (float scaleF, float integral);
  void close();
};

// Declare functions
void readConfigFile (string configFileList, string configFilePath );
void getHistoInfo (LogFile* logTextFile);
void calcCumlEff (int i, TH1D* h, LogFile* logFile1);
void normalizeToOne (TH1* h, int j);
void normalizeHistosToLumi (TH1* h, int j, LogFile* logFile);
//=======================================================

#endif						// NORMALIZEHISTOS_H

