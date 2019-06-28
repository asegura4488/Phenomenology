//*******************************************************************************
//	Filename:	LimitsRootFileMaker.hpp
//	Authors:	Andres Florez 
//	Last update:	Nov. 21, 2012
//
//	This is the header file for the source file LimitsRootFileMaker.cpp.
//*******************************************************************************

// Header guard
#ifndef STACKPLOTTER_H
#define STACKPLOTTER_H

//=======================================================
// Contents of header file
//=======================================================

// Load C++ libraries
#include <iostream>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <string.h>

// Load ROOT library
#include <TError.h>

// Load libraries (automatically included by ROOT)
#include <TFile.h>
#include <TH1.h>
#include <TFrame.h>
#include <TPad.h>
#include <TCollection.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TClass.h>
#include <fstream>
using namespace std;


// Set default style options using global variables
gROOT->Reset();						// Reset all global variables
gStyle->SetOptStat(0);					// No histogram statistics box
gStyle->SetOptTitle(0);					// No histogram title

// Stop printing of error messages. Used to suppress the warning for function
// Sumw2() that says "sums of squares of weights structure already created"
gErrorIgnoreLevel = kError;


// Declare global variables
vector<string>  inputRootFile;			// ROOT file storing cut information
vector<string>  inputDirectory;			// Path containing histograms
vector<string>  inputProcess;			// Name of cut
vector<string>  theHistNameStrings;		// Default name of histogram
vector<string>  inputApplyScaleFactor;          // Apply scale factor?
vector<float>   inputScaleFactor;		// Scale factor
vector<float>   inputScaleFactorError;		// Scale factor error
vector<string>  inputRebinHisto;                // Name of historgams to be rebinned
vector<string>  inputDataDrivenBG;
vector<int>     inputRebinFactor;               // Rebin factor 
vector<int>     inputSetLineStyle;              // Set histo line style
vector<int>     inputSetLineWidth;              // Set histo line wdth
vector<string>  inputLegend;                    // Histo legend
vector<int>     inputColor;                     // Histo fill and line color
vector<string>  inputIsSignal;                  // Is signal or BG?
string          inputHistoName;                 // Name of histo to set style for
vector<double>  inputXrangeMin;                 // Minimum range x axis
vector<double>  inputXrangeMax;                 // Maximum range x axis
vector<double>  inputYrangeMin;                 // Minimum range y axis
vector<double>  inputYrangeMax;                 // Maximum range y axis
vector<string>  inputXtitle;                    // Set x-axis title
vector<string>  inputYtitle;                    // set y-axis title
vector<int>     inputRebin;                     // rebinning factor
string          logFilename;                    // Name of log file


string inputRootFilePath;			// Path to ROOT file storing cut information
string inputType;				// Type of info from input file
string inputStyle;                              // Type of info from style input file
string inputString;				// Content of info from input file
string stackedHistosOutputRootFile;             // Output root file
string use_data;                                // Make plots with collision data
string use_log;

TList* listHistos;		        	// Histograms normalized to cross-section & luminosity

int numHistos;                                  // Number of histograms in input ROOT file
int nHistList;                                 // Total number of histograms

// Declare class

// Declare class
class LogFile
{
  ofstream outFile;                             // Output stream for log file

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
void readConfigFile (string configFileList, string configFilePath);
void readConfigPlotStyle (string configFilePath);
TH1F*    SetStyleh(TH1F* h_j);
TH1F* Rebinh(TH1F* h_j);
void makePlotsForLimitsTool(string);
//=======================================================

#endif						// STACKPLOTTER_H


