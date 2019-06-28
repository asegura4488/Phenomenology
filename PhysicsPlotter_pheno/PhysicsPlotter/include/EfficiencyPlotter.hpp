//*******************************************************************************
//	Filename:	EfficiencyPlotter.hpp
//	Authors:	Andres Florez 
//	Last update:	Dec12, 2012
//
//	This is the header file for the source file EfficiencyPlotter.cpp.
//*******************************************************************************

// Header guard
#ifndef EFFICIENCYPLOTTER_H
#define EFFICIENCYPLOTTER_H

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
vector<string>  inputRootFileNum;		// ROOT file storing cut information
vector<string>  inputRootFileDen;               // ROOT file storing cut information
vector<string>  inputDirectory;			// Path containing histograms
vector<string>  inputProcess;			// Name of cut
vector<string>  inputIsSignal;
vector<string>  theHistNameStrings;		// Default name of histogram
vector<string>  inputLegend;                    // Histo legend
vector<int>     inputColor;                     // Histo fill and line color
vector<int>     inputMarkerStyle;               // Marker style 20, 21, 22..
vector<string>  inputHistoName;                 // Name of histo to set style for
vector<double>  inputXrangeMin;                 // Minimum range x axis
vector<double>  inputXrangeMax;                 // Maximum range x axis
vector<double>  inputYrangeMin;                 // Minimum range y axis
vector<double>  inputYrangeMax;                 // Maximum range y axis
vector<string>  inputXtitle;                    // Set x-axis title
vector<string>  inputYtitle;                    // set y-axis title
vector<int>     inputRebin;                     // rebinning factor
vector<int>     inputdifferentBinSize;
vector<double>  inputBinSize;
vector<int>     inputNbins;
vector<int>     inputNbinMax;

string inputType;				// Type of info from input file
string inputStyle;                              // Type of info from style input file
string inputString;				// Content of info from input file
string HistosOutputRootFile;                    // Output root file

TList* listHistosNum;		        	// Histograms normalized to cross-section & luminosity
TList* listHistosDen;
int numHistos;                                  // Number of histograms in input ROOT file
int nHistList;                                 // Total number of histograms

// Declare functions
void readConfigFile (string configFileList, string configFilePath);
void readConfigPlotStyle (string configFilePath);
TH1F* SetStyleh(TH1F* h_j);
TH1F* Rebinh(TH1F* h_j);
TH1F* normalizeToOne ( TH1F* hp);
TH1F* normalizeToNtrk (TH1F* hp, int nTrks);
void makePlots(string);
//=======================================================

#endif						// EFFICIENCYPLOTTER_H


