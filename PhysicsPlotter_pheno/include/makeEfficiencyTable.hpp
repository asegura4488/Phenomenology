//*******************************************************************************
//	Filename:	makeEfficiencyTables.hpp
//	Authors:	Andres Florez 
//	Last update:	Aug. 23, 2012
//
//	This is the header file for the source file makeEfficiencyTables.cpp.
//*******************************************************************************

// Header guard
#ifndef MAKEEFFICIENCYTABLE_H
#define MAKEEFFICIENCYTABLE_H

//=======================================================
// Contents of header file
//=======================================================

// Load C++ libraries
#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <math.h>

// Load ROOT library
#include <TError.h>

// Load libraries (automatically included by ROOT)
#include <TFile.h>
#include <TClass.h>

using namespace std;

// Declare global variables
vector<float>  inputNumerator;
vector<float>  inputDenominator;
vector<string> inputProcess;
vector<string> inputTheProcess;
vector<string> inputCutName; 
vector<float>  inputEvents;
vector<float>  inputCrossSection;
vector<float>  inputScaleFactor;
vector<float>  inputScaleFactorError;
vector<float>  inputLumi;
vector<float>  inputSkimming;
string inputRootFilePath;			// Path to ROOT file storing cut information
string inputType;				// Type of info from input file
string inputString;				// Content of info from input file
string logFilename;				// Name of log file
string theFile;
// Declare class
class LogFile
{
  ofstream outFile;				// Output stream for log file
  
public:
  LogFile (string filePath);
  void writeCumEff (int numHistos, float cumEff);
  void close();
};

// Declare functions
void readConfigLog (string configFilePath );
//=======================================================

#endif						// MAKEEFFICIENCYTABLE_H

