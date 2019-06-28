//*******************************************************************************
//	Filename:	makeEfficiencyTables.cpp
//	Authors:	Andres Florez 
//	Last update:	Sep. 21, 2012
//*******************************************************************************

//--------------------------------------------------------------------------------------------------------------
// Load supplied header file and source codes
#include "include/makeEfficiencyTable.hpp"
#include "src/calcBayEffError.cpp"
#include "src/calcBinEffError.cpp"
#include "src/LogFile.cpp"
//--------------------------------------------------------------------------------------------------------------

void makeEfficiencyTable (string configFileList, string configFilePath)
{
  // read Config file
  readConfigLog (configFilePath); 
}


void readConfigLog (string configFilePath)
{
  // Input stream for data file
  ifstream inFile;

  string configPlotList=configFilePath; 
 
  // Open input file
  inFile.open (configPlotList.c_str(), ios::in);
  
  // If can't open input file, exit the code
  if (!inFile)
    {
      cerr << "ERROR: Can't open input file: " << configPlotList << endl;
      exit (1);
    }
 
  int counter = 0; 
  int counterProcesses = 0;

  //===================================== 
  // Read the list of log files
  //=====================================

  while (inFile >> theFile)
    {
      // ROOT configuration file
      ifstream inputFile;
      string logOutputFile = "LogFiles/"+theFile;
      counter++;
      inputFile.open(logOutputFile.c_str(), ios::in);

      size_t found = theFile.find("output_NormalizePlotter_");
      string theCut = "output_NormalizePlotter_";

      theFile.erase(found, theCut.size());

      size_t found2 = theFile.find(".log");
      string theCut2 = ".log";

      theFile.erase(found2, theCut2.size()); 
      
      inputCutName.push_back(theFile);

      //============================================= 
      // Read the infirmation from each log file
      // and store it in vectors
      //============================================

      while (inputFile >> inputString >> inputType)
        {
          if (inputString == "EventsAnalysed")
            {
              counterProcesses++;
	      inputDenominator.push_back (atof(inputType.c_str()) );
            }
          else if (inputString == "EventsPassing")
	    inputNumerator.push_back (atof(inputType.c_str()) );
          else if (inputString == "Process")
            inputProcess.push_back (inputType.c_str());
          else if (inputString == "Events")
            inputEvents.push_back(atof(inputType.c_str()));
          else if (inputString == "CrossSection")
            inputCrossSection.push_back(atof(inputType.c_str()));
          else if (inputString == "luminosity")
            inputLumi.push_back(atof(inputType.c_str()));
          else if (inputString == "SkimmingEfficiency")
            inputSkimming.push_back(atof(inputType.c_str()));
          else if (inputString == "ScaleFactor")
            inputScaleFactor.push_back(atof(inputType.c_str()));
          else if (inputString == "ScaleFactorError")
            inputScaleFactorError.push_back(atof(inputType.c_str()));
        }
    }
  inFile.close();

  ifstream inputTheProcesses;
  string logProcesses = "config/LogProcesses.in";
  inputTheProcesses.open(logProcesses.c_str(), ios::in);
  string inputTypeProcess;

  while (inputTheProcesses >> inputTypeProcess)
    {
      inputTheProcess.push_back(inputTypeProcess.c_str());
    }

  cout <<"\\begin{table}[htbp]"<<endl;
  cout <<"\\tiny"<<endl;
  cout <<"\\centering"<<endl;
  cout <<"\\begin{tabular}{lcccc}\\hline\\hline"<<endl;

  for (int i = 0; i < inputTheProcess.size(); i++)
   {
       if(i == 0 ) cout << setw(10)<<" Cut &";
       if (i < (inputTheProcess.size() -1))
         cout << setw(12) <<"$"<<inputTheProcess.at(i) <<"$"<<" & ";
       else 
         cout << setw(12) <<"$"<<inputTheProcess.at(i) <<"$"<<" \\\\ ";
   } 
  cout <<" " <<endl;

  float getElement = (float)counterProcesses/float(counter);

  //cout<<inputCutName.at(0)<<" & ";
  cout<<"0"<<" & ";

  int getIntElement = ceil(getElement);
  int printlast = 0;
 
  for (int j = 0; j < getIntElement; j++)
   {
     for (int k = 0; k < inputTheProcess.size(); k++)
       {
         if ( inputTheProcess.at(k) != inputProcess.at(j)) continue;
         float efficiency = 1000.;
         float efferror = -1000.;
         if (inputDenominator.at(j) != 0)
          {
            efficiency = 100*(inputNumerator.at(j)/inputDenominator.at(j));
            efferror = 100*calcBinEffError(inputNumerator.at(j), inputDenominator.at(j));
          }

         if ((efficiency+efferror >= 1.) || ((efficiency-efferror) <= 0.))
          {
             efferror = 100*calcBayEffError(inputNumerator.at(j), inputDenominator.at(j)); 
          }
         if ((inputNumerator.at(j) == 0.) && (inputDenominator.at(j) > 0.))
          {
            efferror = 100*calcBinEffError(1.0, inputDenominator.at(j));
          }
          printlast++; 
          if ( printlast < inputTheProcess.size() )
            cout << setw(5)<<setprecision(3) <<"$"<<efficiency <<" \\pm " << setw(5)<<setprecision(3) << efferror  <<"$ &";
          else 
            cout << setw(5)<<setprecision(3) <<"$"<<efficiency <<" \\pm " << setw(5)<<setprecision(3) << efferror  <<"$ \\\\";
       }
   }

   cout <<" " <<endl;

   int nextline = 0;
   int theCutLabel = 1; 
   for (int i = 0; i < (inputNumerator.size() - getIntElement); i++)
   {
     for (int p = 0; p < inputTheProcess.size(); p++)
       {
         if ( inputTheProcess.at(p) != inputProcess.at(i)) continue;
         float numerator = inputNumerator.at(i+getIntElement);
         float denominator = inputNumerator.at(i);
         
         float efficiency = 1000;
         float efferror2 = -1000;
         if (denominator != 0.)
          {
            efficiency = 100*(numerator/denominator);
            efferror2 = 100*(calcBinEffError(numerator, denominator));
          }
         if ((efficiency+efferror2 >= 1.) || ((efficiency-efferror2) <= 0.))
          {
             efferror2 = 100*calcBayEffError(numerator, denominator);
          }

         if (numerator == 0.0)
          {
             efferror2 = 100*calcBinEffError(1.0, inputDenominator.at(i));
          }

         if (nextline == 0)
          {
            cout<<theCutLabel<<" & "; 
            cout<<inputCutName.at(theCutLabel)<<" & ";
          }
          
         nextline++;
         if ( nextline < inputTheProcess.size() )
            cout << setw(5) << setprecision(3) <<"$"<<efficiency <<" \\pm " << setw(5)<<setprecision(3) << efferror2 <<"$ & ";
         else 
            cout << setw(5) << setprecision(3) <<"$"<<efficiency <<" \\pm " << setw(5)<<setprecision(3) << efferror2 <<"$ \\\\ ";
         if ( i !=0 && ((nextline%inputTheProcess.size()) == 0) )
          {
            cout <<" " <<endl;
            nextline = 0;
            theCutLabel++;
          }

       }
    }

    cout<<"\\hline"<<endl;

    for (int e = (inputProcess.size() - getIntElement); e < inputProcess.size(); e++)
     {

       for (int g = 0; g < inputTheProcess.size(); g++)
        {
          if (inputTheProcess.at(g) == inputProcess.at(e))
           {
             float EventsError = 999.;
            
             if(inputTheProcess.at(g) == "data") 
               EventsError = sqrt(inputEvents.at(e));
             else {
                  float cumErrError = 999.; 
                  if (inputNumerator.at(e) != 0.0)
                   {
                      cumErrError = calcBinEffError( inputNumerator.at(e), inputDenominator.at(e));
                      float cumEfficiency = inputNumerator.at(e)/inputDenominator.at(e);
                      if ((cumEfficiency+cumErrError) >= 1.0 || (cumEfficiency-cumErrError) <= 0.0)
                       {
                         cumErrError = calcBayEffError( inputNumerator.at(e), inputDenominator.at(e));
                       }
                          EventsError = cumErrError * inputLumi.at(e) * inputCrossSection.at(e) * inputSkimming.at(e) * inputScaleFactor.at(e); 
                          EventsError = sqrt(EventsError*EventsError + inputScaleFactorError.at(e)*inputScaleFactorError.at(e));
                   } 
              }              
             if (e == (inputProcess.size() - getIntElement) ) cout <<"Events "<< setw(10) << " & "; 
             cout << setw(10) << setprecision(5) << "$" << inputEvents.at(e) << " \\pm  "<<EventsError<<"$ &"; 
           }
        }

     }

    cout <<" \\\\ "<<endl;

    cout<<"\\hline \\hline"<<endl;
    cout<<"\\end{tabular}"<<endl;
    cout<<"\\end{table}"<<endl;

}

