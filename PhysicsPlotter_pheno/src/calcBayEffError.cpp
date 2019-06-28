//*******************************************************************************
//	Filename:	calcBayEffError.cpp
//	Authors:	Andres Florez 
//	Last update:	Aug. 8, 2012
//
//	This file contains a function to calculate the Bayesian efficiency error.
//*******************************************************************************

float calcBayEffError (float numerator, float denominator)
{
  float effError;
  
  TH1F* theNumHisto = new TH1F ("theNumHisto", "theNumHisto", 1, 0, 1);
  theNumHisto->SetBinContent (1, numerator);
  theNumHisto->Sumw2();
  
  TH1F* theDenHisto = new TH1F ("theDenHisto", "theDenHisto", 1, 0, 1);
  theDenHisto->SetBinContent (1, denominator);
  theDenHisto->Sumw2();
  
  TGraphAsymmErrors* bayesEff = new TGraphAsymmErrors();
  bayesEff->BayesDivide (theNumHisto, theDenHisto, "b");
  
  if (bayesEff->GetErrorYhigh(0) > bayesEff->GetErrorYlow(0))
    {
      effError = bayesEff->GetErrorYhigh(0);
    }
  else
    {
      effError = bayesEff->GetErrorYlow(0);
    }
  
  delete theNumHisto;
  delete theDenHisto;
  delete bayesEff;
  
  return effError;
}
