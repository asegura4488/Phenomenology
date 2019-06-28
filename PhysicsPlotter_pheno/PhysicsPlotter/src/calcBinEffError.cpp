//*******************************************************************************
//	Filename:	calcBinEffError.cpp
//	Authors:	Andres Florez 
//	Last update:	Sep. 21, 2013
//
//	This file contains a function to calculate the Bayesian efficiency error.
//*******************************************************************************

float calcBinEffError (float numerator, float denominator)
{
  if ( denominator > 0 )
   {
     float efficiency = numerator/denominator;
     float efferror = sqrt(efficiency*(1.0-efficiency)/denominator);
   }
  return efferror;
}
