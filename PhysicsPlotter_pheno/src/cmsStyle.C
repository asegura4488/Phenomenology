//
// CMS Style, based on a style file from BaBar
//

#include <iostream>

#include "../include/cmsStyle.h"

#include "TROOT.h"

void SetcmsStyle ()
{
  static TStyle* cmsStyle = 0;
  std::cout << "\nApplying CMS style settings...\n" << std::endl ;
  if ( cmsStyle==0 ) cmsStyle = cmsStyle();
  gROOT->SetStyle("CMS");
  gROOT->ForceStyle();
}

TStyle* cmsStyle() 
{
  TStyle *cmsStyle = new TStyle("CMS","cms style");

  // use plain black on white colors
  Int_t icol=0; // WHITE
  cmsStyle->SetFrameBorderMode(icol);
  cmsStyle->SetFrameFillColor(icol);
  cmsStyle->SetCanvasBorderMode(icol);
  cmsStyle->SetCanvasColor(icol);
  cmsStyle->SetPadBorderMode(icol);
  cmsStyle->SetPadColor(icol);
  cmsStyle->SetStatColor(icol);
  //cmsStyle->SetFillColor(icol); // don't use: white fill color for *all* objects

  // set the paper & margin sizes
  cmsStyle->SetPaperSize(20,26);

  // set margin sizes
  cmsStyle->SetPadTopMargin(0.05);
  cmsStyle->SetPadRightMargin(0.05);
  cmsStyle->SetPadBottomMargin(0.16);
  cmsStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  cmsStyle->SetTitleXOffset(1.4);
  cmsStyle->SetTitleYOffset(1.4);

  // use large fonts
  //Int_t font=72; // Helvetica italics
  Int_t font=42; // Helvetica
  Double_t tsize=0.03;
  cmsStyle->SetTextFont(font);

  cmsStyle->SetTextSize(tsize);
  cmsStyle->SetLabelFont(font,"x");
  cmsStyle->SetTitleFont(font,"x");
  cmsStyle->SetLabelFont(font,"y");
  cmsStyle->SetTitleFont(font,"y");
  cmsStyle->SetLabelFont(font,"z");
  cmsStyle->SetTitleFont(font,"z");
  
  cmsStyle->SetLabelSize(tsize,"x");
  cmsStyle->SetTitleSize(tsize,"x");
  cmsStyle->SetLabelSize(tsize,"y");
  cmsStyle->SetTitleSize(tsize,"y");
  cmsStyle->SetLabelSize(tsize,"z");
  cmsStyle->SetTitleSize(tsize,"z");

  // use bold lines and markers
/*
  cmsStyle->SetMarkerStyle(20);
  cmsStyle->SetMarkerSize(1.2);
  cmsStyle->SetHistLineWidth(2.);
  cmsStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
*/
  // get rid of X error bars 
  //cmsStyle->SetErrorX(0.001);
  // get rid of error bar caps
  cmsStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  cmsStyle->SetOptTitle(0);
  //cmsStyle->SetOptStat(1111);
  cmsStyle->SetOptStat(0);
  //cmsStyle->SetOptFit(1111);
  cmsStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  cmsStyle->SetPadTickX(1);
  cmsStyle->SetPadTickY(1);

  return cmsStyle;

}

