/**
 * @file MArEXStyle.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-11-13
 */

#include "TColor.h"
#include "TH1.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TStyle.h"

TStyle* MArEXStyle()
{
  printf("Using MArEX default plot style \n");

  TStyle* marexStyle = new TStyle("marexStyle", "MArEX Style");

  // Centre title
  marexStyle->SetTitleAlign(22);
  marexStyle->SetTitleX(.5);
  marexStyle->SetTitleY(.95);
  marexStyle->SetTitleBorderSize(0);

  // No info box
  marexStyle->SetOptStat(0);

  // set the background color to white
//   marexStyle->SetFillColor(10);
//   marexStyle->SetFrameFillColor(10);
//   marexStyle->SetCanvasColor(10);
//   marexStyle->SetPadColor(10);
  marexStyle->SetTitleFillColor(0);
//   marexStyle->SetStatColor(10);

  // set canvas options 
  marexStyle->SetCanvasDefW(1000); //600
  marexStyle->SetCanvasDefH(500); //500 
  marexStyle->SetCanvasColor(0); // canvas...
  marexStyle->SetCanvasBorderMode(0);
  marexStyle->SetCanvasBorderSize(0);   
  marexStyle->SetPadBorderMode(0); 
  marexStyle->SetPadBottomMargin(0.09); //0.09 //0.18 //margins...
  marexStyle->SetPadTopMargin(0.12); //0.12
  marexStyle->SetPadLeftMargin(0.10); //0.18
  marexStyle->SetPadRightMargin(0.09);
  marexStyle->SetPadGridX(1); // grids, tickmarks
  marexStyle->SetPadGridY(1);
  marexStyle->SetPadTickX(1);
  marexStyle->SetPadTickY(1);
  marexStyle->SetFrameBorderMode(0);
  marexStyle->SetPaperSize(20,24); // US letter size 

  // Set the default line color for a fit function to be red
  marexStyle->SetFuncColor(kRed);

  // Marker settings
  marexStyle->SetMarkerStyle(kFullCircle);

  // Legends
  marexStyle->SetLegendBorderSize(1);

  // Scientific notation on axes
  //  TGaxis::SetMaxDigits(3);

  // Axis titles
  marexStyle->SetTitleSize(.06, "xyz");
  marexStyle->SetTitleOffset(1.2, "xyz");
  // More space for y-axis to avoid clashing with big numbers
  marexStyle->SetTitleOffset(1.3, "y");
  // This applies the same settings to the overall plot title
  marexStyle->SetTitleSize(.06, "");
  marexStyle->SetTitleOffset(.9, "");
  // Axis labels (numbering)
  marexStyle->SetLabelSize(.06, "xyz");
  marexStyle->SetLabelOffset(.005, "xyz");

  // Prevent ROOT from occasionally automatically zero-suppressing
  marexStyle->SetHistMinimumZero();

  // Thicker lines
  marexStyle->SetHistLineWidth(2);  //   <- Does not work?
  marexStyle->SetFrameLineWidth(2);
  marexStyle->SetFuncWidth(2);

  // Set the number of tick marks to show
  marexStyle->SetNdivisions(506, "xyz");

  // Set the tick mark style
  marexStyle->SetPadTickX(1);
  marexStyle->SetPadTickY(1);

  // Fonts
  const int kMArEXFont = 42;
  marexStyle->SetStatFont(kMArEXFont);
  marexStyle->SetLabelFont(kMArEXFont, "xyz");
  marexStyle->SetTitleFont(kMArEXFont, "xyz");
  marexStyle->SetTitleFont(kMArEXFont, ""); // Apply same setting to plot titles
  marexStyle->SetTextFont(kMArEXFont);
  marexStyle->SetLegendFont(kMArEXFont);

  // Get moodier colours for colz
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  marexStyle->SetNumberContours(NCont);
  
  return marexStyle;
}

void SetMArEXStyle ()
{
  static TStyle* marexStyle = 0;
  std::cout << "\nApplying MArEX style settings...\n" << std::endl ;
  if ( marexStyle==0 ) marexStyle = MArEXStyle();
  gROOT->SetStyle("marexStyle");
  gROOT->ForceStyle();
}