/**
 * @file hist2DPlots.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-11-15
 */

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>

#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TFrame.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TLine.h"
#include "TAxis.h"
#include "TColor.h"
#include "TLegend.h"
#include "TAttMarker.h"
#include "TRandom3.h"

#include "MArEXStyle.C"

const static Int_t num_hist = 2;
TH2D* h[num_hist];
Int_t histCounter = 0;

void retriveHistograms(const char *file_name, const char *hist_name){
    TFile *hist_file = TFile::Open(file_name, "READ");
    if (!hist_file || hist_file->IsZombie()) {
        cout << "Unable to open " << file_name << " for reading..." <<endl;
        return;
    }

    h[histCounter] = (TH2D*)hist_file->Get(hist_name);
    histCounter++;
}

void hist2DPlots() {

    retriveHistograms("../rootFiles/pkup.root", "delT_PTBC_beam_intensity_hist");
    retriveHistograms("../rootFiles/pkup.root", "delT_FIMG_beam_intensity_hist");

    //Plotting
    SetMArEXStyle();
    gStyle->SetOptStat(1110);

    TCanvas *c[2];

    int i = 0;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();
    h[i]->GetXaxis()->SetTitle("Pulse Intensity (in Num Protons)");
    h[i]->GetYaxis()->SetTitle("#Delta t (in ns)");
    h[i]->SetTitle("Pulse Intensity vs #Delta t for PTBC");
    h[i]->Draw("colz");
    // h[i]->SetMarkerStyle(6);
    // h[i]->SetMarkerSize(0.5);
    gPad->SetLogz();
    gStyle->SetPalette(57);

    c[i]->Print("../plots/delT_PTBC_pulse_intensity_hist.png");

    i++;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();
    h[i]->GetXaxis()->SetTitle("Pulse Intensity (in Num Protons)");
    h[i]->GetYaxis()->SetTitle("#Delta t (in ns)");
    h[i]->SetTitle("Pulse Intensity vs #Delta t for FIMG");
    h[i]->Draw("colz");
    // h[i]->SetMarkerStyle(6);
    // h[i]->SetMarkerSize(0.5);
    gPad->SetLogz();
    gStyle->SetPalette(57);

    c[i]->Print("../plots/delT_FIMG_pulse_intensity_hist.png");
}