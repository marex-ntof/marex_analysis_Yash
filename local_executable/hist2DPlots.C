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
#include "TCutG.h"
#include "TColor.h"
#include "TLegend.h"
#include "TAttMarker.h"
#include "TRandom3.h"

#include "MArEXStyle.C"

const static Int_t num_hist = 2;
TH2D* h[num_hist];
Int_t histCounter = 0;

TCutG* PTBC_tof_amp_cut;

Double_t tof_min = 1e2;
Double_t tof_max = 1e8;
Double_t amp_min = 0.;
Double_t amp_max = 70000.;

//PTBC cuts
Double_t t[6][2];
Double_t a[6][2];

void fillCutsPTBC(){
    t[0][0] = 780.0;
    a[0][0] = 7700.0;

    t[0][1] = 1090.0;
    a[0][1] = 5451.0;

    t[1][0] = 1090.0;
    a[1][0] = 5451.0;

    t[1][1] = 2605.0;
    a[1][1] = 7167.0;

    t[2][0] = 2605.0;
    a[2][0] = 9500.0;

    t[2][1] = 2856.0;
    a[2][1] = 9500.0;

    t[3][0] = 2856.0;
    a[3][0] = 7432.0;

    t[3][1] = 15290.0;
    a[3][1] = 7432.0;

    t[4][0] = 15290.0;
    a[4][0] = 7432.0;

    t[4][1] = 18708.0;
    a[4][1] = 4000.0; //2416

    t[5][0] = 18708.0;
    a[5][0] = 4000.0; //2416

    t[5][1] = 1e8;
    a[5][1] = 4000.0; //3076
}

void fillCutPlot(){
    
    PTBC_tof_amp_cut->SetVarX("x");
    PTBC_tof_amp_cut->SetVarY("y");
    PTBC_tof_amp_cut->SetPoint(0, t[0][0], amp_max);
    PTBC_tof_amp_cut->SetPoint(1, t[0][0], a[0][0]);
    PTBC_tof_amp_cut->SetPoint(2, t[0][1], a[0][1]);
    PTBC_tof_amp_cut->SetPoint(3, t[1][1], a[1][1]);
    PTBC_tof_amp_cut->SetPoint(4, t[2][0], a[2][0]);
    PTBC_tof_amp_cut->SetPoint(5, t[2][1], a[2][1]);
    PTBC_tof_amp_cut->SetPoint(6, t[3][0], a[3][0]);
    PTBC_tof_amp_cut->SetPoint(7, t[3][1], a[3][1]);
    PTBC_tof_amp_cut->SetPoint(8, t[4][1], a[4][1]);
    PTBC_tof_amp_cut->SetPoint(9, t[5][1], a[5][1]);
}

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

    PTBC_tof_amp_cut = new TCutG("PTBC_tof_amp_cut",9);

    fillCutsPTBC();
    fillCutPlot();

    // retriveHistograms("../rootFiles/pkup.root", "delT_PTBC_beam_intensity_hist");
    // retriveHistograms("../rootFiles/pkup.root", "delT_FIMG_beam_intensity_hist");

    retriveHistograms("../rootFiles/cutoffAnalysis_PTBC_ar_bottle_full.root", "PTBC_tof_amp_fOut_total");
    retriveHistograms("../rootFiles/cutoffAnalysis_PTBC_ar_bottle_full.root", "PTBC_tof_amp_fIn_total");

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

    PTBC_tof_amp_cut->SetLineColor(2);
    PTBC_tof_amp_cut->Draw("same");

    // c[i]->Print("../plots/delT_PTBC_pulse_intensity_hist.png");

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

    PTBC_tof_amp_cut->SetLineColor(2);
    PTBC_tof_amp_cut->Draw("same");

    // c[i]->Print("../plots/delT_FIMG_pulse_intensity_hist.png");
}