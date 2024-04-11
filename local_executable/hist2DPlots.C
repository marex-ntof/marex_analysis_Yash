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

const static Int_t num_cuts = 2;
TCutG* cut[num_cuts];
Int_t cutCounter = 0;

TCutG* PTBC_tof_amp_cut_det2;
TCutG* PTBC_tof_amp_cut_det3;
TCutG* PTBC_tof_amp_cut_det4;
TCutG* PTBC_tof_amp_cut_det5;
TCutG* PTBC_tof_amp_cut_det6;
TCutG* PTBC_tof_amp_cut_det7;
TCutG* PTBC_tof_amp_cut_para;

Double_t tof_min = 1e2;
Double_t tof_max = 1e8;
Double_t amp_min = 0.;
Double_t amp_max = 70000.;

//PTBC cuts
Double_t t_det2[4][2];
Double_t a_det2[4][2];
Double_t t_det3to7[5][2][2];
Double_t a_det3to7[5][2][2];

void fillCutsPTBC_det2(){
    t_det2[0][0] = 800.0;
    a_det2[0][0] = 8000.0;

    t_det2[0][1] = 2605.0;
    a_det2[0][1] = 8000.0;

    t_det2[1][0] = 2600.0;
    a_det2[1][0] = 9500.0;

    t_det2[1][1] = 2800.0;
    a_det2[1][1] = 9500.0;

    t_det2[2][0] = 2800.0;
    a_det2[2][0] = 8000.0;

    t_det2[2][1] = 5000.0;
    a_det2[2][1] = 8000.0;

    t_det2[3][0] = 5000.0;
    a_det2[3][0] = 4000.0; 

    t_det2[3][1] = 1e8;
    a_det2[3][1] = 4000.0; 
}

void fillCutsPTBC_det3(){
    t_det3to7[0][0][0] = 800.0;
    a_det3to7[0][0][0] = 5000.0;

    t_det3to7[0][0][1] = 3000.0;
    a_det3to7[0][0][1] = 5000.0;

    t_det3to7[0][1][0] = 3000.0;
    a_det3to7[0][1][0] = 3500.0;

    t_det3to7[0][1][1] = 1e8;
    a_det3to7[0][1][1] = 3500.0;
}

void fillCutsPTBC_det4(){
    t_det3to7[1][0][0] = 800.0;
    a_det3to7[1][0][0] = 6000.0;

    t_det3to7[1][0][1] = 2000.0;
    a_det3to7[1][0][1] = 6000.0;

    t_det3to7[1][1][0] = 2000.0;
    a_det3to7[1][1][0] = 3500.0;

    t_det3to7[1][1][1] = 1e8;
    a_det3to7[1][1][1] = 3500.0;
}

void fillCutsPTBC_det5(){
    t_det3to7[2][0][0] = 800.0;
    a_det3to7[2][0][0] = 7000.0;

    t_det3to7[2][0][1] = 7000.0;
    a_det3to7[2][0][1] = 7000.0;

    t_det3to7[2][1][0] = 7000.0;
    a_det3to7[2][1][0] = 3500.0;

    t_det3to7[2][1][1] = 1e8;
    a_det3to7[2][1][1] = 3500.0;
}

void fillCutsPTBC_det6(){
    t_det3to7[3][0][0] = 800.0;
    a_det3to7[3][0][0] = 6000.0;

    t_det3to7[3][0][1] = 6000.0;
    a_det3to7[3][0][1] = 6000.0;

    t_det3to7[3][1][0] = 6000.0;
    a_det3to7[3][1][0] = 4000.0;

    t_det3to7[3][1][1] = 1e8;
    a_det3to7[3][1][1] = 4000.0;
}

void fillCutsPTBC_det7(){
    t_det3to7[4][0][0] = 800.0;
    a_det3to7[4][0][0] = 4000.0;

    t_det3to7[4][0][1] = 4000.0;
    a_det3to7[4][0][1] = 4000.0;

    t_det3to7[4][1][0] = 4000.0;
    a_det3to7[4][1][0] = 3000.0;

    t_det3to7[4][1][1] = 1e8;
    a_det3to7[4][1][1] = 3000.0;
}

void fillCut_det2(){
    PTBC_tof_amp_cut_det2 = new TCutG("PTBC_tof_amp_cut_det2",8);
    PTBC_tof_amp_cut_det2->SetLineColor(2);
    PTBC_tof_amp_cut_det2->SetLineWidth(2);
    PTBC_tof_amp_cut_det2->SetVarX("x");
    PTBC_tof_amp_cut_det2->SetVarY("y");
    PTBC_tof_amp_cut_det2->SetPoint(0, t_det2[0][0], amp_max);
    PTBC_tof_amp_cut_det2->SetPoint(1, t_det2[0][0], a_det2[0][0]);
    PTBC_tof_amp_cut_det2->SetPoint(2, t_det2[0][1], a_det2[0][1]);
    PTBC_tof_amp_cut_det2->SetPoint(3, t_det2[1][0], a_det2[1][0]);
    PTBC_tof_amp_cut_det2->SetPoint(4, t_det2[1][1], a_det2[1][1]);
    PTBC_tof_amp_cut_det2->SetPoint(5, t_det2[2][0], a_det2[2][0]);
    PTBC_tof_amp_cut_det2->SetPoint(6, t_det2[2][1], a_det2[2][1]);
    PTBC_tof_amp_cut_det2->SetPoint(7, t_det2[3][0], a_det2[3][0]);
    PTBC_tof_amp_cut_det2->SetPoint(8, t_det2[3][1], a_det2[3][1]);
}

void fillCut_det3(){
    PTBC_tof_amp_cut_det3 = new TCutG("PTBC_tof_amp_cut_det3",4);
    PTBC_tof_amp_cut_det3->SetLineColor(2);
    PTBC_tof_amp_cut_det3->SetLineWidth(2);
    PTBC_tof_amp_cut_det3->SetVarX("x");
    PTBC_tof_amp_cut_det3->SetVarY("y");
    PTBC_tof_amp_cut_det3->SetPoint(0, t_det3to7[0][0][0], amp_max);
    PTBC_tof_amp_cut_det3->SetPoint(1, t_det3to7[0][0][0], a_det3to7[0][0][0]);
    PTBC_tof_amp_cut_det3->SetPoint(2, t_det3to7[0][0][1], a_det3to7[0][0][1]);
    PTBC_tof_amp_cut_det3->SetPoint(3, t_det3to7[0][1][0], a_det3to7[0][1][0]);
    PTBC_tof_amp_cut_det3->SetPoint(4, t_det3to7[0][1][1], a_det3to7[0][1][1]);
}

void fillCut_det4(){
    PTBC_tof_amp_cut_det4 = new TCutG("PTBC_tof_amp_cut_det4",4);
    PTBC_tof_amp_cut_det4->SetLineColor(2);
    PTBC_tof_amp_cut_det4->SetLineWidth(2);
    PTBC_tof_amp_cut_det4->SetVarX("x");
    PTBC_tof_amp_cut_det4->SetVarY("y");
    PTBC_tof_amp_cut_det4->SetPoint(0, t_det3to7[1][0][0], amp_max);
    PTBC_tof_amp_cut_det4->SetPoint(1, t_det3to7[1][0][0], a_det3to7[1][0][0]);
    PTBC_tof_amp_cut_det4->SetPoint(2, t_det3to7[1][0][1], a_det3to7[1][0][1]);
    PTBC_tof_amp_cut_det4->SetPoint(3, t_det3to7[1][1][0], a_det3to7[1][1][0]);
    PTBC_tof_amp_cut_det4->SetPoint(4, t_det3to7[1][1][1], a_det3to7[1][1][1]);
}

void fillCut_det5(){
    PTBC_tof_amp_cut_det5 = new TCutG("PTBC_tof_amp_cut_det5",4);
    PTBC_tof_amp_cut_det5->SetLineColor(2);
    PTBC_tof_amp_cut_det5->SetLineWidth(2);
    PTBC_tof_amp_cut_det5->SetVarX("x");
    PTBC_tof_amp_cut_det5->SetVarY("y");
    PTBC_tof_amp_cut_det5->SetPoint(0, t_det3to7[2][0][0], amp_max);
    PTBC_tof_amp_cut_det5->SetPoint(1, t_det3to7[2][0][0], a_det3to7[2][0][0]);
    PTBC_tof_amp_cut_det5->SetPoint(2, t_det3to7[2][0][1], a_det3to7[2][0][1]);
    PTBC_tof_amp_cut_det5->SetPoint(3, t_det3to7[2][1][0], a_det3to7[2][1][0]);
    PTBC_tof_amp_cut_det5->SetPoint(4, t_det3to7[2][1][1], a_det3to7[2][1][1]);
}

void fillCut_det6(){
    PTBC_tof_amp_cut_det6 = new TCutG("PTBC_tof_amp_cut_det6",4);
    PTBC_tof_amp_cut_det6->SetLineColor(2);
    PTBC_tof_amp_cut_det6->SetLineWidth(2);
    PTBC_tof_amp_cut_det6->SetVarX("x");
    PTBC_tof_amp_cut_det6->SetVarY("y");
    PTBC_tof_amp_cut_det6->SetPoint(0, t_det3to7[3][0][0], amp_max);
    PTBC_tof_amp_cut_det6->SetPoint(1, t_det3to7[3][0][0], a_det3to7[3][0][0]);
    PTBC_tof_amp_cut_det6->SetPoint(2, t_det3to7[3][0][1], a_det3to7[3][0][1]);
    PTBC_tof_amp_cut_det6->SetPoint(3, t_det3to7[3][1][0], a_det3to7[3][1][0]);
    PTBC_tof_amp_cut_det6->SetPoint(4, t_det3to7[3][1][1], a_det3to7[3][1][1]);
}

void fillCut_det7(){
    PTBC_tof_amp_cut_det7 = new TCutG("PTBC_tof_amp_cut_det7",4);
    PTBC_tof_amp_cut_det7->SetLineColor(2);
    PTBC_tof_amp_cut_det7->SetLineWidth(2);
    PTBC_tof_amp_cut_det7->SetVarX("x");
    PTBC_tof_amp_cut_det7->SetVarY("y");
    PTBC_tof_amp_cut_det7->SetPoint(0, t_det3to7[4][0][0], amp_max);
    PTBC_tof_amp_cut_det7->SetPoint(1, t_det3to7[4][0][0], a_det3to7[4][0][0]);
    PTBC_tof_amp_cut_det7->SetPoint(2, t_det3to7[4][0][1], a_det3to7[4][0][1]);
    PTBC_tof_amp_cut_det7->SetPoint(3, t_det3to7[4][1][0], a_det3to7[4][1][0]);
    PTBC_tof_amp_cut_det7->SetPoint(4, t_det3to7[4][1][1], a_det3to7[4][1][1]);
}

void fillCut_para(){
    PTBC_tof_amp_cut_para = new TCutG("PTBC_tof_amp_cut_para",4);
    PTBC_tof_amp_cut_para->SetLineColor(2);
    PTBC_tof_amp_cut_para->SetLineWidth(2);
    PTBC_tof_amp_cut_para->SetVarX("x");
    PTBC_tof_amp_cut_para->SetVarY("y");
    PTBC_tof_amp_cut_para->SetPoint(0, 800.0, amp_max);
    PTBC_tof_amp_cut_para->SetPoint(1, 800.0, 5000.0);
    PTBC_tof_amp_cut_para->SetPoint(2, 3000.0, 5000.0);
    PTBC_tof_amp_cut_para->SetPoint(3, 3000.0, 4000.0);
    PTBC_tof_amp_cut_para->SetPoint(4, 1e8, 4000.0);
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

void retriveCuts(const char *file_name, const char *cut_name){
    TFile *hist_file = TFile::Open(file_name, "READ");
    if (!hist_file || hist_file->IsZombie()) {
        cout << "Unable to open " << file_name << " for reading..." <<endl;
        return;
    }

    cut[cutCounter] = (TCutG*)hist_file->Get(cut_name);
    cutCounter++;
}

void hist2DPlots() {

    fillCutsPTBC_det2();
    fillCutsPTBC_det3();
    fillCutsPTBC_det4();
    fillCutsPTBC_det5();
    fillCutsPTBC_det6();
    fillCutsPTBC_det7();
    fillCut_para();
    fillCut_det2();
    fillCut_det3();
    fillCut_det4();
    fillCut_det5();
    fillCut_det6();
    fillCut_det7();

    // retriveHistograms("../rootFiles/pkup.root", "delT_PTBC_beam_intensity_hist");
    // retriveHistograms("../rootFiles/pkup.root", "delT_FIMG_beam_intensity_hist");

    retriveHistograms("../rootFiles/cutoffAnalysis_FIMG_ar_bottle_full.root", "FIMG_tof_amp_dedi_det1");       
    retriveCuts("../rootFiles/cutoffAnalysis_FIMG_ar_bottle_full.root", "FIMG_my_tof_amp_cut_dedi_det1");
    retriveHistograms("../rootFiles/cutoffAnalysis_FIMG_ar_bottle_full.root", "FIMG_tof_amp_dedi_det2");       
    retriveCuts("../rootFiles/cutoffAnalysis_FIMG_ar_bottle_full.root", "FIMG_my_tof_amp_cut_dedi_det2");
    // retriveHistograms("../rootFiles/cutoffAnalysis_PTBC_ar_bottle_full.root", "PTBC_tof_amp_fIn_det3_afterCuts");       //
    // retriveHistograms("../rootFiles/cutoffAnalysis_PTBC_ar_bottle_full.root", "PTBC_tof_amp_fIn_det4_afterCuts");       //
    // retriveHistograms("../rootFiles/cutoffAnalysis_PTBC_ar_bottle_full.root", "PTBC_tof_amp_fIn_det5_afterCuts");       //
    // retriveHistograms("../rootFiles/cutoffAnalysis_PTBC_ar_bottle_full.root", "PTBC_tof_amp_fIn_det6_afterCuts");       //
    // retriveHistograms("../rootFiles/cutoffAnalysis_PTBC_ar_bottle_full.root", "PTBC_tof_amp_fIn_det7_afterCuts");       //
    // retriveHistograms("../rootFiles/cutoffAnalysis_PTBC_ar_bottle_full.root", "PTBC_tof_amp_fOut_parasitic");           //
    // retriveHistograms("../rootFiles/cutoffAnalysis_PTBC_ar_bottle_full.root", "PTBC_tof_amp_fIn_parasitic");            //

    //Plotting
    SetMArEXStyle();
    // gStyle->SetOptStat(1110);
    gStyle->SetPalette(57);

    TCanvas *c[2];

    int i = 0;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();
    h[i]->GetXaxis()->SetTitle("TOF (in ns)");
    h[i]->GetYaxis()->SetTitle("Amplitude (a.u.)");
    // h[i]->GetXaxis()->SetRangeUser(3e2,2e3);
    h[i]->Draw("colz");
    // PTBC_tof_amp_cut_det2->Draw("same");
    cut[i]->Draw("same");
    gPad->SetLogx();
    gPad->SetLogz();

    // i++;

    // c[i] = new TCanvas(Form("c%d", i)," ");
    // c[i]->cd();
    // h[i]->GetXaxis()->SetTitle("TOF (in ns)");
    // h[i]->GetYaxis()->SetTitle("Amplitude (a.u.)");
    // h[i]->Draw("colz");
    // // PTBC_tof_amp_cut_det3->Draw("same");
    // cut[i]->Draw("same");
    // gPad->SetLogx();

    // i++;

    // c[i] = new TCanvas(Form("c%d", i)," ");
    // c[i]->cd();
    // h[i]->GetXaxis()->SetTitle("TOF (in ns)");
    // h[i]->GetYaxis()->SetTitle("Amplitude (a.u.)");
    // h[i]->Draw("colz");
    // PTBC_tof_amp_cut_det4->Draw("same");
    // gPad->SetLogx();

    // i++;

    // c[i] = new TCanvas(Form("c%d", i)," ");
    // c[i]->cd();
    // h[i]->GetXaxis()->SetTitle("TOF (in ns)");
    // h[i]->GetYaxis()->SetTitle("Amplitude (a.u.)");
    // h[i]->Draw("colz");
    // PTBC_tof_amp_cut_det5->Draw("same");
    // gPad->SetLogx();

    // i++;

    // c[i] = new TCanvas(Form("c%d", i)," ");
    // c[i]->cd();
    // h[i]->GetXaxis()->SetTitle("TOF (in ns)");
    // h[i]->GetYaxis()->SetTitle("Amplitude (a.u.)");
    // h[i]->Draw("colz");
    // PTBC_tof_amp_cut_det6->Draw("same");
    // gPad->SetLogx();

    // i++;

    // c[i] = new TCanvas(Form("c%d", i)," ");
    // c[i]->cd();
    // h[i]->GetXaxis()->SetTitle("TOF (in ns)");
    // h[i]->GetYaxis()->SetTitle("Amplitude (a.u.)");
    // h[i]->Draw("colz");
    // PTBC_tof_amp_cut_det7->Draw("same");
    // gPad->SetLogx();

    // i++;

    // c[i] = new TCanvas(Form("c%d", i)," ");
    // c[i]->cd();
    // h[i]->GetXaxis()->SetTitle("TOF (in ns)");
    // h[i]->GetYaxis()->SetTitle("Amplitude (a.u.)");
    // h[i]->Draw("colz");
    // PTBC_tof_amp_cut_para->Draw("same");
    // gPad->SetLogx();

    // i++;

    // c[i] = new TCanvas(Form("c%d", i)," ");
    // c[i]->cd();
    // h[i]->GetXaxis()->SetTitle("TOF (in ns)");
    // h[i]->GetYaxis()->SetTitle("Amplitude (a.u.)");
    // h[i]->Draw("colz");
    // PTBC_tof_amp_cut_para->Draw("same");
    // gPad->SetLogx();

    // c[i] = new TCanvas(Form("c%d", i)," ");
    // c[i]->cd();
    // h[i]->GetXaxis()->SetTitle("Pulse Intensity (in Num Protons)");
    // h[i]->GetYaxis()->SetTitle("#Delta t (in ns)");
    // h[i]->SetTitle("Pulse Intensity vs #Delta t for PTBC");
    // h[i]->Draw("colz");
    // // h[i]->SetMarkerStyle(6);
    // // h[i]->SetMarkerSize(0.5);
    // gPad->SetLogz();
    // gStyle->SetPalette(57);

    // PTBC_tof_amp_cut->SetLineColor(2);
    // PTBC_tof_amp_cut->Draw("same");

    // c[i]->Print("../plots/delT_PTBC_pulse_intensity_hist.png");

    // i++;

    // c[i] = new TCanvas(Form("c%d", i)," ");
    // c[i]->cd();
    // h[i]->GetXaxis()->SetTitle("Pulse Intensity (in Num Protons)");
    // h[i]->GetYaxis()->SetTitle("#Delta t (in ns)");
    // h[i]->SetTitle("Pulse Intensity vs #Delta t for FIMG");
    // h[i]->Draw("colz");
    // // h[i]->SetMarkerStyle(6);
    // // h[i]->SetMarkerSize(0.5);
    // gPad->SetLogz();
    // gStyle->SetPalette(57);

    // PTBC_tof_amp_cut->SetLineColor(2);
    // PTBC_tof_amp_cut->Draw("same");

    // c[i]->Print("../plots/delT_FIMG_pulse_intensity_hist.png");
}