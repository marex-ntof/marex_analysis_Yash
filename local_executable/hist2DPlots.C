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
#include "../remote_executable/parameters.h"

const static Int_t num_hist = 2;
TH2D* h[num_hist];
Int_t histCounter = 0;

const static Int_t num_cuts = 2;
TCutG* cut[num_cuts];
Int_t cutCounter = 0;

TH2D* PTBC_tof_amp_hists[6];
TCutG* PTBC_tof_amp_cuts[6];

// TH2D* PTBC_tof_amp_hist_det2 = 0;
// TH2D* PTBC_tof_amp_hist_det3 = 0;
// TH2D* PTBC_tof_amp_hist_det4 = 0;
// TH2D* PTBC_tof_amp_hist_det5 = 0;
// TH2D* PTBC_tof_amp_hist_det6 = 0;
// TH2D* PTBC_tof_amp_hist_det7 = 0;

// TH1D* PTBC_tof_amp_cut_det2 = 0;
// TH1D* PTBC_tof_amp_cut_det3 = 0;
// TH1D* PTBC_tof_amp_cut_det4 = 0;
// TH1D* PTBC_tof_amp_cut_det5 = 0;
// TH1D* PTBC_tof_amp_cut_det6 = 0;
// TH1D* PTBC_tof_amp_cut_det7 = 0;

TCutG* FIMG_my_tof_amp_cut_dedi_det1 = 0;
TCutG* FIMG_my_tof_amp_cut_dedi_det2 = 0;
TCutG* FIMG_my_tof_amp_cut_para_det1 = 0;
TCutG* FIMG_my_tof_amp_cut_para_det2 = 0;

Double_t tof_min = 1e2;
Double_t tof_max = 1e8;
Double_t amp_min = 0.;
Double_t amp_max = 70000.;

void fillCutGraph_FIMG(){

    FIMG_my_tof_amp_cut_dedi_det1 = new TCutG("FIMG_my_tof_amp_cut_dedi_det1",6);
    FIMG_my_tof_amp_cut_dedi_det1->SetLineColor(2);
    FIMG_my_tof_amp_cut_dedi_det1->SetLineWidth(2);
    FIMG_my_tof_amp_cut_dedi_det1->SetVarX("x");
    FIMG_my_tof_amp_cut_dedi_det1->SetVarY("y");
    FIMG_my_tof_amp_cut_dedi_det1->SetPoint(0, tof_cut_FIMG[0][0][0], amp_max);
    FIMG_my_tof_amp_cut_dedi_det1->SetPoint(1, tof_cut_FIMG[0][0][0], amp_cut_FIMG[0][0][0]);
    FIMG_my_tof_amp_cut_dedi_det1->SetPoint(2, tof_cut_FIMG[0][1][0], amp_cut_FIMG[0][1][0]);
    FIMG_my_tof_amp_cut_dedi_det1->SetPoint(3, tof_cut_FIMG[0][2][0], amp_cut_FIMG[0][2][0]);
    FIMG_my_tof_amp_cut_dedi_det1->SetPoint(4, tof_cut_FIMG[0][3][0], amp_cut_FIMG[0][3][0]);
    FIMG_my_tof_amp_cut_dedi_det1->SetPoint(5, tof_cut_FIMG[0][4][0], amp_cut_FIMG[0][4][0]);
    FIMG_my_tof_amp_cut_dedi_det1->SetPoint(6, tof_cut_FIMG[0][4][1], amp_cut_FIMG[0][4][1]);

    FIMG_my_tof_amp_cut_dedi_det2 = new TCutG("FIMG_my_tof_amp_cut_dedi_det2",6);
    FIMG_my_tof_amp_cut_dedi_det2->SetLineColor(2);
    FIMG_my_tof_amp_cut_dedi_det2->SetLineWidth(2);
    FIMG_my_tof_amp_cut_dedi_det2->SetVarX("x");
    FIMG_my_tof_amp_cut_dedi_det2->SetVarY("y");
    FIMG_my_tof_amp_cut_dedi_det2->SetPoint(0, tof_cut_FIMG[1][0][0], amp_max);
    FIMG_my_tof_amp_cut_dedi_det2->SetPoint(1, tof_cut_FIMG[1][0][0], amp_cut_FIMG[1][0][0]);
    FIMG_my_tof_amp_cut_dedi_det2->SetPoint(2, tof_cut_FIMG[1][1][0], amp_cut_FIMG[1][1][0]);
    FIMG_my_tof_amp_cut_dedi_det2->SetPoint(3, tof_cut_FIMG[1][2][0], amp_cut_FIMG[1][2][0]);
    FIMG_my_tof_amp_cut_dedi_det2->SetPoint(4, tof_cut_FIMG[1][3][0], amp_cut_FIMG[1][3][0]);
    FIMG_my_tof_amp_cut_dedi_det2->SetPoint(5, tof_cut_FIMG[1][4][0], amp_cut_FIMG[1][4][0]);
    FIMG_my_tof_amp_cut_dedi_det2->SetPoint(6, tof_cut_FIMG[1][4][1], amp_cut_FIMG[1][4][1]);

    FIMG_my_tof_amp_cut_para_det1 = new TCutG("FIMG_my_tof_amp_cut_para_det1",6);
    FIMG_my_tof_amp_cut_para_det1->SetLineColor(2);
    FIMG_my_tof_amp_cut_para_det1->SetLineWidth(2);
    FIMG_my_tof_amp_cut_para_det1->SetVarX("x");
    FIMG_my_tof_amp_cut_para_det1->SetVarY("y");
    FIMG_my_tof_amp_cut_para_det1->SetPoint(0, tof_cut_FIMG[0][0][0], amp_max);
    FIMG_my_tof_amp_cut_para_det1->SetPoint(1, tof_cut_FIMG[0][0][0], amp_cut_FIMG[0][0][0]);
    FIMG_my_tof_amp_cut_para_det1->SetPoint(2, tof_cut_FIMG[0][1][0], amp_cut_FIMG[0][1][0]);
    FIMG_my_tof_amp_cut_para_det1->SetPoint(3, tof_cut_FIMG[0][2][0], amp_cut_FIMG[0][2][0]);
    FIMG_my_tof_amp_cut_para_det1->SetPoint(4, tof_cut_FIMG[0][3][0], amp_cut_FIMG[0][3][0]);
    FIMG_my_tof_amp_cut_para_det1->SetPoint(5, tof_cut_FIMG[0][4][0], amp_cut_FIMG[0][4][0]);
    FIMG_my_tof_amp_cut_para_det1->SetPoint(6, tof_cut_FIMG[0][4][1], amp_cut_FIMG[0][4][1]);

    FIMG_my_tof_amp_cut_para_det2 = new TCutG("FIMG_my_tof_amp_cut_para_det2",6);
    FIMG_my_tof_amp_cut_para_det2->SetLineColor(2);
    FIMG_my_tof_amp_cut_para_det2->SetLineWidth(2);
    FIMG_my_tof_amp_cut_para_det2->SetVarX("x");
    FIMG_my_tof_amp_cut_para_det2->SetVarY("y");
    FIMG_my_tof_amp_cut_para_det2->SetPoint(0, tof_cut_FIMG[1][0][0], amp_max);
    FIMG_my_tof_amp_cut_para_det2->SetPoint(1, tof_cut_FIMG[1][0][0], amp_cut_FIMG[1][0][0]);
    FIMG_my_tof_amp_cut_para_det2->SetPoint(2, tof_cut_FIMG[1][1][0], amp_cut_FIMG[1][1][0]);
    FIMG_my_tof_amp_cut_para_det2->SetPoint(3, tof_cut_FIMG[1][2][0], amp_cut_FIMG[1][2][0]);
    FIMG_my_tof_amp_cut_para_det2->SetPoint(4, tof_cut_FIMG[1][3][0], amp_cut_FIMG[1][3][0]);
    FIMG_my_tof_amp_cut_para_det2->SetPoint(5, tof_cut_FIMG[1][4][0], amp_cut_FIMG[1][4][0]);
    FIMG_my_tof_amp_cut_para_det2->SetPoint(6, tof_cut_FIMG[1][4][1], amp_cut_FIMG[1][4][1]);
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


TH1D* GetHist1D(const char* file_name, const char* hist_name) {
    // Open the ROOT file
    TFile *file = TFile::Open(file_name, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: cannot open file " << file_name << std::endl;
        return nullptr;
    }

    // Retrieve the histogram
    TH1D *hist = nullptr;
    file->GetObject(hist_name, hist);
    if (!hist) {
        std::cerr << "Error: cannot find histogram " << hist_name << " in file " << file_name << std::endl;
        file->Close();
        return nullptr;
    }

    return hist;

    // TH1D *histClone = (TH1D*)hist->Clone();
    // file->Close();
    // return histClone;
}

TH2D* GetHist2D(const char* file_name, const char* hist_name) {
    // Open the ROOT file
    TFile *file = TFile::Open(file_name, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: cannot open file " << file_name << std::endl;
        return nullptr;
    }

    // Retrieve the histogram
    TH2D *hist = nullptr;
    file->GetObject(hist_name, hist);
    if (!hist) {
        std::cerr << "Error: cannot find histogram " << hist_name << " in file " << file_name << std::endl;
        file->Close();
        return nullptr;
    }

    return hist;

    // TH2D *histClone = (TH2D*)hist->Clone();
    // file->Close();
    // return histClone;
}

TCutG* retrive_TCutG(const char *file_name, const char *cut_name){
    
    TFile* root_file = TFile::Open(file_name, "READ");
    TCutG* tcutg_new = (TCutG*)root_file->Get(cut_name);

    return tcutg_new;
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

void plot_PTBC_ringing(){

    TH2D* tof_amp_hist = GetHist2D("../rootFiles/cutoffAnalysis_PTBC_none.root", "PTBC_tof_amp_det2");

    //Plotting
    SetMArEXStyle();
    // gStyle->SetOptStat(1110);
    gStyle->SetPalette(57);
    gStyle->SetCanvasDefW(1200); //600
    gStyle->SetCanvasDefH(500); //500 

    TCanvas *ringing_c = new TCanvas("ringing_c"," ");
    ringing_c->cd();
    ringing_c->Draw();

    ringing_c->cd(0);
    TPad* ringing_pad_1 = new TPad("ringing_pad_1", "ringing_pad_1", 0., 0., 0.5, 1.);
    ringing_pad_1->SetFillColor(kWhite);
    ringing_pad_1->SetBorderMode(0);
    ringing_pad_1->SetTopMargin(0.07);
    ringing_pad_1->SetBottomMargin(0.1);
    ringing_pad_1->SetLeftMargin(0.1);
    ringing_pad_1->SetRightMargin(0.11);
    ringing_pad_1->Draw();
    ringing_pad_1->cd();

    tof_amp_hist->SetTitle("");
    // X Axis
    tof_amp_hist->GetXaxis()->SetLabelSize(0.04);
    tof_amp_hist->GetXaxis()->SetTitleSize(0.04);
    tof_amp_hist->GetXaxis()->SetTitleOffset(1.1);
    tof_amp_hist->GetXaxis()->SetTitle("TOF (in ns)");
    tof_amp_hist->GetXaxis()->SetRangeUser(630.957344480193, 10000.0);
    // Y Axis
    tof_amp_hist->GetYaxis()->SetTitle("Amplitude (a.u.)");
    tof_amp_hist->GetYaxis()->SetRangeUser(0., 10000.0);
    tof_amp_hist->GetYaxis()->SetMaxDigits(4);
    tof_amp_hist->Draw("COLZ");
    gPad->SetLogx();
    gPad->SetLogz();

    /////////////////

    ringing_c->cd(0);
    TPad* ringing_pad_2 = new TPad("ringing_pad_2", "ringing_pad_2", 0.5, 0., 1., 1.);
    ringing_pad_2->SetFillColor(kWhite);
    ringing_pad_2->SetBorderMode(0);
    ringing_pad_2->SetTopMargin(0.07);
    ringing_pad_2->SetBottomMargin(0.1);
    ringing_pad_2->SetLeftMargin(0.1);
    ringing_pad_2->SetRightMargin(0.11);
    ringing_pad_2->Draw();
    ringing_pad_2->cd();

    // TH2D* tof_amp_hist_rebin = (TH2D*)(tof_amp_hist->Clone("tof_amp_hist_rebin"));
    TH2D* tof_amp_hist_rebin = (TH2D*)tof_amp_hist->Rebin2D(200, 5, "tof_amp_hist_rebin");
    tof_amp_hist_rebin->SetTitle("");
    // X Axis
    tof_amp_hist_rebin->GetXaxis()->SetLabelSize(0.04);
    tof_amp_hist_rebin->GetXaxis()->SetTitleSize(0.04);
    tof_amp_hist_rebin->GetXaxis()->SetTitleOffset(1.1);
    tof_amp_hist_rebin->GetXaxis()->SetTitle("TOF (in ns)");
    tof_amp_hist_rebin->GetXaxis()->SetRangeUser(630.957344480193, 10000.0);
    // Y Axis
    tof_amp_hist_rebin->GetYaxis()->SetTitle("Amplitude (a.u.)");
    tof_amp_hist_rebin->GetYaxis()->SetRangeUser(0., 10000.0);
    tof_amp_hist_rebin->GetYaxis()->SetMaxDigits(4);
    tof_amp_hist_rebin->Draw("COLZ");
    gPad->SetLogx();
    gPad->SetLogz();
}

void plot_det_cuts_ptbc(){

    for (Int_t i = 0; i < 6; i++)
    {
        PTBC_tof_amp_hists[i] = GetHist2D("../rootFiles/cutoffAnalysis_PTBC_none.root", Form("PTBC_tof_amp_det%i", i+2));
        PTBC_tof_amp_cuts[i] = retrive_TCutG("../inputFiles/PTBC_cuts.root", Form("tof_amp_cut_det%i", i+2));
    }

    //Plotting
    SetMArEXStyle();
    // gStyle->SetOptStat(1110);
    gStyle->SetPalette(57);
    gStyle->SetCanvasDefW(1000); //600
    gStyle->SetCanvasDefH(500); //500 

    TCanvas *c[6];

    int i = 0;

    for (Int_t i = 0; i < 6; i++)
    {
        c[i] = new TCanvas(Form("c%d", i)," ");
        c[i]->cd();
        c[i]->SetBorderMode(0);
        c[i]->SetTopMargin(0.1);
        c[i]->SetBottomMargin(0.1);
        c[i]->SetLeftMargin(0.1);
        c[i]->SetRightMargin(0.09);
        
        PTBC_tof_amp_hists[i]->SetTitle(Form("Detector %i", i+2));
        // X Axis
        PTBC_tof_amp_hists[i]->GetXaxis()->SetLabelSize(0.04);
        PTBC_tof_amp_hists[i]->GetXaxis()->SetTitleSize(0.04);
        PTBC_tof_amp_hists[i]->GetXaxis()->SetTitleOffset(1.1);
        PTBC_tof_amp_hists[i]->GetXaxis()->SetTitle("TOF (in ns)");
        // Y Axis
        PTBC_tof_amp_hists[i]->GetYaxis()->SetLabelSize(0.04);
        PTBC_tof_amp_hists[i]->GetYaxis()->SetTitleSize(0.04);
        PTBC_tof_amp_hists[i]->GetYaxis()->SetTitleOffset(1.2);
        PTBC_tof_amp_hists[i]->GetYaxis()->SetTitle("Amplitude (a.u.)");
        // PTBC_tof_amp_hists[i]->GetXaxis()->SetRangeUser(3e2,2e3);
        PTBC_tof_amp_hists[i]->Draw("COLZ");
        gPad->SetLogx();
        gPad->SetLogz();

        PTBC_tof_amp_cuts[i]->SetLineColor(2);
        PTBC_tof_amp_cuts[i]->SetLineWidth(2);
        PTBC_tof_amp_cuts[i]->Draw("SAME");
        c[i]->Print(Form("/home/yash/Thesis/marex/figures/bkgd_plots/tof_amp_hist_cut_det%i.png", i+2));
    }
}

void plot_det_cuts_fimg(){

    TH2D* FIMG_tof_amp_hists[2];
    TCutG* FIMG_tof_amp_cuts[2];

    for (Int_t i = 0; i < 2; i++)
    {
        FIMG_tof_amp_hists[i] = GetHist2D("../rootFiles/cutoffAnalysis_FIMG_ar_bottle_full.root", Form("FIMG_tof_amp_dedi_det%i", i+1));
        FIMG_tof_amp_cuts[i] = retrive_TCutG("../rootFiles/cutoffAnalysis_FIMG_ar_bottle_full.root", Form("FIMG_my_tof_amp_cut_dedi_det%i", i+1));
    }

    //Plotting
    SetMArEXStyle();
    // gStyle->SetOptStat(1110);
    gStyle->SetPalette(57);
    gStyle->SetCanvasDefW(1000); //600
    gStyle->SetCanvasDefH(500); //500


    TCanvas *c[2];

    for (Int_t i = 0; i < 2; i++)
    {
        c[i] = new TCanvas(Form("c%d", i)," ");
        c[i]->cd();
        c[i]->SetBorderMode(0);
        c[i]->SetTopMargin(0.1);
        c[i]->SetBottomMargin(0.1);
        c[i]->SetLeftMargin(0.1);
        c[i]->SetRightMargin(0.09);
        
        FIMG_tof_amp_hists[i]->SetTitle(Form("Detector %i - Micromegas", i+1));
        // X Axis
        FIMG_tof_amp_hists[i]->GetXaxis()->SetLabelSize(0.04);
        FIMG_tof_amp_hists[i]->GetXaxis()->SetTitleSize(0.04);
        FIMG_tof_amp_hists[i]->GetXaxis()->SetTitleOffset(1.1);
        FIMG_tof_amp_hists[i]->GetXaxis()->SetTitle("TOF (in ns)");
        // Y Axis
        FIMG_tof_amp_hists[i]->GetYaxis()->SetLabelSize(0.04);
        FIMG_tof_amp_hists[i]->GetYaxis()->SetTitleSize(0.04);
        FIMG_tof_amp_hists[i]->GetYaxis()->SetTitleOffset(1.2);
        FIMG_tof_amp_hists[i]->GetYaxis()->SetTitle("Amplitude (a.u.)");
        FIMG_tof_amp_hists[i]->GetYaxis()->SetRangeUser(0,5100);
        FIMG_tof_amp_hists[i]->SetStats(0);
        FIMG_tof_amp_hists[i]->Draw("COLZ");
        gPad->SetLogx();
        gPad->SetLogz();

        FIMG_tof_amp_cuts[i]->SetLineColor(2);
        FIMG_tof_amp_cuts[i]->SetLineWidth(2);
        FIMG_tof_amp_cuts[i]->Draw("SAME");
        c[i]->Print(Form("/home/yash/Thesis/ntofAnalysis/figures/bkgd_plots/tof_amp_hist_cut_det%i_fimg.png", i+1));
    }
}

void hist2DPlots() {

    // plot_PTBC_ringing();

    // plot_det_cuts_ptbc();

    // plot_det_cuts_fimg();

    // fillCutGraph_PTBC();
    // fillCutGraph_FIMG();

    // retriveHistograms("../rootFiles/pkup.root", "delT_PTBC_beam_intensity_hist");
    // retriveHistograms("../rootFiles/pkup.root", "delT_FIMG_beam_intensity_hist");

    

    // c[i] = new TCanvas(Form("c%d", i)," ");
    // c[i]->cd();
    // h[i]->GetXaxis()->SetTitle("TOF (in ns)");
    // h[i]->GetYaxis()->SetTitle("Amplitude (a.u.)");
    // // h[i]->GetXaxis()->SetRangeUser(3e2,2e3);
    // h[i]->Draw("COLZ");
    // PTBC_tof_amp_cut_det2->SetLineColor(2);
    // PTBC_tof_amp_cut_det2->Draw("SAME");
    // // cut[i]->Draw("same");
    // gPad->SetLogx();
    // gPad->SetLogz();

    // i++;

    // c[i] = new TCanvas(Form("c%d", i)," ");
    // c[i]->cd();
    // h[i]->GetXaxis()->SetTitle("TOF (in ns)");
    // h[i]->GetYaxis()->SetTitle("Amplitude (a.u.)");
    // // h[i]->GetXaxis()->SetRangeUser(3e2,2e3);
    // h[i]->Draw("COLZ");
    // PTBC_tof_amp_cut_det3->SetLineColor(2);
    // PTBC_tof_amp_cut_det3->Draw("SAME");
    // // cut[i]->Draw("same");
    // gPad->SetLogx();
    // gPad->SetLogz();
}