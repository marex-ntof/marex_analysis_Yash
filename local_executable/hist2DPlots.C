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
TH1D* PTBC_tof_amp_cuts[6];

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

    // fillCutGraph_PTBC();
    // fillCutGraph_FIMG();

    // retriveHistograms("../rootFiles/pkup.root", "delT_PTBC_beam_intensity_hist");
    // retriveHistograms("../rootFiles/pkup.root", "delT_FIMG_beam_intensity_hist");

    for (Int_t i = 0; i < 6; i++)
    {
        PTBC_tof_amp_hists[i] = GetHist2D("../rootFiles/cutoffAnalysis_PTBC_al5.root", Form("PTBC_tof_amp_fIn_det%i", i+2));
        PTBC_tof_amp_cuts[i] = GetHist1D("../rootFiles/PTBC_cuts.root", Form("PTBC_cuts_det%i", i+2));
    }

    //Plotting
    SetMArEXStyle();
    // gStyle->SetOptStat(1110);
    gStyle->SetPalette(57);

    TCanvas *c[6];

    int i = 0;

    for (Int_t i = 0; i < 6; i++)
    {
        c[i] = new TCanvas(Form("c%d", i)," ");
        c[i]->cd();
        PTBC_tof_amp_hists[i]->GetXaxis()->SetTitle("TOF (in ns)");
        PTBC_tof_amp_hists[i]->GetYaxis()->SetTitle("Amplitude (a.u.)");
        // PTBC_tof_amp_hists[i]->GetXaxis()->SetRangeUser(3e2,2e3);
        PTBC_tof_amp_hists[i]->Draw("COLZ");
        gPad->SetLogx();
        gPad->SetLogz();

        PTBC_tof_amp_cuts[i]->SetLineColor(2);
        PTBC_tof_amp_cuts[i]->SetLineWidth(2);
        PTBC_tof_amp_cuts[i]->Draw("SAME");
    }

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