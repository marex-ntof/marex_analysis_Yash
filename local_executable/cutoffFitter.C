/**
 * @file cutoffFitter.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-11-13
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

// TH1D* h[3];
// TH1D* h_old[3];
// Int_t histCounter = 0;
TH2D* tof_area_hist_fOut_det1 = 0;
TH2D* tof_area_hist_fOut_det2 = 0;
TH2D* tof_amp_hist_fOut_det1 = 0;
TH2D* tof_amp_hist_fOut_det2 = 0;

TH2D* tof_area_hist_fIn_det1 = 0;
TH2D* tof_area_hist_fIn_det2 = 0;
TH2D* tof_amp_hist_fIn_det1 = 0;
TH2D* tof_amp_hist_fIn_det2 = 0;
// TH2D* tof_area_hist_out_cutoff = 0;

TH2D* tof_amp_hist_fOut_det1_cutoff = 0;
TH2D* tof_amp_hist_fOut_det2_cutoff = 0;
TH2D* tof_amp_hist_fIn_det1_cutoff = 0;
TH2D* tof_amp_hist_fIn_det2_cutoff = 0;

void retriveHistograms(const char *fname){
    TKey *key;
    TFile *hist_file = TFile::Open(fname, "READ");
    if (!hist_file || hist_file->IsZombie()) {
        cout << "Unable to open " << fname << " for reading..." <<endl;
        return;
    }

    // tof_hist_filter_in = (TH1D*)hist_file->Get("tof_hist_filter_in");
    // energy_hist_filter_in = (TH1D*)hist_file->Get("energy_hist_filter_in");
    // tof_hist_filter_out = (TH1D*)hist_file->Get("tof_hist_filter_out");
    // energy_hist_filter_out = (TH1D*)hist_file->Get("energy_hist_filter_out");
    tof_area_hist_fOut_det1 = (TH2D*)hist_file->Get("tof_area_hist_fOut_det1");
    tof_area_hist_fOut_det2 = (TH2D*)hist_file->Get("tof_area_hist_fOut_det2");
    tof_amp_hist_fOut_det1 = (TH2D*)hist_file->Get("tof_amp_hist_fOut_det1");
    tof_amp_hist_fOut_det2 = (TH2D*)hist_file->Get("tof_amp_hist_fOut_det2");
    // tof_area_hist_out_cutoff = (TH2D*)hist_file->Get("tof_area_hist_out_cutoff");

    tof_area_hist_fIn_det1 = (TH2D*)hist_file->Get("tof_area_hist_fIn_det1");
    tof_area_hist_fIn_det2 = (TH2D*)hist_file->Get("tof_area_hist_fIn_det2");
    tof_amp_hist_fIn_det1 = (TH2D*)hist_file->Get("tof_amp_hist_fIn_det1");
    tof_amp_hist_fIn_det2 = (TH2D*)hist_file->Get("tof_amp_hist_fIn_det2");

    tof_amp_hist_fOut_det1_cutoff = (TH2D*)hist_file->Get("tof_amp_hist_fOut_det1_cutoff");
    tof_amp_hist_fOut_det2_cutoff = (TH2D*)hist_file->Get("tof_amp_hist_fOut_det2_cutoff");
    tof_amp_hist_fIn_det1_cutoff = (TH2D*)hist_file->Get("tof_amp_hist_fIn_det1_cutoff");
    tof_amp_hist_fIn_det2_cutoff = (TH2D*)hist_file->Get("tof_amp_hist_fIn_det2_cutoff");

    // trans_hist_fOut = (TH1D*)hist_file->Get("trans_hist_fOut");
    // trans_hist_fIn = (TH1D*)hist_file->Get("trans_hist_fIn");
    // trans_hist_fOut_endf = (TH1D*)hist_file->Get("trans_hist_fOut_endf");
    // trans_hist_fIn_endf = (TH1D*)hist_file->Get("trans_hist_fIn_endf");
    // h[histCounter] = (TH1D*)hist_file->Get("transmission_hist_e");
    // histCounter++;
    // cross_section_hist_e = (TH1D*)hist_file->Get("cross_section_hist_e");

    // hist_file->Close();
}

void cutoffFitter() {
    
    retriveHistograms("../rootFiles/cutoffSelector_FIMG.root");

    //Plotting
    SetMArEXStyle();
    
    gStyle->SetStatX(0.27);
    gStyle->SetStatY(0.9);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.17);

    TCanvas *c[2];
    Int_t i = 0;

    // c[i] = new TCanvas(Form("c%d", i)," ");
    // c[i]->cd();
    // tof_area_hist_fOut_det1->GetXaxis()->SetTitle("Time of Flight (in ns)");
    // tof_area_hist_fOut_det1->GetYaxis()->SetTitle("Area (a.u.)");
    // tof_area_hist_fOut_det1->SetTitle("ToF vs Area Hist - FIMG Det 1 - Filter Out");
    // tof_area_hist_fOut_det1->Draw("colz");
    // // tof_area_hist_fOut_det1->SetMarkerStyle(6);
    // // tof_area_hist_fOut_det1->SetMarkerSize(0.5);
    // gPad->SetLogx();
    // gPad->SetLogz();
    // gStyle->SetPalette(57);

    // c[i]->Print("../plots/h_tof_area_fOut_FIMG_det1.png");

    // Double_t x[6] = {61878.8,28827.2,19880,13155.2,9751.95,7456.49};
    // Double_t y[6] = {20380.1,171383,290708,455438,600106,774340};
    // TGraph *cutoff = new TGraph(6,x,y);
    // cutoff->SetMarkerColor(1);
    // cutoff->SetMarkerSize(1);
    // cutoff->SetMarkerStyle(8);
    // // cutoff->Draw("P");
    // TF1 *f = new TF1("f", "[2] * TMath::Log(x) * TMath::Log(x) + [1] * TMath::Log(x) + [0]");
    // //[6] * x * x * x * x * x * x + [5] * x * x * x * x * x + [4] * x * x * x * x + [3] * x * x * x + [2] * x * x + [1] * x + [0]
    // cutoff->Fit(f);
    // cutoff->Draw("PSAME");

    // c[i]->Print("../plots/h_tof_area_FIMG_cuts.png");

    // i++;

    // c[i] = new TCanvas(Form("c%d", i)," ");
    // c[i]->cd();
    // tof_area_hist_fOut_det2->GetXaxis()->SetTitle("Time of Flight (in ns)");
    // tof_area_hist_fOut_det2->GetYaxis()->SetTitle("Area (a.u.)");
    // tof_area_hist_fOut_det2->SetTitle("ToF vs Area Hist - FIMG Det 2 - Filter Out");
    // tof_area_hist_fOut_det2->Draw("colz");
    // // tof_area_hist_fOut_det2->SetMarkerStyle(6);
    // // tof_area_hist_fOut_det2->SetMarkerSize(0.5);
    // gPad->SetLogx();
    // gPad->SetLogz();
    // gStyle->SetPalette(57);

    // c[i]->Print("../plots/h_tof_area_fOut_FIMG_det2.png");

    // i++;

    // c[i] = new TCanvas(Form("c%d", i)," ");
    // c[i]->cd();
    // tof_area_hist_out_cutoff->GetXaxis()->SetTitle("Time of Flight (in ns)");
    // tof_area_hist_out_cutoff->GetYaxis()->SetTitle("Area (a.u.)");
    // tof_area_hist_out_cutoff->SetTitle("Area vs Area Hist - FIMG - Filter Out After Cuts");
    // tof_area_hist_out_cutoff->Draw("colz");
    // // tof_area_hist_out_cutoff->SetMarkerStyle(6);
    // // tof_area_hist_out_cutoff->SetMarkerSize(0.5);
    // gPad->SetLogx();
    // gPad->SetLogz();
    // gStyle->SetPalette(57);

    // c[i]->Print("../plots/h_tof_area_FIMG_afterCuts.png");

    // i++;

    // c[i] = new TCanvas(Form("c%d", i)," ");
    // c[i]->cd();
    // tof_amp_hist_fOut_det1->GetXaxis()->SetTitle("Time of Flight (in ns)");
    // tof_amp_hist_fOut_det1->GetYaxis()->SetTitle("Amplitude (a.u.)");
    // tof_amp_hist_fOut_det1->SetTitle("ToF vs Amplitude Hist - FIMG Det 1 - Filter Out");
    // tof_amp_hist_fOut_det1->Draw("colz");
    // // tof_amp_hist_fOut_det1->SetMarkerStyle(6);
    // // tof_amp_hist_fOut_det1->SetMarkerSize(0.5);
    // gPad->SetLogx();
    // gPad->SetLogz();
    // gStyle->SetPalette(57);

    // cout << "Fitting cut curve for Det 1 Filter Out..." << endl;
    // Double_t x_det1[4] = {44286, 16069, 7889, 7413}; // 117788 // 900000, 100000, 22787, 12392, 8256, 7188
    // Double_t y_det1[4] = {600, 1934, 3914, 4443}; // 250, 600, 1496, 2488, 3499, 4431
    // TGraph *cutoff_det1 = new TGraph(4,x_det1,y_det1);
    // cutoff_det1->SetMarkerColor(2);
    // cutoff_det1->SetMarkerSize(1);
    // cutoff_det1->SetMarkerStyle(8);
    // // cutoff_det1->Draw("P");
    // TF1 *f_det1 = new TF1("f_det1", "[2] * TMath::Log(x) * TMath::Log(x) + [1] * TMath::Log(x) + [0]");
    // //[3] * TMath::Log(x) * TMath::Log(x) * TMath::Log(x) + [2] * TMath::Log(x) * TMath::Log(x) + [1] * TMath::Log(x) + [0]
    // //[6] * x * x * x * x * x * x + [5] * x * x * x * x * x + [4] * x * x * x * x + [3] * x * x * x + [2] * x * x + [1] * x + [0]
    // cutoff_det1->Fit(f_det1);
    // cout << " " << endl;
    // f_det1->SetLineWidth(2);
    // cutoff_det1->Draw("PSAME");

    // // c[i]->Print("../plots/h_tof_amp_fOut_FIMG_det1.png");

    // i++;

    // c[i] = new TCanvas(Form("c%d", i)," ");
    // c[i]->cd();
    // tof_amp_hist_fOut_det2->GetXaxis()->SetTitle("Time of Flight (in ns)");
    // tof_amp_hist_fOut_det2->GetYaxis()->SetTitle("Amplitude (a.u.)");
    // tof_amp_hist_fOut_det2->SetTitle("ToF vs Amplitude Hist - FIMG Det 2 - Filter Out");
    // tof_amp_hist_fOut_det2->Draw("colz");
    // // tof_amp_hist_fOut_det2->SetMarkerStyle(6);
    // // tof_amp_hist_fOut_det2->SetMarkerSize(0.5);
    // gPad->SetLogx();
    // gPad->SetLogz();
    // gStyle->SetPalette(57);

    // cout << "Fitting cut curve for Det 2 Filter Out..." << endl;
    // Double_t x_det2[3] = {43883, 21165, 11196}; //200000, 100000, 23212, 13342, 8971, 7739
    // Double_t y_det2[3] = {400, 1540, 4072}; //250, 400, 1352, 2364, 3130, 4236
    // TGraph *cutoff_det2 = new TGraph(3,x_det2,y_det2);
    // cutoff_det2->SetMarkerColor(2);
    // cutoff_det2->SetMarkerSize(1);
    // cutoff_det2->SetMarkerStyle(8);
    // // cutoff_det2->Draw("P");
    // TF1 *f_det2 = new TF1("f_det2", "[2] * TMath::Log(x) * TMath::Log(x) + [1] * TMath::Log(x) + [0]");
    // //[3] * TMath::Log(x) * TMath::Log(x) * TMath::Log(x) + [2] * TMath::Log(x) * TMath::Log(x) + [1] * TMath::Log(x) + [0]
    // //[6] * x * x * x * x * x * x + [5] * x * x * x * x * x + [4] * x * x * x * x + [3] * x * x * x + [2] * x * x + [1] * x + [0]
    // cutoff_det2->Fit(f_det2);
    // f_det2->SetLineWidth(2);
    // cutoff_det2->Draw("PSAME");

    // // c[i]->Print("../plots/h_tof_amp_fOut_FIMG_det2.png");

    // // TLine* l_zero = new TLine(1e-2,700,1e9,700);
    // // l_zero->SetLineWidth(2);
    // // l_zero->SetLineColor(2);
    // // l_zero->Draw();

    // i++;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();
    tof_amp_hist_fIn_det1->GetXaxis()->SetTitle("Time of Flight (in ns)");
    tof_amp_hist_fIn_det1->GetYaxis()->SetTitle("Amplitude (a.u.)");
    // tof_amp_hist_fIn_det1->SetTitle("ToF vs Amplitude Hist - FIMG Det 1 - Filter In");
    tof_amp_hist_fIn_det1->Draw("colz");
    // tof_amp_hist_fIn_det1->SetMarkerStyle(6);
    // tof_amp_hist_fIn_det1->SetMarkerSize(0.5);
    gPad->SetLogx();
    gPad->SetLogz();
    gStyle->SetPalette(57);

    // cout << "Fitting cut curve for Det 1 Filter In..." << endl;
    // Double_t x_det1_in[4] = {13692, 31869, 73521, 106811}; //7750, 9176, 
    // Double_t y_det1_in[4] = {2406, 1588, 1257, 1160}; //4438, 3248, 
    // TGraph *cutoff_det1_in = new TGraph(4,x_det1_in,y_det1_in);
    // cutoff_det1_in->SetMarkerColor(2);
    // cutoff_det1_in->SetMarkerSize(1);
    // cutoff_det1_in->SetMarkerStyle(8);
    // // cutoff_det1_in->Draw("P");
    // TF1 *f_det1_in = new TF1("f_det1_in", "[1] / ( TMath::Log(x) + [0] )");
    // cutoff_det1_in->Fit(f_det1_in);
    // f_det1_in->SetLineWidth(2);
    // cutoff_det1_in->Draw("PSAME");

    auto tight_cut_det1 = new TF1("tight_cut_det1","(800 / ( TMath::Log(x) - TMath::Log(7750) )) + 600",1e3,1e8);
    tight_cut_det1->SetLineColor(2);
    tight_cut_det1->Draw("SAME");

    auto mid_cut_det1 = new TF1("mid_cut_det1","(1800 / ( TMath::Log(x) - TMath::Log(7750) )) + 600",1e3,1e8);
    mid_cut_det1->SetLineColor(3);
    mid_cut_det1->Draw("SAME");


    auto loose_cut_det1 = new TF1("loose_cut_det1","(2800 / ( TMath::Log(x) - TMath::Log(7750) )) + 600",1e3,1e8);
    loose_cut_det1->SetLineColor(4);
    loose_cut_det1->Draw("SAME");

    auto x_asym_det1 = new TF1("x_asym_det1","600",1e3,1e8);
    x_asym_det1->SetLineColor(1);
    x_asym_det1->SetLineWidth(2);
    x_asym_det1->Draw("SAME");

    auto y_asym_det1 = new TLine(7750.0,0,7750.0,4500);
    y_asym_det1->SetLineColor(1);
    y_asym_det1->SetLineWidth(2);
    y_asym_det1->Draw("SAME");

    i++;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();
    tof_amp_hist_fIn_det2->GetXaxis()->SetTitle("Time of Flight (in ns)");
    tof_amp_hist_fIn_det2->GetYaxis()->SetTitle("Amplitude (a.u.)");
    // tof_amp_hist_fIn_det2->SetTitle("ToF vs Amplitude Hist - FIMG Det 2 - Filter In");
    tof_amp_hist_fIn_det2->Draw("colz");
    // tof_amp_hist_fIn_det2->SetMarkerStyle(6);
    // tof_amp_hist_fIn_det2->SetMarkerSize(0.5);
    gPad->SetLogx();
    gPad->SetLogz();
    gStyle->SetPalette(57);

    // cout << "Fitting cut curve for Det 2 Filter Out..." << endl;
    // Double_t x_det2_in[5] = {84012, 45081, 34525, 15099}; //, 7959
    // Double_t y_det2_in[5] = {1116, 1353, 1530, 2127}; //, 2947
    // TGraph *cutoff_det2_in = new TGraph(5,x_det2_in,y_det2_in);
    // cutoff_det2_in->SetMarkerColor(2);
    // cutoff_det2_in->SetMarkerSize(1);
    // cutoff_det2_in->SetMarkerStyle(8);
    // // cutoff_det2_in->Draw("P");
    // TF1 *f_det2_in = new TF1("f_det2_in", "[1] / ( TMath::Log(x) + [0] )");
    // cutoff_det2_in->Fit(f_det2_in);
    // f_det2_in->SetLineWidth(2);
    // cutoff_det2_in->Draw("PSAME");

    auto tight_cut_det2 = new TF1("tight_cut_det2","(1000 / ( TMath::Log(x) - TMath::Log(7080) )) + 400",1e3,1e8);
    tight_cut_det2->SetLineColor(2);
    tight_cut_det2->Draw("SAME");

    auto mid_cut_det2 = new TF1("mid_cut_det2","(2000 / ( TMath::Log(x) - TMath::Log(7080) )) + 400",1e3,1e8);
    mid_cut_det2->SetLineColor(3);
    mid_cut_det2->Draw("SAME");

    auto loose_cut_det2 = new TF1("loose_cut_det2","(3000 / ( TMath::Log(x) - TMath::Log(7080) )) + 400",1e3,1e8);
    loose_cut_det2->SetLineColor(4);
    loose_cut_det2->Draw("SAME");

    auto x_asym_det2 = new TF1("x_asym_det2","400",1e3,1e8);
    x_asym_det2->SetLineColor(1);
    x_asym_det2->Draw("SAME");

    auto y_asym_det2 = new TLine(7080.0,0,7080.0,4500);
    y_asym_det2->SetLineColor(1);
    y_asym_det2->SetLineWidth(2);
    y_asym_det2->Draw("SAME");

    // i++;

    // c[i] = new TCanvas(Form("c%d", i)," ");
    // c[i]->cd();
    // tof_amp_hist_fOut_det1_cutoff->GetXaxis()->SetTitle("Time of Flight (in ns)");
    // tof_amp_hist_fOut_det1_cutoff->GetYaxis()->SetTitle("Amplitude (a.u.)");
    // // tof_amp_hist_fOut_det1_cutoff->SetTitle("ToF vs Amplitude Hist - FIMG Det 2 - Filter In");
    // tof_amp_hist_fOut_det1_cutoff->Draw("colz");
    // // tof_amp_hist_fOut_det1_cutoff->SetMarkerStyle(6);
    // // tof_amp_hist_fOut_det1_cutoff->SetMarkerSize(0.5);
    // gPad->SetLogx();
    // gPad->SetLogz();
    // gStyle->SetPalette(57);

    // i++;

    // c[i] = new TCanvas(Form("c%d", i)," ");
    // c[i]->cd();
    // tof_amp_hist_fOut_det2_cutoff->GetXaxis()->SetTitle("Time of Flight (in ns)");
    // tof_amp_hist_fOut_det2_cutoff->GetYaxis()->SetTitle("Amplitude (a.u.)");
    // // tof_amp_hist_fOut_det2_cutoff->SetTitle("ToF vs Amplitude Hist - FIMG Det 2 - Filter In");
    // tof_amp_hist_fOut_det2_cutoff->Draw("colz");
    // // tof_amp_hist_fOut_det2_cutoff->SetMarkerStyle(6);
    // // tof_amp_hist_fOut_det2_cutoff->SetMarkerSize(0.5);
    // gPad->SetLogx();
    // gPad->SetLogz();
    // gStyle->SetPalette(57);
}