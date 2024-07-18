/**
 * @file anomalies_plots.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2024-05-14
 */

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <cmath>

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

// #include "parameters.h"
// #include "tools.h"

TH2D* PTBC_al5_anomaly_hists[6];
TH2D* PTBC_filterOut_anomaly_hists[6];
TCutG* PTBC_tof_amp_cuts[6];

// TH2D* PTBC_al5_anomaly_det2 = 0;
// TH2D* PTBC_al5_anomaly_det3 = 0;
// TH2D* PTBC_al5_anomaly_det4 = 0;
// TH2D* PTBC_al5_anomaly_det5 = 0;
// TH2D* PTBC_al5_anomaly_det6 = 0;
// TH2D* PTBC_al5_anomaly_det7 = 0;

// TH2D* FIMG_al5_anomaly_det1 = 0;
// TH2D* FIMG_al5_anomaly_det2 = 0;

// TH2D* PTBC_filterOut_anomaly_det2 = 0;
// TH2D* PTBC_filterOut_anomaly_det3 = 0;
// TH2D* PTBC_filterOut_anomaly_det4 = 0;
// TH2D* PTBC_filterOut_anomaly_det5 = 0;
// TH2D* PTBC_filterOut_anomaly_det6 = 0;
// TH2D* PTBC_filterOut_anomaly_det7 = 0;

// TCutG* PTBC_tof_amp_cut_det2 = 0;
// TCutG* PTBC_tof_amp_cut_det3 = 0;
// TCutG* PTBC_tof_amp_cut_det4 = 0;
// TCutG* PTBC_tof_amp_cut_det5 = 0;
// TCutG* PTBC_tof_amp_cut_det6 = 0;
// TCutG* PTBC_tof_amp_cut_det7 = 0;
// TCutG* FIMG_my_tof_amp_cut_dedi_det1 = 0;
// TCutG* FIMG_my_tof_amp_cut_dedi_det2 = 0;
// TCutG* FIMG_my_tof_amp_cut_para_det1 = 0;
// TCutG* FIMG_my_tof_amp_cut_para_det2 = 0;

TH2D* retrive_TH2D_Histograms(const char *file_name, const char *hist_name){
    
    TFile* hist_file = TFile::Open(file_name, "READ");
    TH2D* hist_new = (TH2D*)hist_file->Get(hist_name);

    return hist_new;
}

TCutG* retrive_TCutG(const char *file_name, const char *cut_name){
    
    TFile* root_file = TFile::Open(file_name, "READ");
    TCutG* tcutg_new = (TCutG*)root_file->Get(cut_name);

    return tcutg_new;
}


// void retrive_all_hists_cuts(){

    // PTBC_al5_anomaly_det2 = retriveHistograms("../rootFiles/anomalies.root", "PTBC_al5_anomaly_det2");
    // PTBC_al5_anomaly_det3 = retriveHistograms("../rootFiles/anomalies.root", "PTBC_al5_anomaly_det3");
    // PTBC_al5_anomaly_det4 = retriveHistograms("../rootFiles/anomalies.root", "PTBC_al5_anomaly_det4");
    // PTBC_al5_anomaly_det5 = retriveHistograms("../rootFiles/anomalies.root", "PTBC_al5_anomaly_det5");
    // PTBC_al5_anomaly_det6 = retriveHistograms("../rootFiles/anomalies.root", "PTBC_al5_anomaly_det6");
    // PTBC_al5_anomaly_det7 = retriveHistograms("../rootFiles/anomalies.root", "PTBC_al5_anomaly_det7");

    // FIMG_al5_anomaly_det1 = retriveHistograms("../rootFiles/anomalies.root", "FIMG_al5_anomaly_det1");
    // FIMG_al5_anomaly_det2 = retriveHistograms("../rootFiles/anomalies.root", "FIMG_al5_anomaly_det2");

    // PTBC_filterOut_anomaly_det2 = retriveHistograms("../rootFiles/anomalies.root", "PTBC_filterOut_anomaly_det2");
    // PTBC_filterOut_anomaly_det3 = retriveHistograms("../rootFiles/anomalies.root", "PTBC_filterOut_anomaly_det3");
    // PTBC_filterOut_anomaly_det4 = retriveHistograms("../rootFiles/anomalies.root", "PTBC_filterOut_anomaly_det4");
    // PTBC_filterOut_anomaly_det5 = retriveHistograms("../rootFiles/anomalies.root", "PTBC_filterOut_anomaly_det5");
    // PTBC_filterOut_anomaly_det6 = retriveHistograms("../rootFiles/anomalies.root", "PTBC_filterOut_anomaly_det6");
    // PTBC_filterOut_anomaly_det7 = retriveHistograms("../rootFiles/anomalies.root", "PTBC_filterOut_anomaly_det7");

    // PTBC_tof_amp_cut_det2 = retrive_TCutG("../inputFiles/PTBC_cuts.root", "PTBC_tof_amp_cut_det2");
    // PTBC_tof_amp_cut_det3 = retrive_TCutG("../inputFiles/PTBC_cuts.root", "PTBC_tof_amp_cut_det3");
    // PTBC_tof_amp_cut_det4 = retrive_TCutG("../inputFiles/PTBC_cuts.root", "PTBC_tof_amp_cut_det4");
    // PTBC_tof_amp_cut_det5 = retrive_TCutG("../inputFiles/PTBC_cuts.root", "PTBC_tof_amp_cut_det5");
    // PTBC_tof_amp_cut_det6 = retrive_TCutG("../inputFiles/PTBC_cuts.root", "PTBC_tof_amp_cut_det6");
    // PTBC_tof_amp_cut_det7 = retrive_TCutG("../inputFiles/PTBC_cuts.root", "PTBC_tof_amp_cut_det7");

    // FIMG_my_tof_amp_cut_dedi_det1 = retrive_TCutG("../rootFiles/anomalies.root", "FIMG_my_tof_amp_cut_dedi_det1");
    // FIMG_my_tof_amp_cut_dedi_det2 = retrive_TCutG("../rootFiles/anomalies.root", "FIMG_my_tof_amp_cut_dedi_det2");
    // FIMG_my_tof_amp_cut_para_det1 = retrive_TCutG("../rootFiles/anomalies.root", "FIMG_my_tof_amp_cut_para_det1");
    // FIMG_my_tof_amp_cut_para_det2 = retrive_TCutG("../rootFiles/anomalies.root", "FIMG_my_tof_amp_cut_para_det2");
// }

void plot_filterOut_anomaly_PTBC(){

    for (Int_t i = 2; i < 8; i++){
        PTBC_filterOut_anomaly_hists[i-2] = retrive_TH2D_Histograms("../rootFiles/anomalies.root", Form("PTBC_filterOut_anomaly_det%i", i));
    }

    for (Int_t i = 2; i < 8; i++){
        PTBC_tof_amp_cuts[i-2] = retrive_TCutG("../inputFiles/PTBC_cuts.root", Form("tof_amp_cut_det%i", i));
    }

    //Plotting
    SetMArEXStyle();
    // gStyle->SetOptStat(1110);
    gStyle->SetPalette(57);

    TCanvas *c[6];

    int i = 0;

    for (Int_t i = 0; i < 6; i++)
    {
        c[i] = new TCanvas(Form("c_ptbc_filterOut%d", i)," ");
        c[i]->cd();
        PTBC_filterOut_anomaly_hists[i]->GetXaxis()->SetTitle("TOF (in ns)");
        PTBC_filterOut_anomaly_hists[i]->GetYaxis()->SetTitle("Amplitude (a.u.)");
        // PTBC_filterOut_anomaly_hists[i]->GetXaxis()->SetRangeUser(3e2,2e3);
        PTBC_filterOut_anomaly_hists[i]->Draw("COLZ");
        gPad->SetLogx();
        gPad->SetLogz();

        PTBC_tof_amp_cuts[i]->SetLineColor(2);
        PTBC_tof_amp_cuts[i]->SetLineWidth(2);
        PTBC_tof_amp_cuts[i]->Draw("SAME");
    }
}

void plot_al5_anomaly_PTBC(){

    for (Int_t i = 2; i < 8; i++){
        PTBC_al5_anomaly_hists[i-2] = retrive_TH2D_Histograms("../rootFiles/anomalies.root", Form("PTBC_al5_anomaly_det%i", i));
    }

    for (Int_t i = 2; i < 8; i++){
        PTBC_tof_amp_cuts[i-2] = retrive_TCutG("../inputFiles/PTBC_cuts.root", Form("tof_amp_cut_det%i", i));
    }

    //Plotting
    SetMArEXStyle();
    // gStyle->SetOptStat(1110);
    gStyle->SetPalette(57);

    TCanvas *c[6];

    int i = 0;

    for (Int_t i = 0; i < 6; i++)
    {
        c[i] = new TCanvas(Form("c_ptbc_al5%d", i)," ");
        c[i]->cd();
        PTBC_al5_anomaly_hists[i]->GetXaxis()->SetTitle("TOF (in ns)");
        PTBC_al5_anomaly_hists[i]->GetYaxis()->SetTitle("Amplitude (a.u.)");
        // PTBC_al5_anomaly_hists[i]->GetXaxis()->SetRangeUser(3e2,2e3);
        PTBC_al5_anomaly_hists[i]->Draw("COLZ");
        gPad->SetLogx();
        gPad->SetLogz();

        PTBC_tof_amp_cuts[i]->SetLineColor(2);
        PTBC_tof_amp_cuts[i]->SetLineWidth(2);
        PTBC_tof_amp_cuts[i]->Draw("SAME");
    }
}

void anomalies_plots(){

    // retrive_all_hists_cuts();

    // plot_al5_anomaly_PTBC();

    plot_filterOut_anomaly_PTBC();
}

