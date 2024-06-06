/**
 * @file cutoffFitter_PTBC.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-11-13
 * @date 2024-05-21 (name change)
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

const std::string target_name("none"); ////bi1, al3, al5, al8, c1p2_ts, al5_ts, bi1p2_ts, cf_bottle, cf_bottle_rot, cf_bottle_rotBack, ar_bottle_full, bi1sep17
//no targets - none, none_ts, ar_bottle_empty
const std::string target_name_title("No Target");
//Bi (1 cm), Target Bi (1.2 cm), Al (3 cm), Al (5 cm), Target Al (5 cm), Al (8 cm), Target C (1.2 cm), Empty Bottle, Empty Bottle Rotated
//Argon Tank, Bi (1 cm) - Sep 17, No Target, SCUBA Tank, Argon Tank Empty

Double_t min_tof = 800.0; //ns
Double_t max_tof = 1e8; //ns
Double_t max_amp = 50000.0;

Int_t bins_per_decade_cuts = 5;

TH2D* PTBC_tof_amp_hists[6];
TH1D* projection_hists[6];
TF1* alphas_fits_total[6];
TH1D* det_cut_hists[6];
Double_t det_alphas_cuts[6];
TCutG* det_cuts[6];

TCanvas *alphas_canvas[6];
TLegend *alphas_legend[6];
Int_t alphas_plot_index = 0;

TCanvas *cuts_canvas[6];
Int_t cuts_plot_index = 0;

Double_t expo_fit_ranges[6][2];
Double_t gaus_fit_ranges[6][2];
Double_t total_fit_ranges[6][2];

void fill_fit_ranges(){

    //[][0] - Lower range value
    //[][1] - Upper range value

    //det2
    expo_fit_ranges[0][0] = 725;
    expo_fit_ranges[0][1] = 1100;
    gaus_fit_ranges[0][0] = 1300;
    gaus_fit_ranges[0][1] = 2800;
    total_fit_ranges[0][0] = 725;
    total_fit_ranges[0][1] = 3000;

    //det3
    expo_fit_ranges[1][0] = 800; //525
    expo_fit_ranges[1][1] = 1200; //650
    gaus_fit_ranges[1][0] = 1300;
    gaus_fit_ranges[1][1] = 2500;
    total_fit_ranges[1][0] = 900; //525
    total_fit_ranges[1][1] = 2800;

    //det4
    expo_fit_ranges[2][0] = 600;
    expo_fit_ranges[2][1] = 1000;
    gaus_fit_ranges[2][0] = 1500;
    gaus_fit_ranges[2][1] = 2500;
    total_fit_ranges[2][0] = 550;
    total_fit_ranges[2][1] = 2800;

    //det5
    expo_fit_ranges[3][0] = 815;
    expo_fit_ranges[3][1] = 950;
    gaus_fit_ranges[3][0] = 1500;
    gaus_fit_ranges[3][1] = 2500;
    total_fit_ranges[3][0] = 815;
    total_fit_ranges[3][1] = 2800;

    //det6
    expo_fit_ranges[4][0] = 815;
    expo_fit_ranges[4][1] = 1100;
    gaus_fit_ranges[4][0] = 1700;
    gaus_fit_ranges[4][1] = 2700;
    total_fit_ranges[4][0] = 815;
    total_fit_ranges[4][1] = 3380;

    //det7
    expo_fit_ranges[5][0] = 500;
    expo_fit_ranges[5][1] = 900;
    gaus_fit_ranges[5][0] = 1300;
    gaus_fit_ranges[5][1] = 2300;
    total_fit_ranges[5][0] = 500;
    total_fit_ranges[5][1] = 2400;

    return;
}

TH2D* retriveHistograms(const char *file_name, const char *hist_name){
    
    TFile* hist_file = TFile::Open(file_name, "READ");
    TH2D* hist_new = (TH2D*)hist_file->Get(hist_name);

    return hist_new;
}

Double_t exponential(Double_t *x, Double_t *par) {
    return TMath::Exp(par[0] + par[1] * x[0]);
}

Double_t gaussian(Double_t *x, Double_t *par) {
    return par[0] * TMath::Exp(-0.5 * ((x[0] - par[1])/par[2]) * ((x[0] - par[1])/par[2]));
    // return par[0] * TMath::Gaus(x[0], par[1], par[2]);
}

Double_t fitFunction(Double_t *x, Double_t *par){
    return exponential(x, par) + gaussian(x, &par[2]);
}

// void convert_hist_to_TCutG(Int_t det_num, TH1D* cut_hist){

//     const Int_t nBins = cut_hist->GetNbinsX();
//     const Int_t nPoints = nBins * 2; // Two points per bin (for the top and bottom of the bin)
    
//     Double_t *x = new Double_t[nPoints];
//     Double_t *y = new Double_t[nPoints];
    
//     Int_t pointIndex = 0;
    
//     for (Int_t i = 1; i <= nBins; ++i) {
//         Double_t binCenter = cut_hist->GetBinCenter(i);
//         Double_t binContent = cut_hist->GetBinContent(i);
        
//         if (binContent > 0) {
//             x[pointIndex] = binCenter;
//             y[pointIndex] = 0; // Bottom of the bin
//             ++pointIndex;

//             x[pointIndex] = binCenter;
//             y[pointIndex] = binContent; // Top of the bin
//             ++pointIndex;
//         }
//     }
// }

// cuts for tof < 10^4
void determine_gamma_flash_cuts(Int_t det_num, TH1D* cut_hist){

    TH2D* PTBC_tof_amp_hist_forCuts;

    PTBC_tof_amp_hist_forCuts = (TH2D*)PTBC_tof_amp_hists[det_num-2]->Rebin2D((Int_t) 1000/bins_per_decade_cuts, 5, Form("PTBC_tof_amp_forCuts_det%i", det_num));

    Int_t num_tof_bins = PTBC_tof_amp_hist_forCuts->GetNbinsX();

    for (Int_t i = 1; i <= num_tof_bins; i++) 
    {
        if (i < 5) //starting from xbin = 5 (around tof = 800ns)
        {
            cut_hist->SetBinContent(i, 50000.);
            continue;
        }

        if (i >= 5 && i <= 10) // For tof below 10^4 ns
        {
            std::string projection_name = "profile_bin_" + std::to_string(i);
            TH1D* proj_hist = (TH1D*)PTBC_tof_amp_hist_forCuts->ProjectionY(projection_name.c_str(),i, i);

            TF1 *gaus_fit_total = new TF1("gaus_fit_total", "gaus", 0, 10000);
            // gaus_fit_total->SetParameters(par);
            proj_hist->Fit(gaus_fit_total, "0R");

            Double_t mean_val = gaus_fit_total->GetParameter(1);
            Double_t std_dev = gaus_fit_total->GetParameter(2);
            
            Double_t cut_val = 0;

            //For detectors 4 and 7, the bins closer to the gamma flash are cut harder (6 sigma)
            if (det_num == 4 || det_num == 7)
            {
                if (i <= 6){
                    cut_val = mean_val + 6*std_dev;
                } else {
                    cut_val = mean_val + 5*std_dev;
                }
            } else {
                cut_val = mean_val + 5*std_dev;
            }

            if (cut_val > det_alphas_cuts[det_num-2])
            {
                cut_hist->SetBinContent(i, cut_val);
            } else {
                cut_hist->SetBinContent(i, det_alphas_cuts[det_num-2]);
            }

            // For bins 5, 6, and 7, we are choosing the cut that is of higher value and setting it equal to all three bins (and bin 8 is also included for detector 6)
            if (i == 8)
            {   
                Double_t bin_7_val = cut_hist->GetBinContent(7);
                Double_t bin_6_val = cut_hist->GetBinContent(6);
                Double_t bin_5_val = cut_hist->GetBinContent(5);

                Double_t max_cut = 0.;

                if (det_num == 6)
                {
                    Double_t bin_8_val = cut_hist->GetBinContent(8);
                    max_cut = std::max({bin_8_val, bin_7_val, bin_6_val, bin_5_val});
                    cut_hist->SetBinContent(5, max_cut);
                    cut_hist->SetBinContent(6, max_cut);
                    cut_hist->SetBinContent(7, max_cut);
                    cut_hist->SetBinContent(8, max_cut);
                } else {
                    max_cut = std::max({bin_7_val, bin_6_val, bin_5_val});
                    cut_hist->SetBinContent(5, max_cut);
                    cut_hist->SetBinContent(6, max_cut);
                    cut_hist->SetBinContent(7, max_cut);
                }
            }   

            continue;
        }

        if (i > 10){ //40 // For tof above 10^4 ns
            cut_hist->SetBinContent(i, det_alphas_cuts[det_num-2]);
            continue;
        }
    }
}

// cuts for tof > 10^4
void determine_alpha_cut(Int_t det_num, TH2D* tof_amp_hist, TH1D* projection_hist, TF1* total_fit) {

    // projection_hist = (TH1D*)tof_amp_hist->ProjectionY(Form("profile_fission_alphas_det%i", det_num),5001, 6000);

    TF1 *expon = new TF1("expon", "expo", expo_fit_ranges[det_num-2][0], expo_fit_ranges[det_num-2][1]);
    TF1 *gauss = new TF1("gauss", "gaus", gaus_fit_ranges[det_num-2][0], gaus_fit_ranges[det_num-2][1]);

    Double_t par[5];

    projection_hist->Fit(expon,"0R");
    projection_hist->Fit(gauss,"0R+");
    
    expon->GetParameters(&par[0]);
    gauss->GetParameters(&par[2]);

    total_fit->SetParameters(par);
    projection_hist->Fit(total_fit, "0R+");

    det_alphas_cuts[det_num-2] = total_fit->GetParameter(3) + 5 * total_fit->GetParameter(4);

    return;
}

void plot_alpha_cuts(Int_t det_num, TH1D* projection_hist, TF1* total_fit) {

    //Plotting
    SetMArEXStyle();
    
    gStyle->SetStatX(0.27);
    gStyle->SetStatY(0.9);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.17);

    alphas_canvas[alphas_plot_index] = new TCanvas(Form("alphas_c_%i", alphas_plot_index)," ");
    alphas_canvas[alphas_plot_index]->cd();

    alphas_legend[alphas_plot_index] = new TLegend(0.7,0.65,0.86,0.85);

    projection_hist->GetXaxis()->SetTitle("Amplitude (a.u.)");
    projection_hist->GetYaxis()->SetTitle("Number of Events");
    projection_hist->SetTitle(Form("Num Events vs Amp within 10^{7} - 10^{8} ToF - Det %i", det_num));
    projection_hist->GetXaxis()->SetRangeUser(400, 3500);
    projection_hist->Draw();
    alphas_legend[alphas_plot_index]->AddEntry(projection_hist,"Data","l");

    ////////// Drawing individual fit
    // Exponential
    Double_t exp_par0 = total_fit->GetParameter(0);
    Double_t exp_par1 = total_fit->GetParameter(1);
    TF1 *expo_fit_plot = new TF1("expo_fit_plot", "exp([0] + [1] * x)", 400, 3500);
    expo_fit_plot->SetParameter(0, exp_par0);
    expo_fit_plot->SetParameter(1, exp_par1);
    expo_fit_plot->SetLineColor(3);
    expo_fit_plot->Draw("SAME");
    alphas_legend[alphas_plot_index]->AddEntry(expo_fit_plot,"Exponential","l");

    // Gaussian
    Double_t constant = total_fit->GetParameter(2);
    Double_t mean = total_fit->GetParameter(3);
    Double_t sigma = total_fit->GetParameter(4);
    TF1 *gaus_fit_plot = new TF1("gaus_fit_plot", "[0]*exp(-0.5*((x-[1])/[2])**2)", 400, 3500);
    gaus_fit_plot->SetParameter(0, constant);
    gaus_fit_plot->SetParameter(1, mean);
    gaus_fit_plot->SetParameter(2, sigma);
    gaus_fit_plot->SetLineColor(2);
    gaus_fit_plot->Draw("SAME");
    alphas_legend[alphas_plot_index]->AddEntry(gaus_fit_plot,"Gaussian","l");

    //Plotting Global Fit
    TF1 *total_fit_plot = new TF1("total_fit_plot", "expo(0)+gaus(2)", 400, 3500);
    total_fit_plot->SetParameter(0, exp_par0);
    total_fit_plot->SetParameter(1, exp_par1);
    total_fit_plot->SetParameter(2, constant);
    total_fit_plot->SetParameter(3, mean);
    total_fit_plot->SetParameter(4, sigma);
    total_fit_plot->SetLineColor(1);
    total_fit_plot->SetLineWidth(3); 
    total_fit_plot->Draw("SAME");
    alphas_legend[alphas_plot_index]->AddEntry(total_fit_plot,"Global Fit","l");

    alphas_legend[alphas_plot_index]->SetMargin(0.4);
    alphas_legend[alphas_plot_index]->Draw();

    alphas_plot_index++;

    return;
}

void plot_det_cuts(Int_t det_num, TH2D* tof_amp_hist, TH1D* cut_hist) {
    
    //Plotting
    SetMArEXStyle();
    
    gStyle->SetStatX(0.27);
    gStyle->SetStatY(0.9);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.17);
    gStyle->SetPalette(57);

    cuts_canvas[cuts_plot_index] = new TCanvas(Form("cuts_c_%i", cuts_plot_index)," ");
    cuts_canvas[cuts_plot_index]->cd();

    tof_amp_hist->GetXaxis()->SetTitle("TOF (in ns)");
    tof_amp_hist->GetYaxis()->SetTitle("Amplitude (a.u.)");
    tof_amp_hist->SetTitle(Form("ToF-Amp Hist - Det %i - %s", det_num, target_name_title.c_str()));
    tof_amp_hist->Draw("COLZ");
    gPad->SetLogx();
    gPad->SetLogz();

    cut_hist->SetLineColor(2);
    cut_hist->SetLineWidth(2);
    cut_hist->Draw("SAME");

    // cuts_canvas[cuts_plot_index]->Print(Form("../plots/cuts_plots/tof_amp_hist_%s_det%i_%ibpd_100Bins.png", target_name.c_str(), det_num, bins_per_decade_cuts));
    
    cuts_plot_index++;
    return;
}

void StoreHist(){
    
    TFile *f = new TFile("../inputFiles/PTBC_cuts.root","recreate");

    for (Int_t i = 0; i < 6; i++)
    {
        det_cut_hists[i]->SetLineColor(2);
        det_cut_hists[i]->SetLineWidth(2);
        det_cut_hists[i]->Write();
    }
    
    f->Close();

    std::cout << "Created output file 'PTBC_cuts.root'" << std::endl;
}

void cutoffFitter_PTBC() {

    fill_fit_ranges();

    //Calculating TOF (x) bin edges FOR CUTS
    Int_t Num_decades = 6;
    Int_t num_bins_tof_cuts = bins_per_decade_cuts * Num_decades;
    Double_t bin_edges_tof_cuts[num_bins_tof_cuts+1];
    Double_t step_tof_cuts = ((Double_t) 1.0/(Double_t) bins_per_decade_cuts);
    for(Int_t i = 0; i < num_bins_tof_cuts+1; i++)
    {
        Double_t base = 10.;
        Double_t exponent = (step_tof_cuts * (Double_t) i) + 2.;
        bin_edges_tof_cuts[i] = (Double_t) std::pow(base, exponent);
    }
    
    PTBC_tof_amp_hists[0] = retriveHistograms("../rootFiles/cutoffAnalysis_PTBC_ar_bottle_empty.root", "PTBC_tof_amp_det2");
    PTBC_tof_amp_hists[1] = retriveHistograms("../rootFiles/cutoffAnalysis_PTBC_ar_bottle_empty.root", "PTBC_tof_amp_det3");
    PTBC_tof_amp_hists[2] = retriveHistograms("../rootFiles/cutoffAnalysis_PTBC_ar_bottle_empty.root", "PTBC_tof_amp_det4");
    PTBC_tof_amp_hists[3] = retriveHistograms("../rootFiles/cutoffAnalysis_PTBC_ar_bottle_empty.root", "PTBC_tof_amp_det5");
    PTBC_tof_amp_hists[4] = retriveHistograms("../rootFiles/cutoffAnalysis_PTBC_ar_bottle_empty.root", "PTBC_tof_amp_det6");
    PTBC_tof_amp_hists[5] = retriveHistograms("../rootFiles/cutoffAnalysis_PTBC_ar_bottle_empty.root", "PTBC_tof_amp_det7");

    for (Int_t i = 0; i < 6; i++)
    {
        projection_hists[i] = (TH1D*)PTBC_tof_amp_hists[i]->ProjectionY(Form("profile_fission_alphas_det%i", i+2), 5001, 6000);
        alphas_fits_total[i] = new TF1(Form("alphas_total_fit_det%i", i+2), "expo(0)+gaus(2)", total_fit_ranges[i][0], total_fit_ranges[i][1]);
        determine_alpha_cut(i+2, PTBC_tof_amp_hists[i], projection_hists[i], alphas_fits_total[i]);
        // plot_alpha_cuts(i+2, projection_hists[i], alphas_fits_total[i]);

        det_cut_hists[i] = new TH1D(Form("PTBC_cuts_det%i", i+2), Form("ToF-Amp cut Hist - PTBC Det %i - %s", i+2, target_name_title.c_str()), num_bins_tof_cuts, bin_edges_tof_cuts);
        determine_gamma_flash_cuts(i+2, det_cut_hists[i]);
        // plot_det_cuts(i+2, PTBC_tof_amp_hists[i], det_cut_hists[i]);
    }

    StoreHist();

}
