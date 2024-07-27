/**
 * @file transfer_function.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * File to take into account the moderator effect and calculated the transfer function between tof and energy
 * @version 0.1
 * @date 2024-04-23
 */

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <array>

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
#include "TKey.h"
#include "TList.h"

#include "MArEXStyle.C"

TH2D* rf_hist = 0;
TH1D* rf_hist_peak = 0;
TH1D* rf_hist_mean = 0;

TH1D* transfer_func_PTBC_noRF = 0;
TH1D* transfer_func_FIMG_noRF = 0;
TH1D* transfer_func_mean_PTBC_5itr = 0;
// TH1D* transfer_func_mean_PTBC_10itr = 0;
// TH1D* transfer_func_mean_PTBC_15itr = 0;
TH1D* transfer_func_mean_FIMG_5itr = 0;
// TH1D* transfer_func_mean_FIMG_10itr = 0;

TH1D* transfer_func_peak_PTBC_5itr = 0;
// TH1D* transfer_func_peak_PTBC_10itr = 0;

Double_t min_tof = 800.0;
Double_t flight_path_length_PTBC = 182.65 - 0.41; //m
Double_t flight_path_length_FIMG = 183.5 - 0.41; //m
Double_t neutron_mass = 939.56542052; //in MeV
Double_t speed_of_light = 299792458.0; //in m/s

Double_t TOFToEnergy(Double_t t, Double_t flight_path_length, Double_t rf_length){ //t is in seconds, rf_length is in m
    Double_t denom_term = (flight_path_length + rf_length)/(speed_of_light * t);
    Double_t denominator = 1.0 - (denom_term * denom_term);
    Double_t factor = std::sqrt(1.0 / denominator) - 1.0;
    Double_t energy = neutron_mass * factor;
    return energy * 1e6;  //e in eV
}

Double_t get_rf_length_mean(Double_t n_energy){ // energy in eV;
    Int_t e_bin_num = rf_hist->GetXaxis()->FindBin(n_energy);
    std::string projection_name = "profile_" + std::to_string(n_energy);
    TH1D* projection = (TH1D*)rf_hist->ProjectionY(
        projection_name.c_str(),
        e_bin_num, e_bin_num
    );
    // Double_t fwhm = FindFWHM(projection); //in cm
    Double_t rf_length = projection->GetMean(1) * 0.01; //converting to m
    
    return rf_length; //in m
}

Double_t get_rf_length_peak(Double_t n_energy){ // energy in eV;
    Int_t e_bin_num = rf_hist->GetXaxis()->FindBin(n_energy);
    std::string projection_name = "profile_" + std::to_string(n_energy);
    TH1D* projection = (TH1D*)rf_hist->ProjectionY(
        projection_name.c_str(),
        e_bin_num, e_bin_num
    );

    Int_t peak_bin_num = projection->GetMaximumBin();
    Double_t rf_length = projection->GetBinCenter(peak_bin_num) * 0.01; //converting to m
    
    return rf_length; //in m
}

void get_rf_hists(){

    //opening the rf hist file
    TFile *rfFile = TFile::Open("../inputFiles/RF.root", "READ");
    rf_hist = (TH2D*)rfFile->Get("histfluka");

    Int_t num_bins_e_rf = rf_hist->GetNbinsX();
    rf_hist_peak = (TH1D*)rf_hist->ProjectionX("RF Hist Peak", 10, 10);
    rf_hist_mean = (TH1D*)rf_hist->ProjectionX("RF Hist Mean", 15, 15);

    // Filling peak and mean histograms
    for(Int_t i = 0; i < num_bins_e_rf; i++){

        TH1D* y_projection = (TH1D*)rf_hist->ProjectionY("Y projection", i+1, i+1);
        rf_hist_mean->SetBinContent(i+1, y_projection->GetMean(1));
        Int_t peak_bin_num = y_projection->GetMaximumBin();
        Double_t peak_value = y_projection->GetBinCenter(peak_bin_num);
        rf_hist_peak->SetBinContent(i+1, peak_value);
    }
}

void calculate_transfer_function(){

    // tof bins - 1000 bins per decade
    // Calculating TOF (x) bin edges
    Int_t bins_per_decade = 1000;
    Int_t Num_decades = 6;
    Int_t num_bins_tof = bins_per_decade * Num_decades;
    Double_t bin_edges_tof[num_bins_tof+1];
    Double_t step_tof = ((Double_t) 1.0/(Double_t) bins_per_decade);
    for(Int_t i = 0; i < num_bins_tof+1; i++)
    {
        Double_t base = 10.;
        Double_t exponent = (step_tof * (Double_t) i) + 2.;
        bin_edges_tof[i] = (Double_t) std::pow(base, exponent);
    }

    get_rf_hists();

    // Initializing histograms
    transfer_func_PTBC_noRF = new TH1D("transfer_func_PTBC_noRF","ToF vs Energy - Transfer Function - PTBC - no RF",num_bins_tof,bin_edges_tof);
    transfer_func_FIMG_noRF = new TH1D("transfer_func_FIMG_noRF","ToF vs Energy - Transfer Function - FIMG - no RF",num_bins_tof,bin_edges_tof);
    transfer_func_mean_PTBC_5itr = new TH1D("transfer_func_mean_PTBC_5itr","ToF vs Energy - Transfer Function (Mean) - PTBC - 5 Itr",num_bins_tof,bin_edges_tof);
    // transfer_func_mean_PTBC_10itr = new TH1D("transfer_func_mean_PTBC_10itr","ToF vs Energy - Transfer Function (Mean) - PTBC - 10 Itr",num_bins_tof,bin_edges_tof);
    // transfer_func_mean_PTBC_15itr = new TH1D("transfer_func_mean_PTBC_15itr","ToF vs Energy - Transfer Function (Mean) - PTBC - 15 Itr",num_bins_tof,bin_edges_tof);
    transfer_func_mean_FIMG_5itr = new TH1D("transfer_func_mean_FIMG_5itr","ToF vs Energy - Transfer Function (Mean) - FIMG - 5 Itr",num_bins_tof,bin_edges_tof);
    // transfer_func_mean_FIMG_10itr = new TH1D("transfer_func_mean_FIMG_10itr","ToF vs Energy - Transfer Function (Mean) - FIMG - 10 Itr",num_bins_tof,bin_edges_tof);
    
    transfer_func_peak_PTBC_5itr = new TH1D("transfer_func_peak_PTBC_5itr","ToF vs Energy - Transfer Function (Peak) - PTBC - 5 Itr",num_bins_tof,bin_edges_tof);
    // transfer_func_peak_PTBC_10itr = new TH1D("transfer_func_peak_PTBC_10itr","ToF vs Energy - Transfer Function (Peak) - PTBC - 10 Itr",num_bins_tof,bin_edges_tof);

    for(Int_t i = 0; i < num_bins_tof; i++)
    {
        Double_t tof = transfer_func_PTBC_noRF->GetBinCenter(i+1);
        if (tof < 800.)
        {
            transfer_func_PTBC_noRF->SetBinContent(i+1, 0);
            transfer_func_FIMG_noRF->SetBinContent(i+1, 0);
            transfer_func_mean_PTBC_5itr->SetBinContent(i+1, 0);
            // transfer_func_mean_PTBC_10itr->SetBinContent(i+1, 0);
            // transfer_func_mean_PTBC_15itr->SetBinContent(i+1, 0);
            transfer_func_mean_FIMG_5itr->SetBinContent(i+1, 0);
            // transfer_func_mean_FIMG_10itr->SetBinContent(i+1, 0);

            transfer_func_peak_PTBC_5itr->SetBinContent(i+1, 0);
            // transfer_func_peak_PTBC_10itr->SetBinContent(i+1, 0);
            continue;
        }

        Double_t rf_len_mean_PTBC = 0.;
        Double_t rf_len_peak_PTBC = 0.;
        Double_t neutron_e_mean_PTBC = TOFToEnergy(tof * 1e-9, flight_path_length_PTBC, rf_len_mean_PTBC);
        Double_t neutron_e_peak_PTBC = neutron_e_mean_PTBC;
        Double_t rf_len_mean_FIMG = 0.;
        Double_t neutron_e_mean_FIMG = TOFToEnergy(tof * 1e-9, flight_path_length_FIMG, rf_len_mean_FIMG);
        transfer_func_PTBC_noRF->SetBinContent(i+1, neutron_e_mean_PTBC);
        transfer_func_FIMG_noRF->SetBinContent(i+1, neutron_e_mean_FIMG);
        for (Int_t j = 0; j < 5; j++)
        {
            
            rf_len_mean_PTBC = get_rf_length_mean(neutron_e_mean_PTBC);
            neutron_e_mean_PTBC = TOFToEnergy(tof * 1e-9, flight_path_length_PTBC, rf_len_mean_PTBC);
            rf_len_peak_PTBC = get_rf_length_peak(neutron_e_peak_PTBC);
            neutron_e_peak_PTBC = TOFToEnergy(tof * 1e-9, flight_path_length_PTBC, rf_len_peak_PTBC);

            rf_len_mean_FIMG = get_rf_length_mean(neutron_e_mean_FIMG);
            neutron_e_mean_FIMG = TOFToEnergy(tof * 1e-9, flight_path_length_FIMG, rf_len_mean_FIMG);

            if (j == 4)
            {
                transfer_func_mean_PTBC_5itr->SetBinContent(i+1, neutron_e_mean_PTBC);
                transfer_func_mean_FIMG_5itr->SetBinContent(i+1, neutron_e_mean_FIMG);
                transfer_func_peak_PTBC_5itr->SetBinContent(i+1, neutron_e_peak_PTBC);
            }

            // if (j == 9)
            // {
            //     transfer_func_mean_PTBC_10itr->SetBinContent(i+1, neutron_e_mean);
            //     transfer_func_peak_PTBC_10itr->SetBinContent(i+1, neutron_e_peak);
            // }
            
            // if (j == 14)
            // {
            //     transfer_func_mean_PTBC_15itr->SetBinContent(i+1, neutron_e_mean);
            // }
        }  
    }
}

void store_hist(){

    // TFile *output_file = new TFile("../inputFiles/transfer_function.root","recreate");
    TFile *output_file = TFile::Open("../inputFiles/transfer_function.root", "recreate");

    output_file->WriteObject(transfer_func_PTBC_noRF, "transfer_func_PTBC_noRF");
    output_file->WriteObject(transfer_func_FIMG_noRF, "transfer_func_FIMG_noRF");
    output_file->WriteObject(transfer_func_mean_PTBC_5itr, "transfer_func_mean_PTBC_5itr");
    // output_file->WriteObject(transfer_func_mean_PTBC_10itr, "transfer_func_mean_PTBC_10itr");
    // output_file->WriteObject(transfer_func_mean_PTBC_15itr, "transfer_func_mean_PTBC_15itr");
    output_file->WriteObject(transfer_func_mean_FIMG_5itr, "transfer_func_mean_FIMG_5itr");
    output_file->WriteObject(transfer_func_peak_PTBC_5itr, "transfer_func_peak_PTBC_5itr");
    // output_file->WriteObject(transfer_func_peak_PTBC_10itr, "transfer_func_peak_PTBC_10itr");
    output_file->Close();

    std::cout << "Created output file 'transfer_function.root'" << std::endl;
}

void plot_hists(){

    TFile *tfFile = TFile::Open("../inputFiles/transfer_function.root", "READ");
    transfer_func_PTBC_noRF = (TH1D*)tfFile->Get("transfer_func_PTBC_noRF");
    transfer_func_FIMG_noRF = (TH1D*)tfFile->Get("transfer_func_FIMG_noRF");
    transfer_func_mean_PTBC_5itr = (TH1D*)tfFile->Get("transfer_func_mean_PTBC_5itr");
    transfer_func_mean_FIMG_5itr = (TH1D*)tfFile->Get("transfer_func_mean_FIMG_5itr");
    transfer_func_peak_PTBC_5itr = (TH1D*)tfFile->Get("transfer_func_peak_PTBC_5itr");

    get_rf_hists();

    //Plotting
    SetMArEXStyle();
    gStyle->SetCanvasDefW(800); //600
    gStyle->SetCanvasDefH(400); //500 

    gStyle->SetStatX(0.27);
    gStyle->SetStatY(0.9);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.17);

    TCanvas *c[3];
    TPad *p[3][2];
    TLegend *l[3];

    int i = 0;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();
    c[i]->Draw();
    // c[i]->SetLeftMargin(0.1);
    // c[i]->SetRightMargin(0.05);
    
    //Pad 1
    c[i]->cd(0);
    p[i][0] = new TPad(Form("p_%d_0", i), Form("p_%d_0", i), 0., 0.35, 1., 1.);
    p[i][0]->SetFillColor(kWhite);
    p[i][0]->SetBottomMargin(0.00001);
    p[i][0]->SetRightMargin(0.05);
    p[i][0]->SetBorderMode(0);
    p[i][0]->Draw();
    p[i][0]->cd();

    // Y Axis
    transfer_func_PTBC_noRF->GetYaxis()->SetTitle("Energy (in eV)");
    transfer_func_PTBC_noRF->GetYaxis()->SetLabelSize(0.05);
    transfer_func_PTBC_noRF->GetYaxis()->SetTitleSize(0.06);
    transfer_func_PTBC_noRF->GetYaxis()->SetTitleOffset(0.65);
    // X Axis
    transfer_func_PTBC_noRF->GetXaxis()->SetLabelOffset(999);
    transfer_func_PTBC_noRF->GetXaxis()->SetLabelSize(0);
    
    transfer_func_PTBC_noRF->SetTitle("Transfer Function - PTBC");
    transfer_func_PTBC_noRF->SetLineWidth(2);
    transfer_func_PTBC_noRF->Draw(); //"HISTE"
    transfer_func_PTBC_noRF->SetStats(0);
    gPad->SetLogx();
    gPad->SetLogy();

    l[i] = new TLegend(0.45,0.65,0.85,0.85); //0.68,0.7,0.86,0.8       ;         0.72,0.8,0.90,0.9
    l[i]->AddEntry(transfer_func_PTBC_noRF,"Without Moderator Effects","l");

    l[i]->AddEntry(transfer_func_mean_PTBC_5itr,"With Moderator Effects","l");
    transfer_func_mean_PTBC_5itr->SetLineColor(6);
    transfer_func_mean_PTBC_5itr->SetLineWidth(2);
    transfer_func_mean_PTBC_5itr->Draw("SAME");

    l[i]->SetMargin(0.4);
    l[i]->Draw();

    //Pad 2
    c[i]->cd(0);
    p[i][1] = new TPad(Form("p_%d_1", i), Form("p_%d_1", i), 0., 0., 1., 0.35);
    p[i][1]->SetFillColor(kWhite);
    p[i][1]->SetTopMargin(0.00001);
    p[i][1]->SetBottomMargin(0.25);
    p[i][1]->SetRightMargin(0.05);
    p[i][1]->SetBorderMode(0);
    p[i][1]->Draw();
    p[i][1]->cd();

    TH1D* residual_plot_PTBC = (TH1D*)(transfer_func_mean_PTBC_5itr->Clone("residual_plot_PTBC"));
    residual_plot_PTBC->Add(transfer_func_PTBC_noRF, -1);
    residual_plot_PTBC->Divide(transfer_func_PTBC_noRF);
    residual_plot_PTBC->Scale(100., "nosw2");
    residual_plot_PTBC->SetTitle("");

    // X Axis
    residual_plot_PTBC->GetXaxis()->SetTitle("ToF (in ns)");
    residual_plot_PTBC->GetXaxis()->SetTitleSize(0.1);
    residual_plot_PTBC->GetXaxis()->SetLabelOffset(0.01);
    residual_plot_PTBC->GetXaxis()->SetTitleOffset(0.85); //decrease to move up
    residual_plot_PTBC->GetXaxis()->SetLabelSize(0.07);
    // Y Axis
    residual_plot_PTBC->GetYaxis()->SetTitle("\% Difference");
    residual_plot_PTBC->GetYaxis()->SetTitleSize(0.1);
    residual_plot_PTBC->GetYaxis()->SetTitleOffset(0.38);
    residual_plot_PTBC->GetYaxis()->SetLabelSize(0.07);
    residual_plot_PTBC->Draw();
    gPad->SetLogx();

    i++;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();
    c[i]->SetTopMargin(0.05);
    c[i]->SetBottomMargin(0.1);
    c[i]->SetRightMargin(0.1);
    gStyle->SetPalette(57);
    rf_hist->GetXaxis()->SetTitle("Neutron Energy (eV)");
    rf_hist->GetYaxis()->SetTitle("Neutron Moderation Distance (cm)");
    rf_hist->SetTitle("");
    rf_hist->GetXaxis()->SetTitleOffset(1.2);
    rf_hist->GetYaxis()->SetRangeUser(-100., 700.);
    rf_hist->Draw("colz");
    gPad->SetLogx();

    l[i] = new TLegend(0.45,0.13,0.85,0.18);

    // l[i]->AddEntry(rf_hist_peak,"Peak","l");
    // rf_hist_peak->SetLineColor(1);
    // rf_hist_peak->SetLineWidth(2);
    // rf_hist_peak->Draw("SAME");

    l[i]->AddEntry(rf_hist_mean,"Mean moderation distance","l");
    rf_hist_mean->SetLineColor(2);
    rf_hist_mean->SetLineWidth(2);
    rf_hist_mean->Draw("SAME");

    l[i]->SetMargin(0.4);
    l[i]->Draw();

    i++;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();
    c[i]->Draw();
    // c[i]->SetLeftMargin(0.1);
    // c[i]->SetRightMargin(0.05);
    
    //Pad 1
    c[i]->cd(0);
    p[i][0] = new TPad(Form("p_%d_0", i), Form("p_%d_0", i), 0., 0.35, 1., 1.);
    p[i][0]->SetFillColor(kWhite);
    p[i][0]->SetBottomMargin(0.00001);
    p[i][0]->SetRightMargin(0.05);
    p[i][0]->SetBorderMode(0);
    p[i][0]->Draw();
    p[i][0]->cd();

    // Y Axis
    transfer_func_FIMG_noRF->GetYaxis()->SetTitle("Energy (in eV)");
    transfer_func_FIMG_noRF->GetYaxis()->SetLabelSize(0.05);
    transfer_func_FIMG_noRF->GetYaxis()->SetTitleSize(0.06);
    transfer_func_FIMG_noRF->GetYaxis()->SetTitleOffset(0.65);
    // X Axis
    transfer_func_FIMG_noRF->GetXaxis()->SetLabelOffset(999);
    transfer_func_FIMG_noRF->GetXaxis()->SetLabelSize(0);
    
    transfer_func_FIMG_noRF->SetTitle("Transfer Function - FIMG");
    transfer_func_FIMG_noRF->SetLineWidth(2);
    transfer_func_FIMG_noRF->Draw(); //"HISTE"
    transfer_func_FIMG_noRF->SetStats(0);
    gPad->SetLogx();
    gPad->SetLogy();

    l[i] = new TLegend(0.45,0.65,0.85,0.85); //0.68,0.7,0.86,0.8       ;         0.72,0.8,0.90,0.9
    l[i]->AddEntry(transfer_func_FIMG_noRF,"Without Moderator Effects","l");

    l[i]->AddEntry(transfer_func_mean_FIMG_5itr,"With Moderator Effects","l");
    transfer_func_mean_FIMG_5itr->SetLineColor(2);
    transfer_func_mean_FIMG_5itr->SetLineWidth(2);
    transfer_func_mean_FIMG_5itr->Draw("SAME");

    l[i]->SetMargin(0.4);
    l[i]->Draw();

    //Pad 2
    c[i]->cd(0);
    p[i][1] = new TPad(Form("p_%d_1", i), Form("p_%d_1", i), 0., 0., 1., 0.35);
    p[i][1]->SetFillColor(kWhite);
    p[i][1]->SetTopMargin(0.00001);
    p[i][1]->SetBottomMargin(0.25);
    p[i][1]->SetRightMargin(0.05);
    p[i][1]->SetBorderMode(0);
    p[i][1]->Draw();
    p[i][1]->cd();

    TH1D* residual_plot_FIMG = (TH1D*)(transfer_func_mean_FIMG_5itr->Clone("residual_plot_FIMG"));
    residual_plot_FIMG->Add(transfer_func_FIMG_noRF, -1);
    residual_plot_FIMG->Divide(transfer_func_FIMG_noRF);
    residual_plot_FIMG->Scale(100., "nosw2");
    residual_plot_FIMG->SetTitle("");

    // X Axis
    residual_plot_FIMG->GetXaxis()->SetTitle("ToF (in ns)");
    residual_plot_FIMG->GetXaxis()->SetTitleSize(0.1);
    residual_plot_FIMG->GetXaxis()->SetLabelOffset(0.01);
    residual_plot_FIMG->GetXaxis()->SetTitleOffset(0.85); //decrease to move up
    residual_plot_FIMG->GetXaxis()->SetLabelSize(0.07);
    // Y Axis
    residual_plot_FIMG->GetYaxis()->SetTitle("\% Difference");
    residual_plot_FIMG->GetYaxis()->SetTitleSize(0.1);
    residual_plot_FIMG->GetYaxis()->SetTitleOffset(0.38);
    residual_plot_FIMG->GetYaxis()->SetLabelSize(0.07);
    residual_plot_FIMG->Draw();
    gPad->SetLogx();
}

void transfer_function(){

    

    // calculate_transfer_function();

    // store_hist();
    plot_hists();
}