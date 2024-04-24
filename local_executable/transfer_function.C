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
TH1D* transfer_func_hist_PTBC_noRF = 0;
TH1D* transfer_func_hist_PTBC_5itr = 0;
TH1D* transfer_func_hist_PTBC_10itr = 0;
TH1D* transfer_func_hist_PTBC_15itr = 0;

Double_t min_tof = 800.0;
Double_t flight_path_length_PTBC = 182.65 - 0.41; //m
Double_t neutron_mass = 939.56542052; //in MeV
Double_t speed_of_light = 299792458.0; //in m/s

Double_t TOFToEnergy(Double_t t, Double_t flight_path_length, Double_t rf_length){ //t is in seconds, rf_length is in m
    Double_t denom_term = (flight_path_length + rf_length)/(speed_of_light * t);
    Double_t denominator = 1.0 - (denom_term * denom_term);
    Double_t factor = std::sqrt(1.0 / denominator) - 1.0;
    Double_t energy = neutron_mass * factor;
    return energy * 1e6;  //e in eV
}

Double_t get_rf_length(Double_t n_energy){ // energy in eV
    Int_t e_bin_num = rf_hist->GetXaxis()->FindBin(n_energy);
    std::string projection_name = "profile_" + std::to_string(n_energy);
    TH1D* projection = (TH1D*)rf_hist->ProjectionY(
        projection_name.c_str(),
        e_bin_num, e_bin_num
    );
    // Double_t fwhm = FindFWHM(projection); //in cm
    Double_t rf_length = projection->GetMean(1) * 0.01; //projection->GetBinCenter( projection->GetMaximumBin() ) * 0.01; //converting to m
    
    return rf_length; //in m
}

void store_hist(){

    // TFile *output_file = new TFile("../inputFiles/transfer_function.root","recreate");
    TFile *output_file = TFile::Open("../inputFiles/transfer_function.root", "recreate");

    output_file->WriteObject(transfer_func_hist_PTBC_noRF, "transfer_func_hist_PTBC_noRF");

    // cout << "---------3----------" << endl;

    output_file->WriteObject(transfer_func_hist_PTBC_5itr, "transfer_func_hist_PTBC_5itr");

    // cout << "---------4----------" << endl;

    output_file->WriteObject(transfer_func_hist_PTBC_10itr, "transfer_func_hist_PTBC_10itr");

    // cout << "---------5----------" << endl;

    output_file->WriteObject(transfer_func_hist_PTBC_15itr, "transfer_func_hist_PTBC_15itr");

    // cout << "---------6----------" << endl;

    output_file->Close();

    std::cout << "Created output file 'transfer_function.root'" << std::endl;
}

void plot_hists(){

    //Plotting
    SetMArEXStyle();

    gStyle->SetStatX(0.27);
    gStyle->SetStatY(0.9);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.17);

    TCanvas *c[1];
    TLegend *l[1];

    int i = 0;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();

    l[i] = new TLegend(0.77,0.7,0.86,0.85); //0.68,0.7,0.86,0.8       ;         0.72,0.8,0.90,0.9
    l[i]->AddEntry(transfer_func_hist_PTBC_noRF,"no RF","l");
    
    transfer_func_hist_PTBC_noRF->GetXaxis()->SetTitle("ToF (in ns)");
    transfer_func_hist_PTBC_noRF->GetYaxis()->SetTitle("Energy (in eV)");
    // transfer_func_hist_PTBC_noRF->SetTitle(Form("Transmission Histogram - %s", filter_name_title.c_str()));
    transfer_func_hist_PTBC_noRF->SetLineWidth(2);
    transfer_func_hist_PTBC_noRF->Draw(); //"HISTE"
    transfer_func_hist_PTBC_noRF->SetStats(0);
    gPad->SetLogx();
    gPad->SetLogy();

    l[i]->AddEntry(transfer_func_hist_PTBC_5itr,"5 Itr","l");
    transfer_func_hist_PTBC_5itr->SetLineColor(6);
    transfer_func_hist_PTBC_5itr->SetLineWidth(2);
    transfer_func_hist_PTBC_5itr->Draw("SAME");

    l[i]->AddEntry(transfer_func_hist_PTBC_10itr,"10 Itr","l");
    transfer_func_hist_PTBC_10itr->SetLineColor(2);
    transfer_func_hist_PTBC_10itr->SetLineWidth(2);
    transfer_func_hist_PTBC_10itr->Draw("SAME");

    l[i]->AddEntry(transfer_func_hist_PTBC_15itr,"15 Itr","l");
    transfer_func_hist_PTBC_15itr->SetLineColor(3);
    transfer_func_hist_PTBC_15itr->SetLineWidth(2);
    transfer_func_hist_PTBC_15itr->Draw("SAME");

    l[i]->SetMargin(0.4);
    l[i]->Draw();
}

void transfer_function(){

    //opening the rf hist file
    
    TFile *rfFile = TFile::Open("../inputFiles/RF.root", "READ");
    rf_hist = (TH2D*)rfFile->Get("histfluka");

    // tof bins - 1000 bins per decade
    // Calculating TOF (x) bin edges
    int bins_per_decade = 1000;
    int Num_decades = 6;
    int num_bins_tof = bins_per_decade * Num_decades;
    Double_t bin_edges_tof[num_bins_tof+1];
    Double_t step_tof = ((Double_t) 1.0/(Double_t) bins_per_decade);
    for(int i = 0; i < num_bins_tof+1; i++)
    {
        Double_t base = 10.;
        Double_t exponent = (step_tof * (Double_t) i) + 2.;
        bin_edges_tof[i] = (Double_t) std::pow(base, exponent);
    }

    transfer_func_hist_PTBC_noRF = new TH1D("transfer_func_hist_PTBC_noRF","ToF vs Energy - Transfer Function - PTBC - no RF",num_bins_tof,bin_edges_tof);
    transfer_func_hist_PTBC_5itr = new TH1D("transfer_func_hist_PTBC_5itr","ToF vs Energy - Transfer Function - PTBC - 5 Itr",num_bins_tof,bin_edges_tof);
    transfer_func_hist_PTBC_10itr = new TH1D("transfer_func_hist_PTBC_10itr","ToF vs Energy - Transfer Function - PTBC - 10 Itr",num_bins_tof,bin_edges_tof);
    transfer_func_hist_PTBC_15itr = new TH1D("transfer_func_hist_PTBC_15itr","ToF vs Energy - Transfer Function - PTBC - 15 Itr",num_bins_tof,bin_edges_tof);

    for(Int_t i = 0; i < num_bins_tof; i++)
    {
        Double_t tof = transfer_func_hist_PTBC_noRF->GetBinCenter(i+1);
        if (tof < 800.)
        {
            transfer_func_hist_PTBC_noRF->SetBinContent(i+1, 0);
            transfer_func_hist_PTBC_5itr->SetBinContent(i+1, 0);
            transfer_func_hist_PTBC_10itr->SetBinContent(i+1, 0);
            transfer_func_hist_PTBC_15itr->SetBinContent(i+1, 0);
            continue;
        }

        Double_t rf_len = 0.;
        Double_t neutron_e = TOFToEnergy(tof * 1e-9, flight_path_length_PTBC, rf_len);
        transfer_func_hist_PTBC_noRF->SetBinContent(i+1, neutron_e);
        for (Int_t j = 0; j < 15; j++)
        {
            if (j == 5)
            {
                transfer_func_hist_PTBC_5itr->SetBinContent(i+1, neutron_e);
                // cout << "(" << transfer_func_hist_PTBC_5itr->GetBinContent(i+1) << ", " << neutron_e << ")";
            }

            if (j == 10)
            {
                transfer_func_hist_PTBC_10itr->SetBinContent(i+1, neutron_e);
            }
            
            rf_len = get_rf_length(neutron_e);
            neutron_e = TOFToEnergy(tof * 1e-9, flight_path_length_PTBC, rf_len);
            // cout << "(" << rf_len << ", " << neutron_e << ")"; 
            if (j == 14)
            {
                transfer_func_hist_PTBC_15itr->SetBinContent(i+1, neutron_e);
            }
        }  
    }

    store_hist();
    plot_hists();
}