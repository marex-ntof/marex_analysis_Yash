/**
 * @file cf_tank_ana.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2024-09-13
 */

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <map>
#include <vector>

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

Int_t bins_per_decade = 1;

TH1D* transmission_hist_e_PTBC[3];
TH1D* transmission_hist_e_FIMG[3];

TH1D* energy_hist_target_in_PTBC[3];
TH1D* energy_hist_target_in_FIMG[3];
TH1D* energy_hist_target_out_PTBC;
TH1D* energy_hist_target_out_FIMG;

Double_t trans_norm_factors[3];

const Double_t flight_path_length_PTBC = 182.65 - 0.41; //m
const Double_t flight_path_length_FIMG = 183.5 - 0.41; //m
const Double_t neutron_mass = 939.56542052; //in MeV
const Double_t speed_of_light = 299792458.0; //in m/s

Int_t FindDecadePower(Double_t num){
    Int_t decadePower = 0;
    Double_t value = num;
    if (num > 1)
    {
        while (value > 1)
        {
            value /= 10;
            decadePower++;
        }   
    }
    if (num < 1)
    {
        while (value < 1)
        {
            value *= 10;
            decadePower--;
        } 
    }
    return decadePower;
}

TH1D* retriveHistograms(const char *fname, const char *hist_name){
    
    TFile* hist_file = TFile::Open(fname, "READ");
    TH1D* hist_new = (TH1D*)hist_file->Get(hist_name);

    return hist_new;
}

void exclude_first_last_bins(TH1D* hist_to_change){
    // excluding the first and the last bin
    Int_t seaching_for_first_last_bin = 1; // 1 - first bin; 2 - last bin
    Int_t num_bins_e = hist_to_change->GetNbinsX();
    for (Int_t i = 1; i <= num_bins_e; i++)
    {
        Double_t new_bin_content = hist_to_change->GetBinContent(i);

        // Searching for the first bin and setting it to zero
        if (seaching_for_first_last_bin == 1)
        {
            if (new_bin_content != 0)
            {
                hist_to_change->SetBinContent(i, 0);
                hist_to_change->SetBinError(i, 0);
                seaching_for_first_last_bin = 2;
                continue;
            }
        }
        
        // Searching for the last bin and setting it to zero
        if (seaching_for_first_last_bin == 2)
        {
            if (new_bin_content == 0)
            {
                hist_to_change->SetBinContent(i-1, 0);
                hist_to_change->SetBinError(i-1, 0);
                seaching_for_first_last_bin = 3;
                break;
            }
        }
    }
}

TH1D* retriveHistogramsChangeBPD(const char *fname, const char *hist_name, const char *hist_name_new, Int_t bpd_old, Int_t bpd_new){

    cout << "Opening File " << fname << endl;
    TFile *hist_file = TFile::Open(fname, "READ");

    cout << "Extracting Histogram " << hist_name << endl;
    TH1D* hist_old = (TH1D*)hist_file->Get(hist_name);
    Int_t num_bins_old = hist_old->GetNbinsX();
    // cout << "Number of old bins = " << num_bins_old << endl;
    Int_t max_sum_bin_count = (bpd_old/bpd_new); // Number of old bins that need to be summed to make a new bin
    Int_t num_bins_new = num_bins_old / max_sum_bin_count;
    // cout << "Number of new bins = " << num_bins_new << endl;

    Double_t bin_edges_new[num_bins_new + 1];
    Double_t bin_content[num_bins_new];
    Double_t bin_error[num_bins_new];
    Int_t sum_bin_count = 0;
    Int_t new_bin_counter = 0;
    for (int i = 1; i < num_bins_old + 1; i++)
    {
        if (i == 1)
        {
            bin_edges_new[new_bin_counter] = hist_old->GetXaxis()->GetBinLowEdge(i);
            bin_content[new_bin_counter] = 0;
            bin_error[new_bin_counter] = 0;
        }
        bin_content[new_bin_counter] += hist_old->GetBinContent(i);
        sum_bin_count++;

        if (i == num_bins_old)
        {
            sum_bin_count = 0;
            bin_error[new_bin_counter] = sqrt(bin_content[new_bin_counter]);
            bin_edges_new[new_bin_counter+1] = hist_old->GetXaxis()->GetBinUpEdge(i);
        }
        else if (sum_bin_count == max_sum_bin_count)
        {
            sum_bin_count = 0;
            bin_error[new_bin_counter] = sqrt(bin_content[new_bin_counter]);
            new_bin_counter++;
            bin_edges_new[new_bin_counter] = hist_old->GetXaxis()->GetBinUpEdge(i);
            bin_content[new_bin_counter] = 0;
            bin_error[new_bin_counter] = 0;
        } 
    }

    auto hist_new = new TH1D(hist_name_new,hist_old->GetTitle(),num_bins_new,bin_edges_new);
    for(int i = 0; i < num_bins_new; i++){
        hist_new->SetBinContent(i+1, bin_content[i]);
        hist_new->SetBinError(i+1, bin_error[i]);
    }

    return hist_new;
}

Double_t get_norm_factor(const char *fname){

    cout << "Extracting Norm Factors " << endl;
    std::vector<Double_t> *norm_factors_temp;
    std::vector<Double_t> norm_factors;
    TFile* hist_file = TFile::Open(fname, "READ");
    hist_file->GetObject("norm_factors", norm_factors_temp);
    norm_factors = *norm_factors_temp;
    return (norm_factors[1] / norm_factors[0]);

}

void fill_trans(TH1D* target_in_hist, TH1D* target_out_hist, TH1D* trans_hist, Double_t Qout_Qin){

    Int_t num_bins_e = target_in_hist->GetNbinsX();

    for (int i = 0; i < num_bins_e; i++)
    {
        //PTBC
        Double_t bin_content_in_PTBC = target_in_hist->GetBinContent(i+1);
        Double_t bin_content_out_PTBC = target_out_hist->GetBinContent(i+1);
        if (bin_content_in_PTBC == 0. || bin_content_out_PTBC == 0.)
        {
            trans_hist->SetBinContent(i+1, 0.);
            trans_hist->SetBinError(i+1, 0.);
        } else {
            Double_t transmission_PTBC = (bin_content_in_PTBC * Qout_Qin)/bin_content_out_PTBC;
            Double_t bin_unc_PTBC = transmission_PTBC * std::sqrt( (1./bin_content_in_PTBC) + (1./bin_content_out_PTBC) );
            trans_hist->SetBinContent(i+1, transmission_PTBC);
            trans_hist->SetBinError(i+1, bin_unc_PTBC);
        }
    }
}

void calc_trans_1ev_10kev(){

    Double_t target_in_counts[3];
    Double_t target_out_counts = 0.;
    Double_t trans_values[3];
    Double_t trans_values_errors[3];

    for (Int_t i = 0; i < 3; i++)
    {
        target_in_counts[i] = 0.;
        for (Int_t j = 3; j < 7; j++)
        {
            target_in_counts[i] = target_in_counts[i] + energy_hist_target_in_PTBC[i]->GetBinContent(j);
        }
    }
    
    for (Int_t j = 3; j < 7; j++)
    {
        target_out_counts = target_out_counts + energy_hist_target_out_PTBC->GetBinContent(j);
    }

    for (Int_t i = 0; i < 3; i++)
    {
        trans_values[i] = (target_in_counts[i] * trans_norm_factors[i])/target_out_counts;
        trans_values_errors[i] = trans_values[i] * std::sqrt( (1./target_in_counts[i]) + (1./target_out_counts) );
    }
    
    cout << "transmission of CF tank = " << trans_values[0] << " +/- " << trans_values_errors[0];
    cout << "; range = " << trans_values[0]-trans_values_errors[0] << " - " << trans_values[0]+trans_values_errors[0] << endl;
    cout << "transmission of CF tank - rotated 90 deg = " << trans_values[1] << " +/- " << trans_values_errors[1];
    cout << "; range = " << trans_values[1]-trans_values_errors[1] << " - " << trans_values[1]+trans_values_errors[1] << endl;
    cout << "transmission of CF tank - rotated back = " << trans_values[2] << " +/- " << trans_values_errors[2];
    cout << "; range = " << trans_values[2]-trans_values_errors[2] << " - " << trans_values[2]+trans_values_errors[2] << endl;

    cout << "percent diff, CF tank and rotated 90 deg = " << ((trans_values[0] - trans_values[1])/trans_values[0])*100 << endl;
    cout << "percent diff, CF tank and rotated back = " << ((trans_values[2] - trans_values[0])/trans_values[0])*100 << endl;
}

void plotting(){
    //Plotting
    SetMArEXStyle();

    gStyle->SetStatX(0.27);
    gStyle->SetStatY(0.9);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.17);

    TCanvas *c[4];
    TLegend *l[4];

    int i = 0;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();

    transmission_hist_e_PTBC[2]->GetXaxis()->SetTitle("Energy (in eV)");
    transmission_hist_e_PTBC[2]->GetYaxis()->SetTitle("Transmission");
    transmission_hist_e_PTBC[2]->SetTitle("Transmission Histogram - CF Tank");
    transmission_hist_e_PTBC[2]->SetLineColor(2);
    transmission_hist_e_PTBC[2]->SetLineWidth(2);
    transmission_hist_e_PTBC[2]->Draw(); //"HISTE"
    transmission_hist_e_PTBC[2]->SetStats(0);
    // transmission_hist_e_PTBC[2]->SetMarkerStyle(6);
    // transmission_hist_e_PTBC[2]->SetMarkerSize(0.5);
    // gPad->SetGrid();
    gPad->SetLogx();
    // gStyle->SetPalette(57);

    transmission_hist_e_PTBC[0]->SetLineColor(3);
    transmission_hist_e_PTBC[0]->SetLineWidth(2);
    // transmission_hist_e_PTBC[0]->GetXaxis()->SetRangeUser(1e-2,1e3);
    transmission_hist_e_PTBC[0]->Draw("SAME");

    transmission_hist_e_PTBC[1]->SetLineColor(4);
    transmission_hist_e_PTBC[1]->SetLineWidth(2);
    // transmission_hist_e_PTBC[1]->GetXaxis()->SetRangeUser(1e-2,1e3);
    transmission_hist_e_PTBC[1]->Draw("SAME");

    l[i] = new TLegend(0.77,0.7,0.86,0.85); //0.68,0.7,0.86,0.8       ;         0.72,0.8,0.90,0.9
    l[i]->AddEntry(transmission_hist_e_PTBC[0],"CF Tank","l");
    l[i]->AddEntry(transmission_hist_e_PTBC[1],"CF Tank - Rotated 90 deg","l");
    l[i]->AddEntry(transmission_hist_e_PTBC[2],"CF Tank - Rotated Back","l");

    l[i]->SetMargin(0.4);
    l[i]->Draw();
}

void cf_tank_ana(){

    energy_hist_target_in_PTBC[0] = retriveHistogramsChangeBPD("../rootFiles/crossSectionAna_cf_bottle.root", "energy_hist_target_in_PTB", "e_hist_target_in_PTB_cfb", 1000, bins_per_decade);
    energy_hist_target_in_PTBC[1] = retriveHistogramsChangeBPD("../rootFiles/crossSectionAna_cf_bottle_rot.root", "energy_hist_target_in_PTB", "e_hist_target_in_PTB_cfbRot", 1000, bins_per_decade);
    energy_hist_target_in_PTBC[2] = retriveHistogramsChangeBPD("../rootFiles/crossSectionAna_cf_bottle_rotBack.root", "energy_hist_target_in_PTB", "e_hist_target_in_PTB_cfbRotBack", 1000, bins_per_decade);

    energy_hist_target_out_PTBC = retriveHistogramsChangeBPD("../rootFiles/crossSectionAna_cf_bottle.root", "energy_hist_target_out_PTB", "e_hist_target_out_PTB_cfb", 1000, bins_per_decade);

    //norm factors
    trans_norm_factors[0] = get_norm_factor("../rootFiles/crossSectionAna_cf_bottle.root");
    trans_norm_factors[1] = get_norm_factor("../rootFiles/crossSectionAna_cf_bottle_rot.root");
    trans_norm_factors[2] = get_norm_factor("../rootFiles/crossSectionAna_cf_bottle_rotBack.root");

    //Transmission
    transmission_hist_e_PTBC[0] = (TH1D*)energy_hist_target_out_PTBC->Clone("tran_hist_PTBC_cfb");
    transmission_hist_e_PTBC[1] = (TH1D*)energy_hist_target_out_PTBC->Clone("tran_hist_PTBC_cfb_cfbRot");
    transmission_hist_e_PTBC[2] = (TH1D*)energy_hist_target_out_PTBC->Clone("tran_hist_PTBC_cfb_cfbRotBack");

    for (Int_t i = 0; i < 3; i++)
    {
        fill_trans(energy_hist_target_in_PTBC[i], energy_hist_target_out_PTBC, transmission_hist_e_PTBC[i], trans_norm_factors[i]);
    }

    calc_trans_1ev_10kev();
    // plotting();
}