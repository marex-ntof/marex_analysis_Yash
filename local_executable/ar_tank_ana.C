/**
 * @file ar_tank_ana.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2024-09-24
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

Int_t bins_per_decade = 20;

bool plotENDF = true;
bool plotJENDL = true;

TH1D* transmission_hist_e_PTBC = 0;
TH1D* transmission_hist_e_FIMG = 0;
TH1D* cross_section_hist_e_PTBC = 0;
TH1D* cross_section_hist_e_FIMG = 0;

TH1D* endf_trans_hist = 0;
TH1D* endf_xsec_hist = 0;
TH1D* jendl_trans_hist = 0;
TH1D* jendl_xsec_hist = 0;

///////////////////////////////////
//Systematics histogram
//backgrounds
// TH1D* trans_bkgd_hist_ptbc = 0;
// TH1D* trans_bkgd_hist_fimg = 0;

//cf bottle thickness
TH1D* trans_cft_hist_ptbc = 0;
TH1D* trans_cft_hist_fimg = 0;
TH1D* xsec_cft_hist_ptbc = 0;
TH1D* xsec_cft_hist_fimg = 0;

//Ar thcikness
TH1D* xsec_art_hist_ptbc = 0;
TH1D* xsec_art_hist_fimg = 0;
///////////////////////////////////

Double_t ar_bottle_pressure = 197.385 * 1e5; // in Pa (SI unit)
Double_t ar_bottle_temp = 293.0; // in Kelvin
Double_t cf_density = 1.55; // g/cc
Double_t n_Ar_bottle = (10 /*cm*/) * ((ar_bottle_pressure)/(8.31446261815324 * ar_bottle_temp * 1e6)) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/); //8.9

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

void fill_evaluation(TH1D* xsec_hist, Double_t min_e, Double_t max_e, const char *eval_file_name){
    //Extracting Cross Section and transmission from evaluations

    Int_t num_bins_e = xsec_hist->GetNbinsX();

    std::ifstream inputFile(Form("../evalData/%s", eval_file_name));

    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the file.\n";
    }

    //Extracting the data from the text file
    Double_t val1, val2; //val1 -> Energy (eV); val2 -> Cross Section (barns)
    std::string line;
    Double_t line_count = 0;

    std::vector<Double_t> endf_e;
    std::vector<Double_t> endf_xsec;

    while (std::getline(inputFile, line)) {

        line_count++;

        if (line_count == 1)
        {
            continue;
        }

        std::istringstream iss(line);

        if (iss >> val1 >> val2) {
            if (val1 < min_e)
            {
                continue;
            }
            if (val1 >= max_e)
            {
                break;
            }
            endf_e.push_back(val1);
            endf_xsec.push_back(val2);
        } else {
            std::cerr << "Invalid data format.\n";
        }
    }

    std::map<Int_t, std::vector<Double_t>> bin_content_map;

    for (Int_t i = 1; i <= num_bins_e; ++i) {
        bin_content_map[i] = std::vector<Double_t>();
    }

    for (Int_t i = 0; i < endf_e.size(); ++i) {
        // if (endf_e[i] < min_e) continue;
        // if (endf_e[i] >= max_e) break;

        Int_t bin_num = xsec_hist->GetXaxis()->FindBin(endf_e[i]);
        if (bin_num > num_bins_e) break;

        bin_content_map[bin_num].push_back(endf_xsec[i]);
    }

    std::vector<Double_t> bin_contents(num_bins_e, 0.0);

    for (const auto& bin : bin_content_map) {
        Int_t bin_n = bin.first;
        const std::vector<Double_t>& xsecs = bin.second;
        if (xsecs.empty()) {
            xsec_hist->SetBinContent(bin_n, 0.0);
            bin_contents[bin_n - 1] = 0.0;
        } else {
            Double_t bin_content = std::accumulate(xsecs.begin(), xsecs.end(), 0.0) / xsecs.size();
            xsec_hist->SetBinContent(bin_n, bin_content);
            bin_contents[bin_n - 1] = bin_content;
        }
    }

    // Find indices of zero bins
    std::vector<Int_t> zero_indices;
    for (Int_t i = 0; i < num_bins_e; ++i) {
        if (bin_contents[i] == 0.0) {
            zero_indices.push_back(i);
        }
    }

    for (Int_t idx : zero_indices) {
        Int_t left_idx = idx - 1;
        Int_t right_idx = idx + 1;

        // Find the nearest non-zero bins to the left
        while (left_idx >= 0 && bin_contents[left_idx] == 0.0) {
            --left_idx;
        }

        // Find the nearest non-zero bins to the right
        while (right_idx < num_bins_e && bin_contents[right_idx] == 0.0) {
            ++right_idx;
        }

        if (right_idx == num_bins_e)
        {
            break;
        }

        if (left_idx >= 0 && right_idx < num_bins_e) {
            Double_t left_value = bin_contents[left_idx];
            Double_t left_e_val = xsec_hist->GetBinCenter(left_idx + 1);
            Double_t right_value = bin_contents[right_idx];
            Double_t right_e_val = xsec_hist->GetBinCenter(right_idx + 1);

            // Linear Interpolation
            Double_t slope = (right_value - left_value) / (right_e_val - left_e_val);
            Int_t num_bins_to_update = right_idx - left_idx - 1;
            for (Int_t j = 1; j <= num_bins_to_update; ++j) {
                Double_t bin_center = xsec_hist->GetBinCenter(left_idx + j + 1);
                bin_contents[left_idx + j] = slope * (bin_center - left_e_val) + left_value;
            }
        }
    }

    // Update histogram with interpolated values
    for (int i = 1; i <= num_bins_e; ++i) {
        xsec_hist->SetBinContent(i, bin_contents[i - 1]);
    }
}

void fill_cf_thickness_sys(){

    trans_cft_hist_ptbc = (TH1D*)jendl_xsec_hist->Clone("trans_cft_hist_ptbc");
    trans_cft_hist_fimg = (TH1D*)jendl_xsec_hist->Clone("trans_cft_hist_fimg");
    xsec_cft_hist_ptbc = (TH1D*)jendl_xsec_hist->Clone("xsec_cft_hist_ptbc");
    xsec_cft_hist_fimg = (TH1D*)jendl_xsec_hist->Clone("xsec_cft_hist_fimg");

    //calculatnig CF transission hist
    TH1D* c_xsec_hist = (TH1D*)jendl_xsec_hist->Clone("c_xsec_hist");
    TH1D* h_xsec_hist = (TH1D*)jendl_xsec_hist->Clone("h_xsec_hist");
    fill_evaluation(c_xsec_hist, 1e-2, 1e8, "../evalData/JENDL_C_tot_xsec.txt");
    fill_evaluation(h_xsec_hist, 1e-2, 1e8, "../evalData/JENDL_H_tot_xsec.txt");

    Double_t c_const = cf_density * 0.9 * (6.02214076e23) * (1e-24) / (12.011);
    Double_t h_const = cf_density * 0.1 * (6.02214076e23) * (1e-24) / (1.008);
    Double_t del_cf_t = 0.004;
    Int_t num_bins = h_xsec_hist->GetNbinsX();
    TH1D* cf_trans_hist = (TH1D*)h_xsec_hist->Clone("cf_trans_hist");

    for (Int_t i = 1; i < num_bins+1; i++)
    {
        Double_t cf_const = c_const*c_xsec_hist->GetBinContent(i) + h_const*h_xsec_hist->GetBinContent(i);
        cf_trans_hist->SetBinContent(i, std::exp(- del_cf_t * cf_const));
    }
    ///////////////////

    for (Int_t i = 1; i < num_bins+1; i++)
    {
        //PTBC
        Double_t trans_bin_lim_low = transmission_hist_e_PTBC->GetBinContent(i) * cf_trans_hist->GetBinContent(i);
        Double_t trans_bin_lim_up = transmission_hist_e_PTBC->GetBinContent(i) / cf_trans_hist->GetBinContent(i);
        Double_t trans_bin_content = (trans_bin_lim_up+trans_bin_lim_low)/2;
        Double_t trans_bin_err = (trans_bin_lim_up-trans_bin_lim_low)/2;
        trans_cft_hist_ptbc->SetBinContent(i,trans_bin_content);
        trans_cft_hist_ptbc->SetBinError(i,trans_bin_err);

        Double_t xsec_bin_lim_low = 0;
        Double_t xsec_bin_lim_up = 0;

        if (trans_bin_lim_up != 0)
        {
            xsec_bin_lim_low = - std::log(trans_bin_lim_up) / n_Ar_bottle;
        }

        if (trans_bin_lim_low != 0)
        {
            xsec_bin_lim_up = - std::log(trans_bin_lim_low) / n_Ar_bottle;
        }

        Double_t xsec_bin_content = (xsec_bin_lim_up+xsec_bin_lim_low)/2;
        Double_t xsec_bin_err = (xsec_bin_lim_up-xsec_bin_lim_low)/2;
        xsec_cft_hist_ptbc->SetBinContent(i,xsec_bin_content);
        xsec_cft_hist_ptbc->SetBinError(i,xsec_bin_err);
    }

    for (Int_t i = 1; i < num_bins+1; i++)
    {
        //FIMG
        Double_t trans_bin_lim_low = transmission_hist_e_FIMG->GetBinContent(i) * cf_trans_hist->GetBinContent(i);
        Double_t trans_bin_lim_up = transmission_hist_e_FIMG->GetBinContent(i) / cf_trans_hist->GetBinContent(i);
        Double_t trans_bin_content = (trans_bin_lim_up+trans_bin_lim_low)/2;
        Double_t trans_bin_err = (trans_bin_lim_up-trans_bin_lim_low)/2;
        trans_cft_hist_fimg->SetBinContent(i,trans_bin_content);
        trans_cft_hist_fimg->SetBinError(i,trans_bin_err);

        Double_t xsec_bin_lim_low = 0;
        Double_t xsec_bin_lim_up = 0;

        if (trans_bin_lim_up != 0)
        {
            xsec_bin_lim_low = - std::log(trans_bin_lim_up) / n_Ar_bottle;
        }

        if (trans_bin_lim_low != 0)
        {
            xsec_bin_lim_up = - std::log(trans_bin_lim_low) / n_Ar_bottle;
        }

        Double_t xsec_bin_content = (xsec_bin_lim_up+xsec_bin_lim_low)/2;
        Double_t xsec_bin_err = (xsec_bin_lim_up-xsec_bin_lim_low)/2;
        xsec_cft_hist_fimg->SetBinContent(i,xsec_bin_content);
        xsec_cft_hist_fimg->SetBinError(i,xsec_bin_err);
    }
    
    // SetMArEXStyle();
    // TCanvas* c_1 = new TCanvas("C_1"," ");
    // c_1->cd();
    // transmission_hist_e_FIMG->GetXaxis()->SetRangeUser(1e-2,1e8);
    // transmission_hist_e_FIMG->Draw();
    // trans_cft_hist_fimg->SetFillColor(kRed);
    // trans_cft_hist_fimg->SetFillStyle(3001);
    // trans_cft_hist_fimg->Draw("e2same");
    // gPad->SetLogx();
}

void fill_ar_thickness_sys(){

    xsec_art_hist_ptbc = (TH1D*)jendl_xsec_hist->Clone("xsec_art_hist_ptbc");
    xsec_art_hist_fimg = (TH1D*)jendl_xsec_hist->Clone("xsec_art_hist_fimg");

    Int_t num_bins = jendl_xsec_hist->GetNbinsX();
    Double_t n_Ar_bottle_more = (11 /*cm*/) * ((ar_bottle_pressure)/(8.31446261815324 * ar_bottle_temp * 1e6)) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/); //8.9
    Double_t n_Ar_bottle_less = (9 /*cm*/) * ((ar_bottle_pressure)/(8.31446261815324 * ar_bottle_temp * 1e6)) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/); //8.9

    for (Int_t i = 1; i < num_bins+1; i++)
    {
        //PTBC
        Double_t xsec_bin_lim_low = 0;
        Double_t xsec_bin_lim_up = 0;

        if (transmission_hist_e_PTBC->GetBinContent(i) != 0)
        {
            xsec_bin_lim_low = - std::log(transmission_hist_e_PTBC->GetBinContent(i)) / n_Ar_bottle_more;
            xsec_bin_lim_up = - std::log(transmission_hist_e_PTBC->GetBinContent(i)) / n_Ar_bottle_less;
        }

        Double_t xsec_bin_content = (xsec_bin_lim_up+xsec_bin_lim_low)/2;
        Double_t xsec_bin_err = (xsec_bin_lim_up-xsec_bin_lim_low)/2;
        xsec_art_hist_ptbc->SetBinContent(i,xsec_bin_content);
        xsec_art_hist_ptbc->SetBinError(i,xsec_bin_err);
    }

    for (Int_t i = 1; i < num_bins+1; i++)
    {
        //FIMG
        Double_t xsec_bin_lim_low = 0;
        Double_t xsec_bin_lim_up = 0;

        if (transmission_hist_e_FIMG->GetBinContent(i) != 0)
        {
            xsec_bin_lim_low = - std::log(transmission_hist_e_FIMG->GetBinContent(i)) / n_Ar_bottle_more;
            xsec_bin_lim_up = - std::log(transmission_hist_e_FIMG->GetBinContent(i)) / n_Ar_bottle_less;
        }

        Double_t xsec_bin_content = (xsec_bin_lim_up+xsec_bin_lim_low)/2;
        Double_t xsec_bin_err = (xsec_bin_lim_up-xsec_bin_lim_low)/2;
        xsec_art_hist_fimg->SetBinContent(i,xsec_bin_content);
        xsec_art_hist_fimg->SetBinError(i,xsec_bin_err);
    }

    // SetMArEXStyle();
    // TCanvas* c_1 = new TCanvas("C_1"," ");
    // c_1->cd();
    // cross_section_hist_e_FIMG->GetXaxis()->SetRangeUser(1e-2,1e8);
    // cross_section_hist_e_FIMG->Draw();
    // xsec_art_hist_fimg->SetFillColor(kRed);
    // xsec_art_hist_fimg->SetFillStyle(3001);
    // xsec_art_hist_fimg->Draw("e2same");
    // gPad->SetLogx();
    
}

void ar_tank_ana(){

    transmission_hist_e_PTBC = retriveHistograms("../rootFiles/trans_xsec_hists_ar_bottle_full_20bpd.root", "transmission_hist_e_PTBC");
    transmission_hist_e_FIMG = retriveHistograms("../rootFiles/trans_xsec_hists_ar_bottle_full_20bpd.root", "transmission_hist_e_FIMG");
    cross_section_hist_e_PTBC = retriveHistograms("../rootFiles/trans_xsec_hists_ar_bottle_full_20bpd.root", "cross_section_hist_e_PTBC");
    cross_section_hist_e_FIMG = retriveHistograms("../rootFiles/trans_xsec_hists_ar_bottle_full_20bpd.root", "cross_section_hist_e_FIMG");
    // endf_trans_hist = retriveHistograms("../rootFiles/trans_xsec_hists_ar_bottle_full_20bpd.root", "natAr_trans_hist_20bpd_endf");
    // endf_xsec_hist = retriveHistograms("../rootFiles/trans_xsec_hists_ar_bottle_full_20bpd.root", "natAr_xsec_hist_20bpd_endf");
    jendl_trans_hist = retriveHistograms("../rootFiles/trans_xsec_hists_ar_bottle_full_20bpd.root", "natAr_trans_hist_20bpd_jendl");
    jendl_xsec_hist = retriveHistograms("../rootFiles/trans_xsec_hists_ar_bottle_full_20bpd.root", "natAr_xsec_hist_20bpd_jendl");

    fill_cf_thickness_sys();
    fill_ar_thickness_sys();
}