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
TH1D* trans_bkgd_hist_ptbc = 0;
TH1D* trans_bkgd_hist_fimg = 0;
TH1D* xsec_bkgd_hist_ptbc = 0;
TH1D* xsec_bkgd_hist_fimg = 0;

//cf bottle thickness
TH1D* trans_cft_hist_ptbc = 0;
TH1D* trans_cft_hist_fimg = 0;
TH1D* xsec_cft_hist_ptbc = 0;
TH1D* xsec_cft_hist_fimg = 0;

//Ar thcikness
TH1D* xsec_art_hist_ptbc = 0;
TH1D* xsec_art_hist_fimg = 0;

//combined
TH1D* trans_combined_sys_hist_ptbc = 0;
TH1D* trans_combined_sys_hist_fimg = 0;
TH1D* xsec_combined_sys_hist_ptbc = 0;
TH1D* xsec_combined_sys_hist_fimg = 0;
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

void fill_cf_thickness_sys(bool plot_sys){

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
    
    if (plot_sys)
    {
        SetMArEXStyle();
        gStyle->SetCanvasDefW(1000); //600
        gStyle->SetCanvasDefH(500); //500 
        gStyle->SetPadRightMargin(0.04);
        gStyle->SetTitleH(0.1);
        gStyle->SetTitleAlign(33);
        gStyle->SetTitleX(.76);
        // gStyle->SetStats(0);
        
        ////////////////////// PTBC
        TCanvas* c_cf_ptbc = new TCanvas("c_cf_ptbc"," ");
        TPad* p_cf_ptbc[2];
        TLegend* l_cf_ptbc[2];
        //Pad - 1 - CF
        c_cf_ptbc->cd(0);
        p_cf_ptbc[0] = new TPad("p_cf_ptbc_0", "p_cf_ptbc_0", 0., 0.5, 1., 1.);
        p_cf_ptbc[0]->SetFillColor(kWhite);
        p_cf_ptbc[0]->SetBottomMargin(0.00001);
        p_cf_ptbc[0]->SetTopMargin(0.2);
        p_cf_ptbc[0]->SetBorderMode(0);
        p_cf_ptbc[0]->Draw();
        p_cf_ptbc[0]->cd();
        trans_cft_hist_ptbc->SetTitle("Fission Chamber - CF Tank Systematic");
        trans_cft_hist_ptbc->SetTitleOffset(0.5);
        // X-Axis
        trans_cft_hist_ptbc->GetXaxis()->SetLabelOffset(999);
        trans_cft_hist_ptbc->GetXaxis()->SetLabelSize(0);
        trans_cft_hist_ptbc->GetXaxis()->SetTitle("");
        // Y-Axis
        trans_cft_hist_ptbc->GetYaxis()->SetTitle("Transmission");
        trans_cft_hist_ptbc->GetYaxis()->SetLabelSize(0.09);
        trans_cft_hist_ptbc->GetYaxis()->SetTitleSize(0.09);
        trans_cft_hist_ptbc->GetYaxis()->SetTitleOffset(0.45);
        trans_cft_hist_ptbc->SetFillColor(kRed);
        // trans_cft_hist_ptbc->SetFillStyle(3001);
        trans_cft_hist_ptbc->GetXaxis()->SetRangeUser(1e-1,1e8);
        trans_cft_hist_ptbc->GetYaxis()->SetRangeUser(0.22,1.25);
        trans_cft_hist_ptbc->Draw("e2");

        transmission_hist_e_PTBC->SetTitle("");
        transmission_hist_e_PTBC->GetXaxis()->SetRangeUser(1e-1,1e8);
        // transmission_hist_e_PTBC->SetLineWidth(2);
        transmission_hist_e_PTBC->Draw("same");
        gPad->SetLogx();
        gPad->SetGrid();

        l_cf_ptbc[0] = new TLegend(0.15, 0.10, 0.45, 0.30);
        l_cf_ptbc[0]->AddEntry(trans_cft_hist_ptbc, "CF Tank thickness systematic", "f");
        l_cf_ptbc[0]->AddEntry(transmission_hist_e_PTBC, "Transmission w/ stat errors", "l");
        l_cf_ptbc[0]->Draw();

        //Pad - 2 - Ar
        c_cf_ptbc->cd(0);
        p_cf_ptbc[1] = new TPad("p_cf_ptbc_1", "p_cf_ptbc_1", 0., 0., 1., 0.5);
        p_cf_ptbc[1]->SetFillColor(kWhite);
        p_cf_ptbc[1]->SetTopMargin(0.00001);
        p_cf_ptbc[1]->SetBottomMargin(0.2);
        p_cf_ptbc[1]->SetBorderMode(0);
        p_cf_ptbc[1]->Draw();
        p_cf_ptbc[1]->cd();

        xsec_cft_hist_ptbc->SetTitle("");
        // X-Axis
        xsec_cft_hist_ptbc->GetXaxis()->SetTitle("Energy (in eV)");
        xsec_cft_hist_ptbc->GetXaxis()->SetTitleOffset(1.05); //increase to move down
        xsec_cft_hist_ptbc->GetXaxis()->SetLabelSize(0.09);
        xsec_cft_hist_ptbc->GetXaxis()->SetTitleSize(0.09);
        // Y-Axis
        xsec_cft_hist_ptbc->GetYaxis()->SetTitle("Cross Section (in b)");
        xsec_cft_hist_ptbc->GetYaxis()->SetLabelSize(0.09);
        xsec_cft_hist_ptbc->GetYaxis()->SetTitleSize(0.09);
        xsec_cft_hist_ptbc->GetYaxis()->SetTitleOffset(0.45);
        xsec_cft_hist_ptbc->SetFillColor(kRed);
        // xsec_cft_hist_ptbc->SetFillStyle(3001);
        xsec_cft_hist_ptbc->GetXaxis()->SetRangeUser(1e-1,1e8);
        xsec_cft_hist_ptbc->GetYaxis()->SetRangeUser(-4.8,23.);
        xsec_cft_hist_ptbc->Draw("e2");

        cross_section_hist_e_PTBC->SetTitle("");
        cross_section_hist_e_PTBC->GetXaxis()->SetRangeUser(1e-1,1e8);
        // cross_section_hist_e_PTBC->SetLineWidth(2);
        cross_section_hist_e_PTBC->Draw("same");
        gPad->SetLogx();

        l_cf_ptbc[1] = new TLegend(0.15, 0.70, 0.45, 0.90);
        l_cf_ptbc[1]->AddEntry(xsec_cft_hist_ptbc, "CF Tank thickness systematic", "f");
        l_cf_ptbc[1]->AddEntry(cross_section_hist_e_PTBC, "Cross section w/ stat errors", "l");
        l_cf_ptbc[1]->Draw("same");

        ////////////////////// FIMG
        TCanvas* c_cf_fimg = new TCanvas("c_cf_fimg"," ");
        TPad* p_cf_fimg[2];
        TLegend* l_cf_fimg[2];
        //Pad - 1 - CF
        c_cf_fimg->cd(0);
        p_cf_fimg[0] = new TPad("p_cf_fimg_0", "p_cf_fimg_0", 0., 0.5, 1., 1.);
        p_cf_fimg[0]->SetFillColor(kWhite);
        p_cf_fimg[0]->SetBottomMargin(0.00001);
        p_cf_fimg[0]->SetTopMargin(0.2);
        p_cf_fimg[0]->SetBorderMode(0);
        p_cf_fimg[0]->Draw();
        p_cf_fimg[0]->cd();
        trans_cft_hist_fimg->SetTitle("Micromegas - CF Tank Systematic");
        trans_cft_hist_fimg->SetTitleOffset(0.5);
        // X-Axis
        trans_cft_hist_fimg->GetXaxis()->SetLabelOffset(999);
        trans_cft_hist_fimg->GetXaxis()->SetLabelSize(0);
        trans_cft_hist_fimg->GetXaxis()->SetTitle("");
        // Y-Axis
        trans_cft_hist_fimg->GetYaxis()->SetTitle("Transmission");
        trans_cft_hist_fimg->GetYaxis()->SetLabelSize(0.09);
        trans_cft_hist_fimg->GetYaxis()->SetTitleSize(0.09);
        trans_cft_hist_fimg->GetYaxis()->SetTitleOffset(0.45);
        trans_cft_hist_fimg->SetFillColor(kRed);
        // trans_cft_hist_fimg->SetFillStyle(3001);
        trans_cft_hist_fimg->GetXaxis()->SetRangeUser(1e-1,1e6);
        trans_cft_hist_fimg->GetYaxis()->SetRangeUser(0.42,1.25);
        trans_cft_hist_fimg->Draw("e2");

        transmission_hist_e_FIMG->SetTitle("");
        transmission_hist_e_FIMG->GetXaxis()->SetRangeUser(1e-1,1e6);
        // transmission_hist_e_FIMG->SetLineWidth(2);
        transmission_hist_e_FIMG->Draw("same");
        gPad->SetLogx();
        gPad->SetGrid();

        l_cf_fimg[0] = new TLegend(0.15, 0.10, 0.45, 0.30);
        l_cf_fimg[0]->AddEntry(trans_cft_hist_fimg, "CF Tank thickness systematic", "f");
        l_cf_fimg[0]->AddEntry(transmission_hist_e_FIMG, "Transmission w/ stat errors", "l");
        l_cf_fimg[0]->Draw();

        //Pad - 2 - Ar
        c_cf_fimg->cd(0);
        p_cf_fimg[1] = new TPad("p_cf_fimg_1", "p_cf_fimg_1", 0., 0., 1., 0.5);
        p_cf_fimg[1]->SetFillColor(kWhite);
        p_cf_fimg[1]->SetTopMargin(0.00001);
        p_cf_fimg[1]->SetBottomMargin(0.2);
        p_cf_fimg[1]->SetBorderMode(0);
        p_cf_fimg[1]->Draw();
        p_cf_fimg[1]->cd();

        xsec_cft_hist_fimg->SetTitle("");
        // X-Axis
        xsec_cft_hist_fimg->GetXaxis()->SetTitle("Energy (in eV)");
        xsec_cft_hist_fimg->GetXaxis()->SetTitleOffset(1.05); //increase to move down
        xsec_cft_hist_fimg->GetXaxis()->SetLabelSize(0.09);
        xsec_cft_hist_fimg->GetXaxis()->SetTitleSize(0.09);
        // Y-Axis
        xsec_cft_hist_fimg->GetYaxis()->SetTitle("Cross Section (in b)");
        xsec_cft_hist_fimg->GetYaxis()->SetLabelSize(0.09);
        xsec_cft_hist_fimg->GetYaxis()->SetTitleSize(0.09);
        xsec_cft_hist_fimg->GetYaxis()->SetTitleOffset(0.45);
        xsec_cft_hist_fimg->SetFillColor(kRed);
        // xsec_cft_hist_fimg->SetFillStyle(3001);
        xsec_cft_hist_fimg->GetXaxis()->SetRangeUser(1e-1,1e6);
        xsec_cft_hist_fimg->GetYaxis()->SetRangeUser(-4.8,19.);
        xsec_cft_hist_fimg->Draw("e2");

        cross_section_hist_e_FIMG->SetTitle("");
        cross_section_hist_e_FIMG->GetXaxis()->SetRangeUser(1e-1,1e6);
        // cross_section_hist_e_FIMG->SetLineWidth(2);
        cross_section_hist_e_FIMG->Draw("same");
        gPad->SetLogx();

        l_cf_fimg[1] = new TLegend(0.15, 0.70, 0.45, 0.90);
        l_cf_fimg[1]->AddEntry(xsec_cft_hist_fimg, "CF Tank thickness systematic", "f");
        l_cf_fimg[1]->AddEntry(cross_section_hist_e_FIMG, "Cross section w/ stat errors", "l");
        l_cf_fimg[1]->Draw("same");
    }
    
}

void fill_ar_thickness_sys(bool plot_sys){

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

    if (plot_sys)
    {
        SetMArEXStyle();
        gStyle->SetPadRightMargin(0.04);
        // gStyle->SetPadBottomMargin(0.17);
        gStyle->SetCanvasDefW(800); //600
        gStyle->SetCanvasDefH(400); //500 
        // gStyle->SetTitleH(0.1);
        // gStyle->SetTitleAlign(33);
        // gStyle->SetTitleX(.76);
        // gStyle->SetStats(0);
        
        TCanvas* c_ar[2];
        c_ar[0] = new TCanvas("c_ar_0"," ");
        TPad* p_ar[2];
        TLegend* l_ar[2];
        //canvas - 1 - PTBC
        c_ar[0]->cd(0);
        p_ar[0] = new TPad("p_ar_0", "p_ar_0", 0., 0., 1., 1.);
        p_ar[0]->SetFillColor(kWhite);
        p_ar[0]->SetBottomMargin(0.15);
        p_ar[0]->SetTopMargin(0.1);
        p_ar[0]->SetBorderMode(0);
        p_ar[0]->Draw();
        p_ar[0]->cd();
        xsec_art_hist_ptbc->SetTitle("Fission Chamber");
        // X-Axis
        xsec_art_hist_ptbc->GetXaxis()->SetTitle("Energy (in eV)");
        xsec_art_hist_ptbc->GetXaxis()->SetTitleOffset(1.2); //increase to move down
        xsec_art_hist_ptbc->GetXaxis()->SetLabelSize(0.06);
        xsec_art_hist_ptbc->GetXaxis()->SetTitleSize(0.06);
        // Y-Axis
        xsec_art_hist_ptbc->GetYaxis()->SetTitle("Cross Section (in b)");
        xsec_art_hist_ptbc->GetYaxis()->SetLabelSize(0.06);
        xsec_art_hist_ptbc->GetYaxis()->SetTitleSize(0.06);
        xsec_art_hist_ptbc->GetYaxis()->SetTitleOffset(0.55);
        xsec_art_hist_ptbc->SetFillColor(kRed);
        // xsec_art_hist_ptbc->SetFillStyle(3001);
        xsec_art_hist_ptbc->GetXaxis()->SetRangeUser(1e-1,1e8);
        xsec_art_hist_ptbc->GetYaxis()->SetRangeUser(-4.8,23.);
        xsec_art_hist_ptbc->Draw("e2");

        cross_section_hist_e_PTBC->SetTitle("");
        cross_section_hist_e_PTBC->GetXaxis()->SetRangeUser(1e-1,1e8);
        // cross_section_hist_e_PTBC->SetLineWidth(2);
        cross_section_hist_e_PTBC->Draw("same");
        gPad->SetLogx();
        gPad->SetGrid();

        l_ar[0] = new TLegend(0.15, 0.60, 0.45, 0.80);
        l_ar[0]->AddEntry(xsec_art_hist_ptbc, "Argon thickness systematic", "f");
        l_ar[0]->AddEntry(cross_section_hist_e_PTBC, "Cross Section w/ stat errors", "l");
        l_ar[0]->Draw();

        //canvas - 2 - FIMG
        c_ar[1] = new TCanvas("c_ar_1"," ");
        c_ar[1]->cd(0);
        p_ar[1] = new TPad("p_ar_1", "p_ar_1", 0., 0., 1., 1.);
        p_ar[1]->SetFillColor(kWhite);
        p_ar[1]->SetTopMargin(0.1);
        p_ar[1]->SetBottomMargin(0.15);
        p_ar[1]->SetBorderMode(0);
        p_ar[1]->Draw();
        p_ar[1]->cd();

        xsec_art_hist_fimg->SetTitle("Micromegas");
        // X-Axis
        xsec_art_hist_fimg->GetXaxis()->SetTitle("Energy (in eV)");
        xsec_art_hist_fimg->GetXaxis()->SetTitleOffset(1.2); //increase to move up
        xsec_art_hist_fimg->GetXaxis()->SetLabelSize(0.06);
        xsec_art_hist_fimg->GetXaxis()->SetTitleSize(0.06);
        // Y-Axis
        xsec_art_hist_fimg->GetYaxis()->SetTitle("Cross Section (in b)");
        xsec_art_hist_fimg->GetYaxis()->SetLabelSize(0.06);
        xsec_art_hist_fimg->GetYaxis()->SetTitleSize(0.06);
        xsec_art_hist_fimg->GetYaxis()->SetTitleOffset(0.55);
        xsec_art_hist_fimg->SetFillColor(kRed);
        // xsec_art_hist_fimg->SetFillStyle(3001);
        xsec_art_hist_fimg->GetXaxis()->SetRangeUser(1e-1,1e6);
        xsec_art_hist_fimg->GetYaxis()->SetRangeUser(-4.8,19.);
        xsec_art_hist_fimg->Draw("e2");

        cross_section_hist_e_FIMG->SetTitle("");
        cross_section_hist_e_FIMG->GetXaxis()->SetRangeUser(1e-1,1e6);
        // cross_section_hist_e_FIMG->SetLineWidth(2);
        cross_section_hist_e_FIMG->Draw("same");
        gPad->SetLogx();

        l_ar[1] = new TLegend(0.15, 0.60, 0.45, 0.80);
        l_ar[1]->AddEntry(xsec_art_hist_fimg, "Argon thickness systematic", "f");
        l_ar[1]->AddEntry(cross_section_hist_e_FIMG, "Cross section w/ stat errors", "l");
        l_ar[1]->Draw("same");
    }
    
}

void fill_bkgd_sys(bool plot_sys){

    //////////////////////////////// PTBC
    Double_t bin_edges_ptbc[] = {1.00000000e-01, 1.12201845e-01, 1.25892541e-01, 1.41253754e-01,
       1.58489319e-01, 1.77827941e-01, 1.99526231e-01, 2.23872114e-01,
       2.51188643e-01, 2.81838293e-01, 3.16227766e-01, 3.54813389e-01,
       3.98107171e-01, 4.46683592e-01, 5.01187234e-01, 5.62341325e-01,
       6.30957344e-01, 7.07945784e-01, 7.94328235e-01, 8.91250938e-01,
       1.00000000e+00, 1.12201845e+00, 1.25892541e+00, 1.41253754e+00,
       1.58489319e+00, 1.77827941e+00, 1.99526231e+00, 2.23872114e+00,
       2.51188643e+00, 2.81838293e+00, 3.16227766e+00, 3.54813389e+00,
       3.98107171e+00, 4.46683592e+00, 5.01187234e+00, 5.62341325e+00,
       6.30957344e+00, 7.07945784e+00, 7.94328235e+00, 8.91250938e+00,
       1.00000000e+01, 1.12201845e+01, 1.25892541e+01, 1.41253754e+01,
       1.58489319e+01, 1.77827941e+01, 1.99526231e+01, 2.23872114e+01,
       2.51188643e+01, 2.81838293e+01, 3.16227766e+01, 3.54813389e+01,
       3.98107171e+01, 4.46683592e+01, 5.01187234e+01, 5.62341325e+01,
       6.30957344e+01, 7.07945784e+01, 7.94328235e+01, 8.91250938e+01,
       1.00000000e+02, 1.12201845e+02, 1.25892541e+02, 1.41253754e+02,
       1.58489319e+02, 1.77827941e+02, 1.99526231e+02, 2.23872114e+02,
       2.51188643e+02, 2.81838293e+02, 3.16227766e+02, 3.54813389e+02,
       3.98107171e+02, 4.46683592e+02, 5.01187234e+02, 5.62341325e+02,
       6.30957344e+02, 7.07945784e+02, 7.94328235e+02, 8.91250938e+02,
       1.00000000e+03, 1.12201845e+03, 1.25892541e+03, 1.41253754e+03,
       1.58489319e+03, 1.77827941e+03, 1.99526231e+03, 2.23872114e+03,
       2.51188643e+03, 2.81838293e+03, 3.16227766e+03, 3.54813389e+03,
       3.98107171e+03, 4.46683592e+03, 5.01187234e+03, 5.62341325e+03,
       6.30957344e+03, 7.07945784e+03, 7.94328235e+03, 8.91250938e+03,
       1.00000000e+04, 1.12201845e+04, 1.25892541e+04, 1.41253754e+04,
       1.58489319e+04, 1.77827941e+04, 1.99526231e+04, 2.23872114e+04,
       2.51188643e+04, 2.81838293e+04, 3.16227766e+04, 3.54813389e+04,
       3.98107171e+04, 4.46683592e+04, 5.01187234e+04, 5.62341325e+04,
       6.30957344e+04, 7.07945784e+04, 7.94328235e+04, 8.91250938e+04,
       1.00000000e+05, 1.12201845e+05, 1.25892541e+05, 1.41253754e+05,
       1.58489319e+05, 1.77827941e+05, 1.99526231e+05, 2.23872114e+05,
       2.51188643e+05, 2.81838293e+05, 3.16227766e+05, 3.54813389e+05,
       3.98107171e+05, 4.46683592e+05, 5.01187234e+05, 5.62341325e+05,
       6.30957344e+05, 7.07945784e+05, 7.94328235e+05, 8.91250938e+05,
       1.00000000e+06, 1.12201845e+06, 1.25892541e+06, 1.41253754e+06,
       1.58489319e+06, 1.77827941e+06, 1.99526231e+06, 2.23872114e+06,
       2.51188643e+06, 2.81838293e+06, 3.16227766e+06, 3.54813389e+06,
       3.98107171e+06, 4.46683592e+06, 5.01187234e+06, 5.62341325e+06,
       6.30957344e+06, 7.07945784e+06, 7.94328235e+06, 8.91250938e+06,
       1.00000000e+07, 1.12201845e+07, 1.25892541e+07, 1.41253754e+07,
       1.58489319e+07, 1.77827941e+07, 1.99526231e+07, 2.23872114e+07,
       2.51188643e+07, 2.81838293e+07, 3.16227766e+07, 3.54813389e+07,
       3.98107171e+07, 4.46683592e+07, 5.01187234e+07, 5.62341325e+07,
       6.30957344e+07, 7.07945784e+07, 7.94328235e+07, 8.91250938e+07,
       1.00000000e+08};

    //    , 1.12201845e+08, 1.25892541e+08, 1.41253754e+08,
    //    1.58489319e+08, 1.77827941e+08, 1.99526231e+08, 2.23872114e+08,
    //    2.51188643e+08, 2.81838293e+08, 3.16227766e+08, 3.54813389e+08,
    //    3.98107171e+08, 4.46683592e+08, 5.01187234e+08, 5.62341325e+08,
    //    6.30957344e+08, 7.07945784e+08, 7.94328235e+08, 8.91250938e+08,
    //    1.00000000e+09

    Double_t bin_err_up_ptbc[] = {1.74975385e-04, 1.77789899e-03, 9.71639382e-03, 5.56957747e-03,
       2.69431208e-03, 4.12744606e-03, 1.00495401e-03, 2.45248720e-03,
       2.81277881e-03, 1.34403154e-03, 2.50337706e-03, 1.48388191e-03,
       2.36810348e-03, 1.80907618e-03, 4.28997504e-03, 2.18613267e-03,
       4.72117255e-03, 1.58047925e-03, 3.86693302e-03, 2.68686275e-03,
       1.32161731e-03, 2.46367309e-03, 8.21935907e-03, 1.23130993e-03,
       1.05014331e-02, 1.70658057e-03, 7.27017805e-03, 7.21772294e-03,
       6.69200856e-04, 1.55958526e-03, 1.20836149e-03, 1.65827630e-03,
       1.14569325e-01, 2.90015383e-03, 3.29929218e-03, 3.59862445e-04,
       9.87573958e-04, 1.00189834e-03, 1.95575975e-04, 2.02291816e-04,
       9.63478151e-04, 2.34843618e-04, 9.58569894e-04, 2.86313183e-04,
       9.22101350e-04, 2.06553465e-04, 3.74874447e-04, 5.16159652e-05,
       1.65406810e-04, 9.54722483e-04, 8.91320235e-05, 3.08045127e-04,
       2.17570088e-04, 2.31220078e-04, 1.71652459e-04, 3.07543684e-04,
       4.35111551e-04, 3.60897765e-04, 2.12452386e-04, 5.66901038e-04,
       3.64017238e-04, 1.84380388e-04, 1.45841756e-04, 2.12486645e-04,
       2.31705287e-04, 4.94828979e-04, 7.08063187e-05, 1.93933111e-04,
       1.97142964e-04, 8.29116471e-05, 2.01391949e-04, 6.42524690e-05,
       1.11514450e-04, 2.01948063e-04, 9.98550605e-05, 1.72170732e-04,
       1.31828851e-04, 1.88894054e-04, 1.08855910e-04, 1.07843913e-04,
       2.26549752e-04, 3.33293011e-04, 1.61659234e-04, 2.16744763e-04,
       9.97418349e-05, 6.09398578e-05, 1.25834089e-04, 2.90760530e-04,
       2.70687684e-04, 1.35615236e-04, 9.33381869e-05, 7.29358156e-05,
       2.12653672e-04, 1.02215951e-05, 3.06298044e-04, 5.21523311e-04,
       9.87779370e-05, 4.28191845e-07, 2.50528742e-04, 3.40154235e-04,
       1.41756145e-04, 1.23459159e-04, 4.31182615e-06, 1.85610971e-04,
       1.69432923e-04, 1.82373695e-04, 3.24564091e-05, 1.29538178e-05,
       7.74220279e-06, 1.58485569e-04, 1.90026553e-04, 3.85466059e-04,
       1.22308758e-04, 1.84279145e-04, 5.56658375e-05, 2.06122899e-05,
       3.49506471e-05, 6.64841479e-04, 7.49939882e-04, 1.50512069e-04,
       1.09866898e-04, 1.19606539e-04, 4.40337759e-05, 7.54972331e-05,
       1.79830442e-04, 7.02792755e-05, 6.42393174e-05, 3.12599532e-05,
       5.34998318e-05, 1.90533985e-05, 2.33736622e-05, 1.88211916e-05,
       2.49424727e-05, 7.20867117e-06, 3.96395096e-06, 1.26810714e-05,
       5.64645717e-06, 7.35283219e-06, 1.09353243e-05, 6.19190233e-06,
       8.71933360e-06, 5.72417175e-06, 9.55355899e-06, 7.27569405e-06,
       4.92171027e-06, 7.87184964e-06, 7.44361760e-06, 8.29569112e-06,
       9.16604449e-06, 1.09247109e-05, 1.20869299e-05, 1.44878017e-05,
       1.35292453e-05, 1.28972390e-05, 1.77776500e-05, 1.69761568e-05,
       1.25994188e-05, 1.30911494e-05, 8.97797850e-06, 4.02501867e-06,
       9.35157876e-06, 7.38719985e-06, 6.56560022e-06, 4.13202023e-06,
       1.01763062e-06, 4.43522854e-06, 5.75938723e-07, 6.11262235e-06,
       3.52198404e-06, 4.69817423e-06, 4.69913888e-06, 2.56850110e-06,
       2.12042572e-06, 2.17356540e-06, 2.97497399e-06, 2.17695811e-06,
       3.29274020e-06, 9.21747967e-07, 7.70255985e-07, 1.60965758e-06,
       1.11168331e-06};

    //    , 9.53959687e-07, 8.81867841e-07, 2.99813557e-08,
    //    4.14260140e-07, 1.35251259e-07, 8.19599534e-07, 1.03032699e-06,
    //    5.68604086e-07, 2.89348316e-07, 1.58732444e-06, 1.06449694e-06,
    //    7.68250228e-07, 8.05180113e-07, 8.58082643e-07, 7.19911744e-07,
    //    6.33704798e-07, 1.20902356e-07, 1.07916729e-06, 6.76301102e-08

    Double_t bin_err_low_ptbc[] = {1.20182835e-03, 1.08698661e-02, 5.18117960e-02, 2.74876590e-02,
       1.26554658e-02, 1.83645266e-02, 4.21287600e-03, 9.40441549e-03,
       1.00497835e-02, 4.64644591e-03, 8.67181713e-03, 5.29992799e-03,
       8.72923857e-03, 6.81444422e-03, 1.63956043e-02, 8.49694951e-03,
       1.81442889e-02, 6.04656025e-03, 1.43503067e-02, 9.46519853e-03,
       4.22415558e-03, 8.09713328e-03, 4.10339929e-02, 7.98108564e-03,
       6.74228015e-02, 1.02140496e-02, 3.45399305e-02, 5.90529740e-02,
       5.79941265e-03, 5.27717182e-03, 4.10156228e-03, 5.51545920e-03,
       7.38866274e-01, 1.79333599e-02, 1.39088395e-02, 1.13492572e-03,
       2.99525086e-03, 3.83448891e-03, 5.53181265e-04, 5.85989578e-04,
       3.10987163e-03, 6.71934333e-04, 2.79832956e-03, 8.56647380e-04,
       2.75792274e-03, 5.82680392e-04, 1.09181972e-03, 1.47802969e-04,
       4.74638372e-04, 2.95163067e-03, 2.50677911e-04, 8.86361597e-04,
       6.24257782e-04, 6.61545122e-04, 4.83407044e-04, 8.72657634e-04,
       1.25470223e-03, 1.03134863e-03, 6.04954654e-04, 1.62164536e-03,
       1.05032493e-03, 5.24449784e-04, 4.14219266e-04, 6.06003198e-04,
       6.57166557e-04, 1.41170027e-03, 2.00627829e-04, 5.48968026e-04,
       5.56069911e-04, 2.37502370e-04, 5.74288957e-04, 1.84086222e-04,
       3.16720543e-04, 5.71792682e-04, 2.82181569e-04, 4.86577942e-04,
       3.73688705e-04, 5.34381494e-04, 3.09685392e-04, 3.06406573e-04,
       6.43520016e-04, 9.42353051e-04, 4.57864755e-04, 6.15174379e-04,
       2.82434543e-04, 1.72630920e-04, 3.57488839e-04, 8.24607198e-04,
       7.67728994e-04, 3.84282202e-04, 2.64226767e-04, 2.06530572e-04,
       6.00752223e-04, 2.89240879e-05, 8.65187084e-04, 1.47535079e-03,
       2.79308736e-04, 1.21130919e-06, 7.07265110e-04, 9.59612555e-04,
       4.00126307e-04, 3.48200484e-04, 1.21441838e-05, 5.22547386e-04,
       4.76910542e-04, 5.12222447e-04, 9.12638357e-05, 3.63752288e-05,
       2.17191260e-05, 4.44145312e-04, 5.35644569e-04, 1.08634359e-03,
       3.42732960e-04, 5.15849346e-04, 1.55742577e-04, 5.76303572e-05,
       9.76811793e-05, 1.85674796e-03, 2.09922577e-03, 4.20905656e-04,
       3.06586298e-04, 3.33693558e-04, 1.22758673e-04, 2.10657917e-04,
       5.01277210e-04, 1.95808186e-04, 1.78981457e-04, 8.70558478e-05,
       1.48964214e-04, 5.30513548e-05, 6.50614177e-05, 5.23866175e-05,
       6.94370051e-05, 2.00604505e-05, 1.10296353e-05, 3.52829167e-05,
       1.57093540e-05, 2.04560894e-05, 3.04220717e-05, 1.72255767e-05,
       2.42556340e-05, 1.59229465e-05, 2.65747374e-05, 2.02379294e-05,
       1.36901727e-05, 2.18958336e-05, 2.07049261e-05, 2.30749702e-05,
       2.54968108e-05, 3.03898128e-05, 3.36250053e-05, 4.03068003e-05,
       3.76405217e-05, 3.58831215e-05, 4.94622866e-05, 4.72324561e-05,
       3.50518061e-05, 3.64183139e-05, 2.49758519e-05, 1.11973811e-05,
       2.60153912e-05, 2.05508896e-05, 1.82645431e-05, 1.14943341e-05,
       2.83085002e-06, 1.23373978e-05, 1.60205724e-06, 1.70025563e-05,
       9.79657610e-06, 1.30679261e-05, 1.30705150e-05, 7.14413771e-06,
       5.89779252e-06, 6.04561278e-06, 8.27457286e-06, 6.05497964e-06,
       9.15831739e-06, 2.56370772e-06, 2.14235068e-06, 4.47698096e-06,
       3.09197540e-06};

    //    , 2.65326521e-06, 2.45273749e-06, 8.33875616e-08,
    //    1.15218442e-06, 3.76174938e-07, 2.27954448e-06, 2.86563006e-06,
    //    1.58145199e-06, 8.04761240e-07, 4.41478266e-06, 2.96064260e-06,
    //    2.13670267e-06, 2.23941077e-06, 2.38653693e-06, 2.00224687e-06,
    //    1.76248027e-06, 3.36257369e-07, 3.00141680e-06, 1.88096616e-07

    Int_t num_bins_ptbc = static_cast<Int_t>(sizeof(bin_edges_ptbc)/sizeof(*bin_edges_ptbc) - 1);
    trans_bkgd_hist_ptbc = new TH1D("trans_bkgd_hist_ptbc", "Background errors transmission - ptbc", num_bins_ptbc, bin_edges_ptbc);
    xsec_bkgd_hist_ptbc = new TH1D("xsec_bkgd_hist_ptbc", "Background errors cross section - ptbc", num_bins_ptbc, bin_edges_ptbc);

    for (Int_t i = 1; i < num_bins_ptbc+1; i++)
    {
        Double_t trans_bin_lim_low = transmission_hist_e_PTBC->GetBinContent(20+i) + bin_err_up_ptbc[i-1];
        Double_t trans_bin_lim_up = transmission_hist_e_PTBC->GetBinContent(20+i) - bin_err_low_ptbc[i-1];
        Double_t trans_bin_content = (trans_bin_lim_low + trans_bin_lim_up)/2;
        Double_t trans_bin_error = (trans_bin_lim_low - trans_bin_lim_up)/2;
        trans_bkgd_hist_ptbc->SetBinContent(i, trans_bin_content);
        trans_bkgd_hist_ptbc->SetBinError(i, trans_bin_error);

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
        Double_t xsec_bin_content = (xsec_bin_lim_up + xsec_bin_lim_low)/2;
        Double_t xsec_bin_err = (xsec_bin_lim_up - xsec_bin_lim_low)/2;
        xsec_bkgd_hist_ptbc->SetBinContent(i, xsec_bin_content);
        xsec_bkgd_hist_ptbc->SetBinError(i, xsec_bin_err);
    }

    //////////////////////////////// FIMG
    Double_t bin_edges_fimg[] = {1.00000000e-01, 1.12201845e-01, 1.25892541e-01, 1.41253754e-01,
       1.58489319e-01, 1.77827941e-01, 1.99526231e-01, 2.23872114e-01,
       2.51188643e-01, 2.81838293e-01, 3.16227766e-01, 3.54813389e-01,
       3.98107171e-01, 4.46683592e-01, 5.01187234e-01, 5.62341325e-01,
       6.30957344e-01, 7.07945784e-01, 7.94328235e-01, 8.91250938e-01,
       1.00000000e+00, 1.12201845e+00, 1.25892541e+00, 1.41253754e+00,
       1.58489319e+00, 1.77827941e+00, 1.99526231e+00, 2.23872114e+00,
       2.51188643e+00, 2.81838293e+00, 3.16227766e+00, 3.54813389e+00,
       3.98107171e+00, 4.46683592e+00, 5.01187234e+00, 5.62341325e+00,
       6.30957344e+00, 7.07945784e+00, 7.94328235e+00, 8.91250938e+00,
       1.00000000e+01, 1.12201845e+01, 1.25892541e+01, 1.41253754e+01,
       1.58489319e+01, 1.77827941e+01, 1.99526231e+01, 2.23872114e+01,
       2.51188643e+01, 2.81838293e+01, 3.16227766e+01, 3.54813389e+01,
       3.98107171e+01, 4.46683592e+01, 5.01187234e+01, 5.62341325e+01,
       6.30957344e+01, 7.07945784e+01, 7.94328235e+01, 8.91250938e+01,
       1.00000000e+02, 1.12201845e+02, 1.25892541e+02, 1.41253754e+02,
       1.58489319e+02, 1.77827941e+02, 1.99526231e+02, 2.23872114e+02,
       2.51188643e+02, 2.81838293e+02, 3.16227766e+02, 3.54813389e+02,
       3.98107171e+02, 4.46683592e+02, 5.01187234e+02, 5.62341325e+02,
       6.30957344e+02, 7.07945784e+02, 7.94328235e+02, 8.91250938e+02,
       1.00000000e+03, 1.12201845e+03, 1.25892541e+03, 1.41253754e+03,
       1.58489319e+03, 1.77827941e+03, 1.99526231e+03, 2.23872114e+03,
       2.51188643e+03, 2.81838293e+03, 3.16227766e+03, 3.54813389e+03,
       3.98107171e+03, 4.46683592e+03, 5.01187234e+03, 5.62341325e+03,
       6.30957344e+03, 7.07945784e+03, 7.94328235e+03, 8.91250938e+03,
       1.00000000e+04, 1.12201845e+04, 1.25892541e+04, 1.41253754e+04,
       1.58489319e+04, 1.77827941e+04, 1.99526231e+04, 2.23872114e+04,
       2.51188643e+04, 2.81838293e+04, 3.16227766e+04, 3.54813389e+04,
       3.98107171e+04, 4.46683592e+04, 5.01187234e+04, 5.62341325e+04,
       6.30957344e+04, 7.07945784e+04, 7.94328235e+04, 8.91250938e+04,
       1.00000000e+05, 1.12201845e+05, 1.25892541e+05, 1.41253754e+05,
       1.58489319e+05, 1.77827941e+05, 1.99526231e+05, 2.23872114e+05,
       2.51188643e+05, 2.81838293e+05, 3.16227766e+05, 3.54813389e+05,
       3.98107171e+05, 4.46683592e+05, 5.01187234e+05, 5.62341325e+05,
       6.30957344e+05, 7.07945784e+05, 7.94328235e+05, 8.91250938e+05,
       1.00000000e+06};

    Double_t bin_err_low_fimg[] = {7.44330791e-03, 5.62878727e-03, 5.26444842e-03, 6.75438194e-03,
       3.12791401e-03, 2.90086201e-03, 4.48315669e-03, 3.91170255e-03,
       2.14444987e-03, 2.30664008e-03, 1.13053623e-03, 1.64090160e-03,
       2.05720995e-03, 2.37283602e-03, 2.18010796e-03, 1.93648164e-03,
       1.38560311e-03, 1.35092112e-03, 1.17317283e-03, 1.54036094e-03,
       1.04081453e-03, 1.20651728e-03, 6.39178853e-04, 1.89387109e-05,
       1.09842063e-03, 1.47494875e-03, 8.74458268e-04, 1.23845742e-03,
       5.47558594e-04, 3.84481131e-04, 1.34714309e-03, 8.30113399e-04,
       1.13366202e-03, 9.35199597e-04, 9.66706323e-04, 5.71441755e-04,
       9.54208129e-04, 5.68995270e-04, 8.47901032e-04, 6.16907593e-04,
       1.23996441e-03, 8.67421216e-04, 1.03594350e-03, 2.73848297e-04,
       6.29219978e-04, 1.24356944e-03, 7.70711124e-04, 6.47578149e-04,
       1.23030498e-03, 5.73690005e-04, 6.05492436e-04, 9.83491158e-04,
       1.07551893e-03, 8.16745635e-04, 9.64682134e-04, 9.76865872e-04,
       9.73238967e-04, 3.76589005e-04, 3.29487648e-04, 7.31992660e-04,
       8.61220272e-04, 2.87126833e-04, 8.23604322e-04, 1.36343406e-04,
       8.80332073e-04, 6.43982274e-04, 1.00180520e-03, 4.70811949e-04,
       8.31391457e-04, 3.12675253e-04, 1.72144642e-04, 8.69067034e-04,
       1.33347248e-03, 4.17241134e-05, 5.42569418e-05, 7.48816671e-04,
       1.56404135e-04, 5.04689944e-04, 2.43135019e-05, 1.14578072e-03,
       6.92957834e-05, 5.01273617e-04, 8.90840463e-04, 2.20257346e-04,
       1.73387700e-04, 2.79023653e-04, 5.27487403e-04, 1.05466487e-04,
       8.39149235e-04, 2.10129028e-04, 1.06399849e-03, 6.63144122e-04,
       9.09787300e-04, 2.50709881e-04, 8.60593813e-04, 1.74378632e-04,
       1.93325981e-03, 1.32451436e-03, 7.30178226e-04, 5.04132632e-05,
       1.60487341e-05, 1.02149804e-03, 3.42765165e-04, 9.28974046e-04,
       4.04032425e-04, 2.68252863e-04, 4.93231507e-04, 2.79072213e-04,
       1.02956254e-04, 8.22532754e-04, 1.07982202e-03, 8.80689818e-04,
       5.19607262e-04, 1.31751188e-03, 1.86395227e-04, 8.33640188e-04,
       3.25817010e-04, 3.48080933e-03, 4.47408307e-03, 1.73881958e-03,
       6.04278650e-04, 3.67749595e-04, 1.37849805e-04, 5.90002759e-05,
       6.90638312e-04, 5.08368916e-04, 1.22614701e-04, 2.35397438e-04,
       1.88284630e-04, 3.07398958e-04, 3.13007855e-04, 4.55336069e-05,
       7.13025457e-05, 1.45768354e-04, 1.38653151e-04, 2.26687761e-04,
       1.52132690e-04, 3.35091105e-04, 7.77289033e-06, 3.90340452e-04,
       3.61172261e-04, 5.25950513e-04, 2.33374148e-04, 2.08962209e-04};

    Int_t num_bins_fimg = static_cast<Int_t>(sizeof(bin_edges_fimg)/sizeof(*bin_edges_fimg) - 1);
    trans_bkgd_hist_fimg = new TH1D("trans_bkgd_hist_fimg", "Background errors transmission - fimg", num_bins_fimg, bin_edges_fimg);
    xsec_bkgd_hist_fimg = new TH1D("xsec_bkgd_hist_fimg", "Background errors cross section - fimg", num_bins_fimg, bin_edges_fimg);

    for (Int_t i = 1; i < num_bins_fimg+1; i++)
    {
        Double_t trans_bin_lim_low = transmission_hist_e_FIMG->GetBinContent(20+i);
        Double_t trans_bin_lim_up = transmission_hist_e_FIMG->GetBinContent(20+i) - bin_err_low_fimg[i-1];
        Double_t trans_bin_content = (trans_bin_lim_low + trans_bin_lim_up)/2;
        Double_t trans_bin_error = (trans_bin_lim_low - trans_bin_lim_up)/2;
        trans_bkgd_hist_fimg->SetBinContent(i, trans_bin_content);
        trans_bkgd_hist_fimg->SetBinError(i, trans_bin_error);

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
        Double_t xsec_bin_content = (xsec_bin_lim_up + xsec_bin_lim_low)/2;
        Double_t xsec_bin_err = (xsec_bin_lim_up - xsec_bin_lim_low)/2;
        xsec_bkgd_hist_fimg->SetBinContent(i, xsec_bin_content);
        xsec_bkgd_hist_fimg->SetBinError(i, xsec_bin_err);
    }

    if (plot_sys)
    {
        SetMArEXStyle();
        gStyle->SetCanvasDefW(1000); //600
        gStyle->SetCanvasDefH(500); //500 
        gStyle->SetPadRightMargin(0.04);
        gStyle->SetTitleH(0.1);
        gStyle->SetTitleAlign(33);
        gStyle->SetTitleX(.76);
        // gStyle->SetStats(0);
        
        ////////////////////// PTBC
        TCanvas* c_bkgd_ptbc = new TCanvas("c_bkgd_ptbc"," ");
        TPad* p_bkgd_ptbc[2];
        TLegend* l_bkgd_ptbc[2];
        //Pad - 1 - CF
        c_bkgd_ptbc->cd(0);
        p_bkgd_ptbc[0] = new TPad("p_bkgd_ptbc_0", "p_bkgd_ptbc_0", 0., 0.5, 1., 1.);
        p_bkgd_ptbc[0]->SetFillColor(kWhite);
        p_bkgd_ptbc[0]->SetBottomMargin(0.00001);
        p_bkgd_ptbc[0]->SetTopMargin(0.2);
        p_bkgd_ptbc[0]->SetBorderMode(0);
        p_bkgd_ptbc[0]->Draw();
        p_bkgd_ptbc[0]->cd();
        trans_bkgd_hist_ptbc->SetTitle("Fission Chamber Background Systematic");
        trans_bkgd_hist_ptbc->SetTitleOffset(0.5);
        // X-Axis
        trans_bkgd_hist_ptbc->GetXaxis()->SetLabelOffset(999);
        trans_bkgd_hist_ptbc->GetXaxis()->SetLabelSize(0);
        trans_bkgd_hist_ptbc->GetXaxis()->SetTitle("");
        // Y-Axis
        trans_bkgd_hist_ptbc->GetYaxis()->SetTitle("Transmission");
        trans_bkgd_hist_ptbc->GetYaxis()->SetLabelSize(0.09);
        trans_bkgd_hist_ptbc->GetYaxis()->SetTitleSize(0.09);
        trans_bkgd_hist_ptbc->GetYaxis()->SetTitleOffset(0.45);
        trans_bkgd_hist_ptbc->SetFillColor(kRed);
        // trans_bkgd_hist_ptbc->SetFillStyle(3001);
        trans_bkgd_hist_ptbc->GetXaxis()->SetRangeUser(1e-1,1e8);
        trans_bkgd_hist_ptbc->GetYaxis()->SetRangeUser(0.22,1.25);
        trans_bkgd_hist_ptbc->Draw("e2");

        transmission_hist_e_PTBC->SetTitle("");
        transmission_hist_e_PTBC->GetXaxis()->SetRangeUser(1e-1,1e8);
        // transmission_hist_e_PTBC->SetLineWidth(2);
        transmission_hist_e_PTBC->Draw("same");
        gPad->SetLogx();
        gPad->SetGrid();

        l_bkgd_ptbc[0] = new TLegend(0.15, 0.10, 0.45, 0.30);
        l_bkgd_ptbc[0]->AddEntry(trans_bkgd_hist_ptbc, "Background systematic", "f");
        l_bkgd_ptbc[0]->AddEntry(transmission_hist_e_PTBC, "Transmission w/ stat errors", "l");
        l_bkgd_ptbc[0]->Draw();

        //Pad - 2 - Ar
        c_bkgd_ptbc->cd(0);
        p_bkgd_ptbc[1] = new TPad("p_bkgd_ptbc_1", "p_bkgd_ptbc_1", 0., 0., 1., 0.5);
        p_bkgd_ptbc[1]->SetFillColor(kWhite);
        p_bkgd_ptbc[1]->SetTopMargin(0.00001);
        p_bkgd_ptbc[1]->SetBottomMargin(0.2);
        p_bkgd_ptbc[1]->SetBorderMode(0);
        p_bkgd_ptbc[1]->Draw();
        p_bkgd_ptbc[1]->cd();

        xsec_bkgd_hist_ptbc->SetTitle("");
        // X-Axis
        xsec_bkgd_hist_ptbc->GetXaxis()->SetTitle("Energy (in eV)");
        xsec_bkgd_hist_ptbc->GetXaxis()->SetTitleOffset(1.05); //increase to move down
        xsec_bkgd_hist_ptbc->GetXaxis()->SetLabelSize(0.09);
        xsec_bkgd_hist_ptbc->GetXaxis()->SetTitleSize(0.09);
        // Y-Axis
        xsec_bkgd_hist_ptbc->GetYaxis()->SetTitle("Cross Section (in b)");
        xsec_bkgd_hist_ptbc->GetYaxis()->SetLabelSize(0.09);
        xsec_bkgd_hist_ptbc->GetYaxis()->SetTitleSize(0.09);
        xsec_bkgd_hist_ptbc->GetYaxis()->SetTitleOffset(0.45);
        xsec_bkgd_hist_ptbc->SetFillColor(kRed);
        // xsec_bkgd_hist_ptbc->SetFillStyle(3001);
        xsec_bkgd_hist_ptbc->GetXaxis()->SetRangeUser(1e-1,1e8);
        xsec_bkgd_hist_ptbc->GetYaxis()->SetRangeUser(-4.8,23.);
        xsec_bkgd_hist_ptbc->Draw("e2");

        cross_section_hist_e_PTBC->SetTitle("");
        cross_section_hist_e_PTBC->GetXaxis()->SetRangeUser(1e-1,1e8);
        // cross_section_hist_e_PTBC->SetLineWidth(2);
        cross_section_hist_e_PTBC->Draw("same");
        gPad->SetLogx();

        l_bkgd_ptbc[1] = new TLegend(0.15, 0.70, 0.45, 0.90);
        l_bkgd_ptbc[1]->AddEntry(xsec_bkgd_hist_ptbc, "Background systematic", "f");
        l_bkgd_ptbc[1]->AddEntry(cross_section_hist_e_PTBC, "Cross section w/ stat errors", "l");
        l_bkgd_ptbc[1]->Draw("same");

        ////////////////////// FIMG
        TCanvas* c_bkgd_fimg = new TCanvas("c_bkgd_fimg"," ");
        TPad* p_bkgd_fimg[2];
        TLegend* l_bkgd_fimg[2];
        //Pad - 1 - CF
        c_bkgd_fimg->cd(0);
        p_bkgd_fimg[0] = new TPad("p_bkgd_fimg_0", "p_bkgd_fimg_0", 0., 0.5, 1., 1.);
        p_bkgd_fimg[0]->SetFillColor(kWhite);
        p_bkgd_fimg[0]->SetBottomMargin(0.00001);
        p_bkgd_fimg[0]->SetTopMargin(0.2);
        p_bkgd_fimg[0]->SetBorderMode(0);
        p_bkgd_fimg[0]->Draw();
        p_bkgd_fimg[0]->cd();
        trans_bkgd_hist_fimg->SetTitle("Micromegas Background Systematic");
        trans_bkgd_hist_fimg->SetTitleOffset(0.5);
        // X-Axis
        trans_bkgd_hist_fimg->GetXaxis()->SetLabelOffset(999);
        trans_bkgd_hist_fimg->GetXaxis()->SetLabelSize(0);
        trans_bkgd_hist_fimg->GetXaxis()->SetTitle("");
        // Y-Axis
        trans_bkgd_hist_fimg->GetYaxis()->SetTitle("Transmission");
        trans_bkgd_hist_fimg->GetYaxis()->SetLabelSize(0.09);
        trans_bkgd_hist_fimg->GetYaxis()->SetTitleSize(0.09);
        trans_bkgd_hist_fimg->GetYaxis()->SetTitleOffset(0.45);
        trans_bkgd_hist_fimg->SetFillColor(kRed);
        // trans_bkgd_hist_fimg->SetFillStyle(3001);
        trans_bkgd_hist_fimg->GetXaxis()->SetRangeUser(1e-1,1e6);
        trans_bkgd_hist_fimg->GetYaxis()->SetRangeUser(0.42,1.25);
        trans_bkgd_hist_fimg->Draw("e2");

        transmission_hist_e_FIMG->SetTitle("");
        transmission_hist_e_FIMG->GetXaxis()->SetRangeUser(1e-1,1e6);
        // transmission_hist_e_FIMG->SetLineWidth(2);
        transmission_hist_e_FIMG->Draw("same");
        gPad->SetLogx();
        gPad->SetGrid();

        l_bkgd_fimg[0] = new TLegend(0.15, 0.10, 0.45, 0.30);
        l_bkgd_fimg[0]->AddEntry(trans_bkgd_hist_fimg, "Background systematic", "f");
        l_bkgd_fimg[0]->AddEntry(transmission_hist_e_FIMG, "Transmission w/ stat errors", "l");
        l_bkgd_fimg[0]->Draw();

        //Pad - 2 - Ar
        c_bkgd_fimg->cd(0);
        p_bkgd_fimg[1] = new TPad("p_bkgd_fimg_1", "p_bkgd_fimg_1", 0., 0., 1., 0.5);
        p_bkgd_fimg[1]->SetFillColor(kWhite);
        p_bkgd_fimg[1]->SetTopMargin(0.00001);
        p_bkgd_fimg[1]->SetBottomMargin(0.2);
        p_bkgd_fimg[1]->SetBorderMode(0);
        p_bkgd_fimg[1]->Draw();
        p_bkgd_fimg[1]->cd();

        xsec_bkgd_hist_fimg->SetTitle("");
        // X-Axis
        xsec_bkgd_hist_fimg->GetXaxis()->SetTitle("Energy (in eV)");
        xsec_bkgd_hist_fimg->GetXaxis()->SetTitleOffset(1.05); //increase to move down
        xsec_bkgd_hist_fimg->GetXaxis()->SetLabelSize(0.09);
        xsec_bkgd_hist_fimg->GetXaxis()->SetTitleSize(0.09);
        // Y-Axis
        xsec_bkgd_hist_fimg->GetYaxis()->SetTitle("Cross Section (in b)");
        xsec_bkgd_hist_fimg->GetYaxis()->SetLabelSize(0.09);
        xsec_bkgd_hist_fimg->GetYaxis()->SetTitleSize(0.09);
        xsec_bkgd_hist_fimg->GetYaxis()->SetTitleOffset(0.45);
        xsec_bkgd_hist_fimg->SetFillColor(kRed);
        // xsec_bkgd_hist_fimg->SetFillStyle(3001);
        xsec_bkgd_hist_fimg->GetXaxis()->SetRangeUser(1e-1,1e6);
        xsec_bkgd_hist_fimg->GetYaxis()->SetRangeUser(-4.8,19.);
        xsec_bkgd_hist_fimg->Draw("e2");

        cross_section_hist_e_FIMG->SetTitle("");
        cross_section_hist_e_FIMG->GetXaxis()->SetRangeUser(1e-1,1e6);
        // cross_section_hist_e_FIMG->SetLineWidth(2);
        cross_section_hist_e_FIMG->Draw("same");
        gPad->SetLogx();

        l_bkgd_fimg[1] = new TLegend(0.15, 0.70, 0.45, 0.90);
        l_bkgd_fimg[1]->AddEntry(xsec_bkgd_hist_fimg, "Background systematic", "f");
        l_bkgd_fimg[1]->AddEntry(cross_section_hist_e_FIMG, "Cross section w/ stat errors", "l");
        l_bkgd_fimg[1]->Draw("same");
    }
}

Double_t calc_avg_xsec_evaluations(TH1D* xsec_hist, Int_t min_bin, Int_t max_bin){

    Double_t sum_xsec = 0;
    for (Int_t i = min_bin; i < max_bin+1; i++)
    {
        sum_xsec += xsec_hist->GetBinContent(i);
    }
    Double_t xsec_vec = sum_xsec/(max_bin-min_bin+1);
    return xsec_vec;
}

std::vector<Double_t> calc_avg_xsec_det(TH1D* xsec_hist, TH1D* sys_hist, Int_t min_bin, Int_t max_bin){

    Double_t sum_xsec = 0;
    Double_t sum_xsec_stat_err = 0;
    Double_t sum_xsec_sys_err_up = 0;
    Double_t sum_xsec_sys_err_low = 0;

    for (Int_t i = min_bin; i < max_bin+1; i++)
    {
        sum_xsec += xsec_hist->GetBinContent(20+i);
        sum_xsec_stat_err += xsec_hist->GetBinError(20+i) * xsec_hist->GetBinError(20+i);
        Double_t sys_err_up = sys_hist->GetBinError(i) + sys_hist->GetBinContent(i) - xsec_hist->GetBinContent(20+i);
        sum_xsec_sys_err_up += sys_err_up * sys_err_up;
        Double_t sys_err_low = sys_hist->GetBinError(i) - sys_hist->GetBinContent(i) + xsec_hist->GetBinContent(20+i);
        sum_xsec_sys_err_low += sys_err_low * sys_err_low;
    }

    std::vector<Double_t> xsec_vec;
    xsec_vec.push_back(sum_xsec/(max_bin-min_bin+1));
    xsec_vec.push_back(std::sqrt(sum_xsec_stat_err)/(max_bin-min_bin+1));
    xsec_vec.push_back(std::sqrt(sum_xsec_sys_err_up)/(max_bin-min_bin+1));
    xsec_vec.push_back(std::sqrt(sum_xsec_sys_err_low)/(max_bin-min_bin+1));

    return xsec_vec;
}

void combine_sys(bool plot_hists){
      
    trans_combined_sys_hist_ptbc = (TH1D*)trans_bkgd_hist_ptbc->Clone("trans_combined_sys_hist_ptbc");
    trans_combined_sys_hist_fimg = (TH1D*)xsec_bkgd_hist_fimg->Clone("trans_combined_sys_hist_fimg");
    xsec_combined_sys_hist_ptbc = (TH1D*)trans_bkgd_hist_ptbc->Clone("xsec_combined_sys_hist_ptbc");
    xsec_combined_sys_hist_fimg = (TH1D*)xsec_bkgd_hist_fimg->Clone("xsec_combined_sys_hist_fimg");

    Int_t num_bins_ptbc = trans_bkgd_hist_ptbc->GetNbinsX();
    Int_t num_bins_fimg = xsec_bkgd_hist_fimg->GetNbinsX();

    for (Int_t i = 1; i < num_bins_ptbc+1; i++)
    {
        Double_t trans_val = transmission_hist_e_PTBC->GetBinContent(20+i);
        Double_t trans_bkgd_err_up = trans_bkgd_hist_ptbc->GetBinError(i) + trans_bkgd_hist_ptbc->GetBinContent(i) - trans_val;
        Double_t trans_bkgd_err_low = trans_bkgd_hist_ptbc->GetBinError(i) - trans_bkgd_hist_ptbc->GetBinContent(i) + trans_val;
        Double_t trans_cft_err_up = trans_cft_hist_ptbc->GetBinError(20+i) + trans_cft_hist_ptbc->GetBinContent(20+i) - trans_val;
        Double_t trans_cft_err_low = trans_cft_hist_ptbc->GetBinError(20+i) - trans_cft_hist_ptbc->GetBinContent(20+i) + trans_val;

        Double_t xsec_val = cross_section_hist_e_PTBC->GetBinContent(20+i);
        Double_t xsec_bkgd_err_up = xsec_bkgd_hist_ptbc->GetBinError(i) + xsec_bkgd_hist_ptbc->GetBinContent(i) - xsec_val;
        Double_t xsec_bkgd_err_low = xsec_bkgd_hist_ptbc->GetBinError(i) - xsec_bkgd_hist_ptbc->GetBinContent(i) + xsec_val;
        Double_t xsec_cft_err_up = xsec_cft_hist_ptbc->GetBinError(20+i) + xsec_cft_hist_ptbc->GetBinContent(20+i) - xsec_val;
        Double_t xsec_cft_err_low = xsec_cft_hist_ptbc->GetBinError(20+i) - xsec_cft_hist_ptbc->GetBinContent(20+i) + xsec_val;
        Double_t xsec_art_err_up = xsec_art_hist_ptbc->GetBinError(20+i) + xsec_art_hist_ptbc->GetBinContent(20+i) - xsec_val;
        Double_t xsec_art_err_low = xsec_art_hist_ptbc->GetBinError(20+i) - xsec_art_hist_ptbc->GetBinContent(20+i) + xsec_val;

        Double_t trans_tot_err_up = std::sqrt( trans_bkgd_err_up*trans_bkgd_err_up + trans_cft_err_up*trans_cft_err_up );
        Double_t trans_tot_err_low = std::sqrt( trans_bkgd_err_low*trans_bkgd_err_low + trans_cft_err_low*trans_cft_err_low );
        Double_t xsec_tot_err_up = std::sqrt( xsec_bkgd_err_up*xsec_bkgd_err_up + xsec_cft_err_up*xsec_cft_err_up + xsec_art_err_up*xsec_art_err_up );
        Double_t xsec_tot_err_low = std::sqrt( xsec_bkgd_err_low*xsec_bkgd_err_low + xsec_cft_err_low*xsec_cft_err_low + xsec_art_err_low*xsec_art_err_low );

        trans_combined_sys_hist_ptbc->SetBinContent(i, trans_val + (trans_tot_err_up - trans_tot_err_low)/2);
        trans_combined_sys_hist_ptbc->SetBinError(i, (trans_tot_err_up + trans_tot_err_low)/2);

        xsec_combined_sys_hist_ptbc->SetBinContent(i, xsec_val + (xsec_tot_err_up - xsec_tot_err_low)/2);
        xsec_combined_sys_hist_ptbc->SetBinError(i, (xsec_tot_err_up + xsec_tot_err_low)/2);

    }

    for (Int_t i = 1; i < num_bins_fimg+1; i++)
    {
        Double_t trans_val = transmission_hist_e_FIMG->GetBinContent(20+i);
        Double_t trans_bkgd_err_up = trans_bkgd_hist_fimg->GetBinError(i) + trans_bkgd_hist_fimg->GetBinContent(i) - trans_val;
        Double_t trans_bkgd_err_low = trans_bkgd_hist_fimg->GetBinError(i) - trans_bkgd_hist_fimg->GetBinContent(i) + trans_val;
        Double_t trans_cft_err_up = trans_cft_hist_fimg->GetBinError(20+i) + trans_cft_hist_fimg->GetBinContent(20+i) - trans_val;
        Double_t trans_cft_err_low = trans_cft_hist_fimg->GetBinError(20+i) - trans_cft_hist_fimg->GetBinContent(20+i) + trans_val;

        Double_t xsec_val = cross_section_hist_e_FIMG->GetBinContent(20+i);
        Double_t xsec_bkgd_err_up = xsec_bkgd_hist_fimg->GetBinError(i) + xsec_bkgd_hist_fimg->GetBinContent(i) - xsec_val;
        Double_t xsec_bkgd_err_low = xsec_bkgd_hist_fimg->GetBinError(i) - xsec_bkgd_hist_fimg->GetBinContent(i) + xsec_val;
        Double_t xsec_cft_err_up = xsec_cft_hist_fimg->GetBinError(20+i) + xsec_cft_hist_fimg->GetBinContent(20+i) - xsec_val;
        Double_t xsec_cft_err_low = xsec_cft_hist_fimg->GetBinError(20+i) - xsec_cft_hist_fimg->GetBinContent(20+i) + xsec_val;
        Double_t xsec_art_err_up = xsec_art_hist_fimg->GetBinError(20+i) + xsec_art_hist_fimg->GetBinContent(20+i) - xsec_val;
        Double_t xsec_art_err_low = xsec_art_hist_fimg->GetBinError(20+i) - xsec_art_hist_fimg->GetBinContent(20+i) + xsec_val;

        Double_t trans_tot_err_up = std::sqrt( trans_bkgd_err_up*trans_bkgd_err_up + trans_cft_err_up*trans_cft_err_up );
        Double_t trans_tot_err_low = std::sqrt( trans_bkgd_err_low*trans_bkgd_err_low + trans_cft_err_low*trans_cft_err_low );
        Double_t xsec_tot_err_up = std::sqrt( xsec_bkgd_err_up*xsec_bkgd_err_up + xsec_cft_err_up*xsec_cft_err_up + xsec_art_err_up*xsec_art_err_up );
        Double_t xsec_tot_err_low = std::sqrt( xsec_bkgd_err_low*xsec_bkgd_err_low + xsec_cft_err_low*xsec_cft_err_low + xsec_art_err_low*xsec_art_err_low );

        trans_combined_sys_hist_fimg->SetBinContent(i, trans_val + (trans_tot_err_up - trans_tot_err_low)/2);
        trans_combined_sys_hist_fimg->SetBinError(i, (trans_tot_err_up + trans_tot_err_low)/2);

        xsec_combined_sys_hist_fimg->SetBinContent(i, xsec_val + (xsec_tot_err_up - xsec_tot_err_low)/2);
        xsec_combined_sys_hist_fimg->SetBinError(i, (xsec_tot_err_up + xsec_tot_err_low)/2);
    }

    Double_t avg_xsec_endf = calc_avg_xsec_evaluations(endf_xsec_hist, 61, 80);
    Double_t avg_xsec_jendl = calc_avg_xsec_evaluations(jendl_xsec_hist, 61, 80);
    std::vector<Double_t> avg_xsec_ptbc = calc_avg_xsec_det(cross_section_hist_e_PTBC, xsec_combined_sys_hist_ptbc, 41, 60);
    std::vector<Double_t> avg_xsec_fimg = calc_avg_xsec_det(cross_section_hist_e_FIMG, xsec_combined_sys_hist_fimg, 41, 60);

    cout << "Avg xsec 10eV - 100eV fimg = " << avg_xsec_fimg.at(0) << " +- " << avg_xsec_fimg.at(1) << " +" << avg_xsec_fimg.at(2) << " -" << avg_xsec_fimg.at(3) << endl;
    cout << "Avg xsec 10eV - 100eV ptbc = " << avg_xsec_ptbc.at(0) << " +- " << avg_xsec_ptbc.at(1) << " +" << avg_xsec_ptbc.at(2) << " -" << avg_xsec_ptbc.at(3) << endl;
    cout << "Avg xsec 10eV - 100eV endf = " << avg_xsec_endf << endl;
    cout << "Avg xsec 10eV - 100eV jendl = " << avg_xsec_jendl << endl;

    if (plot_hists)
    {
        SetMArEXStyle();
        gStyle->SetPadRightMargin(0.04);
        // gStyle->SetPadBottomMargin(0.17);
        gStyle->SetCanvasDefW(800); //600
        gStyle->SetCanvasDefH(400); //500 
        // gStyle->SetTitleH(0.1);
        // gStyle->SetTitleAlign(33);
        // gStyle->SetTitleX(.76);
        // gStyle->SetStats(0);
        
        TCanvas* c_tot[4];
        TPad* p_tot[4];
        TLegend* l_tot[4];

        //canvas - 1 - PTBC - xsec
        c_tot[0] = new TCanvas("c_tot_0"," ");
        c_tot[0]->cd(0);
        p_tot[0] = new TPad("p_tot_0", "p_tot_0", 0., 0., 1., 1.);
        p_tot[0]->SetFillColor(kWhite);
        p_tot[0]->SetBottomMargin(0.15);
        p_tot[0]->SetTopMargin(0.1);
        p_tot[0]->SetBorderMode(0);
        p_tot[0]->Draw();
        p_tot[0]->cd();
        xsec_combined_sys_hist_ptbc->SetTitle("Argon - Fission Chamber");
        // X-Axis
        xsec_combined_sys_hist_ptbc->GetXaxis()->SetTitle("Energy (in eV)");
        xsec_combined_sys_hist_ptbc->GetXaxis()->SetTitleOffset(1.2); //increase to move down
        xsec_combined_sys_hist_ptbc->GetXaxis()->SetLabelSize(0.06);
        xsec_combined_sys_hist_ptbc->GetXaxis()->SetTitleSize(0.06);
        // Y-Axis
        xsec_combined_sys_hist_ptbc->GetYaxis()->SetTitle("Cross Section (in b)");
        xsec_combined_sys_hist_ptbc->GetYaxis()->SetLabelSize(0.06);
        xsec_combined_sys_hist_ptbc->GetYaxis()->SetTitleSize(0.06);
        xsec_combined_sys_hist_ptbc->GetYaxis()->SetTitleOffset(0.55);
        xsec_combined_sys_hist_ptbc->SetFillColor(kCyan);
        // xsec_combined_sys_hist_ptbc->SetFillStyle(3001);
        xsec_combined_sys_hist_ptbc->GetXaxis()->SetRangeUser(1e-1,1e8);
        xsec_combined_sys_hist_ptbc->GetYaxis()->SetRangeUser(-4.8,23.);
        xsec_combined_sys_hist_ptbc->Draw("e2");

        cross_section_hist_e_PTBC->SetTitle("");
        cross_section_hist_e_PTBC->GetXaxis()->SetRangeUser(1e-1,1e8);
        // cross_section_hist_e_PTBC->SetLineWidth(2);
        cross_section_hist_e_PTBC->Draw("same");

        endf_xsec_hist->SetLineColor(2);
        endf_xsec_hist->SetLineWidth(2);
        endf_xsec_hist->GetXaxis()->SetRangeUser(1e-1,1e8);
        endf_xsec_hist->Draw("][SAME");

        jendl_xsec_hist->SetLineColor(1);
        jendl_xsec_hist->SetLineWidth(2);
        jendl_xsec_hist->GetXaxis()->SetRangeUser(1e-1,1e8);
        jendl_xsec_hist->Draw("][SAME");

        gPad->SetLogx();
        gPad->SetGrid();

        l_tot[0] = new TLegend(0.15, 0.60, 0.45, 0.80);
        l_tot[0]->AddEntry(xsec_combined_sys_hist_ptbc, "Total systematic error", "f");
        l_tot[0]->AddEntry(cross_section_hist_e_PTBC, "Cross Section w/ stat errors", "l");
        l_tot[0]->AddEntry(endf_xsec_hist,"ENDF","l");
        l_tot[0]->AddEntry(jendl_xsec_hist,"JENDL-5","l");
        l_tot[0]->Draw();

        // c_tot[0]->Print("../plots/results_plots/xsec_ar_tank_ptbc.png");

        //canvas - 2 - FIMG - xsec
        c_tot[1] = new TCanvas("c_tot_1"," ");
        c_tot[1]->cd(0);
        p_tot[1] = new TPad("p_tot_1", "p_tot_1", 0., 0., 1., 1.);
        p_tot[1]->SetFillColor(kWhite);
        p_tot[1]->SetTopMargin(0.1);
        p_tot[1]->SetBottomMargin(0.15);
        p_tot[1]->SetBorderMode(0);
        p_tot[1]->Draw();
        p_tot[1]->cd();

        xsec_combined_sys_hist_fimg->SetTitle("Argon - Micromegas");
        // X-Axis
        xsec_combined_sys_hist_fimg->GetXaxis()->SetTitle("Energy (in eV)");
        xsec_combined_sys_hist_fimg->GetXaxis()->SetTitleOffset(1.2); //increase to move up
        xsec_combined_sys_hist_fimg->GetXaxis()->SetLabelSize(0.06);
        xsec_combined_sys_hist_fimg->GetXaxis()->SetTitleSize(0.06);
        // Y-Axis
        xsec_combined_sys_hist_fimg->GetYaxis()->SetTitle("Cross Section (in b)");
        xsec_combined_sys_hist_fimg->GetYaxis()->SetLabelSize(0.06);
        xsec_combined_sys_hist_fimg->GetYaxis()->SetTitleSize(0.06);
        xsec_combined_sys_hist_fimg->GetYaxis()->SetTitleOffset(0.55);
        xsec_combined_sys_hist_fimg->SetFillColor(kCyan);
        // xsec_combined_sys_hist_fimg->SetFillStyle(3001);
        xsec_combined_sys_hist_fimg->GetXaxis()->SetRangeUser(1e-1,1e6);
        xsec_combined_sys_hist_fimg->GetYaxis()->SetRangeUser(-4.8,23.);
        xsec_combined_sys_hist_fimg->Draw("e2");

        cross_section_hist_e_FIMG->SetTitle("");
        cross_section_hist_e_FIMG->GetXaxis()->SetRangeUser(1e-1,1e6);
        // cross_section_hist_e_FIMG->SetLineWidth(2);
        cross_section_hist_e_FIMG->Draw("same");

        endf_xsec_hist->SetLineColor(2);
        endf_xsec_hist->SetLineWidth(2);
        endf_xsec_hist->GetXaxis()->SetRangeUser(1e-1,1e6);
        endf_xsec_hist->Draw("][SAME");

        jendl_xsec_hist->SetLineColor(1);
        jendl_xsec_hist->SetLineWidth(2);
        jendl_xsec_hist->GetXaxis()->SetRangeUser(1e-1,1e6);
        jendl_xsec_hist->Draw("][SAME");

        gPad->SetLogx();
        gPad->SetGrid();

        l_tot[1] = new TLegend(0.15, 0.60, 0.45, 0.80);
        l_tot[1]->AddEntry(xsec_combined_sys_hist_fimg, "Total systematic error", "f");
        l_tot[1]->AddEntry(cross_section_hist_e_FIMG, "Cross section w/ stat errors", "l");
        l_tot[1]->AddEntry(endf_xsec_hist,"ENDF","l");
        l_tot[1]->AddEntry(jendl_xsec_hist,"JENDL-5","l");
        l_tot[1]->Draw("same");

        // c_tot[1]->Print("../plots/results_plots/xsec_ar_tank_fimg.png");

        ////////////////////// transmission

        //canvas - 3 - PTBC - trans
        c_tot[2] = new TCanvas("c_tot_2"," ");
        c_tot[2]->cd(0);
        p_tot[2] = new TPad("p_tot_2", "p_tot_2", 0., 0., 1., 1.);
        p_tot[2]->SetFillColor(kWhite);
        p_tot[2]->SetBottomMargin(0.15);
        p_tot[2]->SetTopMargin(0.1);
        p_tot[2]->SetBorderMode(0);
        p_tot[2]->Draw();
        p_tot[2]->cd();
        trans_combined_sys_hist_ptbc->SetTitle("Argon Tank - Fission Chamber");
        // X-Axis
        trans_combined_sys_hist_ptbc->GetXaxis()->SetTitle("Energy (in eV)");
        trans_combined_sys_hist_ptbc->GetXaxis()->SetTitleOffset(1.2); //increase to move down
        trans_combined_sys_hist_ptbc->GetXaxis()->SetLabelSize(0.06);
        trans_combined_sys_hist_ptbc->GetXaxis()->SetTitleSize(0.06);
        // Y-Axis
        trans_combined_sys_hist_ptbc->GetYaxis()->SetTitle("Transmission");
        trans_combined_sys_hist_ptbc->GetYaxis()->SetLabelSize(0.06);
        trans_combined_sys_hist_ptbc->GetYaxis()->SetTitleSize(0.06);
        trans_combined_sys_hist_ptbc->GetYaxis()->SetTitleOffset(0.65);
        trans_combined_sys_hist_ptbc->SetFillColor(kCyan);
        // trans_combined_sys_hist_ptbc->SetFillStyle(3001);
        trans_combined_sys_hist_ptbc->GetXaxis()->SetRangeUser(1e-1,1e8);
        trans_combined_sys_hist_ptbc->GetYaxis()->SetRangeUser(0.3,1.25);
        trans_combined_sys_hist_ptbc->Draw("e2");

        transmission_hist_e_PTBC->SetTitle("");
        transmission_hist_e_PTBC->GetXaxis()->SetRangeUser(1e-1,1e8);
        // transmission_hist_e_PTBC->SetLineWidth(2);
        transmission_hist_e_PTBC->Draw("same");

        endf_trans_hist->SetLineColor(2);
        endf_trans_hist->SetLineWidth(2);
        endf_trans_hist->GetXaxis()->SetRangeUser(1e-1,1e8);
        endf_trans_hist->Draw("][SAME");

        jendl_trans_hist->SetLineColor(1);
        jendl_trans_hist->SetLineWidth(2);
        jendl_trans_hist->GetXaxis()->SetRangeUser(1e-1,1e8);
        jendl_trans_hist->Draw("][SAME");
        
        gPad->SetLogx();
        gPad->SetGrid();

        l_tot[2] = new TLegend(0.15, 0.20, 0.45, 0.40);
        l_tot[2]->AddEntry(trans_combined_sys_hist_ptbc, "Total systematic error", "f");
        l_tot[2]->AddEntry(transmission_hist_e_PTBC, "Transmission w/ stat errors", "l");
        l_tot[2]->AddEntry(endf_trans_hist,"ENDF","l");
        l_tot[2]->AddEntry(jendl_trans_hist,"JENDL-5","l");
        l_tot[2]->Draw();

        // c_tot[2]->Print("../plots/results_plots/transmission_ar_tank_ptbc.png");

        //canvas - 4 - FIMG - trans
        c_tot[3] = new TCanvas("c_tot_3"," ");
        c_tot[3]->cd(0);
        p_tot[3] = new TPad("p_tot_3", "p_tot_3", 0., 0., 1., 1.);
        p_tot[3]->SetFillColor(kWhite);
        p_tot[3]->SetTopMargin(0.1);
        p_tot[3]->SetBottomMargin(0.15);
        p_tot[3]->SetBorderMode(0);
        p_tot[3]->Draw();
        p_tot[3]->cd();

        trans_combined_sys_hist_fimg->SetTitle("Argon Tank - Micromegas");
        // X-Axis
        trans_combined_sys_hist_fimg->GetXaxis()->SetTitle("Energy (in eV)");
        trans_combined_sys_hist_fimg->GetXaxis()->SetTitleOffset(1.2); //increase to move up
        trans_combined_sys_hist_fimg->GetXaxis()->SetLabelSize(0.06);
        trans_combined_sys_hist_fimg->GetXaxis()->SetTitleSize(0.06);
        // Y-Axis
        trans_combined_sys_hist_fimg->GetYaxis()->SetTitle("Transmission");
        trans_combined_sys_hist_fimg->GetYaxis()->SetLabelSize(0.06);
        trans_combined_sys_hist_fimg->GetYaxis()->SetTitleSize(0.06);
        trans_combined_sys_hist_fimg->GetYaxis()->SetTitleOffset(0.65);
        trans_combined_sys_hist_fimg->SetFillColor(kCyan);
        // trans_combined_sys_hist_fimg->SetFillStyle(3001);
        trans_combined_sys_hist_fimg->GetXaxis()->SetRangeUser(1e-1,1e6);
        trans_combined_sys_hist_fimg->GetYaxis()->SetRangeUser(0.3,1.25);
        trans_combined_sys_hist_fimg->Draw("e2");

        transmission_hist_e_FIMG->SetTitle("");
        transmission_hist_e_FIMG->GetXaxis()->SetRangeUser(1e-1,1e6);
        // transmission_hist_e_FIMG->SetLineWidth(2);
        transmission_hist_e_FIMG->Draw("same");

        endf_trans_hist->SetLineColor(2);
        endf_trans_hist->SetLineWidth(2);
        endf_trans_hist->GetXaxis()->SetRangeUser(1e-1,1e6);
        endf_trans_hist->Draw("][SAME");

        jendl_trans_hist->SetLineColor(1);
        jendl_trans_hist->SetLineWidth(2);
        jendl_trans_hist->GetXaxis()->SetRangeUser(1e-1,1e6);
        jendl_trans_hist->Draw("][SAME");

        gPad->SetLogx();
        gPad->SetGrid();

        l_tot[3] = new TLegend(0.15, 0.20, 0.45, 0.40);
        l_tot[3]->AddEntry(trans_combined_sys_hist_fimg, "Total systematic error", "f");
        l_tot[3]->AddEntry(transmission_hist_e_FIMG, "Transmission w/ stat errors", "l");
        l_tot[3]->AddEntry(endf_trans_hist,"ENDF","l");
        l_tot[3]->AddEntry(jendl_trans_hist,"JENDL-5","l");
        l_tot[3]->Draw("same");

        // c_tot[3]->Print("../plots/results_plots/transmission_ar_tank_fimg.png");
    }
     
}

    
// }

void ar_tank_ana(){

    transmission_hist_e_PTBC = retriveHistograms("../rootFiles/trans_xsec_hists_ar_bottle_full_20bpd.root", "transmission_hist_e_PTBC");
    transmission_hist_e_FIMG = retriveHistograms("../rootFiles/trans_xsec_hists_ar_bottle_full_20bpd.root", "transmission_hist_e_FIMG");
    cross_section_hist_e_PTBC = retriveHistograms("../rootFiles/trans_xsec_hists_ar_bottle_full_20bpd.root", "cross_section_hist_e_PTBC");
    cross_section_hist_e_FIMG = retriveHistograms("../rootFiles/trans_xsec_hists_ar_bottle_full_20bpd.root", "cross_section_hist_e_FIMG");
    endf_trans_hist = retriveHistograms("../rootFiles/trans_xsec_hists_ar_bottle_full_20bpd.root", "natAr_trans_hist_20bpd_endf");
    endf_xsec_hist = retriveHistograms("../rootFiles/trans_xsec_hists_ar_bottle_full_20bpd.root", "natAr_xsec_hist_20bpd_endf");
    jendl_trans_hist = retriveHistograms("../rootFiles/trans_xsec_hists_ar_bottle_full_20bpd.root", "natAr_trans_hist_20bpd_jendl");
    jendl_xsec_hist = retriveHistograms("../rootFiles/trans_xsec_hists_ar_bottle_full_20bpd.root", "natAr_xsec_hist_20bpd_jendl");

    fill_cf_thickness_sys(false);
    fill_ar_thickness_sys(false);
    fill_bkgd_sys(false);
    combine_sys(true);
}