/**
 * @file plotsFromRootFile.C
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

// const static Int_t num_hist = 2
// TH1D* h[num_hist];
// TH1D* h_old[num_hist];
// Int_t histCounter = 0;
Int_t plot_index=0;

TCanvas *c_subtractionPlots[2];
TLegend *l_subtractionPlots[2];

// void retriveHistograms(const char *fname){
//     TKey *key;
//     TFile *hist_file = TFile::Open(fname, "READ");
//     if (!hist_file || hist_file->IsZombie()) {
//         cout << "Unable to open " << fname << " for reading..." <<endl;
//         return;
//     }

//     // tof_hist_filter_in = (TH1D*)hist_file->Get("tof_hist_filter_in");
//     // energy_hist_filter_in = (TH1D*)hist_file->Get("energy_hist_filter_in");
//     // tof_hist_filter_out = (TH1D*)hist_file->Get("tof_hist_filter_out");
//     // energy_hist_filter_out = (TH1D*)hist_file->Get("energy_hist_filter_out");

//     // trans_hist_fOut = (TH1D*)hist_file->Get("trans_hist_fOut");
//     // trans_hist_fIn = (TH1D*)hist_file->Get("trans_hist_fIn");
//     // trans_hist_fOut_endf = (TH1D*)hist_file->Get("trans_hist_fOut_endf");
//     // trans_hist_fIn_endf = (TH1D*)hist_file->Get("trans_hist_fIn_endf");
//     h[histCounter] = (TH1D*)hist_file->Get("transmission_hist_e");
//     histCounter++;
//     // cross_section_hist_e = (TH1D*)hist_file->Get("cross_section_hist_e");

//     // hist_file->Close();
// }

TH1D* Get_1D_hist(const char *fname, const char *hist_name){
    
    TFile* hist_file = TFile::Open(fname, "READ");
    TH1D* hist_new = (TH1D*)hist_file->Get(hist_name);

    return hist_new;
}

// void retriveHistogramsChangeBPD(const char *fname, const char *hist_name, Int_t bpd_old, Int_t bpd_new){

//     // if(bpd_old < bpd_new) {

//     // }

//     // TKey *key;
//     cout << "Opening File " << fname << endl;
//     TFile *hist_file = TFile::Open(fname, "READ");
//     if (!hist_file || hist_file->IsZombie()) {
//         cout << "Unable to open " << fname << " for reading..." << endl;
//         return;
//     }

//     h_old[histCounter] = (TH1D*)hist_file->Get(hist_name);
//     Int_t num_bins_old = h_old[histCounter]->GetNbinsX();
//     cout << "Number of old bins = " << num_bins_old << endl;
//     Int_t max_sum_bin_count = (bpd_old/bpd_new); // Number of old bins that need to be summed to make a new bin
//     Int_t num_bins_new = num_bins_old / max_sum_bin_count;
//     cout << "Number of new bins = " << num_bins_new << endl;

//     Double_t bin_edges_new[num_bins_new + 1];
//     Double_t bin_content[num_bins_new];
//     Double_t bin_error[num_bins_new];
//     Int_t sum_bin_count = 0;
//     Int_t new_bin_counter = 0;
//     for (int i = 1; i < num_bins_old + 1; i++)
//     {
//         if (i == 1)
//         {
//             bin_edges_new[new_bin_counter] = h_old[histCounter]->GetXaxis()->GetBinLowEdge(i);
//             bin_content[new_bin_counter] = 0;
//             bin_error[new_bin_counter] = 0;
//         }
//         bin_content[new_bin_counter] += h_old[histCounter]->GetBinContent(i);
//         bin_error[new_bin_counter] += h_old[histCounter]->GetBinError(i) * h_old[histCounter]->GetBinError(i);
//         sum_bin_count++;

//         if (i == num_bins_old)
//         {
//             // sum_bin_count = 0;
//             bin_error[new_bin_counter] = sqrt(bin_error[new_bin_counter])/max_sum_bin_count;
//             bin_content[new_bin_counter] = bin_content[new_bin_counter]/max_sum_bin_count;
//             bin_edges_new[new_bin_counter+1] = h_old[histCounter]->GetXaxis()->GetBinUpEdge(i);
//         }
//         else if (sum_bin_count == max_sum_bin_count)
//         {
//             sum_bin_count = 0;
//             bin_error[new_bin_counter] = sqrt(bin_error[new_bin_counter])/max_sum_bin_count;
//             bin_content[new_bin_counter] = bin_content[new_bin_counter]/max_sum_bin_count;
//             new_bin_counter++;
//             bin_edges_new[new_bin_counter] = h_old[histCounter]->GetXaxis()->GetBinUpEdge(i);
//             bin_content[new_bin_counter] = 0;
//             bin_error[new_bin_counter] = 0;
//         } 
//     }
    
//     h[histCounter] = new TH1D(Form("h_%i", histCounter),"Transmission Hist",num_bins_new,bin_edges_new);
//     for(int i = 0; i < num_bins_new; i++){
//         h[histCounter]->SetBinContent(i+1, bin_content[i]);
//         h[histCounter]->SetBinError(i+1, bin_error[i]);
//     }
//     histCounter++;
// }

TH1D* subtractHists(TH1D* h1, TH1D* h2){
    TH1D *h_sub = (TH1D*)h1->Clone();
    h_sub->Add(h2, -1);
    Int_t num_bins = h1->GetNbinsX();
    for (int i = 0; i < num_bins; i++){
        // Double_t h1_bin_content = h1->GetBinContent(i+1);
        // Double_t h2_bin_content = h2->GetBinContent(i+1);
        Double_t h1_err = h1->GetBinError(i+1);
        Double_t h2_err = h2->GetBinError(i+1);
        Double_t hist_bin_error = sqrt( h1_err * h1_err + h2_err * h2_err);
        // h_sub->SetBinContent(i+1, h1_bin_content - h2_bin_content);
        h_sub->SetBinError(i+1, hist_bin_error);
    }
    return h_sub;
}

TH1D* subtractHists_theoryEval(TH1D* h1, TH1D* h2){
    TH1D *h_sub = (TH1D*)h1->Clone();
    Int_t num_bins = h1->GetNbinsX();
    for (int i = 0; i < num_bins; i++){
        Double_t h1_bin_content = h1->GetBinContent(i+1);
        Double_t h2_bin_content = h2->GetBinContent(i+1);
        if (h1_bin_content == 0 || h2_bin_content == 0)
        {
            h_sub->SetBinContent(i+1, 0);
            h_sub->SetBinError(i+1, 0);
            continue;
        }
        Double_t h1_err = h1->GetBinError(i+1);
        h_sub->SetBinContent(i+1, h1_bin_content - h2_bin_content);
        h_sub->SetBinError(i+1, h1_err);
    }
    return h_sub;
}

void plot_subtraction_plots(TH1D* hist_1, const char* x_label, const char* legend_1, TH1D* hist_2, const char* legend_2, TH1D* subtraction_plot, const char* plot_title, const char* output_file_name, Double_t legend_coord[4], Double_t max_e){

    TLine* zero_line = new TLine(1e-1, 0., max_e, 0.);
    zero_line->SetLineWidth(2);
    zero_line->SetLineColor(1);
    zero_line->SetLineStyle(2);

    //Plotting
    SetMArEXStyle();
    
    gStyle->SetStatX(0.27);
    gStyle->SetStatY(0.9);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.17);

    gStyle->SetCanvasDefW(1000); //600
    gStyle->SetCanvasDefH(500); //500
    gStyle->SetPadRightMargin(0.05);

    TCanvas *canv = new TCanvas(Form("c%d", plot_index)," ");
    canv->cd();
    canv->Draw();

    TPad *p_upper = new TPad(Form("p_upper_%d", plot_index), Form("p_upper_%d", plot_index), 0., 0.4, 1., 1.);
    p_upper->SetFillColor(kWhite);
    p_upper->SetBottomMargin(0.00001);
    p_upper->SetBorderMode(0);
    p_upper->Draw();
    p_upper->cd();

    hist_1->GetYaxis()->SetTitle(x_label);
    hist_1->GetYaxis()->SetLabelSize(0.08);
    hist_1->GetYaxis()->SetTitleSize(0.08);
    hist_1->GetYaxis()->SetTitleOffset(0.50);
    // X Axis
    hist_1->GetXaxis()->SetLabelOffset(999);
    hist_1->GetXaxis()->SetLabelSize(0);
    hist_1->SetTitle(plot_title);
    hist_1->SetLineWidth(1);
    // hist_1->GetXaxis()->SetRangeUser(1e-1, max_e);
    hist_1->Draw();
    hist_1->SetStats(0);
    gPad->SetGrid();
    gPad->SetLogx();
    // gPad->SetLogy();

    hist_2->SetLineWidth(1);
    hist_2->SetLineColor(2);
    // hist_2->GetXaxis()->SetRangeUser(1e-1, max_e);
    hist_2->Draw("][SAME");

    TLegend *sub_legend = new TLegend(legend_coord[0], legend_coord[1], legend_coord[2], legend_coord[3]);
    sub_legend->AddEntry(hist_1, legend_1, "l");
    sub_legend->AddEntry(hist_2, legend_2, "l");
    sub_legend->Draw();

    canv->cd(0);
    TPad *p_lower = new TPad(Form("p_upper_%d", plot_index), Form("p_upper_%d", plot_index), 0., 0., 1., 0.4);
    p_lower->SetFillColor(kWhite);
    p_lower->SetTopMargin(0.00001);
    p_lower->SetBottomMargin(0.3);
    p_lower->SetBorderMode(0);
    p_lower->Draw();
    p_lower->cd();

    subtraction_plot->SetTitle("");
    // X Axis
    subtraction_plot->GetXaxis()->SetTitle("Energy (in eV)");
    subtraction_plot->GetXaxis()->SetTitleSize(0.12);
    subtraction_plot->GetXaxis()->SetLabelOffset(0.01);
    subtraction_plot->GetXaxis()->SetTitleOffset(1.13); //decrease to move up
    subtraction_plot->GetXaxis()->SetLabelSize(0.12);
    // Y Axis
    subtraction_plot->GetYaxis()->SetTitle("Difference");
    subtraction_plot->GetYaxis()->SetTitleSize(0.12);
    subtraction_plot->GetYaxis()->SetTitleOffset(0.32);
    subtraction_plot->GetYaxis()->SetLabelSize(0.12);
    subtraction_plot->GetYaxis()->SetNdivisions(5);
    // TGaxis::SetExponentOffset(-0.06, -0.8, "y"); // X and Y offset for Y axis
    // subtraction_plot->GetXaxis()->SetRangeUser(1e-1, max_e);
    subtraction_plot->Draw();
    // gPad->SetGrid();
    gPad->SetLogx();
    zero_line->Draw("SAME");

    canv->Print(output_file_name);

    plot_index++;
    return;

}

std::vector<Double_t> calculate_avg_xsec(TH1D* xsec_hist, Int_t min_bin, Int_t max_bin, bool is_endf){

    Double_t sum_xsec = 0;
    Double_t sum_xsec_err = 0;

    for (Int_t i = min_bin; i < max_bin+1; i++)
    {
        sum_xsec += xsec_hist->GetBinContent(i);
        if (!is_endf)
        {
            sum_xsec_err += xsec_hist->GetBinError(i) * xsec_hist->GetBinError(i);   
        }
    }

    std::vector<Double_t> xsec_vec;
    xsec_vec.push_back(sum_xsec/(max_bin-min_bin+1));
    xsec_vec.push_back(std::sqrt(sum_xsec_err)/(max_bin-min_bin+1));

    return xsec_vec;
    
}

void plot_transmission_station_effects(){

    TH1D* trans_hist_al5_ptbc = Get_1D_hist("../rootFiles/trans_xsec_hists_al5_20bpd.root", "transmission_hist_e_PTBC");
    TH1D* trans_hist_al5_fimg = Get_1D_hist("../rootFiles/trans_xsec_hists_al5_20bpd.root", "transmission_hist_e_FIMG");

    TH1D* trans_hist_al5_ts_ptbc = Get_1D_hist("../rootFiles/trans_xsec_hists_al5_ts_20bpd.root", "transmission_hist_e_PTBC");
    TH1D* trans_hist_al5_ts_fimg = Get_1D_hist("../rootFiles/trans_xsec_hists_al5_ts_20bpd.root", "transmission_hist_e_FIMG");

    TH1D* subtraction_hist_ptbc = subtractHists(trans_hist_al5_ptbc, trans_hist_al5_ts_ptbc);
    TH1D* subtraction_hist_fimg = subtractHists(trans_hist_al5_fimg, trans_hist_al5_ts_fimg);

    trans_hist_al5_ptbc->GetYaxis()->SetRangeUser(-0.1, 1.1);
    trans_hist_al5_fimg->GetYaxis()->SetRangeUser(-0.1, 1.1);
    trans_hist_al5_ts_ptbc->GetYaxis()->SetRangeUser(-0.1, 1.1);
    trans_hist_al5_ts_fimg->GetYaxis()->SetRangeUser(-0.1, 1.1);

    trans_hist_al5_ptbc->GetXaxis()->SetRangeUser(1e-1, 1e8);
    trans_hist_al5_fimg->GetXaxis()->SetRangeUser(1e-1, 1e6);
    trans_hist_al5_ts_ptbc->GetXaxis()->SetRangeUser(1e-1, 1e8);
    trans_hist_al5_ts_fimg->GetXaxis()->SetRangeUser(1e-1, 1e6);

    Double_t trans_legend_coord[4] = {0.15,0.1,0.5,0.3};

    plot_subtraction_plots(trans_hist_al5_ptbc, "Transmission", "Al - Filter Station", trans_hist_al5_ts_ptbc, "Al - Transmission Station", subtraction_hist_ptbc, "Transmission - Al 5cm - Fission Chamber", "../plots/stability_plots/al5_filter_ts_ptbc.png", trans_legend_coord, 1e8);
    plot_subtraction_plots(trans_hist_al5_fimg, "Transmission", "Al - Filter Station", trans_hist_al5_ts_fimg, "Al - Transmission Station", subtraction_hist_fimg, "Transmission - Al 5cm - Micromegas", "../plots/stability_plots/al5_filter_ts_fimg.png", trans_legend_coord, 1e6);
}

void plot_al5_plots(){

    TH1D* trans_hist_al5_ptbc = Get_1D_hist("../rootFiles/trans_xsec_hists_al5_50bpd.root", "transmission_hist_e_PTBC");
    TH1D* trans_hist_al5_fimg = Get_1D_hist("../rootFiles/trans_xsec_hists_al5_50bpd.root", "transmission_hist_e_FIMG");
    TH1D* xsec_hist_al5_ptbc = Get_1D_hist("../rootFiles/trans_xsec_hists_al5_50bpd.root", "cross_section_hist_e_PTBC");
    TH1D* xsec_hist_al5_fimg = Get_1D_hist("../rootFiles/trans_xsec_hists_al5_50bpd.root", "cross_section_hist_e_FIMG");
    TH1D* endf_trans_hist = Get_1D_hist("../rootFiles/trans_xsec_hists_al5_50bpd.root", "endf_trans_hist");
    TH1D* endf_xsec_hist = Get_1D_hist("../rootFiles/trans_xsec_hists_al5_50bpd.root", "endf_xsec_hist");

    TH1D* trans_subtraction_hist_ptbc = subtractHists_theoryEval(trans_hist_al5_ptbc, endf_trans_hist);
    TH1D* trans_subtraction_hist_fimg = subtractHists_theoryEval(trans_hist_al5_fimg, endf_trans_hist);
    TH1D* xsec_subtraction_hist_ptbc = subtractHists_theoryEval(xsec_hist_al5_ptbc, endf_xsec_hist);
    TH1D* xsec_subtraction_hist_fimg = subtractHists_theoryEval(xsec_hist_al5_fimg, endf_xsec_hist);

    trans_hist_al5_fimg->GetXaxis()->SetRangeUser(1e-1, 1e6);
    trans_hist_al5_ptbc->GetXaxis()->SetRangeUser(1e-1, 1e8);
    xsec_hist_al5_fimg->GetXaxis()->SetRangeUser(1e-1, 1e6);
    xsec_hist_al5_ptbc->GetXaxis()->SetRangeUser(1e-1, 1e8);

    trans_subtraction_hist_ptbc->GetXaxis()->SetRangeUser(1e-1, 1e8);
    trans_subtraction_hist_fimg->GetXaxis()->SetRangeUser(1e-1, 1e6);
    xsec_subtraction_hist_ptbc->GetXaxis()->SetRangeUser(1e-1, 1e8);
    xsec_subtraction_hist_fimg->GetXaxis()->SetRangeUser(1e-1, 1e6);

    trans_hist_al5_ptbc->GetYaxis()->SetRangeUser(-0.1,1.05);
    trans_hist_al5_fimg->GetYaxis()->SetRangeUser(-0.1,1.05);
    xsec_hist_al5_ptbc->GetYaxis()->SetRangeUser(-1,24);
    xsec_hist_al5_fimg->GetYaxis()->SetRangeUser(-1,24);

    trans_subtraction_hist_ptbc->GetYaxis()->SetRangeUser(-0.24,0.24);
    trans_subtraction_hist_fimg->GetYaxis()->SetRangeUser(-0.24,0.24);
    xsec_subtraction_hist_ptbc->GetYaxis()->SetRangeUser(-1.7,1.7);
    xsec_subtraction_hist_fimg->GetYaxis()->SetRangeUser(-1.7,1.7);

    Double_t trans_legend_coord[4] = {0.2,0.1,0.35,0.3};
    Double_t xsec_legend_coord[4] = {0.2,0.55,0.35,0.75};

    plot_subtraction_plots(trans_hist_al5_fimg, "Transmission", "5 cm Al", endf_trans_hist, "ENDF-VIII", trans_subtraction_hist_fimg, "Transmission - Al 5cm - Micromegas", "../plots/results_plots/trans_al5_fimg_50bpd.png", trans_legend_coord, 1e6);
    plot_subtraction_plots(trans_hist_al5_ptbc, "Transmission", "5 cm Al", endf_trans_hist, "ENDF-VIII", trans_subtraction_hist_ptbc, "Transmission - Al 5cm - Fission Chamber", "../plots/results_plots/trans_al5_ptbc_50bpd.png", trans_legend_coord, 1e8);
    plot_subtraction_plots(xsec_hist_al5_fimg, "Cross Section", "5 cm Al", endf_xsec_hist, "ENDF-VIII", xsec_subtraction_hist_fimg, "Cross Section - Al 5cm - Micromegas", "../plots/results_plots/xsec_al5_fimg_50bpd.png", xsec_legend_coord, 1e6);
    plot_subtraction_plots(xsec_hist_al5_ptbc, "Cross Section", "5 cm Al", endf_xsec_hist, "ENDF-VIII", xsec_subtraction_hist_ptbc, "Cross Section - Al 5cm - Fission Chamber", "../plots/results_plots/xsec_al5_ptbc_50bpd.png", xsec_legend_coord, 1e8);

    std::vector<Double_t> avg_xsec_fimg = calculate_avg_xsec(xsec_hist_al5_fimg, 151, 250, false);
    std::vector<Double_t> avg_xsec_ptbc = calculate_avg_xsec(xsec_hist_al5_ptbc, 151, 250, false);
    std::vector<Double_t> avg_xsec_endf = calculate_avg_xsec(endf_xsec_hist, 151, 250, true);

    cout << "Avg xsec 10eV - 1keV fimg = " << avg_xsec_fimg.at(0) << " +- " << avg_xsec_fimg.at(1) << endl;
    cout << "Avg xsec 10eV - 1keV ptbc = " << avg_xsec_ptbc.at(0) << " +- " << avg_xsec_ptbc.at(1) << endl;
    cout << "Avg xsec 10eV - 1keV endf = " << avg_xsec_endf.at(0) << " +- " << avg_xsec_endf.at(1) << endl;
    
}

void plot_bi1_plots(){

    TH1D* trans_hist_bi1_ptbc = Get_1D_hist("../rootFiles/trans_xsec_hists_bi1_20bpd.root", "transmission_hist_e_PTBC");
    TH1D* trans_hist_bi1_fimg = Get_1D_hist("../rootFiles/trans_xsec_hists_bi1_20bpd.root", "transmission_hist_e_FIMG");
    TH1D* xsec_hist_bi1_ptbc = Get_1D_hist("../rootFiles/trans_xsec_hists_bi1_20bpd.root", "cross_section_hist_e_PTBC");
    TH1D* xsec_hist_bi1_fimg = Get_1D_hist("../rootFiles/trans_xsec_hists_bi1_20bpd.root", "cross_section_hist_e_FIMG");
    TH1D* endf_trans_hist = Get_1D_hist("../rootFiles/trans_xsec_hists_bi1_20bpd.root", "endf_trans_hist");
    TH1D* endf_xsec_hist = Get_1D_hist("../rootFiles/trans_xsec_hists_bi1_20bpd.root", "endf_xsec_hist");

    TH1D* trans_subtraction_hist_ptbc = subtractHists_theoryEval(trans_hist_bi1_ptbc, endf_trans_hist);
    TH1D* trans_subtraction_hist_fimg = subtractHists_theoryEval(trans_hist_bi1_fimg, endf_trans_hist);
    TH1D* xsec_subtraction_hist_ptbc = subtractHists_theoryEval(xsec_hist_bi1_ptbc, endf_xsec_hist);
    TH1D* xsec_subtraction_hist_fimg = subtractHists_theoryEval(xsec_hist_bi1_fimg, endf_xsec_hist);

    trans_hist_bi1_fimg->GetXaxis()->SetRangeUser(1e-1, 1e6);
    trans_hist_bi1_ptbc->GetXaxis()->SetRangeUser(1e-1, 1e8);
    xsec_hist_bi1_fimg->GetXaxis()->SetRangeUser(1e-1, 1e6);
    xsec_hist_bi1_ptbc->GetXaxis()->SetRangeUser(1e-1, 1e8);

    trans_subtraction_hist_ptbc->GetXaxis()->SetRangeUser(1e-1, 1e8);
    trans_subtraction_hist_fimg->GetXaxis()->SetRangeUser(1e-1, 1e6);
    xsec_subtraction_hist_ptbc->GetXaxis()->SetRangeUser(1e-1, 1e8);
    xsec_subtraction_hist_fimg->GetXaxis()->SetRangeUser(1e-1, 1e6);

    trans_hist_bi1_ptbc->GetYaxis()->SetRangeUser(-0.1,1.05);
    trans_hist_bi1_fimg->GetYaxis()->SetRangeUser(-0.1,1.05);
    xsec_hist_bi1_ptbc->GetYaxis()->SetRangeUser(-1,24);
    xsec_hist_bi1_fimg->GetYaxis()->SetRangeUser(-1,24);

    trans_subtraction_hist_ptbc->GetYaxis()->SetRangeUser(-0.24,0.24);
    trans_subtraction_hist_fimg->GetYaxis()->SetRangeUser(-0.24,0.24);
    xsec_subtraction_hist_ptbc->GetYaxis()->SetRangeUser(-5,5);
    xsec_subtraction_hist_fimg->GetYaxis()->SetRangeUser(-5,5);

    Double_t trans_legend_coord[4] = {0.15,0.1,0.3,0.3};
    Double_t xsec_legend_coord[4] = {0.15,0.55,0.3,0.75};

    plot_subtraction_plots(trans_hist_bi1_fimg, "Transmission", "1 cm Bi", endf_trans_hist, "ENDF-VIII", trans_subtraction_hist_fimg, "Transmission - Bi 1cm - Micromegas", "../plots/results_plots/trans_bi1_fimg_20bpd.png", trans_legend_coord, 1e6);
    plot_subtraction_plots(trans_hist_bi1_ptbc, "Transmission", "1 cm Bi", endf_trans_hist, "ENDF-VIII", trans_subtraction_hist_ptbc, "Transmission - Bi 1cm - Fission Chamber", "../plots/results_plots/trans_bi1_ptbc_20bpd.png", trans_legend_coord, 1e8);
    plot_subtraction_plots(xsec_hist_bi1_fimg, "Cross Section", "1 cm Bi", endf_xsec_hist, "ENDF-VIII", xsec_subtraction_hist_fimg, "Cross Section - Bi 1cm - Micromegas", "../plots/results_plots/xsec_bi1_fimg_20bpd.png", xsec_legend_coord, 1e6);
    plot_subtraction_plots(xsec_hist_bi1_ptbc, "Cross Section", "1 cm Bi", endf_xsec_hist, "ENDF-VIII", xsec_subtraction_hist_ptbc, "Cross Section - Bi 1cm - Fission Chamber", "../plots/results_plots/xsec_bi1_ptbc_20bpd.png", xsec_legend_coord, 1e8);
    
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

void plot_c1p2_plots(){

    TH1D* trans_hist_c1p2_ts_ptbc = Get_1D_hist("../rootFiles/trans_xsec_hists_c1p2_ts_20bpd.root", "transmission_hist_e_PTBC");
    TH1D* trans_hist_c1p2_ts_fimg = Get_1D_hist("../rootFiles/trans_xsec_hists_c1p2_ts_20bpd.root", "transmission_hist_e_FIMG");
    TH1D* xsec_hist_c1p2_ts_ptbc = Get_1D_hist("../rootFiles/trans_xsec_hists_c1p2_ts_20bpd.root", "cross_section_hist_e_PTBC");
    TH1D* xsec_hist_c1p2_ts_fimg = Get_1D_hist("../rootFiles/trans_xsec_hists_c1p2_ts_20bpd.root", "cross_section_hist_e_FIMG");
    TH1D* endf_trans_hist = Get_1D_hist("../rootFiles/trans_xsec_hists_c1p2_ts_20bpd.root", "endf_trans_hist");
    TH1D* endf_xsec_hist = Get_1D_hist("../rootFiles/trans_xsec_hists_c1p2_ts_20bpd.root", "endf_xsec_hist");

    TH1D* trans_subtraction_hist_ptbc = subtractHists_theoryEval(trans_hist_c1p2_ts_ptbc, endf_trans_hist);
    TH1D* trans_subtraction_hist_fimg = subtractHists_theoryEval(trans_hist_c1p2_ts_fimg, endf_trans_hist);
    TH1D* xsec_subtraction_hist_ptbc = subtractHists_theoryEval(xsec_hist_c1p2_ts_ptbc, endf_xsec_hist);
    TH1D* xsec_subtraction_hist_fimg = subtractHists_theoryEval(xsec_hist_c1p2_ts_fimg, endf_xsec_hist);

    trans_hist_c1p2_ts_fimg->GetXaxis()->SetRangeUser(1e-1, 1e6);
    trans_hist_c1p2_ts_ptbc->GetXaxis()->SetRangeUser(1e-1, 1e8);
    xsec_hist_c1p2_ts_fimg->GetXaxis()->SetRangeUser(1e-1, 1e6);
    xsec_hist_c1p2_ts_ptbc->GetXaxis()->SetRangeUser(1e-1, 1e8);

    trans_subtraction_hist_ptbc->GetXaxis()->SetRangeUser(1e-1, 1e8);
    trans_subtraction_hist_fimg->GetXaxis()->SetRangeUser(1e-1, 1e6);
    xsec_subtraction_hist_ptbc->GetXaxis()->SetRangeUser(1e-1, 1e8);
    xsec_subtraction_hist_fimg->GetXaxis()->SetRangeUser(1e-1, 1e6);

    trans_hist_c1p2_ts_ptbc->GetYaxis()->SetRangeUser(0.42,1.05);
    trans_hist_c1p2_ts_fimg->GetYaxis()->SetRangeUser(0.42,0.95);
    xsec_hist_c1p2_ts_ptbc->GetYaxis()->SetRangeUser(-0.6,7.1);
    xsec_hist_c1p2_ts_fimg->GetYaxis()->SetRangeUser(0.1,7.1);

    trans_subtraction_hist_ptbc->GetYaxis()->SetRangeUser(-0.19,0.19);
    trans_subtraction_hist_fimg->GetYaxis()->SetRangeUser(-0.19,0.19);
    xsec_subtraction_hist_ptbc->GetYaxis()->SetRangeUser(-1.9,1.9);
    xsec_subtraction_hist_fimg->GetYaxis()->SetRangeUser(-1.9,1.9);

    Double_t xsec_legend_coord[4] = {0.2,0.1,0.35,0.3};
    Double_t trans_legend_coord[4] = {0.2,0.55,0.35,0.75};

    plot_subtraction_plots(trans_hist_c1p2_ts_fimg, "Transmission", "1.2cm C", endf_trans_hist, "ENDF-VIII", trans_subtraction_hist_fimg, "Transmission - C 1.2cm - Micromegas", "../plots/results_plots/trans_c1p2_ts_fimg_20bpd.png", trans_legend_coord, 1e6);
    plot_subtraction_plots(trans_hist_c1p2_ts_ptbc, "Transmission", "1.2cm C", endf_trans_hist, "ENDF-VIII", trans_subtraction_hist_ptbc, "Transmission - C 1.2cm - Fission Chamber", "../plots/results_plots/trans_c1p2_ts_ptbc_20bpd.png", trans_legend_coord, 1e8);
    plot_subtraction_plots(xsec_hist_c1p2_ts_fimg, "Cross Section", "1.2cm C", endf_xsec_hist, "ENDF-VIII", xsec_subtraction_hist_fimg, "Cross Section - C 1.2cm - Micromegas", "../plots/results_plots/xsec_c1p2_ts_fimg_20bpd.png", xsec_legend_coord, 1e6);
    plot_subtraction_plots(xsec_hist_c1p2_ts_ptbc, "Cross Section", "1.2cm C", endf_xsec_hist, "ENDF-VIII", xsec_subtraction_hist_ptbc, "Cross Section - C 1.2cm - Fission Chamber", "../plots/results_plots/xsec_c1p2_ts_ptbc_20bpd.png", xsec_legend_coord, 1e8);
    
    // std::vector<Double_t> avg_xsec_endf = calculate_avg_xsec(endf_xsec_hist, 151, 250, true);

    std::vector<Double_t> avg_xsec_fimg = calculate_avg_xsec(xsec_hist_c1p2_ts_fimg, 61, 100, false);
    std::vector<Double_t> avg_xsec_ptbc = calculate_avg_xsec(xsec_hist_c1p2_ts_ptbc, 61, 100, false);
    std::vector<Double_t> avg_xsec_endf = calculate_avg_xsec(endf_xsec_hist, 61, 100, true);

    cout << "Avg xsec 10eV - 1keV fimg = " << avg_xsec_fimg.at(0) << " +- " << avg_xsec_fimg.at(1) << endl;
    cout << "Avg xsec 10eV - 1keV ptbc = " << avg_xsec_ptbc.at(0) << " +- " << avg_xsec_ptbc.at(1) << endl;
    cout << "Avg xsec 10eV - 1keV endf = " << avg_xsec_endf.at(0) << " +- " << avg_xsec_endf.at(1) << endl;
    
}

void plotsFromRootFile() {

    // plot_c1p2_plots();
    plot_al5_plots();
    // plot_bi1_plots();
    // plot_transmission_station_effects();

    // Int_t newBPD = 50;
    // Int_t oldBPD = 100;

    // retriveHistogramsChangeBPD("../rootFiles/crossSectionAna_cf_bottle.root", oldBPD, newBPD);
    // retriveHistogramsChangeBPD("../rootFiles/crossSectionAna_cf_bottle_rot.root", oldBPD, newBPD);
    // retriveHistogramsChangeBPD("../rootFiles/crossSectionAna_cf_bottle_rotBack.root", oldBPD, newBPD);

    // retriveHistogramsChangeBPD("../rootFiles/crossSectionAna_ar_bottle_full.root", oldBPD, newBPD);

    // TH1D* h_initial_rotBack = subtractHists(h[0], h[2]);
    // TH1D* h_initial_rot = subtractHists(h[0], h[1]);

    //Plotting
    // SetMArEXStyle();
    
    // gStyle->SetStatX(0.27);
    // gStyle->SetStatY(0.9);
    // gStyle->SetStatH(0.1);
    // gStyle->SetStatW(0.17);

    // TCanvas *c[2];
    // TLegend *l[2];

    // int i = 0;

    // c[i] = new TCanvas(Form("c%d", i)," ");
    // c[i]->cd();

    // l[i] = new TLegend(0.72,0.2,0.90,0.3);
    // l[i]->AddEntry(h[0],"Empty Bottle","l");
    
    // h[0]->GetXaxis()->SetTitle("Energy (in eV)");
    // h[0]->GetYaxis()->SetTitle("Transmission");
    // h[0]->SetTitle("Transmission Histogram - Empty Bottle runs - PTBC");
    // h[0]->SetLineWidth(2);
    // h[0]->Draw(); //"HISTE"
    // h[0]->SetStats(0);
    // // h[0]->SetMarkerStyle(6);
    // // h[0]->SetMarkerSize(2);
    // gPad->SetGrid();
    // gPad->SetLogx();

    // l[i]->AddEntry(h[1],"Empty Bottle - Rotated","l");
    // h[1]->SetLineWidth(2);
    // h[1]->SetLineColor(2);
    // h[1]->Draw("SAME");

    // l[i]->AddEntry(h[2],"Empty Bottle - Rotated Back","l");
    // h[2]->SetLineWidth(2);
    // h[2]->SetLineColor(3);
    // h[2]->Draw("SAME");

    // l[i]->Draw();

    // //////////////////////////////////

    // i++;

    // c[i] = new TCanvas(Form("c%d", i)," ");
    // c[i]->cd();

    // l[i] = new TLegend(0.72,0.2,0.90,0.3);
    // l[i]->AddEntry(h_initial_rotBack,"Initial - Rotated Back","l");
    // h_initial_rotBack->GetXaxis()->SetTitle("Energy (in eV)");
    // h_initial_rotBack->GetYaxis()->SetTitle("Difference in Transmission");
    // h_initial_rotBack->SetTitle("Difference in Transmission Hist - Empty Bottle runs - PTBC");
    // h_initial_rotBack->SetLineWidth(2);
    // h_initial_rotBack->Draw(); //"HISTE"
    // h_initial_rotBack->SetStats(0);
    // // h_initial_rotBack->SetMarkerStyle(6);
    // // h_initial_rotBack->SetMarkerSize(2);
    // gPad->SetGrid();
    // gPad->SetLogx();

    // l[i]->AddEntry(h_initial_rot,"Initial - Rotated","l");
    // h_initial_rot->SetLineWidth(2);
    // h_initial_rot->SetLineColor(2);
    // h_initial_rot->Draw("SAME");

    // TLine* l_zero = new TLine(1e-2,0,1e9,0);
    // l_zero->SetLineWidth(2);
    // l_zero->SetLineColor(1);
    // l_zero->Draw();

    // l[i]->Draw();

}