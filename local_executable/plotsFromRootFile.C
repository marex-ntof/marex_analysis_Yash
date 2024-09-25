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

void plot_subtraction_plots(TH1D* hist_1, TH1D* hist_2, TH1D* subtraction_plot, const char* plot_title, const char* output_file_name, Int_t max_e){

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

    gStyle->SetCanvasDefW(800); //600
    gStyle->SetCanvasDefH(400); //500
    gStyle->SetPadRightMargin(0.05);

    TCanvas *canv = new TCanvas(Form("c%d", plot_index)," ");
    canv->cd();
    canv->Draw();

    TPad *p_upper = new TPad(Form("p_upper_%d", plot_index), Form("p_upper_%d", plot_index), 0., 0.35, 1., 1.);
    p_upper->SetFillColor(kWhite);
    p_upper->SetBottomMargin(0.00001);
    p_upper->SetBorderMode(0);
    p_upper->Draw();
    p_upper->cd();

    hist_1->GetYaxis()->SetTitle("Transmission");
    hist_1->GetYaxis()->SetLabelSize(0.05);
    hist_1->GetYaxis()->SetTitleSize(0.06);
    hist_1->GetYaxis()->SetTitleOffset(0.65);
    // X Axis
    hist_1->GetXaxis()->SetLabelOffset(999);
    hist_1->GetXaxis()->SetLabelSize(0);
    hist_1->SetTitle(plot_title);
    hist_1->SetLineWidth(2);
    hist_1->GetXaxis()->SetRangeUser(1e-1, max_e);
    hist_1->Draw();
    hist_1->SetStats(0);
    gPad->SetGrid();
    gPad->SetLogx();
    // gPad->SetLogy();

    hist_2->SetLineWidth(2);
    hist_2->SetLineColor(2);
    hist_2->GetXaxis()->SetRangeUser(1e-1, max_e);
    hist_2->Draw("SAME");

    TLegend *sub_legend = new TLegend(0.15,0.1,0.4,0.3);
    sub_legend->AddEntry(hist_1, "Al - Filter Station", "l");
    sub_legend->AddEntry(hist_2, "Al - Transmission Station", "l");
    sub_legend->Draw();

    canv->cd(0);
    TPad *p_lower = new TPad(Form("p_upper_%d", plot_index), Form("p_upper_%d", plot_index), 0., 0., 1., 0.35);
    p_lower->SetFillColor(kWhite);
    p_lower->SetTopMargin(0.00001);
    p_lower->SetBottomMargin(0.2);
    p_lower->SetBorderMode(0);
    p_lower->Draw();
    p_lower->cd();

    subtraction_plot->SetTitle("");
    // X Axis
    subtraction_plot->GetXaxis()->SetTitle("Energy (in eV)");
    subtraction_plot->GetXaxis()->SetTitleSize(0.1);
    subtraction_plot->GetXaxis()->SetLabelOffset(0.01);
    subtraction_plot->GetXaxis()->SetTitleOffset(0.85); //decrease to move up
    subtraction_plot->GetXaxis()->SetLabelSize(0.07);
    // Y Axis
    subtraction_plot->GetYaxis()->SetTitle("Difference");
    subtraction_plot->GetYaxis()->SetTitleSize(0.1);
    subtraction_plot->GetYaxis()->SetTitleOffset(0.38);
    subtraction_plot->GetYaxis()->SetLabelSize(0.07);
    // TGaxis::SetExponentOffset(-0.06, -0.8, "y"); // X and Y offset for Y axis
    subtraction_plot->GetXaxis()->SetRangeUser(1e-1, max_e);
    subtraction_plot->Draw();
    gPad->SetGrid();
    gPad->SetLogx();
    canv->Print( Form("../plots/stability_plots/%s.png", output_file_name) );

    zero_line->Draw("SAME");

    plot_index++;
    return;

}

void plotsFromRootFile() {

    // Int_t newBPD = 50;
    // Int_t oldBPD = 100;

    // retriveHistogramsChangeBPD("../rootFiles/crossSectionAna_cf_bottle.root", oldBPD, newBPD);
    // retriveHistogramsChangeBPD("../rootFiles/crossSectionAna_cf_bottle_rot.root", oldBPD, newBPD);
    // retriveHistogramsChangeBPD("../rootFiles/crossSectionAna_cf_bottle_rotBack.root", oldBPD, newBPD);

    // retriveHistogramsChangeBPD("../rootFiles/crossSectionAna_ar_bottle_full.root", oldBPD, newBPD);

    // TH1D* h_initial_rotBack = subtractHists(h[0], h[2]);
    // TH1D* h_initial_rot = subtractHists(h[0], h[1]);

    TH1D* trans_hist_al5_ptbc = Get_1D_hist("../rootFiles/trans_xsec_hists_al5_50bpd.root", "transmission_hist_e_PTBC");
    TH1D* trans_hist_al5_fimg = Get_1D_hist("../rootFiles/trans_xsec_hists_al5_50bpd.root", "transmission_hist_e_FIMG");

    TH1D* trans_hist_al5_ts_ptbc = Get_1D_hist("../rootFiles/trans_xsec_hists_al5_ts_50bpd.root", "transmission_hist_e_PTBC");
    TH1D* trans_hist_al5_ts_fimg = Get_1D_hist("../rootFiles/trans_xsec_hists_al5_ts_50bpd.root", "transmission_hist_e_FIMG");

    TH1D* subtraction_hist_ptbc = subtractHists(trans_hist_al5_ptbc, trans_hist_al5_ts_ptbc);
    TH1D* subtraction_hist_fimg = subtractHists(trans_hist_al5_fimg, trans_hist_al5_ts_fimg);

    plot_subtraction_plots(trans_hist_al5_ptbc, trans_hist_al5_ts_ptbc, subtraction_hist_ptbc, "Transmission - Al 5cm - Fission Chamber", "al5_filter_ts_ptbc", 1e8);
    plot_subtraction_plots(trans_hist_al5_fimg, trans_hist_al5_ts_fimg, subtraction_hist_fimg, "Transmission - Al 5cm - Micromegas", "al5_filter_ts_fimg", 1e6);

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