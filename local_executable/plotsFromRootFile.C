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

TH1D* h[2];
TH1D* h_old[2];
Int_t histCounter = 0;

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

    // trans_hist_fOut = (TH1D*)hist_file->Get("trans_hist_fOut");
    // trans_hist_fIn = (TH1D*)hist_file->Get("trans_hist_fIn");
    // trans_hist_fOut_endf = (TH1D*)hist_file->Get("trans_hist_fOut_endf");
    // trans_hist_fIn_endf = (TH1D*)hist_file->Get("trans_hist_fIn_endf");
    h[histCounter] = (TH1D*)hist_file->Get("transmission_hist_e");
    histCounter++;
    // cross_section_hist_e = (TH1D*)hist_file->Get("cross_section_hist_e");

    // hist_file->Close();
}

void retriveHistogramsChangeBPD(const char *fname, const char *hist_name, Int_t bpd_old, Int_t bpd_new){

    // if(bpd_old < bpd_new) {

    // }

    // TKey *key;
    cout << "Opening File " << fname << endl;
    TFile *hist_file = TFile::Open(fname, "READ");
    if (!hist_file || hist_file->IsZombie()) {
        cout << "Unable to open " << fname << " for reading..." << endl;
        return;
    }

    h_old[histCounter] = (TH1D*)hist_file->Get(hist_name);
    Int_t num_bins_old = h_old[histCounter]->GetNbinsX();
    cout << "Number of old bins = " << num_bins_old << endl;
    Int_t max_sum_bin_count = (bpd_old/bpd_new); // Number of old bins that need to be summed to make a new bin
    Int_t num_bins_new = num_bins_old / max_sum_bin_count;
    cout << "Number of new bins = " << num_bins_new << endl;

    Double_t bin_edges_new[num_bins_new + 1];
    Double_t bin_content[num_bins_new];
    Double_t bin_error[num_bins_new];
    Int_t sum_bin_count = 0;
    Int_t new_bin_counter = 0;
    for (int i = 1; i < num_bins_old + 1; i++)
    {
        if (i == 1)
        {
            bin_edges_new[new_bin_counter] = h_old[histCounter]->GetXaxis()->GetBinLowEdge(i);
            bin_content[new_bin_counter] = 0;
            bin_error[new_bin_counter] = 0;
        }
        bin_content[new_bin_counter] += h_old[histCounter]->GetBinContent(i);
        bin_error[new_bin_counter] += h_old[histCounter]->GetBinError(i) * h_old[histCounter]->GetBinError(i);
        sum_bin_count++;

        if (i == num_bins_old)
        {
            // sum_bin_count = 0;
            bin_error[new_bin_counter] = sqrt(bin_error[new_bin_counter])/max_sum_bin_count;
            bin_content[new_bin_counter] = bin_content[new_bin_counter]/max_sum_bin_count;
            bin_edges_new[new_bin_counter+1] = h_old[histCounter]->GetXaxis()->GetBinUpEdge(i);
        }
        else if (sum_bin_count == max_sum_bin_count)
        {
            sum_bin_count = 0;
            bin_error[new_bin_counter] = sqrt(bin_error[new_bin_counter])/max_sum_bin_count;
            bin_content[new_bin_counter] = bin_content[new_bin_counter]/max_sum_bin_count;
            new_bin_counter++;
            bin_edges_new[new_bin_counter] = h_old[histCounter]->GetXaxis()->GetBinUpEdge(i);
            bin_content[new_bin_counter] = 0;
            bin_error[new_bin_counter] = 0;
        } 
    }
    
    h[histCounter] = new TH1D(Form("h_%i", histCounter),"Transmission Hist",num_bins_new,bin_edges_new);
    for(int i = 0; i < num_bins_new; i++){
        h[histCounter]->SetBinContent(i+1, bin_content[i]);
        h[histCounter]->SetBinError(i+1, bin_error[i]);
    }
    histCounter++;
}

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

void plotsFromRootFile() {

    Int_t newBPD = 50;
    Int_t oldBPD = 100;

    // retriveHistogramsChangeBPD("../rootFiles/crossSectionAna_cf_bottle.root", oldBPD, newBPD);
    // retriveHistogramsChangeBPD("../rootFiles/crossSectionAna_cf_bottle_rot.root", oldBPD, newBPD);
    // retriveHistogramsChangeBPD("../rootFiles/crossSectionAna_cf_bottle_rotBack.root", oldBPD, newBPD);

    retriveHistogramsChangeBPD("../rootFiles/crossSectionAna_ar_bottle_full.root", oldBPD, newBPD);

    // TH1D* h_initial_rotBack = subtractHists(h[0], h[2]);
    // TH1D* h_initial_rot = subtractHists(h[0], h[1]);

    //Plotting
    SetMArEXStyle();
    
    gStyle->SetStatX(0.27);
    gStyle->SetStatY(0.9);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.17);

    TCanvas *c[2];
    TLegend *l[2];

    int i = 0;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();

    l[i] = new TLegend(0.72,0.2,0.90,0.3);
    l[i]->AddEntry(h[0],"Empty Bottle","l");
    
    h[0]->GetXaxis()->SetTitle("Energy (in eV)");
    h[0]->GetYaxis()->SetTitle("Transmission");
    h[0]->SetTitle("Transmission Histogram - Empty Bottle runs - PTBC");
    h[0]->SetLineWidth(2);
    h[0]->Draw(); //"HISTE"
    h[0]->SetStats(0);
    // h[0]->SetMarkerStyle(6);
    // h[0]->SetMarkerSize(2);
    gPad->SetGrid();
    gPad->SetLogx();

    l[i]->AddEntry(h[1],"Empty Bottle - Rotated","l");
    h[1]->SetLineWidth(2);
    h[1]->SetLineColor(2);
    h[1]->Draw("SAME");

    l[i]->AddEntry(h[2],"Empty Bottle - Rotated Back","l");
    h[2]->SetLineWidth(2);
    h[2]->SetLineColor(3);
    h[2]->Draw("SAME");

    l[i]->Draw();

    //////////////////////////////////

    i++;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();

    l[i] = new TLegend(0.72,0.2,0.90,0.3);
    l[i]->AddEntry(h_initial_rotBack,"Initial - Rotated Back","l");
    h_initial_rotBack->GetXaxis()->SetTitle("Energy (in eV)");
    h_initial_rotBack->GetYaxis()->SetTitle("Difference in Transmission");
    h_initial_rotBack->SetTitle("Difference in Transmission Hist - Empty Bottle runs - PTBC");
    h_initial_rotBack->SetLineWidth(2);
    h_initial_rotBack->Draw(); //"HISTE"
    h_initial_rotBack->SetStats(0);
    // h_initial_rotBack->SetMarkerStyle(6);
    // h_initial_rotBack->SetMarkerSize(2);
    gPad->SetGrid();
    gPad->SetLogx();

    l[i]->AddEntry(h_initial_rot,"Initial - Rotated","l");
    h_initial_rot->SetLineWidth(2);
    h_initial_rot->SetLineColor(2);
    h_initial_rot->Draw("SAME");

    TLine* l_zero = new TLine(1e-2,0,1e9,0);
    l_zero->SetLineWidth(2);
    l_zero->SetLineColor(1);
    l_zero->Draw();

    l[i]->Draw();

}