/**
 * @file fimgCutsAna.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-11-29
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

TH1D* trans_loose_cut_det_1 = 0;
TH1D* trans_mid_cut_det_1 = 0;
TH1D* trans_tight_cut_det_1 = 0;
TH1D* trans_loose_cut_det_2 = 0;
TH1D* trans_mid_cut_det_2 = 0;
TH1D* trans_tight_cut_det_2 = 0;

TH1D* endf_trans_hist = 0;

Double_t n_C_1p2cm = 0.105;// (1.2 /*cm*/) * (2.267 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (12.011 /*g/mole*/);

void endf(Double_t n, Double_t energy_bin_edges[]){
    //Extracting ENDF Cross Section and transmission

    std::ifstream inputFile("../evalData/C_tot_xsec.txt"); // endf_file_name.c_str()

    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the file.\n";
    }

    //Extracting the data from the text file
    Double_t val1, val2; //val1 -> Energy (eV); val2 -> Cross Section (barns)
    std::string line;
    Double_t line_count = 0;
    std::vector<Double_t> energy;
    std::vector<Double_t> xsec;

    while (std::getline(inputFile, line)) {

        line_count++;

        if (line_count == 1)
        {
            continue;
        }

        if (line_count == 1430) //Bi - 27019; Al - 9745
        {
            break;
        }

        std::istringstream iss(line);

        if (iss >> val1 >> val2) {
            energy.push_back(val1);
            xsec.push_back(val2);
        } else {
            std::cerr << "Invalid data format.\n";
        }
    }

    //Filling the histograms
    Double_t trans_sum = 0;
    Int_t sum_counter = 0;
    Int_t bin_counter = 1;

    for (int i = 0; i < energy.size(); i++)
    {        
        if (energy[i] > energy_bin_edges[bin_counter])
        {
            if (sum_counter == 0)
            {
                endf_trans_hist->SetBinContent(bin_counter, 0);
            } else {
                endf_trans_hist->SetBinContent(bin_counter, trans_sum/sum_counter);
            }
            
            trans_sum = std::exp(- n * xsec[i]);
            sum_counter = 1;
            bin_counter++;
            i--;
            continue;
        }

        if (energy[i] < energy_bin_edges[bin_counter])
        {
            trans_sum += std::exp(- n * xsec[i]);
            sum_counter++;
        }
    }
}

TH1D* retriveHistogramsChangeBPD(const char *fname, const char *hist_name, Int_t bpd_old, Int_t bpd_new){

    // if(bpd_old < bpd_new) {

    // }

    cout << "Opening File " << fname << endl;
    TFile *hist_file = TFile::Open(fname, "READ");
    // if (!hist_file || hist_file->IsZombie()) {
    //     cout << "Unable to open " << fname << " for reading..." << endl;
    //     return;
    // }

    // for (auto&& keyAsObj : *hist_file->GetListOfKeys()){
    //     auto key = (TKey*) keyAsObj;
    //     cout << key->GetName() << " " << key->GetClassName() << endl;
    // }

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
        bin_error[new_bin_counter] += hist_old->GetBinError(i) * hist_old->GetBinError(i);
        sum_bin_count++;

        if (i == num_bins_old)
        {
            // sum_bin_count = 0;
            bin_error[new_bin_counter] = sqrt(bin_error[new_bin_counter])/max_sum_bin_count;
            bin_content[new_bin_counter] = bin_content[new_bin_counter]/max_sum_bin_count;
            bin_edges_new[new_bin_counter+1] = hist_old->GetXaxis()->GetBinUpEdge(i);
        }
        else if (sum_bin_count == max_sum_bin_count)
        {
            sum_bin_count = 0;
            bin_error[new_bin_counter] = sqrt(bin_error[new_bin_counter])/max_sum_bin_count;
            bin_content[new_bin_counter] = bin_content[new_bin_counter]/max_sum_bin_count;
            new_bin_counter++;
            bin_edges_new[new_bin_counter] = hist_old->GetXaxis()->GetBinUpEdge(i);
            bin_content[new_bin_counter] = 0;
            bin_error[new_bin_counter] = 0;
        } 
    }

    auto hist_new = new TH1D(hist_name,hist_old->GetTitle(),num_bins_new,bin_edges_new);
    for(int i = 0; i < num_bins_new; i++){
        hist_new->SetBinContent(i+1, bin_content[i]);
        hist_new->SetBinError(i+1, bin_error[i]);
    }

    return hist_new;
}

void fimgCutsAna() {

    Int_t newBPD = 10;
    Int_t oldBPD = 1000;

    trans_loose_cut_det_1 = retriveHistogramsChangeBPD("../rootFiles/cutoffAnalysis_FIMG_c1p2_ts.root", "trans_loose_cut_det_1", oldBPD, newBPD);
    trans_mid_cut_det_1 = retriveHistogramsChangeBPD("../rootFiles/cutoffAnalysis_FIMG_c1p2_ts.root", "trans_mid_cut_det_1", oldBPD, newBPD);
    trans_tight_cut_det_1 = retriveHistogramsChangeBPD("../rootFiles/cutoffAnalysis_FIMG_c1p2_ts.root", "trans_tight_cut_det_1", oldBPD, newBPD);
    trans_loose_cut_det_2 = retriveHistogramsChangeBPD("../rootFiles/cutoffAnalysis_FIMG_c1p2_ts.root", "trans_loose_cut_det_2", oldBPD, newBPD);
    trans_mid_cut_det_2 = retriveHistogramsChangeBPD("../rootFiles/cutoffAnalysis_FIMG_c1p2_ts.root", "trans_mid_cut_det_2", oldBPD, newBPD);
    trans_tight_cut_det_2 = retriveHistogramsChangeBPD("../rootFiles/cutoffAnalysis_FIMG_c1p2_ts.root", "trans_tight_cut_det_2", oldBPD, newBPD);

    //Getting energy bin edges
    Int_t num_bins_e = trans_loose_cut_det_1->GetNbinsX();
    Double_t bin_edges_e[num_bins_e + 1];
    for (int i = 0; i < num_bins_e; i++)
    {
        bin_edges_e[i] = trans_loose_cut_det_1->GetXaxis()->GetBinLowEdge(i+1);

        if (i == num_bins_e - 1)
        {
            bin_edges_e[i+1] = trans_loose_cut_det_1->GetXaxis()->GetBinUpEdge(i+1);
        }
    }

    //ENDF Hists
    endf_trans_hist = new TH1D("endf_trans_hist","ENDF Transmission Hist",num_bins_e,bin_edges_e);
    //Filling ENDF hist
    endf(n_C_1p2cm, bin_edges_e);
    
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

    l[i] = new TLegend(0.77,0.2,0.86,0.3); //0.68,0.7,0.86,0.8       ;         0.72,0.8,0.90,0.9

    l[i]->AddEntry(trans_loose_cut_det_1,"Loose Cut (A = 800)","l");    
    trans_loose_cut_det_1->GetXaxis()->SetTitle("Energy (in eV)");
    trans_loose_cut_det_1->GetYaxis()->SetTitle("Transmission");
    trans_loose_cut_det_1->SetTitle("Transmission Histogram - Det 1 - C (1.2 cm)");
    trans_loose_cut_det_1->SetLineWidth(2);
    trans_loose_cut_det_1->Draw(); //"HISTE"
    trans_loose_cut_det_1->SetStats(0);
    // trans_loose_cut_det_1->SetMarkerStyle(6);
    // trans_loose_cut_det_1->SetMarkerSize(0.5);
    // gPad->SetGrid();
    gPad->SetLogx();
    // gStyle->SetPalette(57);

    l[i]->AddEntry(trans_mid_cut_det_1,"Mid Cut (A = 1800)","l");
    trans_mid_cut_det_1->SetLineColor(2);
    trans_mid_cut_det_1->SetLineWidth(2);
    trans_mid_cut_det_1->Draw("SAME");

    l[i]->AddEntry(trans_tight_cut_det_1,"Tight Cut (A = 2800)","l");
    trans_tight_cut_det_1->SetLineColor(3);
    trans_tight_cut_det_1->SetLineWidth(2);
    trans_tight_cut_det_1->Draw("SAME");
    
    l[i]->AddEntry(endf_trans_hist,"ENDF","l");
    endf_trans_hist->SetLineColor(1);
    endf_trans_hist->SetLineWidth(2);
    // endf_trans_hist->GetXaxis()->SetRange(1e-2,2e7);
    endf_trans_hist->Draw("SAME");    

    l[i]->SetMargin(0.4);
    l[i]->Draw();
    // c[i]->Print(Form("../plots/h_trans_e_%s_%iBPD.png", filter_name.c_str(), bins_per_decade));

    i++;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();

    l[i] = new TLegend(0.77,0.2,0.86,0.3); //0.68,0.7,0.86,0.8       ;         0.72,0.8,0.90,0.9

    l[i]->AddEntry(trans_loose_cut_det_2,"Loose Cut (A = 1000)","l");    
    trans_loose_cut_det_2->GetXaxis()->SetTitle("Energy (in eV)");
    trans_loose_cut_det_2->GetYaxis()->SetTitle("Transmission");
    trans_loose_cut_det_2->SetTitle("Transmission Histogram - Det 2 - C (1.2 cm)");
    trans_loose_cut_det_2->SetLineWidth(2);
    trans_loose_cut_det_2->Draw(); //"HISTE"
    trans_loose_cut_det_2->SetStats(0);
    // trans_loose_cut_det_2->SetMarkerStyle(6);
    // trans_loose_cut_det_2->SetMarkerSize(0.5);
    // gPad->SetGrid();
    gPad->SetLogx();
    // gStyle->SetPalette(57);

    l[i]->AddEntry(trans_mid_cut_det_2,"Mid Cut (A = 2000)","l");
    trans_mid_cut_det_2->SetLineColor(2);
    trans_mid_cut_det_2->SetLineWidth(2);
    trans_mid_cut_det_2->Draw("SAME");

    l[i]->AddEntry(trans_tight_cut_det_2,"Tight Cut (A = 3000)","l");
    trans_tight_cut_det_2->SetLineColor(3);
    trans_tight_cut_det_2->SetLineWidth(2);
    trans_tight_cut_det_2->Draw("SAME");
    
    l[i]->AddEntry(endf_trans_hist,"ENDF","l");
    endf_trans_hist->SetLineColor(1);
    endf_trans_hist->SetLineWidth(2);
    // endf_trans_hist->GetXaxis()->SetRange(1e-2,2e7);
    endf_trans_hist->Draw("SAME");    

    l[i]->SetMargin(0.4);
    l[i]->Draw();
    // c[i]->Print(Form("../plots/h_trans_e_%s_%iBPD.png", filter_name.c_str(), bins_per_decade));

}