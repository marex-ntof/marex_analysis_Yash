/**
 * @file data_stability_plots.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2024-02-26
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

Int_t histCounter = 1;

//////////// PTBC - counts plots
TH1D* norm_counts_empty_sep19_PTBC = 0;
TH1D* norm_counts_empty_oct01_PTBC = 0;

TH1D* norm_counts_Bi_sep17_PTBC = 0;
TH1D* norm_counts_Bi_sep23_PTBC = 0;

TH1D* norm_counts_Al5_sep25_PTBC = 0;
TH1D* norm_counts_Al5_oct02_PTBC = 0;

TH1D* norm_counts_emptyTank_oct12_PTBC = 0;
TH1D* norm_counts_emptyTank_oct16_PTBC = 0;

TH1D* norm_counts_ArgonFull_oct18_PTBC = 0;
TH1D* norm_counts_ArgonFull_oct22_PTBC = 0;

//////////// FIMG - counts plots
TH1D* norm_counts_empty_sep19_FIMG = 0;
TH1D* norm_counts_empty_oct01_FIMG = 0;

TH1D* norm_counts_Bi_sep17_FIMG = 0;
TH1D* norm_counts_Bi_sep23_FIMG = 0;

TH1D* norm_counts_Al5_sep25_FIMG = 0;
TH1D* norm_counts_Al5_oct02_FIMG = 0;

TH1D* norm_counts_emptyTank_oct12_FIMG = 0;
TH1D* norm_counts_emptyTank_oct16_FIMG = 0;

TH1D* norm_counts_ArgonFull_oct18_FIMG = 0;
TH1D* norm_counts_ArgonFull_oct22_FIMG = 0;

////////////                     
TCanvas *c[10];
TLegend *l[10];
TPad *p[10][2];
TH1D *residual_plot[10];
Int_t plot_index = 0;

void plot_norm_counts_plots(TH1D* norm_count_plot_1, const char* date_1, TH1D* norm_count_plot_2, const char* date_2, const char* plot_title, const char* output_file_name){

    c[plot_index] = new TCanvas(Form("c%d", plot_index)," ");
    c[plot_index]->cd();
    c[plot_index]->Draw();

    p[plot_index][0] = new TPad(Form("p_0_%d", plot_index), Form("p_0_%d", plot_index), 0., 0.35, 1., 1.);
    p[plot_index][0]->SetFillColor(kWhite);
    p[plot_index][0]->SetBottomMargin(0.00001);
    p[plot_index][0]->SetBorderMode(0);
    p[plot_index][0]->Draw();
    p[plot_index][0]->cd();
    // norm_count_plot_1->GetXaxis()->SetTitle("Energy (in eV)");
    // Y Axis
    norm_count_plot_1->GetYaxis()->SetTitle("Normalized Counts");
    norm_count_plot_1->GetYaxis()->SetLabelSize(0.05);
    norm_count_plot_1->GetYaxis()->SetTitleSize(0.06);
    norm_count_plot_1->GetYaxis()->SetTitleOffset(0.65);
    // X Axis
    norm_count_plot_1->GetXaxis()->SetLabelOffset(999);
    norm_count_plot_1->GetXaxis()->SetLabelSize(0);
    norm_count_plot_1->SetTitle(plot_title);
    norm_count_plot_1->SetLineWidth(2);
    norm_count_plot_1->Draw();
    norm_count_plot_1->SetStats(0);
    gPad->SetGrid();
    gPad->SetLogx();
    gPad->SetLogy();

    norm_count_plot_2->SetLineWidth(2);
    norm_count_plot_2->SetLineColor(2);
    norm_count_plot_2->Draw("SAME");

    l[plot_index] = new TLegend(0.75,0.7,0.85,0.8);
    l[plot_index]->AddEntry(norm_count_plot_1, date_1, "l");
    l[plot_index]->AddEntry(norm_count_plot_2, date_2, "l");
    l[plot_index]->Draw();

    c[plot_index]->cd(0);
    p[plot_index][1] = new TPad(Form("p_1_%d", plot_index), Form("p_1_%d", plot_index), 0., 0., 1., 0.35);
    p[plot_index][1]->SetFillColor(kWhite);
    p[plot_index][1]->SetTopMargin(0.00001);
    p[plot_index][1]->SetBottomMargin(0.2);
    p[plot_index][1]->SetBorderMode(0);
    p[plot_index][1]->Draw();
    p[plot_index][1]->cd();

    //Computing the residual plots
    // TH1D* subtraction_plot = (TH1D*)(norm_count_plot_1->Clone(Form("residual_plot_%d", plot_index)));
    // subtraction_plot->Add(norm_count_plot_2, -1);
    TH1D* addition_plot = (TH1D*)(norm_count_plot_1->Clone(Form("residual_plot_%d", plot_index)));
    addition_plot->Add(norm_count_plot_2, 1);
    addition_plot->Scale(0.5);
    residual_plot[plot_index] = (TH1D*)(norm_count_plot_1->Clone(Form("residual_plot_%d", plot_index))); //subtraction plot
    residual_plot[plot_index]->Add(norm_count_plot_2, -1);
    residual_plot[plot_index]->Divide(addition_plot);
    residual_plot[plot_index]->SetTitle("");

    // X Axis
    residual_plot[plot_index]->GetXaxis()->SetTitle("Energy (in eV)");
    residual_plot[plot_index]->GetXaxis()->SetTitleSize(0.1);
    residual_plot[plot_index]->GetXaxis()->SetLabelOffset(0.01);
    residual_plot[plot_index]->GetXaxis()->SetTitleOffset(0.85); //decrease to move up
    residual_plot[plot_index]->GetXaxis()->SetLabelSize(0.07);
    // Y Axis
    residual_plot[plot_index]->GetYaxis()->SetTitle("Fractional Residuals");
    residual_plot[plot_index]->GetYaxis()->SetTitleSize(0.1);
    residual_plot[plot_index]->GetYaxis()->SetTitleOffset(0.38);
    residual_plot[plot_index]->GetYaxis()->SetLabelSize(0.07);
    TGaxis::SetExponentOffset(-0.06, -0.8, "y"); // X and Y offset for Y axis
    residual_plot[plot_index]->Draw();
    gPad->SetGrid();
    gPad->SetLogx();
    c[plot_index]->Print( Form("../plots/stability_plots/%s.png", output_file_name) );

    plot_index++;
    return;
}

TH1D* retriveHistograms(const char *file_name, const char *hist_name){
    
    TFile* hist_file = TFile::Open(file_name, "READ");
    TH1D* hist_new = (TH1D*)hist_file->Get(hist_name);
    // hist_file->Close();

    return hist_new;
}

TH1D* retriveHistogramsChangeBPD(const char *file_name, const char *hist_name, Int_t bpd_old, Int_t bpd_new){
    
    TFile* hist_file = TFile::Open(file_name, "READ");
    TH1D* hist_old = (TH1D*)hist_file->Get(hist_name);
    
    Int_t num_bins_old = hist_old->GetNbinsX();
    Int_t max_sum_bin_count = (bpd_old/bpd_new); // Number of old bins that need to be summed to make a new bin
    Int_t num_bins_new = num_bins_old / max_sum_bin_count;

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
            bin_error[new_bin_counter] = sqrt(bin_error[new_bin_counter]);
            bin_edges_new[new_bin_counter+1] = hist_old->GetXaxis()->GetBinUpEdge(i);
        }
        else if (sum_bin_count == max_sum_bin_count)
        {
            sum_bin_count = 0;
            bin_error[new_bin_counter] = sqrt(bin_error[new_bin_counter]);
            new_bin_counter++;
            bin_edges_new[new_bin_counter] = hist_old->GetXaxis()->GetBinUpEdge(i);
            bin_content[new_bin_counter] = 0;
            bin_error[new_bin_counter] = 0;
        } 
    }

    TH1D* hist_new = new TH1D(Form("h_%i", histCounter),"Normalized Counts",num_bins_new,bin_edges_new);
    for(int i = 0; i < num_bins_new; i++){
        hist_new->SetBinContent(i+1, bin_content[i]);
        hist_new->SetBinError(i+1, bin_error[i]);
    }
    histCounter++;

    // excluding the first and the last bin
    Int_t seaching_for_first_last_bin = 1; // 1 - first bin; 2 - last bin
    for (Int_t i = 1; i <= num_bins_new; i++)
    {
        Double_t new_bin_content = hist_new->GetBinContent(i);

        // Searching for the first bin and setting it to zero
        if (seaching_for_first_last_bin == 1)
        {
            if (new_bin_content != 0)
            {
                hist_new->SetBinContent(i, 0);
                seaching_for_first_last_bin = 2;
                continue;
            }
        }
        
        // Searching for the last bin and setting it to zero
        if (seaching_for_first_last_bin == 2)
        {
            if (new_bin_content == 0)
            {
                hist_new->SetBinContent(i-1, 0);
                seaching_for_first_last_bin = 3;
                break;
            }
        }
    }
    
    return hist_new;
}

void data_stability_plots(){

    norm_counts_empty_sep19_PTBC = retriveHistogramsChangeBPD("../rootFiles/data_stability.root", "norm_counts_empty_sep19_PTBC", 100, 20);
    norm_counts_empty_oct01_PTBC = retriveHistogramsChangeBPD("../rootFiles/data_stability.root", "norm_counts_empty_oct01_PTBC", 100, 20);
    norm_counts_Bi_sep17_PTBC = retriveHistogramsChangeBPD("../rootFiles/data_stability.root", "norm_counts_Bi_sep17_PTBC", 100, 20);
    norm_counts_Bi_sep23_PTBC = retriveHistogramsChangeBPD("../rootFiles/data_stability.root", "norm_counts_Bi_sep23_PTBC", 100, 20);
    norm_counts_Al5_sep25_PTBC = retriveHistogramsChangeBPD("../rootFiles/data_stability.root", "norm_counts_Al5_sep25_PTBC", 100, 20);
    norm_counts_Al5_oct02_PTBC = retriveHistogramsChangeBPD("../rootFiles/data_stability.root", "norm_counts_Al5_oct02_PTBC", 100, 20);
    norm_counts_emptyTank_oct12_PTBC = retriveHistogramsChangeBPD("../rootFiles/data_stability.root", "norm_counts_emptyTank_oct12_PTBC", 100, 20);
    norm_counts_emptyTank_oct16_PTBC = retriveHistogramsChangeBPD("../rootFiles/data_stability.root", "norm_counts_emptyTank_oct16_PTBC", 100, 20);
    norm_counts_ArgonFull_oct18_PTBC = retriveHistogramsChangeBPD("../rootFiles/data_stability.root", "norm_counts_ArgonFull_oct18_PTBC", 100, 20);
    norm_counts_ArgonFull_oct22_PTBC = retriveHistogramsChangeBPD("../rootFiles/data_stability.root", "norm_counts_ArgonFull_oct22_PTBC", 100, 20);

    norm_counts_empty_sep19_FIMG = retriveHistogramsChangeBPD("../rootFiles/data_stability.root", "norm_counts_empty_sep19_FIMG", 100, 20);
    norm_counts_empty_oct01_FIMG = retriveHistogramsChangeBPD("../rootFiles/data_stability.root", "norm_counts_empty_oct01_FIMG", 100, 20);
    norm_counts_Bi_sep17_FIMG = retriveHistogramsChangeBPD("../rootFiles/data_stability.root", "norm_counts_Bi_sep17_FIMG", 100, 20);
    norm_counts_Bi_sep23_FIMG = retriveHistogramsChangeBPD("../rootFiles/data_stability.root", "norm_counts_Bi_sep23_FIMG", 100, 20);
    norm_counts_Al5_sep25_FIMG = retriveHistogramsChangeBPD("../rootFiles/data_stability.root", "norm_counts_Al5_sep25_FIMG", 100, 20);
    norm_counts_Al5_oct02_FIMG = retriveHistogramsChangeBPD("../rootFiles/data_stability.root", "norm_counts_Al5_oct02_FIMG", 100, 20);
    norm_counts_emptyTank_oct12_FIMG = retriveHistogramsChangeBPD("../rootFiles/data_stability.root", "norm_counts_emptyTank_oct12_FIMG", 100, 20);
    norm_counts_emptyTank_oct16_FIMG = retriveHistogramsChangeBPD("../rootFiles/data_stability.root", "norm_counts_emptyTank_oct16_FIMG", 100, 20);
    norm_counts_ArgonFull_oct18_FIMG = retriveHistogramsChangeBPD("../rootFiles/data_stability.root", "norm_counts_ArgonFull_oct18_FIMG", 100, 20);
    norm_counts_ArgonFull_oct22_FIMG = retriveHistogramsChangeBPD("../rootFiles/data_stability.root", "norm_counts_ArgonFull_oct22_FIMG", 100, 20);

    //Plotting
    SetMArEXStyle();
    
    gStyle->SetStatX(0.27);
    gStyle->SetStatY(0.9);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.17);

    // Plotting PTBC Normalized counts plots
    plot_norm_counts_plots(norm_counts_empty_sep19_PTBC, "Sep 19", norm_counts_empty_oct01_PTBC, "Oct 01", "Normalized Counts - No Filter - PTBC", "norm_counts_empty_Sep19_Oct01_PTBC");
    plot_norm_counts_plots(norm_counts_Bi_sep17_PTBC, "Sep 17", norm_counts_Bi_sep23_PTBC, "Sep 23", "Normalized Counts - Bi (1 cm) - PTBC", "norm_counts_bi_Sep17_Sep23_PTBC");
    plot_norm_counts_plots(norm_counts_Al5_sep25_PTBC, "Sep 25", norm_counts_Al5_oct02_PTBC, "Oct 02", "Normalized Counts - Al (5 cm) - PTBC", "norm_counts_al5_Sep25_Oct02_PTBC");
    plot_norm_counts_plots(norm_counts_emptyTank_oct12_PTBC, "Oct 12", norm_counts_emptyTank_oct16_PTBC, "Oct 16", "Normalized Counts - Empty CF Tank - PTBC", "norm_counts_emptyTank_Oct12_Oct16_PTBC");
    plot_norm_counts_plots(norm_counts_ArgonFull_oct18_PTBC, "Oct 18", norm_counts_ArgonFull_oct22_PTBC, "Oct 22", "Normalized Counts - Argon Tank - PTBC", "norm_counts_argonTank_Oct18_Oct22_PTBC");

    // Plotting FIMG Normalized counts plots
    plot_norm_counts_plots(norm_counts_empty_sep19_FIMG, "Sep 19", norm_counts_empty_oct01_FIMG, "Oct 01", "Normalized Counts - No Filter - FIMG", "norm_counts_empty_Sep19_Oct01_FIMG");
    plot_norm_counts_plots(norm_counts_Bi_sep17_FIMG, "Sep 17", norm_counts_Bi_sep23_FIMG, "Sep 23", "Normalized Counts - Bi (1 cm) - FIMG", "norm_counts_bi_Sep17_Sep23_FIMG");
    plot_norm_counts_plots(norm_counts_Al5_sep25_FIMG, "Sep 25", norm_counts_Al5_oct02_FIMG, "Oct 02", "Normalized Counts - Al (5 cm) - FIMG", "norm_counts_al5_Sep25_Oct02_FIMG");
    plot_norm_counts_plots(norm_counts_emptyTank_oct12_FIMG, "Oct 12", norm_counts_emptyTank_oct16_FIMG, "Oct 16", "Normalized Counts - Empty CF Tank - FIMG", "norm_counts_emptyTank_Oct12_Oct16_FIMG");
    plot_norm_counts_plots(norm_counts_ArgonFull_oct18_FIMG, "Oct 18", norm_counts_ArgonFull_oct22_FIMG, "Oct 22", "Normalized Counts - Argon Tank - FIMG", "norm_counts_argonTank_Oct18_Oct22_FIMG");

}