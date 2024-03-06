/**
 * @file day_night_plots.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2024-03-05
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

//////////// PTBC - day-night plots
TH1D* day_night_Bi_sep18_PTBC = 0;
TH1D* day_night_Al5_sep27_PTBC = 0;
TH1D* day_night_emptyTS_oct05_PTBC = 0;
TH1D* day_night_emptyTank_oct15_PTBC = 0;
TH1D* day_night_Argon_oct22_PTBC = 0;

//////////// FIMG - day-night plots
TH1D* day_night_Bi_sep18_FIMG = 0;
TH1D* day_night_Al5_sep27_FIMG = 0;
TH1D* day_night_emptyTS_oct05_FIMG = 0;
TH1D* day_night_emptyTank_oct15_FIMG = 0;
TH1D* day_night_Argon_oct22_FIMG = 0;

TF1* dn_fits[10];
TH1D* dn_resi_plots[10];
TH1D* dn_pull_plots[10];
Int_t dn_plot_index = 0;

TH1D* retriveHistograms(const char *file_name, const char *hist_name){
    
    TFile* hist_file = TFile::Open(file_name, "READ");
    TH1D* hist_new = (TH1D*)hist_file->Get(hist_name);
    // hist_file->Close();

    return hist_new;
}

void add_plot_to_canvas(TCanvas *canvas, TH1D* hist, TLegend* legend, const char* legend_entry, Int_t line_color_num){

    canvas->cd();
    hist->SetLineWidth(2);
    hist->SetLineColor(line_color_num);
    hist->Draw("SAME");
    legend->AddEntry(hist,legend_entry,"l");

    return;
}

void add_fit_to_canvas(TCanvas *canvas, TLine* fit_line, Double_t fit_value, Int_t line_color_num){

    canvas->cd();
    fit_line = new TLine(0., fit_value, 86400., fit_value);
    fit_line->SetLineColor(line_color_num);
    fit_line->SetLineWidth(2);
    fit_line->Draw("SAME");

    return;
}

void compute_residuals(TH1D* dn_hist){

    std::cout << "Residual plot number " << dn_plot_index+1 << std::endl;
    std::cout << "Fitting " << dn_hist->GetName() << " plot..." << std::endl;
    dn_fits[dn_plot_index] = new TF1(Form("dn_fit_%d", dn_plot_index), "[0]");
    dn_fits[dn_plot_index]->SetLineColor(dn_hist->GetLineColor());
    dn_hist->Fit(dn_fits[dn_plot_index], "0");
    std::cout << "" << std::endl;

    dn_resi_plots[dn_plot_index] = (TH1D*)(dn_hist->Clone(Form("dn_resi_plot_%d", dn_plot_index)));
    Int_t num_bins = dn_hist->GetNbinsX();
    for (Int_t i = 0; i < num_bins; i++)
    {
        Double_t bin_content = dn_resi_plots[dn_plot_index]->GetBinContent(i+1);
        if (bin_content != 0)
        {
            Double_t bin_error = dn_resi_plots[dn_plot_index]->GetBinError(i+1);
            Double_t fit_para = dn_fits[dn_plot_index]->GetParameter(0);
            Double_t fit_para_error = dn_fits[dn_plot_index]->GetParError(0);

            dn_resi_plots[dn_plot_index]->SetBinContent(i+1, bin_content - fit_para);
            Double_t new_error = std::sqrt(bin_error*bin_error + fit_para_error*fit_para_error);
            dn_resi_plots[dn_plot_index]->SetBinError(i+1, new_error);
        }
    }
    
    Int_t num_pull_bins = 11;
    dn_pull_plots[dn_plot_index] = new TH1D(Form("dn_pull_plots_%i", dn_plot_index), "", num_pull_bins, -5, 5);
    for (Int_t i = 0; i < num_bins; i++){
        Double_t fill_content = dn_resi_plots[dn_plot_index]->GetBinContent(i+1) / dn_hist->GetBinError(i+1);
        dn_pull_plots[dn_plot_index]->Fill(fill_content);
    }
    dn_pull_plots[dn_plot_index]->SetLineColor(dn_hist->GetLineColor());
    
    dn_plot_index++;
}

void pull_distributions(){
    
    Int_t tot_bins = day_night_Bi_sep18_PTBC->GetNbinsX();

    for (Int_t i = 0; i < tot_bins; i++)
    {
        ////////////////////////////////////////////// PTBC /////////////////////////////////////////////////
        if (day_night_Bi_sep18_PTBC->GetBinContent(i+1) < 3e-12 || day_night_Bi_sep18_PTBC->GetBinContent(i+1) > 4e-12)
        {
            day_night_Bi_sep18_PTBC->SetBinContent(i+1, 0);
            day_night_Bi_sep18_PTBC->SetBinError(i+1, 0);
        }

        if (day_night_Al5_sep27_PTBC->GetBinContent(i+1) < 2.5e-12 || day_night_Al5_sep27_PTBC->GetBinContent(i+1) > 3.5e-12)
        {
            day_night_Al5_sep27_PTBC->SetBinContent(i+1, 0);
            day_night_Al5_sep27_PTBC->SetBinError(i+1, 0);
        }

        if (day_night_emptyTS_oct05_PTBC->GetBinContent(i+1) < 3.5e-12 || day_night_emptyTS_oct05_PTBC->GetBinContent(i+1) > 4.5e-12)
        {
            day_night_emptyTS_oct05_PTBC->SetBinContent(i+1, 0);
            day_night_emptyTS_oct05_PTBC->SetBinError(i+1, 0);
        }

        if (day_night_emptyTank_oct15_PTBC->GetBinContent(i+1) < 1e-12 || day_night_emptyTank_oct15_PTBC->GetBinContent(i+1) > 2e-12)
        {
            day_night_emptyTank_oct15_PTBC->SetBinContent(i+1, 0);
            day_night_emptyTank_oct15_PTBC->SetBinError(i+1, 0);
        }

        if (day_night_Argon_oct22_PTBC->GetBinContent(i+1) < 1e-12 || day_night_Argon_oct22_PTBC->GetBinContent(i+1) > 2e-12)
        {
            day_night_Argon_oct22_PTBC->SetBinContent(i+1, 0);
            day_night_Argon_oct22_PTBC->SetBinError(i+1, 0);
        }

        ////////////////////////////////////////////// FIMG /////////////////////////////////////////////////
        if (day_night_Bi_sep18_FIMG->GetBinContent(i+1) < 5e-12 || day_night_Bi_sep18_FIMG->GetBinContent(i+1) > 7e-12)
        {
            day_night_Bi_sep18_FIMG->SetBinContent(i+1, 0);
            day_night_Bi_sep18_FIMG->SetBinError(i+1, 0);
        }

        if (day_night_Al5_sep27_FIMG->GetBinContent(i+1) < 4e-12 || day_night_Al5_sep27_FIMG->GetBinContent(i+1) > 6e-12)
        {
            day_night_Al5_sep27_FIMG->SetBinContent(i+1, 0);
            day_night_Al5_sep27_FIMG->SetBinError(i+1, 0);
        }

        if (day_night_emptyTS_oct05_FIMG->GetBinContent(i+1) < 6e-12 || day_night_emptyTS_oct05_FIMG->GetBinContent(i+1) > 8e-12)
        {
            day_night_emptyTS_oct05_FIMG->SetBinContent(i+1, 0);
            day_night_emptyTS_oct05_FIMG->SetBinError(i+1, 0);
        }

        if (day_night_emptyTank_oct15_FIMG->GetBinContent(i+1) < 1e-12 || day_night_emptyTank_oct15_FIMG->GetBinContent(i+1) > 3e-12)
        {
            day_night_emptyTank_oct15_FIMG->SetBinContent(i+1, 0);
            day_night_emptyTank_oct15_FIMG->SetBinError(i+1, 0);
        }

        if (day_night_Argon_oct22_FIMG->GetBinContent(i+1) < 1e-12 || day_night_Argon_oct22_FIMG->GetBinContent(i+1) > 3e-12)
        {
            day_night_Argon_oct22_FIMG->SetBinContent(i+1, 0);
            day_night_Argon_oct22_FIMG->SetBinError(i+1, 0);
        }
    }

    compute_residuals(day_night_Bi_sep18_PTBC);
    compute_residuals(day_night_Al5_sep27_PTBC);
    compute_residuals(day_night_emptyTS_oct05_PTBC);
    compute_residuals(day_night_emptyTank_oct15_PTBC);
    compute_residuals(day_night_Argon_oct22_PTBC);
    compute_residuals(day_night_Bi_sep18_FIMG);
    compute_residuals(day_night_Al5_sep27_FIMG);
    compute_residuals(day_night_emptyTS_oct05_FIMG);
    compute_residuals(day_night_emptyTank_oct15_FIMG);
    compute_residuals(day_night_Argon_oct22_FIMG);

    //////// Residual Plots
    TCanvas *c_dn_resi[2];
    TLegend *l_dn_resi[2];
    Int_t j = 0;

    c_dn_resi[0] = new TCanvas(Form("c_dn_resi_%d", 0)," ");
    c_dn_resi[0]->cd();
    c_dn_resi[0]->Draw();
    // TGaxis::SetExponentOffset(0, 0, "y"); // X and Y offset for Y axis

    dn_resi_plots[j]->GetXaxis()->SetTitle("Seconds in a day");
    dn_resi_plots[j]->GetYaxis()->SetTitle("Residuals from the fit");
    dn_resi_plots[j]->SetTitle("Day-Night Residuals after fit - PTBC");
    dn_resi_plots[j]->SetLineWidth(2);
    dn_resi_plots[j]->GetYaxis()->SetRangeUser(-1e-12,1e-12);
    dn_resi_plots[j]->Draw("E");
    dn_resi_plots[j]->SetStats(0);
    gPad->SetGrid();
    // gPad->SetLogy();

    l_dn_resi[0] = new TLegend(0.72,0.30,0.90,0.45);
    l_dn_resi[0]->AddEntry(dn_resi_plots[j],"Bi (1cm) Sep 18","l");
    j++;

    add_plot_to_canvas(c_dn_resi[0], dn_resi_plots[j], l_dn_resi[0], "Al (5cm) Sep 27", 40);
    j++;
    add_plot_to_canvas(c_dn_resi[0], dn_resi_plots[j], l_dn_resi[0], "Empty (TS) Oct 05", 2);
    j++;
    add_plot_to_canvas(c_dn_resi[0], dn_resi_plots[j], l_dn_resi[0], "Empty Tank Oct 15", 3);
    j++;
    add_plot_to_canvas(c_dn_resi[0], dn_resi_plots[j], l_dn_resi[0], "Argon Tank Oct 22", 7);
    j++;

    l_dn_resi[0]->Draw();
    // c_dn_resi[0]->Print("../plots/stability_plots/day_night_residuals_PTBC.png");

    ///////////////////////////////////////////////////////////////////////////////////////////////
    c_dn_resi[1] = new TCanvas(Form("c_dn_resi_%d", 1)," ");
    c_dn_resi[1]->cd();
    c_dn_resi[1]->Draw();

    dn_resi_plots[j]->GetXaxis()->SetTitle("Seconds in a day");
    dn_resi_plots[j]->GetYaxis()->SetTitle("Residuals from the fit");
    dn_resi_plots[j]->SetTitle("Day-Night Residuals after fit - FIMG");
    dn_resi_plots[j]->SetLineWidth(2);
    dn_resi_plots[j]->GetYaxis()->SetRangeUser(-1e-12,1e-12);
    dn_resi_plots[j]->Draw("E");
    dn_resi_plots[j]->SetStats(0);
    gPad->SetGrid();
    // gPad->SetLogy();

    l_dn_resi[1] = new TLegend(0.72,0.30,0.90,0.45);
    l_dn_resi[1]->AddEntry(dn_resi_plots[j],"Bi (1cm) Sep 18","l");
    j++;

    add_plot_to_canvas(c_dn_resi[1], dn_resi_plots[j], l_dn_resi[1], "Al (5cm) Sep 27", 40);
    j++;
    add_plot_to_canvas(c_dn_resi[1], dn_resi_plots[j], l_dn_resi[1], "Empty (TS) Oct 05", 2);
    j++;
    add_plot_to_canvas(c_dn_resi[1], dn_resi_plots[j], l_dn_resi[1], "Empty Tank Oct 15", 3);
    j++;
    add_plot_to_canvas(c_dn_resi[1], dn_resi_plots[j], l_dn_resi[1], "Argon Tank Oct 22", 7);
    j++;

    l_dn_resi[1]->Draw();
    // c_dn_resi[1]->Print("../plots/stability_plots/day_night_residuals_FIMG.png");

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////// Pull Plots
    TCanvas *c_dn_pull[2];
    TLegend *l_dn_pull[2];
    j = 0;

    c_dn_pull[0] = new TCanvas(Form("c_dn_pull_%d", 0)," ");
    c_dn_pull[0]->cd();
    c_dn_pull[0]->Draw();
    // TGaxis::SetExponentOffset(0, 0, "y"); // X and Y offset for Y axis

    dn_pull_plots[j]->GetXaxis()->SetTitle("(data - fit) / $$/delta data$$");
    dn_pull_plots[j]->GetYaxis()->SetTitle("counts");
    dn_pull_plots[j]->SetTitle("Day-Night Pull Distributions - PTBC");
    dn_pull_plots[j]->SetLineWidth(2);
    // dn_pull_plots[j]->GetYaxis()->SetRangeUser(-1e-12,1e-12);
    dn_pull_plots[j]->Draw();
    dn_pull_plots[j]->SetStats(0);
    gPad->SetGrid();
    // gPad->SetLogy();

    l_dn_pull[0] = new TLegend(0.72,0.30,0.90,0.45);
    l_dn_pull[0]->AddEntry(dn_pull_plots[j],"Bi (1cm) Sep 18","l");
    j++;

    add_plot_to_canvas(c_dn_pull[0], dn_pull_plots[j], l_dn_pull[0], "Al (5cm) Sep 27", 40);
    j++;
    add_plot_to_canvas(c_dn_pull[0], dn_pull_plots[j], l_dn_pull[0], "Empty (TS) Oct 05", 2);
    j++;
    add_plot_to_canvas(c_dn_pull[0], dn_pull_plots[j], l_dn_pull[0], "Empty Tank Oct 15", 3);
    j++;
    add_plot_to_canvas(c_dn_pull[0], dn_pull_plots[j], l_dn_pull[0], "Argon Tank Oct 22", 7);
    j++;

    l_dn_pull[0]->Draw();
    // c_dn_pull[0]->Print("../plots/stability_plots/day_night_pull_PTBC.png");

    ///////////////////////////////////////////////////////////////////////////////////////////////
    c_dn_pull[1] = new TCanvas(Form("c_dn_pull_%d", 1)," ");
    c_dn_pull[1]->cd();
    c_dn_pull[1]->Draw();
    // TGaxis::SetExponentOffset(0, 0, "y"); // X and Y offset for Y axis

    dn_pull_plots[j]->GetXaxis()->SetTitle("(data - fit) / $$/delta data$$");
    dn_pull_plots[j]->GetYaxis()->SetTitle("counts");
    dn_pull_plots[j]->SetTitle("Day-Night Pull Distributions - FIMG");
    dn_pull_plots[j]->SetLineWidth(2);
    // dn_pull_plots[j]->GetYaxis()->SetRangeUser(-1e-12,1e-12);
    dn_pull_plots[j]->Draw();
    dn_pull_plots[j]->SetStats(0);
    gPad->SetGrid();
    // gPad->SetLogy();

    l_dn_pull[1] = new TLegend(0.72,0.30,0.90,0.45);
    l_dn_pull[1]->AddEntry(dn_pull_plots[j],"Bi (1cm) Sep 18","l");
    j++;

    add_plot_to_canvas(c_dn_pull[1], dn_pull_plots[j], l_dn_pull[1], "Al (5cm) Sep 27", 40);
    j++;
    add_plot_to_canvas(c_dn_pull[1], dn_pull_plots[j], l_dn_pull[1], "Empty (TS) Oct 05", 2);
    j++;
    add_plot_to_canvas(c_dn_pull[1], dn_pull_plots[j], l_dn_pull[1], "Empty Tank Oct 15", 3);
    j++;
    add_plot_to_canvas(c_dn_pull[1], dn_pull_plots[j], l_dn_pull[1], "Argon Tank Oct 22", 7);
    j++;

    l_dn_pull[1]->Draw();
    // c_dn_pull[1]->Print("../plots/stability_plots/day_night_pull_FIMG.png");

    return;
}

void day_night_plots(){

    day_night_Bi_sep18_PTBC = retriveHistograms("../rootFiles/data_stability.root", "day_night_Bi_sep18_PTBC");
    day_night_Al5_sep27_PTBC = retriveHistograms("../rootFiles/data_stability.root", "day_night_Al5_sep27_PTBC");
    day_night_emptyTS_oct05_PTBC = retriveHistograms("../rootFiles/data_stability.root", "day_night_emptyTS_oct05_PTBC");
    day_night_emptyTank_oct15_PTBC = retriveHistograms("../rootFiles/data_stability.root", "day_night_emptyTank_oct15_PTBC");
    day_night_Argon_oct22_PTBC = retriveHistograms("../rootFiles/data_stability.root", "day_night_Argon_oct22_PTBC");

    day_night_Bi_sep18_FIMG = retriveHistograms("../rootFiles/data_stability.root", "day_night_Bi_sep18_FIMG");
    day_night_Al5_sep27_FIMG = retriveHistograms("../rootFiles/data_stability.root", "day_night_Al5_sep27_FIMG");
    day_night_emptyTS_oct05_FIMG = retriveHistograms("../rootFiles/data_stability.root", "day_night_emptyTS_oct05_FIMG");
    day_night_emptyTank_oct15_FIMG = retriveHistograms("../rootFiles/data_stability.root", "day_night_emptyTank_oct15_FIMG");
    day_night_Argon_oct22_FIMG = retriveHistograms("../rootFiles/data_stability.root", "day_night_Argon_oct22_FIMG");

    // day_night_Bi_sep18_PTBC->SetLineColor();
    day_night_Al5_sep27_PTBC->SetLineColor(40);
    day_night_emptyTS_oct05_PTBC->SetLineColor(2);
    day_night_emptyTank_oct15_PTBC->SetLineColor(3);
    day_night_Argon_oct22_PTBC->SetLineColor(7);
    // day_night_Bi_sep18_FIMG->SetLineColor();
    day_night_Al5_sep27_FIMG->SetLineColor(40);
    day_night_emptyTS_oct05_FIMG->SetLineColor(2);
    day_night_emptyTank_oct15_FIMG->SetLineColor(3);
    day_night_Argon_oct22_FIMG->SetLineColor(7);


    //Plotting
    SetMArEXStyle();
    
    gStyle->SetStatX(0.27);
    gStyle->SetStatY(0.9);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.17);

    pull_distributions();

    ///////////////////////////////////////////////////////////////////////////////////////////////
    TCanvas *c_dn[2];
    TLegend *l_dn[2];
    TLine *fit_lines[10];
    Int_t j = 0;
    Int_t fit_index = 0;

    c_dn[j] = new TCanvas(Form("c_dn%d", j)," ");
    c_dn[j]->cd();
    c_dn[j]->Draw();
    // TGaxis::SetExponentOffset(0, 0, "y"); // X and Y offset for Y axis

    day_night_Bi_sep18_PTBC->GetXaxis()->SetTitle("Seconds in a day");
    day_night_Bi_sep18_PTBC->GetYaxis()->SetTitle("Normalized Counts");
    day_night_Bi_sep18_PTBC->SetTitle("Day-Night Variation in Counts - PTBC");
    day_night_Bi_sep18_PTBC->SetLineWidth(2);
    day_night_Bi_sep18_PTBC->GetYaxis()->SetRangeUser(5e-13,5e-12);
    day_night_Bi_sep18_PTBC->Draw();
    day_night_Bi_sep18_PTBC->SetStats(0);
    gPad->SetGrid();
    // gPad->SetLogy();

    l_dn[j] = new TLegend(0.72,0.30,0.90,0.45);
    l_dn[j]->AddEntry(day_night_Bi_sep18_PTBC,"Bi (1cm) Sep 18","l");

    add_plot_to_canvas(c_dn[j], day_night_Al5_sep27_PTBC, l_dn[j], "Al (5cm) Sep 27", 40);
    add_plot_to_canvas(c_dn[j], day_night_emptyTS_oct05_PTBC, l_dn[j], "Empty (TS) Oct 05", 2);
    add_plot_to_canvas(c_dn[j], day_night_emptyTank_oct15_PTBC, l_dn[j], "Empty Tank Oct 15", 3);
    add_plot_to_canvas(c_dn[j], day_night_Argon_oct22_PTBC, l_dn[j], "Argon Tank Oct 22", 7);

    add_fit_to_canvas(c_dn[j], fit_lines[fit_index], dn_fits[fit_index]->GetParameter(0), 4);
    fit_index++;
    add_fit_to_canvas(c_dn[j], fit_lines[fit_index], dn_fits[fit_index]->GetParameter(0), 40);
    fit_index++;
    add_fit_to_canvas(c_dn[j], fit_lines[fit_index], dn_fits[fit_index]->GetParameter(0), 2);
    fit_index++;
    add_fit_to_canvas(c_dn[j], fit_lines[fit_index], dn_fits[fit_index]->GetParameter(0), 3);
    fit_index++;
    add_fit_to_canvas(c_dn[j], fit_lines[fit_index], dn_fits[fit_index]->GetParameter(0), 7);
    fit_index++;

    l_dn[j]->Draw();
    // c_dn[j]->Print("../plots/stability_plots/day_night_vaiation_PTBC.png");

    ///////////////////////////////////////////////////////////////////////////////////////////////
    j++;
    c_dn[j] = new TCanvas(Form("c_dn%d", j)," ");
    c_dn[j]->cd();
    c_dn[j]->Draw();

    day_night_Bi_sep18_FIMG->GetXaxis()->SetTitle("Seconds in a day");
    day_night_Bi_sep18_FIMG->GetYaxis()->SetTitle("Normalized Counts");
    day_night_Bi_sep18_FIMG->SetTitle("Day-Night Variation in Counts - FIMG");
    day_night_Bi_sep18_FIMG->SetLineWidth(2);
    day_night_Bi_sep18_FIMG->GetYaxis()->SetRangeUser(1e-12,9e-12);
    day_night_Bi_sep18_FIMG->Draw();
    day_night_Bi_sep18_FIMG->SetStats(0);
    gPad->SetGrid();
    // gPad->SetLogy();

    l_dn[j] = new TLegend(0.72,0.30,0.90,0.45);
    l_dn[j]->AddEntry(day_night_Bi_sep18_FIMG,"Bi (1cm) Sep 18","l");

    add_plot_to_canvas(c_dn[j], day_night_Al5_sep27_FIMG, l_dn[j], "Al (5cm) Sep 27", 40);
    add_plot_to_canvas(c_dn[j], day_night_emptyTS_oct05_FIMG, l_dn[j], "Empty (TS) Oct 05", 2);
    add_plot_to_canvas(c_dn[j], day_night_emptyTank_oct15_FIMG, l_dn[j], "Empty Tank Oct 15", 3);
    add_plot_to_canvas(c_dn[j], day_night_Argon_oct22_FIMG, l_dn[j], "Argon Tank Oct 22", 7);

    add_fit_to_canvas(c_dn[j], fit_lines[fit_index], dn_fits[fit_index]->GetParameter(0), 4);
    fit_index++;
    add_fit_to_canvas(c_dn[j], fit_lines[fit_index], dn_fits[fit_index]->GetParameter(0), 40);
    fit_index++;
    add_fit_to_canvas(c_dn[j], fit_lines[fit_index], dn_fits[fit_index]->GetParameter(0), 2);
    fit_index++;
    add_fit_to_canvas(c_dn[j], fit_lines[fit_index], dn_fits[fit_index]->GetParameter(0), 3);
    fit_index++;
    add_fit_to_canvas(c_dn[j], fit_lines[fit_index], dn_fits[fit_index]->GetParameter(0), 7);
    fit_index++;

    l_dn[j]->Draw();
    // c_dn[j]->Print("../plots/stability_plots/day_night_vaiation_FIMG.png");
}

// resi_Bi_sep18_PTBC
// resi_Al5_sep27_PTBC
// resi_emptyTS_oct05_PTBC
// resi_emptyTank_oct15_PTBC
// resi_Argon_oct22_PTBC
// resi_Bi_sep18_FIMG
// resi_Al5_sep27_FIMG
// resi_emptyTS_oct05_FIMG
// resi_emptyTank_oct15_FIMG
// resi_Argon_oct22_FIMG

// fit_Bi_sep18_PTBC
// fit_Al5_sep27_PTBC
// fit_emptyTS_oct05_PTBC
// fit_emptyTank_oct15_PTBC
// fit_Argon_oct22_PTBC
// fit_Bi_sep18_FIMG
// fit_Al5_sep27_FIMG
// fit_emptyTS_oct05_FIMG
// fit_emptyTank_oct15_FIMG
// fit_Argon_oct22_FIMG

// day_night_Bi_sep18_PTBC
// day_night_Al5_sep27_PTBC
// day_night_emptyTS_oct05_PTBC
// day_night_emptyTank_oct15_PTBC
// day_night_Argon_oct22_PTBC
// day_night_Bi_sep18_FIMG
// day_night_Al5_sep27_FIMG
// day_night_emptyTS_oct05_FIMG
// day_night_emptyTank_oct15_FIMG
// day_night_Argon_oct22_FIMG