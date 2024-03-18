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
TLine* dn_fit_lines[10];
TCanvas* dn_canvas[10];
TLegend* dn_legends[10];
TPad* dn_pads[10][3];
TLine* zero_line;
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
    
    Int_t num_pull_bins = 15;
    dn_pull_plots[dn_plot_index] = new TH1D(Form("dn_pull_plots_%i", dn_plot_index), "", num_pull_bins, -5, 5);
    for (Int_t i = 0; i < num_bins; i++){
        Double_t fill_content = dn_resi_plots[dn_plot_index]->GetBinContent(i+1) / dn_resi_plots[dn_plot_index]->GetBinError(i+1);
        dn_pull_plots[dn_plot_index]->Fill(fill_content);
    }
    dn_pull_plots[dn_plot_index]->SetLineColor(dn_hist->GetLineColor());
    
    dn_plot_index++;
}

void plot_dn_plots(TH1D* dn_hist, TH1D* resi_hist, TH1D* pull_hist, const char* plot_title_1, const char* output_file_name){

    dn_canvas[dn_plot_index] = new TCanvas(Form("dn_canvas_%d", dn_plot_index)," ");
    dn_canvas[dn_plot_index]->cd();
    dn_canvas[dn_plot_index]->Draw();

    gStyle->SetTitleH(0.1);

    // Pad 1
    dn_canvas[dn_plot_index]->cd(0);
    dn_pads[dn_plot_index][0] = new TPad(Form("dn_pads_0_%d", dn_plot_index), Form("dn_pads_0_%d", dn_plot_index), 0., 0.65, 1., 1.);
    dn_pads[dn_plot_index][0]->SetFillColor(kWhite);
    dn_pads[dn_plot_index][0]->SetBottomMargin(0.00001);
    dn_pads[dn_plot_index][0]->SetBorderMode(0);
    dn_pads[dn_plot_index][0]->Draw();
    dn_pads[dn_plot_index][0]->cd();
    // X-Axis
    dn_hist->GetXaxis()->SetLabelOffset(999);
    dn_hist->GetXaxis()->SetLabelSize(0);
    dn_hist->GetXaxis()->SetTitle("");
    // Y-Axis
    dn_hist->GetYaxis()->SetTitle("Normalized Counts");
    dn_hist->GetYaxis()->SetLabelSize(0.07);
    dn_hist->GetYaxis()->SetTitleSize(0.09);
    dn_hist->GetYaxis()->SetTitleOffset(0.4);
    // Setting the y axis range
    Double_t y_range_min = dn_fits[dn_plot_index]->GetParameter(0) - 0.3e-12;
    Double_t y_range_max = dn_fits[dn_plot_index]->GetParameter(0) + 0.3e-12;
    dn_hist->SetMaximum(y_range_max);
    dn_hist->SetMinimum(y_range_min);
    TGaxis::SetExponentOffset(-0.06, -0.8, "y"); // X and Y offset for Y axis
    dn_hist->SetTitle(plot_title_1);
    dn_hist->SetLineWidth(2);
    dn_hist->SetLineColor(2);
    // dn_hist->GetYaxis()->SetRangeUser(3e-12, 4e-12);
    dn_hist->Draw();
    dn_hist->SetStats(0);
    gPad->SetGrid();

    dn_fit_lines[dn_plot_index] = new TLine(0., dn_fits[dn_plot_index]->GetParameter(0), 86400., dn_fits[dn_plot_index]->GetParameter(0));
    dn_fit_lines[dn_plot_index]->SetLineWidth(2);
    dn_fit_lines[dn_plot_index]->SetLineColor(4);
    dn_fit_lines[dn_plot_index]->SetLineStyle(9);
    dn_fit_lines[dn_plot_index]->Draw("SAME");

    dn_legends[dn_plot_index] = new TLegend(0.80, 0.60, 0.85, 0.80);
    dn_legends[dn_plot_index]->AddEntry(dn_hist, " data", "l");
    dn_legends[dn_plot_index]->AddEntry(dn_fit_lines[dn_plot_index], " fit", "l");
    dn_legends[dn_plot_index]->Draw();

    // Pad 2
    dn_canvas[dn_plot_index]->cd(0);
    dn_pads[dn_plot_index][1] = new TPad(Form("dn_pads_1_%d", dn_plot_index), Form("dn_pads_1_%d", dn_plot_index), 0., 0.3, 1., 0.65);
    dn_pads[dn_plot_index][1]->SetFillColor(kWhite);
    dn_pads[dn_plot_index][1]->SetTopMargin(0.00001);
    dn_pads[dn_plot_index][1]->SetBottomMargin(0.15);
    dn_pads[dn_plot_index][1]->SetBorderMode(0);
    dn_pads[dn_plot_index][1]->Draw();
    dn_pads[dn_plot_index][1]->cd();
    
    resi_hist->SetTitle("");
    // X-Axis
    resi_hist->GetXaxis()->SetTitle("Seconds in a day");
    resi_hist->GetXaxis()->SetTitleOffset(0.75);
    resi_hist->GetXaxis()->SetLabelSize(0.07);
    resi_hist->GetXaxis()->SetTitleSize(0.09);
    // Y-Axis
    resi_hist->GetYaxis()->SetTitle("Residuals");
    resi_hist->GetYaxis()->SetLabelSize(0.07);
    resi_hist->GetYaxis()->SetTitleSize(0.09);
    resi_hist->GetYaxis()->SetTitleOffset(0.4);
    resi_hist->SetLineColor(2);
    resi_hist->SetLineWidth(2);
    resi_hist->GetYaxis()->SetRangeUser(-0.3e-12, 0.3e-12);
    resi_hist->Draw("E");
    resi_hist->SetStats(0);

    zero_line->Draw("SAME");

    // Pad 3
    dn_canvas[dn_plot_index]->cd(0);
    dn_pads[dn_plot_index][2] = new TPad(Form("dn_pads_2_%d", dn_plot_index), Form("dn_pads_2_%d", dn_plot_index), 0., 0., 1., 0.3);
    dn_pads[dn_plot_index][2]->SetFillColor(kWhite);
    dn_pads[dn_plot_index][2]->SetBottomMargin(0.2);
    dn_pads[dn_plot_index][2]->SetBorderMode(0);
    dn_pads[dn_plot_index][2]->Draw();
    dn_pads[dn_plot_index][2]->cd();
    // X-Axis
    pull_hist->GetXaxis()->SetTitle("(data - fit) / #Delta data");
    pull_hist->GetXaxis()->SetTitleSize(0.1);
    // pull_hist->GetXaxis()->SetLabelOffset(0.01);
    pull_hist->GetXaxis()->SetTitleOffset(0.8);
    pull_hist->GetXaxis()->SetLabelSize(0.09);
    // Y-Axis
    pull_hist->GetYaxis()->SetTitle("Counts");
    pull_hist->GetYaxis()->SetTitleSize(0.1);
    pull_hist->GetYaxis()->SetTitleOffset(0.36);
    pull_hist->GetYaxis()->SetLabelSize(0.09);

    pull_hist->SetTitle("Pull Distribution");
    pull_hist->SetLineColor(2);
    pull_hist->SetLineWidth(2);
    pull_hist->Draw();
    pull_hist->SetStats(0);
    gPad->SetGrid();
    
    dn_canvas[dn_plot_index]->Print( Form("../plots/stability_plots/%s.png", output_file_name) );

    dn_plot_index++;
    return;
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

    dn_plot_index = 0;
    plot_dn_plots(day_night_Bi_sep18_PTBC, dn_resi_plots[dn_plot_index], dn_pull_plots[dn_plot_index], "Day-Night Variation in Counts - Bi (1cm) Sep 18 - PTBC", "dn_vaiation_bi_sep18_PTBC");
    plot_dn_plots(day_night_Al5_sep27_PTBC, dn_resi_plots[dn_plot_index], dn_pull_plots[dn_plot_index], "Day-Night Variation in Counts - Al (5cm) Sep 27 - PTBC", "dn_vaiation_al5_sep27_PTBC");
    plot_dn_plots(day_night_emptyTS_oct05_PTBC, dn_resi_plots[dn_plot_index], dn_pull_plots[dn_plot_index], "Day-Night Variation in Counts - Empty (TS) Oct 05 - PTBC", "dn_vaiation_emptyTS_oct05_PTBC");
    plot_dn_plots(day_night_emptyTank_oct15_PTBC, dn_resi_plots[dn_plot_index], dn_pull_plots[dn_plot_index], "Day-Night Variation in Counts - Empty Tank Oct 15 - PTBC", "dn_vaiation_emptyTank_oct15_PTBC");
    plot_dn_plots(day_night_Argon_oct22_PTBC, dn_resi_plots[dn_plot_index], dn_pull_plots[dn_plot_index], "Day-Night Variation in Counts - Argon Tank Oct 22 - PTBC", "dn_vaiation_argon_oct22_PTBC");
    plot_dn_plots(day_night_Bi_sep18_FIMG, dn_resi_plots[dn_plot_index], dn_pull_plots[dn_plot_index], "Day-Night Variation in Counts - Bi (1cm) Sep 18 - FIMG", "dn_vaiation_bi_sep18_FIMG");
    plot_dn_plots(day_night_Al5_sep27_FIMG, dn_resi_plots[dn_plot_index], dn_pull_plots[dn_plot_index], "Day-Night Variation in Counts - Al (5cm) Sep 27 - FIMG", "dn_vaiation_al5_sep27_FIMG");
    plot_dn_plots(day_night_emptyTS_oct05_FIMG, dn_resi_plots[dn_plot_index], dn_pull_plots[dn_plot_index], "Day-Night Variation in Counts - Empty (TS) Oct 05 - FIMG", "dn_vaiation_emptyTS_oct05_FIMG");
    plot_dn_plots(day_night_emptyTank_oct15_FIMG, dn_resi_plots[dn_plot_index], dn_pull_plots[dn_plot_index], "Day-Night Variation in Counts - Empty Tank Oct 15 - FIMG", "dn_vaiation_emptyTank_oct15_FIMG");
    plot_dn_plots(day_night_Argon_oct22_FIMG, dn_resi_plots[dn_plot_index], dn_pull_plots[dn_plot_index], "Day-Night Variation in Counts - Argon Tank Oct 22 - FIMG", "dn_vaiation_argon_oct22_FIMG");


    return;
}

void plot_all_dn_plots(){
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
    return;
}

void day_night_plots(){

    day_night_Bi_sep18_PTBC = retriveHistograms("../rootFiles/data_stability_1hr_bins.root", "day_night_Bi_sep18_PTBC");
    day_night_Al5_sep27_PTBC = retriveHistograms("../rootFiles/data_stability_1hr_bins.root", "day_night_Al5_sep27_PTBC");
    day_night_emptyTS_oct05_PTBC = retriveHistograms("../rootFiles/data_stability_1hr_bins.root", "day_night_emptyTS_oct05_PTBC");
    day_night_emptyTank_oct15_PTBC = retriveHistograms("../rootFiles/data_stability_1hr_bins.root", "day_night_emptyTank_oct15_PTBC");
    day_night_Argon_oct22_PTBC = retriveHistograms("../rootFiles/data_stability_1hr_bins.root", "day_night_Argon_oct22_PTBC");

    day_night_Bi_sep18_FIMG = retriveHistograms("../rootFiles/data_stability_1hr_bins.root", "day_night_Bi_sep18_FIMG");
    day_night_Al5_sep27_FIMG = retriveHistograms("../rootFiles/data_stability_1hr_bins.root", "day_night_Al5_sep27_FIMG");
    day_night_emptyTS_oct05_FIMG = retriveHistograms("../rootFiles/data_stability_1hr_bins.root", "day_night_emptyTS_oct05_FIMG");
    day_night_emptyTank_oct15_FIMG = retriveHistograms("../rootFiles/data_stability_1hr_bins.root", "day_night_emptyTank_oct15_FIMG");
    day_night_Argon_oct22_FIMG = retriveHistograms("../rootFiles/data_stability_1hr_bins.root", "day_night_Argon_oct22_FIMG");

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

    zero_line = new TLine(0., 0., 86400., 0.);
    zero_line->SetLineWidth(2);
    zero_line->SetLineColor(1);
    zero_line->SetLineStyle(2);

    //Plotting
    SetMArEXStyle();
    
    gStyle->SetStatX(0.27);
    gStyle->SetStatY(0.9);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.17);

    pull_distributions();

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