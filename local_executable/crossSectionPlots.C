/**
 * @file crossSectionPlots.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-11-13
 */

#include "crossSectionPlots.h"

TH1D* retriveHistograms(const char *fname, const char *hist_name){
    
    TFile* hist_file = TFile::Open(fname, "READ");
    // if (!hist_file || hist_file->IsZombie()) {
    //     cout << "Unable to open " << fname << " for reading..." <<endl;
    //     return;
    // }

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

void endf(Double_t n, Int_t num_bins_e, Double_t e_bin_low_edge, Double_t e_bin_up_edge){
    //Extracting ENDF Cross Section and transmission

    fillMaxLineCount();

    std::ifstream inputFile(Form("../evalData/%s", eval_file_name_map[filter_name].c_str())); // endf_file_name.c_str()

    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the file.\n";
    }

    //Extracting the data from the text file
    Double_t val1, val2; //val1 -> Energy (eV); val2 -> Cross Section (barns)
    std::string line;
    Double_t line_count = 0;

    while (std::getline(inputFile, line)) {

        line_count++;

        if (line_count == 1)
        {
            continue;
        }

        if (line_count == evalFile_max_line_num) //Bi - 27019; Al - 9745
        {
            break;
        }

        std::istringstream iss(line);

        if (iss >> val1 >> val2) {
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
        if (endf_e[i] < e_bin_low_edge) continue;
        if (endf_e[i] >= e_bin_up_edge) break;

        Int_t bin_num = endf_xsec_hist->GetXaxis()->FindBin(endf_e[i]);
        if (bin_num > num_bins_e) break;

        bin_content_map[bin_num].push_back(endf_xsec[i]);
    }

    std::vector<Double_t> bin_contents(num_bins_e, 0.0);

    for (const auto& bin : bin_content_map) {
        Int_t bin_n = bin.first;
        const std::vector<Double_t>& xsecs = bin.second;
        if (xsecs.empty()) {
            endf_xsec_hist->SetBinContent(bin_n, 0.0);
            bin_contents[bin_n - 1] = 0.0;
        } else {
            Double_t bin_content = std::accumulate(xsecs.begin(), xsecs.end(), 0.0) / xsecs.size();
            endf_xsec_hist->SetBinContent(bin_n, bin_content);
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
            Double_t left_e_val = endf_xsec_hist->GetBinCenter(left_idx + 1);
            Double_t right_value = bin_contents[right_idx];
            Double_t right_e_val = endf_xsec_hist->GetBinCenter(right_idx + 1);

            // Linear Interpolation
            Double_t slope = (right_value - left_value) / (right_e_val - left_e_val);
            Int_t num_bins_to_update = right_idx - left_idx - 1;
            for (Int_t j = 1; j <= num_bins_to_update; ++j) {
                Double_t bin_center = endf_xsec_hist->GetBinCenter(left_idx + j + 1);
                bin_contents[left_idx + j] = slope * (bin_center - left_e_val) + left_value;
            }
        }
    }

    // Update histogram with interpolated values
    for (int i = 1; i <= num_bins_e; ++i) {
        endf_xsec_hist->SetBinContent(i, bin_contents[i - 1]);
        if (bin_contents[i - 1] == 0)
        {
            endf_trans_hist->SetBinContent(i, 0);
            continue;
        }
        Double_t trans_val = std::exp(- n * bin_contents[i - 1]);
        endf_trans_hist->SetBinContent(i, trans_val);
    }
}

void endf_argon(Double_t n, Int_t bpd){
    //Extracting ENDF Cross Section and transmission

    cout << "Opening ENDF Argon hist File ../inputFiles/natAr_xsec_hists.root" << endl;
    TFile *hist_file = TFile::Open("../inputFiles/natAr_xsec_hists.root", "READ");

    if (bpd == 20)
    {
        endf_xsec_hist = (TH1D*)hist_file->Get("natAr_xsec_hist_20bpd_endf");
        endf_trans_hist = (TH1D*)endf_xsec_hist->Clone("natAr_trans_hist_20bpd_endf");
    }

    if (bpd == 50)
    {
        endf_xsec_hist = (TH1D*)hist_file->Get("natAr_xsec_hist_50bpd_endf");
        endf_trans_hist = (TH1D*)endf_xsec_hist->Clone("natAr_trans_hist_50bpd_endf");
    }

    if (bpd == 100)
    {
        endf_xsec_hist = (TH1D*)hist_file->Get("natAr_xsec_hist_100bpd_endf");
        endf_trans_hist = (TH1D*)endf_xsec_hist->Clone("natAr_trans_hist_100bpd_endf");
    }
    
    Int_t num_bins_hist = endf_trans_hist->GetNbinsX();

    for (Int_t i = 1; i < num_bins_hist+1; i++)
    {
        Double_t xsec_val = endf_trans_hist->GetBinContent(i);
        Double_t trans_val = std::exp(- n * xsec_val);
        endf_trans_hist->SetBinContent(i, trans_val);
    }
}

void jendl_argon(Double_t n, Int_t bpd){
    //Extracting JENDL Cross Section and transmission

    cout << "Opening JENDL Argon hist File ../inputFiles/natAr_xsec_hists.root" << endl;
    TFile *hist_file = TFile::Open("../inputFiles/natAr_xsec_hists.root", "READ");

    if (bpd == 20)
    {
        jendl_xsec_hist = (TH1D*)hist_file->Get("natAr_xsec_hist_20bpd_jendl");
        jendl_trans_hist = (TH1D*)jendl_xsec_hist->Clone("natAr_trans_hist_20bpd_jendl");
    }

    if (bpd == 50)
    {
        jendl_xsec_hist = (TH1D*)hist_file->Get("natAr_xsec_hist_50bpd_jendl");
        jendl_trans_hist = (TH1D*)jendl_xsec_hist->Clone("natAr_trans_hist_50bpd_jendl");
    }

    if (bpd == 100)
    {
        jendl_xsec_hist = (TH1D*)hist_file->Get("natAr_xsec_hist_100bpd_jendl");
        jendl_trans_hist = (TH1D*)jendl_xsec_hist->Clone("natAr_trans_hist_100bpd_jendl");
    }

    Int_t num_bins_hist = jendl_trans_hist->GetNbinsX();

    for (Int_t i = 1; i < num_bins_hist+1; i++)
    {
        Double_t xsec_val = jendl_trans_hist->GetBinContent(i);
        Double_t trans_val = std::exp(- n * xsec_val);
        jendl_trans_hist->SetBinContent(i, trans_val);
    }
}


////// Do not use for transmission calculation
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

void changeHistBPD(TH1D* hist_new, TH1D* hist_old, Int_t num_bins_new, Double_t bin_edges_new[]){

    Int_t num_bins_old = hist_old->GetNbinsX();
    Int_t max_sum_bin_count = (num_bins_old/num_bins_new); // Number of old bins that need to be summed to make a new bin

    // Double_t bin_edges_new[num_bins_new + 1];
    Double_t bin_content[num_bins_new];
    Double_t bin_error[num_bins_new];
    Int_t sum_bin_count = 0;
    Int_t new_bin_counter = 0;
    for (int i = 1; i < num_bins_old + 1; i++)
    {
        if (i == 1)
        {
            // bin_edges_new[new_bin_counter] = hist_old->GetXaxis()->GetBinLowEdge(i);
            bin_content[new_bin_counter] = 0;
            bin_error[new_bin_counter] = 0;
        }
        bin_content[new_bin_counter] += hist_old->GetBinContent(i);
        sum_bin_count++;

        if (i == num_bins_old)
        {
            sum_bin_count = 0;
            bin_error[new_bin_counter] = sqrt(bin_content[new_bin_counter]);
            // bin_edges_new[new_bin_counter+1] = hist_old->GetXaxis()->GetBinUpEdge(i);
        }
        else if (sum_bin_count == max_sum_bin_count)
        {
            sum_bin_count = 0;
            bin_error[new_bin_counter] = sqrt(bin_content[new_bin_counter]);
            new_bin_counter++;
            // bin_edges_new[new_bin_counter] = hist_old->GetXaxis()->GetBinUpEdge(i);
            bin_content[new_bin_counter] = 0;
            bin_error[new_bin_counter] = 0;
        } 
    }

    for(int i = 0; i < num_bins_new; i++){
        hist_new->SetBinContent(i+1, bin_content[i]);
        hist_new->SetBinError(i+1, bin_error[i]);
    }

    return;
}

void calc_xsec(const char *fname, Int_t num_bins_e, Double_t bin_edges_e[]){
    cout << "Opening File " << fname << endl;
    TFile *hist_file = TFile::Open(fname, "READ");

    cout << "Extracting Histograms " << endl;
    energy_hist_target_in_PTBC = (TH1D*)hist_file->Get("energy_hist_target_in_PTB");
    energy_hist_target_in_FIMG = (TH1D*)hist_file->Get("energy_hist_target_in_FIMG");
    energy_hist_target_out_PTBC = (TH1D*)hist_file->Get("energy_hist_target_out_PTB");
    energy_hist_target_out_FIMG = (TH1D*)hist_file->Get("energy_hist_target_out_FIMG");

    cout << "Extracting Norm Factors " << endl;
    // norm_factor_target_in = (Double_t)hist_file->Get("norm_factor_target_in");
    // norm_factor_target_out = (Double_t)hist_file->Get("norm_factor_target_out");
    std::vector<Double_t> *norm_factors_temp;
    hist_file->GetObject("norm_factors", norm_factors_temp);
    norm_factors = *norm_factors_temp;
    Double_t Qout_Qin = norm_factors[1]/norm_factors[0]; //norm_factor_target_out / norm_factor_target_in;

    TH1D* energy_hist_target_in_PTBC_new = 0;
    TH1D* energy_hist_target_in_FIMG_new = 0;
    TH1D* energy_hist_target_out_PTBC_new = 0;
    TH1D* energy_hist_target_out_FIMG_new = 0;

    //Target In Hists
    energy_hist_target_in_PTBC_new = new TH1D("energy_hist_target_in_PTBC_new","Energy Hist Target In - PTBC", num_bins_e, bin_edges_e);
    energy_hist_target_in_FIMG_new = new TH1D("energy_hist_target_in_FIMG_new","Energy Hist Target In - FIMG", num_bins_e, bin_edges_e);
    //Target Out Hists
    energy_hist_target_out_PTBC_new = new TH1D("energy_hist_target_out_PTBC_new","Energy Hist Target Out - PTBC", num_bins_e, bin_edges_e);
    energy_hist_target_out_FIMG_new = new TH1D("energy_hist_target_out_FIMG_new","Energy Hist Target Out - FIMG", num_bins_e, bin_edges_e);

    changeHistBPD(energy_hist_target_in_PTBC_new, energy_hist_target_in_PTBC, num_bins_e, bin_edges_e);
    changeHistBPD(energy_hist_target_in_FIMG_new, energy_hist_target_in_FIMG, num_bins_e, bin_edges_e);
    changeHistBPD(energy_hist_target_out_PTBC_new, energy_hist_target_out_PTBC, num_bins_e, bin_edges_e);
    changeHistBPD(energy_hist_target_out_FIMG_new, energy_hist_target_out_FIMG, num_bins_e, bin_edges_e);

    //transmission histogram
    transmission_hist_e_PTBC = new TH1D("transmission_hist_e_PTBC","Transmission Hist - PTBC",num_bins_e,bin_edges_e);
    transmission_hist_e_FIMG = new TH1D("transmission_hist_e_FIMG","Transmission Hist - FIMG",num_bins_e,bin_edges_e);
    //cross section histogram
    cross_section_hist_e_PTBC = new TH1D("cross_section_hist_e_PTBC","Cross Section Hist - PTBC",num_bins_e,bin_edges_e);
    cross_section_hist_e_FIMG = new TH1D("cross_section_hist_e_FIMG","Cross Section Hist - FIMG",num_bins_e,bin_edges_e);

    //Transmission
    for (int i = 0; i < num_bins_e; i++)
    {
        //PTBC
        Double_t bin_content_in_PTBC = energy_hist_target_in_PTBC_new->GetBinContent(i+1);
        Double_t bin_content_out_PTBC = energy_hist_target_out_PTBC_new->GetBinContent(i+1);
        if (bin_content_in_PTBC == 0. || bin_content_out_PTBC == 0.)
        {
            transmission_hist_e_PTBC->SetBinContent(i+1, 0.);
            transmission_hist_e_PTBC->SetBinError(i+1, 0.);
        } else {
            Double_t transmission_PTBC = (bin_content_in_PTBC * Qout_Qin)/bin_content_out_PTBC;
            Double_t bin_unc_PTBC = transmission_PTBC * std::sqrt( (1./bin_content_in_PTBC) + (1./bin_content_out_PTBC) );
            transmission_hist_e_PTBC->SetBinContent(i+1, transmission_PTBC);
            transmission_hist_e_PTBC->SetBinError(i+1, bin_unc_PTBC);
        }

        //FIMG
        Double_t bin_content_in_FIMG = energy_hist_target_in_FIMG_new->GetBinContent(i+1);
        Double_t bin_content_out_FIMG = energy_hist_target_out_FIMG_new->GetBinContent(i+1);
        if (bin_content_in_FIMG == 0. || bin_content_out_FIMG == 0.)
        {
            transmission_hist_e_FIMG->SetBinContent(i+1, 0.);
            transmission_hist_e_FIMG->SetBinError(i+1, 0.);
        } else {
            Double_t transmission_FIMG = (bin_content_in_FIMG * Qout_Qin)/bin_content_out_FIMG;
            Double_t bin_unc_FIMG = transmission_FIMG * std::sqrt( (1./bin_content_in_FIMG) + (1./bin_content_out_FIMG) );
            transmission_hist_e_FIMG->SetBinContent(i+1, transmission_FIMG);
            transmission_hist_e_FIMG->SetBinError(i+1, bin_unc_FIMG);
        }
    }

    //Cross Section
    Double_t num_density = num_density_map[filter_name];
    cout << Form("Number density of %s = ", filter_name.c_str()) << num_density << endl;
    Double_t n_inverse = ( (Double_t) 1.0/num_density);
    cout << Form("1/n of %s = ", filter_name.c_str()) << n_inverse << endl;

    for (int i = 0; i < num_bins_e; i++)
    {
        //PTBC
        Double_t trans_bin_content_PTBC = transmission_hist_e_PTBC->GetBinContent(i+1);
        Double_t trans_bin_error_PTBC = transmission_hist_e_PTBC->GetBinError(i+1);
        if (trans_bin_content_PTBC == 0)
        {
            cross_section_hist_e_PTBC->SetBinContent(i+1, 0);
            cross_section_hist_e_PTBC->SetBinError(i+1, 0);
        } else {
            Double_t cross_section_PTBC = - n_inverse * std::log(trans_bin_content_PTBC);
            Double_t bin_unc_PTBC = n_inverse * (1./trans_bin_content_PTBC) * trans_bin_error_PTBC;
            cross_section_hist_e_PTBC->SetBinContent(i+1, cross_section_PTBC);
            cross_section_hist_e_PTBC->SetBinError(i+1, bin_unc_PTBC);
        }

        //FIMG
        Double_t trans_bin_content_FIMG = transmission_hist_e_FIMG->GetBinContent(i+1);
        Double_t trans_bin_error_FIMG = transmission_hist_e_FIMG->GetBinError(i+1);
        if (trans_bin_content_FIMG == 0)
        {
            cross_section_hist_e_FIMG->SetBinContent(i+1, 0);
            cross_section_hist_e_FIMG->SetBinError(i+1, 0);
        } else {
            Double_t cross_section_FIMG = - n_inverse * std::log(trans_bin_content_FIMG);
            Double_t bin_unc_FIMG = n_inverse * (1./trans_bin_content_FIMG) * trans_bin_error_FIMG;
            cross_section_hist_e_FIMG->SetBinContent(i+1, cross_section_FIMG);
            cross_section_hist_e_FIMG->SetBinError(i+1, bin_unc_FIMG);
        }
    }

    exclude_first_last_bins(transmission_hist_e_PTBC);
    exclude_first_last_bins(transmission_hist_e_FIMG);
    exclude_first_last_bins(cross_section_hist_e_PTBC);
    exclude_first_last_bins(cross_section_hist_e_FIMG);

    cout << Form("Total Protons %s = ", filter_name.c_str()) << norm_factors[0] << endl;
    cout << "Total Protons Target Out = " << norm_factors[1] << endl;
}

void plot_hists(){
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

    l[i] = new TLegend(0.77,0.7,0.86,0.85); //0.68,0.7,0.86,0.8       ;         0.72,0.8,0.90,0.9
    l[i]->AddEntry(transmission_hist_e_PTBC,"PTBC","l");
    
    transmission_hist_e_PTBC->GetXaxis()->SetTitle("Energy (in eV)");
    transmission_hist_e_PTBC->GetYaxis()->SetTitle("Transmission");
    transmission_hist_e_PTBC->SetTitle(Form("Transmission Histogram - %s", filter_name_title.c_str()));
    transmission_hist_e_PTBC->SetLineWidth(1);
    transmission_hist_e_PTBC->GetXaxis()->SetRangeUser(1e-1,1e8);
    transmission_hist_e_PTBC->Draw(); //"HISTE"
    transmission_hist_e_PTBC->SetStats(0);
    // transmission_hist_e_PTBC->SetMarkerStyle(6);
    // transmission_hist_e_PTBC->SetMarkerSize(0.5);
    // gPad->SetGrid();
    gPad->SetLogx();
    // gStyle->SetPalette(57);

    l[i]->AddEntry(transmission_hist_e_FIMG,"FIMG","l");
    transmission_hist_e_FIMG->SetLineColor(6);
    transmission_hist_e_FIMG->SetLineWidth(1);
    transmission_hist_e_FIMG->GetXaxis()->SetRangeUser(1e-1,1e6);
    transmission_hist_e_FIMG->Draw("SAME");

    if (fillENDF){
        l[i]->AddEntry(endf_trans_hist,"ENDF","l");
        endf_trans_hist->SetLineColor(2);
        endf_trans_hist->SetLineWidth(2);
        endf_trans_hist->GetXaxis()->SetRangeUser(1e-1,1e8);
        endf_trans_hist->Draw("SAME");
    }

    if (fillJENDL)
    {
        l[i]->AddEntry(jendl_trans_hist,"JENDL-5","l");
        jendl_trans_hist->SetLineColor(1);
        jendl_trans_hist->SetLineWidth(2);
        // jendl_trans_hist->GetXaxis()->SetRangeUser(1e-2,2e7);
        jendl_trans_hist->Draw("SAME");
    }

    // auto endf_xsec_graph = new TGraph();
    // auto jendl_xsec_graph = new TGraph();

    // if (fillENDF){
    //     auto endf_trans_graph = new TGraph();
    //     for (Int_t i = 0; i < endf_e.size(); i++)
    //     {
    //         endf_trans_graph->SetPoint(i+1, endf_e.at(i), endf_trans.at(i));
    //     }
    //     l[i]->AddEntry(endf_trans_graph,"ENDF","l");
    //     endf_trans_graph->SetLineColor(2);
    //     endf_trans_graph->SetLineWidth(1);
    //     // endf_trans_graph->GetXaxis()->SetRangeUser(1e-2,2e7);
    //     endf_trans_graph->Draw("SAME");
    // }

    // if (fillJENDL)
    // {
    //     auto jendl_trans_graph = new TGraph();
    //     for (Int_t i = 0; i < jendl_e.size(); i++)
    //     {
    //         jendl_trans_graph->SetPoint(i+1, jendl_e.at(i), jendl_trans.at(i));
    //     }
    //     l[i]->AddEntry(jendl_trans_graph,"JENDL-5","l");
    //     jendl_trans_graph->SetLineColor(1);
    //     jendl_trans_graph->SetLineWidth(1);
    //     // jendl_trans_graph->GetXaxis()->SetRangeUser(1e-2,2e7);
    //     jendl_trans_graph->Draw("SAME");
    // }
    

    l[i]->SetMargin(0.4);
    l[i]->Draw();
    // c[i]->Print(Form("../plots/h_trans_e_%s_%iBPD.png", filter_name.c_str(), bins_per_decade));

    i++;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();

    l[i] = new TLegend(0.77,0.7,0.86,0.85);
    l[i]->AddEntry(cross_section_hist_e_PTBC,"PTBC","l");

    cross_section_hist_e_PTBC->GetXaxis()->SetTitle("Energy (in eV)");
    cross_section_hist_e_PTBC->GetYaxis()->SetTitle("Cross Section (in barns)");
    cross_section_hist_e_PTBC->SetTitle(Form("Cross Section Histogram - %s", filter_name_title.c_str()));
    cross_section_hist_e_PTBC->SetLineWidth(1);
    cross_section_hist_e_PTBC->GetXaxis()->SetRangeUser(1e-1,1e8);
    cross_section_hist_e_PTBC->Draw(); //"HISTE"
    cross_section_hist_e_PTBC->SetStats(0);
    // cross_section_hist_e_PTBC->SetMarkerStyle(6);
    // cross_section_hist_e_PTBC->SetMarkerSize(0.5);
    // gPad->SetGrid();
    gPad->SetLogx();
    // gStyle->SetPalette(57);

    l[i]->AddEntry(cross_section_hist_e_FIMG,"FIMG","l");
    cross_section_hist_e_FIMG->SetLineColor(6);
    cross_section_hist_e_FIMG->SetLineWidth(1);
    cross_section_hist_e_FIMG->GetXaxis()->SetRangeUser(1e-1,1e6);
    cross_section_hist_e_FIMG->Draw("SAME");

    if (fillENDF){
        l[i]->AddEntry(endf_xsec_hist,"ENDF","l");
        endf_xsec_hist->SetLineColor(2);
        endf_xsec_hist->SetLineWidth(2);
        endf_xsec_hist->GetXaxis()->SetRangeUser(1e-1,1e8);
        endf_xsec_hist->Draw("SAME");
    }

    if (fillJENDL)
    {
        l[i]->AddEntry(jendl_xsec_hist,"JENDL-5","l");
        jendl_xsec_hist->SetLineColor(1);
        jendl_xsec_hist->SetLineWidth(2);
        // jendl_xsec_hist->GetXaxis()->SetRange(1e-2,2e7);
        jendl_xsec_hist->Draw("SAME");
    }

    l[i]->SetMargin(0.4);
    l[i]->Draw();
    // c[i]->Print(Form("../plots/h_xsec_e_%s_%iBPD.png", filter_name.c_str(), bins_per_decade));
} 

void store_hists(){
    
    TFile *f = new TFile(Form("../rootFiles/trans_xsec_hists_%s_%ibpd.root", filter_name.c_str(), bins_per_decade),"recreate");

    transmission_hist_e_PTBC->Write();
    transmission_hist_e_FIMG->Write();
    cross_section_hist_e_PTBC->Write();
    cross_section_hist_e_FIMG->Write();

    if (fillENDF)
    {
        endf_trans_hist->Write();
        endf_xsec_hist->Write();
    }

    if (fillJENDL)
    {
        jendl_trans_hist->Write();
        jendl_xsec_hist->Write();
    }
    
    f->Close();

    std::cout << Form("Created output file 'trans_xsec_hists_%s_%ibpd.root'", filter_name.c_str(), bins_per_decade) << std::endl;
}

void crossSectionPlots(){

    fillNumDensityMap();
    fillEValFileNameMap();

    // //Getting energy bin edges
    // Int_t num_bins_e = transmission_hist_e_PTBC->GetNbinsX();
    // Double_t bin_edges_e[num_bins_e + 1];
    // for (int i = 0; i < num_bins_e; i++)
    // {
    //     bin_edges_e[i] = transmission_hist_e_PTBC->GetXaxis()->GetBinLowEdge(i+1);

    //     if (i == num_bins_e - 1)
    //     {
    //         bin_edges_e[i+1] = transmission_hist_e_PTBC->GetXaxis()->GetBinUpEdge(i+1);
    //     }
    // }

    //Calculating Energy bin edges
    Double_t tof_min = 1e3; //ns
    Double_t tof_max = 1e8; //ns
    Double_t e_min = TOFToEnergy(tof_max * 1e-9, flight_path_length_FIMG); //converting into seconds
    Double_t e_max = TOFToEnergy(tof_min * 1e-9, flight_path_length_FIMG); //converting into seconds
    Int_t min_power = FindDecadePower(e_min);
    Int_t max_power = FindDecadePower(e_max);
    Int_t num_decades_e = max_power - min_power;
    Int_t num_bins_e = bins_per_decade * num_decades_e;
    Double_t bin_edges_e[num_bins_e+1];
    Double_t step_e = ((Double_t) 1.0/(Double_t) bins_per_decade);
    for(Int_t i = 0; i < num_bins_e+1; i++)
    {
        Double_t base = 10.;
        Double_t exponent = (step_e * (Double_t) i) + (Double_t) min_power;
        bin_edges_e[i] = (Double_t) std::pow(base, exponent);
    }

    cout << "Number of e bins = " << num_bins_e << endl;

    calc_xsec(Form("../rootFiles/crossSectionAna_%s.root", filter_name.c_str()), num_bins_e, bin_edges_e);

    Double_t num_density = 0.;
    if (fillENDF || fillJENDL)
    {
        num_density = num_density_map[filter_name];
    }

    if (!filter_name.compare("ar_bottle_full"))
    {

        //ENDF Hists
        if (fillENDF){
            // endf_trans_hist = new TH1D("endf_trans_hist","ENDF Transmission Hist",num_bins_e,bin_edges_e);
            // endf_xsec_hist = new TH1D("endf_xsec_hist","ENDF Cross Section Hist",num_bins_e,bin_edges_e);
            endf_argon(num_density, bins_per_decade);
        }

        if(fillJENDL){
            // jendl_trans_hist = new TH1D("jendl_trans_hist","JENDL Transmission Hist",num_bins_e,bin_edges_e);
            // jendl_xsec_hist = new TH1D("jendl_xsec_hist","JENDL Cross Section Hist",num_bins_e,bin_edges_e);  
            jendl_argon(num_density, bins_per_decade);
        }
        
    } else {

        //ENDF Hists
        if (fillENDF){
            endf_trans_hist = new TH1D("endf_trans_hist","ENDF Transmission Hist",num_bins_e,bin_edges_e);
            endf_xsec_hist = new TH1D("endf_xsec_hist","ENDF Cross Section Hist",num_bins_e,bin_edges_e);
            endf(num_density, num_bins_e, bin_edges_e[0], bin_edges_e[num_bins_e]);
        }  

    }

    if (storeHists)
    {
        store_hists();
    }
    
    if (plotHists)
    {
        plot_hists();
    }
}