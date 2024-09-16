/**
 * @file ptbcCutsAna.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2024-09-06
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

//////// Run variables
// std::string mode = "run"; //"run", "test"
const std::string filter_name("al5"); //bi1, al3, al5, al8, c1p2_ts, al5_ts, al5_c_ts, bi1p2_ts, cf_bottle, cf_bottle_rot, ar_bottle_full
const std::string filter_name_title("Al (5 cm)");
Int_t bins_per_decade = 50;
//Bi (1 cm), Target Bi (1.2 cm), Al (3 cm), Al (5 cm), Target Al (5 cm), Al (8 cm), Target C (1.2 cm), Empty Tank, Empty Tank Rotated
//Argon Tank, Al (5 cm) & C (1.2 cm)

bool fillENDF = true;
bool fillJENDL = false;

TH1D* energy_hist_target_in[5];
TH1D* energy_hist_target_out[5];

TH1D* transmission_hist[5];

std::vector<Double_t> norm_factors;

TH1D* endf_trans_hist = 0;
TH1D* endf_xsec_hist = 0;

std::vector<Double_t> endf_e;
std::vector<Double_t> endf_xsec;
std::vector<Double_t> endf_trans;

std::vector<Double_t> jendl_e;
std::vector<Double_t> jendl_xsec;
std::vector<Double_t> jendl_trans;

Int_t evalFile_max_line_num = 0;
const Double_t flight_path_length_PTBC = 182.65 - 0.41; //m
const Double_t flight_path_length_FIMG = 183.5 - 0.41; //m
const Double_t neutron_mass = 939.56542052; //in MeV
const Double_t speed_of_light = 299792458.0; //in m/s
Double_t ar_bottle_pressure = 197.385 * 1e5; // in Pa (SI unit)
Double_t ar_bottle_temp = 293.0; // in Kelvin

Double_t n_Bi_1cm = (1.0 /*cm*/) * (9.78 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (208.9804 /*g/mole*/);
Double_t n_Bi_1p2cm = (1.2 /*cm*/) * (9.78 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (208.9804 /*g/mole*/);
Double_t n_Al_3cm = (3.0 /*cm*/) * (2.70 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (26.9815 /*g/mole*/);
Double_t n_Al_5cm = (5.0 /*cm*/) * (2.70 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (26.9815 /*g/mole*/);
Double_t n_Al_8cm = (8.0 /*cm*/) * (2.70 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (26.9815 /*g/mole*/);
Double_t n_C_1p2cm = 0.105;// (1.2 /*cm*/) * (2.267 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (12.011 /*g/mole*/);
Double_t n_CFib_1cm = (1.0 /*cm*/) * (1.8 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (12.011 /*g/mole*/);
Double_t n_Ar_bottle = (8.9 /*cm*/) * ((ar_bottle_pressure)/(8.31446261815324 * ar_bottle_temp * 1e6)) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/);

////////////////////////
std::map<std::string, Double_t> num_density_map;
void fillNumDensityMap(){
    num_density_map.emplace("bi1", n_Bi_1cm);
    num_density_map.emplace("bi1p2_ts", n_Bi_1p2cm);
    num_density_map.emplace("al3", n_Al_3cm);
    num_density_map.emplace("al5", n_Al_5cm);
    num_density_map.emplace("al5_ts", n_Al_5cm);
    num_density_map.emplace("al8", n_Al_8cm);
    num_density_map.emplace("c1p2_ts", n_C_1p2cm);
    num_density_map.emplace("cf_bottle", n_CFib_1cm);
    num_density_map.emplace("cf_bottle_rot", n_CFib_1cm);
    num_density_map.emplace("cf_bottle_rotBack", n_CFib_1cm);
    num_density_map.emplace("ar_bottle_full", n_Ar_bottle);
}

std::map<std::string, std::string> eval_file_name_map;
void fillEValFileNameMap(){

    eval_file_name_map.emplace("bi1", "Bi_tot_xsec.txt");
    eval_file_name_map.emplace("bi1p2_ts", "Bi_tot_xsec.txt");
    eval_file_name_map.emplace("al3", "Al_tot_xsec.txt");
    eval_file_name_map.emplace("al5", "Al_tot_xsec.txt");
    eval_file_name_map.emplace("al5_ts", "Al_tot_xsec.txt");
    eval_file_name_map.emplace("al8", "Al_tot_xsec.txt");
    eval_file_name_map.emplace("c1p2_ts", "C_tot_xsec.txt");
    eval_file_name_map.emplace("cf_bottle", "C_tot_xsec.txt");
    eval_file_name_map.emplace("cf_bottle_rot", "C_tot_xsec.txt");
    eval_file_name_map.emplace("cf_bottle_rotBack", "C_tot_xsec.txt");
    eval_file_name_map.emplace("ar_bottle_full", "Ar40_tot_xsec.txt");
}
////////////////

void fillMaxLineCount(){

    if (!filter_name.compare("bi1") || !filter_name.compare("bi1p2_ts")){
        evalFile_max_line_num = 27019;
    }

    if (!filter_name.compare("al3") || !filter_name.compare("al5") || !filter_name.compare("al5_ts") || !filter_name.compare("al8")){
        evalFile_max_line_num = 9745;
    }

    if (!filter_name.compare("c1p2_ts") || !filter_name.compare("cf_bottle") || !filter_name.compare("cf_bottle_rot") || !filter_name.compare("cf_bottle_rotBack")){
        evalFile_max_line_num = 1430;
    }

    if (!filter_name.compare("ar_bottle_full")){
        evalFile_max_line_num = 51037;
    }

}

Double_t EnergyToTOF(Double_t e, Double_t flight_path_length){ //e is in eV, flight_path_length is in m
    Double_t KE_M = (e * 1e-6)/neutron_mass;
    Double_t denominator = 1.0 - 1.0/((KE_M + 1.0)*(KE_M + 1.0));
    Double_t correction_factor = std::sqrt(1.0 / denominator);
    Double_t TOF = (flight_path_length/speed_of_light) * correction_factor;
    return TOF; //tof in seconds
}
    
Double_t TOFToEnergy(Double_t t, Double_t flight_path_length){ //t is in seconds
    Double_t denom_term = (flight_path_length)/(speed_of_light * t);
    Double_t denominator = 1.0 - (denom_term * denom_term);
    Double_t factor = std::sqrt(1.0 / denominator) - 1.0;
    Double_t energy = neutron_mass * factor;
    return energy * 1e6;  //e in eV
}

Double_t TOFToEnergy(Double_t t, Double_t flight_path_length, Double_t rf_length){ //t is in seconds, rf_length is in m
    Double_t denom_term = (flight_path_length + rf_length)/(speed_of_light * t);
    Double_t denominator = 1.0 - (denom_term * denom_term);
    Double_t factor = std::sqrt(1.0 / denominator) - 1.0;
    Double_t energy = neutron_mass * factor;
    return energy * 1e6;  //e in eV
}

Double_t EnergyToVelocity(Double_t e){ //e is in eV
    Double_t KE_M = (e * 1e-6)/neutron_mass;
    Double_t factor = std::sqrt(1.0 - 1.0/((KE_M + 1.0)*(KE_M + 1.0)));
    return speed_of_light * factor; //velocity in m/s
}

// Double_t xsecToTransmission(Double_t n, Double_t xsec){
//     return std::exp(- n * xsec);
// }

Int_t FindDecadePower(Double_t num){
    Int_t decadePower = 0;
    Double_t value = num;
    if (num > 1)
    {
        while (value > 1)
        {
            value /= 10;
            decadePower++;
        }   
    }
    if (num < 1)
    {
        while (value < 1)
        {
            value *= 10;
            decadePower--;
        } 
    }
    return decadePower;
}

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

void fill_trans(TH1D* target_in, TH1D* target_out, TH1D* trans_hist, Int_t num_bins_e, Double_t norm_factor){

    for (int i = 0; i < num_bins_e; i++)
    {
        //PTBC
        Double_t bin_content_in_PTBC = target_in->GetBinContent(i+1);
        Double_t bin_content_out_PTBC = target_out->GetBinContent(i+1);
        if (bin_content_in_PTBC == 0. || bin_content_out_PTBC == 0.)
        {
            trans_hist->SetBinContent(i+1, 0.);
            trans_hist->SetBinError(i+1, 0.);
        } else {
            Double_t transmission_PTBC = (bin_content_in_PTBC * norm_factor)/bin_content_out_PTBC;
            Double_t bin_unc_PTBC = transmission_PTBC * std::sqrt( (1./bin_content_in_PTBC) + (1./bin_content_out_PTBC) );
            trans_hist->SetBinContent(i+1, transmission_PTBC);
            trans_hist->SetBinError(i+1, bin_unc_PTBC);
        }
    }
}

void calc_trans(const char *fname, Int_t num_bins_e, Double_t bin_edges_e[]){
    cout << "Opening File " << fname << endl;
    TFile *hist_file = TFile::Open(fname, "READ");

    cout << "Extracting Histograms " << endl;
    energy_hist_target_in[0] = (TH1D*)hist_file->Get("energy_hist_target_in_thecut");
    energy_hist_target_in[1] = (TH1D*)hist_file->Get("energy_hist_target_in_20more");
    energy_hist_target_in[2] = (TH1D*)hist_file->Get("energy_hist_target_in_50more");
    energy_hist_target_in[3] = (TH1D*)hist_file->Get("energy_hist_target_in_20less");
    energy_hist_target_in[4] = (TH1D*)hist_file->Get("energy_hist_target_in_50less");

    energy_hist_target_out[0] = (TH1D*)hist_file->Get("energy_hist_target_out_thecut");
    energy_hist_target_out[1] = (TH1D*)hist_file->Get("energy_hist_target_out_20more");
    energy_hist_target_out[2] = (TH1D*)hist_file->Get("energy_hist_target_out_50more");
    energy_hist_target_out[3] = (TH1D*)hist_file->Get("energy_hist_target_out_20less");
    energy_hist_target_out[4] = (TH1D*)hist_file->Get("energy_hist_target_out_50less");

    cout << "Extracting Norm Factors " << endl;
    std::vector<Double_t> *norm_factors_temp;
    hist_file->GetObject("norm_factors", norm_factors_temp);
    norm_factors = *norm_factors_temp;
    Double_t Qout_Qin = norm_factors[1]/norm_factors[0]; //norm_factor_target_out / norm_factor_target_in;

    //transmission histogram
    transmission_hist[0] = new TH1D("transmission_hist_thecut","Transmission Hist - The Cut",num_bins_e,bin_edges_e);
    transmission_hist[1] = new TH1D("transmission_hist_20more","Transmission Hist - 20%% More",num_bins_e,bin_edges_e);
    transmission_hist[2] = new TH1D("transmission_hist_50more","Transmission Hist - 50%% More",num_bins_e,bin_edges_e);
    transmission_hist[3] = new TH1D("transmission_hist_20less","Transmission Hist - 20%% Less",num_bins_e,bin_edges_e);
    transmission_hist[4] = new TH1D("transmission_hist_50less","Transmission Hist - 50%% Less",num_bins_e,bin_edges_e);
    
    //Transmission
    for (Int_t i = 0; i < 5; i++)
    {
        fill_trans(energy_hist_target_in[i], energy_hist_target_out[i], transmission_hist[i], num_bins_e, Qout_Qin);
    }

    for (Int_t i = 0; i < 5; i++)
    {
        exclude_first_last_bins(transmission_hist[i]);
    }
    
}

void ptbcCutsAna(){

    fillNumDensityMap();
    fillEValFileNameMap();

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

    calc_trans(Form("../rootFiles/PTBC_cuts_ana_%s.root", filter_name.c_str()), num_bins_e, bin_edges_e);

    Double_t num_density = 0.;
    if (fillENDF || fillJENDL)
    {
        num_density = num_density_map[filter_name];
        endf_trans_hist = new TH1D("endf_trans_hist","ENDF Transmission Hist",num_bins_e,bin_edges_e);
        endf_xsec_hist = new TH1D("endf_xsec_hist","ENDF Cross Section Hist",num_bins_e,bin_edges_e);
        endf(num_density, num_bins_e, bin_edges_e[0], bin_edges_e[num_bins_e]);
    }

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

    l[i] = new TLegend(0.77,0.7,0.86,0.85); //0.68,0.7,0.86,0.8

    transmission_hist[0]->GetXaxis()->SetTitle("Energy (in eV)");
    transmission_hist[0]->GetYaxis()->SetTitle("Transmission");
    transmission_hist[0]->SetTitle(Form("Transmission Histogram - %s", filter_name_title.c_str()));
    // transmission_hist[0]->SetLineWidth(2);
    transmission_hist[0]->SetLineColor(2);
    transmission_hist[0]->Draw(); //"HISTE"
    transmission_hist[0]->SetStats(0);
    gPad->SetLogx();
    l[i]->AddEntry(transmission_hist[0],"Original Cut","l");

    transmission_hist[1]->SetLineColor(3);
    transmission_hist[1]->Draw("SAME");
    l[i]->AddEntry(transmission_hist[1],"Cut 20%% More","l");

    transmission_hist[2]->SetLineColor(4);
    transmission_hist[2]->Draw("SAME");
    l[i]->AddEntry(transmission_hist[2],"Cut 50%% More","l");

    transmission_hist[3]->SetLineColor(6);
    transmission_hist[3]->Draw("SAME");
    l[i]->AddEntry(transmission_hist[3],"Cut 20%% Less","l");

    transmission_hist[4]->SetLineColor(7);
    transmission_hist[4]->Draw("SAME");
    l[i]->AddEntry(transmission_hist[4],"Cut 50%% Less","l");

    if (fillENDF){
        l[i]->AddEntry(endf_trans_hist,"ENDF","l");
        endf_trans_hist->SetLineColor(1);
        endf_trans_hist->SetLineWidth(2);
        // endf_trans_hist->GetXaxis()->SetRange(1e-2,2e7);
        endf_trans_hist->Draw("SAME");
    }

    l[i]->SetMargin(0.4);
    l[i]->Draw();

    i++;

    
}  
