/**
 * @file transmissionAna.h
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

// std::string mode = "run"; //"run", "test"

TH1D* tof_hist_filter_al5 = 0;
TH1D* energy_hist_filter_al5 = 0;
TH1D* tof_hist_filter_al8 = 0;
TH1D* energy_hist_filter_al8 = 0;
// TH1D* normalized_events_hist_e = 0;
// TH1D* cross_section_hist_e = 0;

// auto endf_trans = new TGraph();
// auto endf_xsec = new TGraph();
// auto endf_rf_trans = new TGraph();
// auto endf_rf_xsec = new TGraph();

TH1D* endf_trans_hist = 0;
TH1D* endf_xsec_hist = 0;
TH1D* endf_rf_trans_hist = 0;
TH1D* endf_rf_xsec_hist = 0;

Int_t bins_per_decade = 100;
Double_t flight_path_length_PTB = 182.65 - 0.41; //m
Double_t neutron_mass = 939.56542052; //in MeV
Double_t speed_of_light = 299792458.0; //in m/s

Double_t n_Bi_1cm = (1.0 /*cm*/) * (9.78 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (208.9804 /*g/mole*/);
Double_t n_Al_3cm = (3.0 /*cm*/) * (2.70 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (26.9815 /*g/mole*/);
Double_t n_Al_5cm = (5.0 /*cm*/) * (2.70 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (26.9815 /*g/mole*/);
Double_t n_Al_8cm = (8.0 /*cm*/) * (2.70 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (26.9815 /*g/mole*/);

Double_t t_gamma_PTB = (flight_path_length_PTB / speed_of_light) * 1e9; //converting into ns

Double_t t11 = 780.0;
Double_t a11 = 7700.0;

Double_t t12 = 1090.0;
Double_t a12 = 5451.0;

Double_t t21 = 1090.0;
Double_t a21 = 5451.0;

Double_t t22 = 2605.0;
Double_t a22 = 7167.0;

Double_t t31 = 2605.0;
Double_t a31 = 7960.0;

Double_t t32 = 2856.0;
Double_t a32 = 7960.0;

Double_t t41 = 2856.0;
Double_t a41 = 7432.0;

Double_t t42 = 15290.0;
Double_t a42 = 7432.0;

Double_t t51 = 15290.0;
Double_t a51 = 7432.0;

Double_t t52 = 18708.0;
Double_t a52 = 3600.0; //2416

Double_t t61 = 18708.0;
Double_t a61 = 3600.0; //2416

Double_t t62 = 1e8;
Double_t a62 = 3600.0; //3076

//Filter In runs
//Bi (1 cm) Filter
std::vector<Int_t> c_au_f_bi_t_out = {117369, 117370, 117371, 117372, 117373, 117374, 117375, 117376, 117377, 117378, 117379, 117380, 117381, 117382, 117383, 117384, 117385, 117386, 117387, 117388, 117389, 117390, 117436};
std::vector<Int_t> c_out_f_bi_t_out = {117391, 117392, 117393, 117394, 117395, 117396, 117397};
std::vector<Int_t> c_pb_f_bi_t_out = {117422, 117423, 117424, 117425, 117426, 117427, 117428};
std::vector<Int_t> c_ta_f_bi_t_out = {117437, 117439, 117440, 117441, 117442, 117443};

//Al (8 cm) Filter
std::vector<Int_t> c_ta_f_al8_t_out = {117449, 117450, 117451, 117453, 117454, 117455, 117456, 117457, 117458, 117459, 117460};

//Al (5 cm) Filter
std::vector<Int_t> c_ta_f_al5_t_out = {117461, 117462, 117463, 117464, 117465, 117466, 117467, 117476, 117477, 117478, 117479, 117480, 117481, 117482}; //, 117483, 117484
std::vector<Int_t> c_out_f_al5_t_out = {117470, 117471, 117472, 117473, 117474, 117475};

//Filter Out runs
std::vector<Int_t> c_au_f_out_t_out = {117357, 117358, 117359, 117362, 117363, 117364, 117365, 117366, 117367, 117368};
std::vector<Int_t> c_out_f_out_t_out = {117405, 117406, 117408, 117409, 117410, 117411, 117412};
std::vector<Int_t> c_pb_f_out_t_out = {117429, 117430, 117431, 117432, 117433, 117429, 117430, 117431, 117432, 117433, 117434, 117435};
std::vector<Int_t> c_ta_f_out_t_out = {117444, 117445, 117446, 117447, 117448, 117452};

std::vector<Int_t> filter_al5_runs;
std::vector<Int_t> filter_al8_runs;

Double_t yOnTheCutLine(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t x3){
    return ((y2 - y1)*(x3 - x1)/(x2 - x1) + y1);
}

Double_t EnergyToTOF(Double_t e){ //e is in eV
    Double_t KE_M = (e * 1e-6)/neutron_mass;
    Double_t denominator = 1.0 - 1.0/((KE_M + 1.0)*(KE_M + 1.0));
    Double_t correction_factor = std::sqrt(1.0 / denominator);
    Double_t TOF = (flight_path_length_PTB/speed_of_light) * correction_factor;
    return TOF; //tof in seconds
}
    
Double_t TOFToEnergy(Double_t t){ //t is in seconds
    Double_t denom_term = (flight_path_length_PTB)/(speed_of_light * t);
    Double_t denominator = 1.0 - (denom_term * denom_term);
    Double_t factor = std::sqrt(1.0 / denominator) - 1.0;
    Double_t energy = neutron_mass * factor;
    return energy * 1e6;  //e in eV
}

Double_t TOFToEnergy(Double_t t, Double_t rf_length){ //t is in seconds, rf_length is in m
    Double_t denom_term = (flight_path_length_PTB + rf_length)/(speed_of_light * t);
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

Double_t FindFWHM(TH1D* projection_hist){
    
    Int_t max_bin_num = projection_hist->GetMaximumBin();
    Double_t max_value = projection_hist->GetBinContent(max_bin_num);
    Int_t tot_bins = projection_hist->GetNbinsX();
    Double_t left_edge = 0;
    Double_t right_edge = 0;
    //Find left edge
    for (Int_t i = max_bin_num-1; i > 0; i--)
    {
        Double_t bin_value = projection_hist->GetBinContent(i);
        if (bin_value > max_value/2.0)
        {
            continue;
        } else {
            left_edge = projection_hist->GetXaxis()->GetBinUpEdge(i);
            break;
        }
    }

    //Find right edge
    for (Int_t i = max_bin_num+1; i < tot_bins+1; i++)
    {
        Double_t bin_value = projection_hist->GetBinContent(i);
        if (bin_value > max_value/2.0)
        {
            continue;
        } else {
            right_edge = projection_hist->GetXaxis()->GetBinLowEdge(i);
            break;
        }
    }
    return (right_edge - left_edge);
}

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