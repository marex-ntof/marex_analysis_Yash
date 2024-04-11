/**
 * @file crossSectionPlots.h
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
#include "TKey.h"
#include "TList.h"

#include "MArEXStyle.C"

//////// Run variables
// std::string mode = "run"; //"run", "test"
const std::string filter_name("al5"); //bi1, al3, al5, al8, c1p2_ts, al5_ts, bi1p2_ts, cf_bottle, cf_bottle_rot, ar_bottle_full
const std::string filter_name_title("Al (5 cm)");
Int_t bins_per_decade = 20;
//Bi (1 cm), Target Bi (1.2 cm), Al (3 cm), Al (5 cm), Target Al (5 cm), Al (8 cm), Target C (1.2 cm), Empty Bottle, Empty Bottle Rotated
//Argon Tank

// const std::string root_file_name("rootFiles/crossSectionAna_al_8cm.root"); //Al_tot_xsec.txt, Bi_tot_xsec.txt
// const std::string endf_file_name("evalData/Bi_tot_xsec.txt"); //Al_tot_xsec.txt, Bi_tot_xsec.txt

// TH1D* tof_hist_filter_in = 0;
// TH1D* energy_hist_filter_in = 0;
// TH1D* tof_hist_filter_out = 0;
// TH1D* energy_hist_filter_out = 0;

TH1D* transmission_hist_e_PTBC = 0;
TH1D* transmission_hist_e_FIMG = 0;
TH1D* cross_section_hist_e_PTBC = 0;
TH1D* cross_section_hist_e_FIMG = 0;

TH1D* energy_hist_target_in_PTBC = 0;
TH1D* energy_hist_target_in_FIMG = 0;
TH1D* energy_hist_target_out_PTBC = 0;
TH1D* energy_hist_target_out_FIMG = 0;

TH1D* transmission_hist_e_PTBC_nTOF_Cuts = 0;

std::vector<Double_t> norm_factors;

// Double_t norm_factor_target_in = 0;
// Double_t norm_factor_target_out = 0;

// TH1D* trans_hist_fOut = 0;
// TH1D* trans_hist_fIn = 0;
// TH1D* trans_hist_fOut_endf = 0;
// TH1D* trans_hist_fIn_endf = 0;

// auto endf_trans = new TGraph();
// auto endf_xsec = new TGraph();
// auto endf_rf_trans = new TGraph();
// auto endf_rf_xsec = new TGraph();

TH2D* rf_hist = 0;
TH1D* endf_trans_hist = 0;
TH1D* endf_xsec_hist = 0;
TH1D* jendl_trans_hist = 0;
TH1D* jendl_xsec_hist = 0;
TH1D* endf_rf_trans_hist = 0;
TH1D* endf_rf_xsec_hist = 0;

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
Double_t n_Ar_bottle = (11.0 /*cm*/) * ((ar_bottle_pressure)/(8.31446261815324 * ar_bottle_temp * 1e6)) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/);

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
    eval_file_name_map.emplace("ar_bottle_full", "Ar_tot_xsec.txt");
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
