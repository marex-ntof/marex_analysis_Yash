/**
 * @file pkup.h
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

//////// Run variables
bool fill_filterOut_runs = true;
bool fill_bi1_runs = false;
bool fill_al3_runs = false;
bool fill_al5_runs = false;
bool fill_al8_runs = false;

TH2D* tpkup_beam_intensity_hist = 0;
TH2D* delT_PTBC_beam_intensity_hist = 0;
TH2D* delT_FIMG_beam_intensity_hist = 0;

Int_t bins_per_decade = 100;
Double_t flight_path_length_PTB = 182.65 - 0.41; //m
Double_t flight_path_length_FIMG = 183.65 - 0.41; //m
Double_t neutron_mass = 939.56542052; //in MeV
Double_t speed_of_light = 299792458.0; //in m/s

///////// FIMG stats from 117386
//Filter In runs
std::vector<Int_t> list_of_runs;
//Bi (1 cm) Filter
std::vector<Int_t> c_au_f_bi_t_out = {117369, 117370, 117371, 117372, 117373, 117374, 117375, 117376, 117377, 117378, 117379, 117380, 117381, 117382, 117383, 117384, 117385, 117386, 117387, 117388, 117389, 117390, 117436};
std::vector<Int_t> c_out_f_bi_t_out = {117391, 117392, 117393, 117394, 117395, 117396, 117397};
std::vector<Int_t> c_pb_f_bi_t_out = {117422, 117423, 117424, 117425, 117426, 117427, 117428};
std::vector<Int_t> c_ta_f_bi_t_out = {117437, 117439, 117440, 117441, 117442, 117443};

//Al (8 cm) Filter
std::vector<Int_t> c_ta_f_al8_t_out = {117449, 117450, 117451, 117453, 117454, 117455, 117456, 117457, 117458, 117459, 117460};

//Al (5 cm) Filter
std::vector<Int_t> c_ta_f_al5_t_out = {117461, 117462, 117463, 117464, 117465, 117466, 117467, 117476, 117477, 117478, 117479, 117480, 117481, 117482, 117483, 117484};
std::vector<Int_t> c_out_f_al5_t_out = {117470, 117471, 117472, 117473, 117474, 117475};
std::vector<Int_t> c_pb_f_al5_t_out = {117497, 117498, 117499, 117500, 117501, 117502};
std::vector<Int_t> c_c_f_al5_t_out = {117503, 117504, 117505, 117506, 117507, 117508, 117509, 117510};

//Al (3 cm) Filter
std::vector<Int_t> c_ta_f_al3_t_out = {117485, 117486, 117487, 117488, 117489, 117490, 117491, 117492, 117493};
std::vector<Int_t> c_pb_f_al3_t_out = {117496};

//Filter Out runs
// std::vector<Int_t> c_au_f_out_t_out = {117350, 117357, 117358, 117359, 117362, 117363, 117364, 117365, 117366, 117367, 117368};
// std::vector<Int_t> c_atic_f_out_t_out = {117351, 117355, 117356}; // 117352, 117353,
std::vector<Int_t> c_out_f_out_t_out = {117398, 117405, 117406, 117408, 117409, 117410, 117411, 117412, 117398}; //
std::vector<Int_t> c_pb_f_out_t_out = {117429, 117430, 117431, 117432, 117433, 117434, 117435};
std::vector<Int_t> c_ta_f_out_t_out = {117444, 117445, 117446, 117447, 117448}; //117452
std::vector<Int_t> c_c_f_out_t_out = {117511, 117512, 117513, 117514, 117515, 117516, 117517, 117518};

void fillRuns(){

    // Int_t total_runs = 0;
    // total_runs = c_au_f_out_t_out.size() + c_out_f_out_t_out.size() + c_pb_f_out_t_out.size() + c_ta_f_out_t_out.size();
    // total_runs += c_au_f_bi_t_out.size() + c_out_f_bi_t_out.size() + c_pb_f_bi_t_out.size() + c_ta_f_bi_t_out.size();
    // total_runs += c_ta_f_al8_t_out.size();
    // total_runs += c_ta_f_al5_t_out.size() + c_out_f_al5_t_out.size() + c_pb_f_al5_t_out.size() + c_c_f_al5_t_out.size();
    // total_runs += c_ta_f_al3_t_out.size() + c_pb_f_al3_t_out.size();

    // list_of_runs.reserve( total_runs );

    // target out runs
    if (fill_filterOut_runs == true){
        // list_of_runs.insert( list_of_runs.end(), c_au_f_out_t_out.begin(), c_au_f_out_t_out.end() );
        // list_of_runs.insert( list_of_runs.end(), c_atic_f_out_t_out.begin(), c_atic_f_out_t_out.end() );
        list_of_runs.insert( list_of_runs.end(), c_out_f_out_t_out.begin(), c_out_f_out_t_out.end() );
        list_of_runs.insert( list_of_runs.end(), c_pb_f_out_t_out.begin(), c_pb_f_out_t_out.end() );
        list_of_runs.insert( list_of_runs.end(), c_ta_f_out_t_out.begin(), c_ta_f_out_t_out.end() );
        list_of_runs.insert( list_of_runs.end(), c_c_f_out_t_out.begin(), c_c_f_out_t_out.end() ); 
    }

    // Bi 1 cm filter
    if (fill_bi1_runs == true){
        list_of_runs.insert( list_of_runs.end(), c_au_f_bi_t_out.begin(), c_au_f_bi_t_out.end() );
        list_of_runs.insert( list_of_runs.end(), c_out_f_bi_t_out.begin(), c_out_f_bi_t_out.end() );
        list_of_runs.insert( list_of_runs.end(), c_pb_f_bi_t_out.begin(), c_pb_f_bi_t_out.end() );
        list_of_runs.insert( list_of_runs.end(), c_ta_f_bi_t_out.begin(), c_ta_f_bi_t_out.end() );
    }

    // Al 8 cm filter
    if (fill_al8_runs == true){
        list_of_runs.insert( list_of_runs.end(), c_ta_f_al8_t_out.begin(), c_ta_f_al8_t_out.end() );
    }

    // Al 5 cm filter
    if (fill_al5_runs == true){
        list_of_runs.insert( list_of_runs.end(), c_ta_f_al5_t_out.begin(), c_ta_f_al5_t_out.end() );
        list_of_runs.insert( list_of_runs.end(), c_out_f_al5_t_out.begin(), c_out_f_al5_t_out.end() );
        list_of_runs.insert( list_of_runs.end(), c_pb_f_al5_t_out.begin(), c_pb_f_al5_t_out.end() );
        list_of_runs.insert( list_of_runs.end(), c_c_f_al5_t_out.begin(), c_c_f_al5_t_out.end() );
    }

    // Al 3 cm filter
    if (fill_al3_runs == true){
        list_of_runs.insert( list_of_runs.end(), c_ta_f_al3_t_out.begin(), c_ta_f_al3_t_out.end() );
        list_of_runs.insert( list_of_runs.end(), c_pb_f_al3_t_out.begin(), c_pb_f_al3_t_out.end() );
    }
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