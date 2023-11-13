/**
 * @file ptbc.h
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
const std::string filter_name("al8"); //bi1, al3, al5, al8, c2
const std::string endf_file_name("Al_tot_xsec.txt"); //Al_tot_xsec.txt, Bi_tot_xsec.txt
const std::string mode("run");

//Root file
TFile *outputRootFile = 0;

TH1D* tof_hist_filter_in = 0;
TH1D* energy_hist_filter_in = 0;
TH1D* tof_hist_filter_out = 0;
TH1D* energy_hist_filter_out = 0;
TH1D* transmission_hist_e = 0;
TH1D* cross_section_hist_e = 0;

Int_t bins_per_decade = 100;
Double_t flight_path_length_PTB = 182.65 - 0.41; //m
Double_t neutron_mass = 939.56542052; //in MeV
Double_t speed_of_light = 299792458.0; //in m/s

Double_t n_Bi_1cm = (1.0 /*cm*/) * (9.78 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (208.9804 /*g/mole*/);
Double_t n_Al_3cm = (3.0 /*cm*/) * (2.70 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (26.9815 /*g/mole*/);
Double_t n_Al_5cm = (5.0 /*cm*/) * (2.70 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (26.9815 /*g/mole*/);
Double_t n_Al_8cm = (8.0 /*cm*/) * (2.70 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (26.9815 /*g/mole*/);
std::map<std::string, Double_t> num_density_map;
void fillNumDensityMap(){
    num_density_map.emplace("bi1", n_Bi_1cm);
    num_density_map.emplace("al3", n_Al_3cm);
    num_density_map.emplace("al5", n_Al_5cm);
    num_density_map.emplace("al8", n_Al_8cm);
}

Double_t t_gamma_PTB = (flight_path_length_PTB / speed_of_light) * 1e9; //converting into ns

Double_t t[6][2];
Double_t a[6][2];

//// vectors to store run numbers
std::vector<Int_t> filter_in_runs;
std::vector<Int_t> filter_out_runs;

//Filter In runs
//Bi (1 cm) Filter
std::vector<Int_t> c_au_f_bi_t_out = {117369, 117370, 117371, 117372, 117373, 117374, 117375, 117376, 117377, 117378, 117379, 117380, 117381, 117382, 117383, 117384, 117385, 117386, 117387, 117388, 117389, 117390, 117436};
std::vector<Int_t> c_out_f_bi_t_out = {117391, 117392, 117393, 117394, 117395, 117396, 117397};
std::vector<Int_t> c_pb_f_bi_t_out = {117422, 117423, 117424, 117425, 117426, 117427, 117428};
std::vector<Int_t> c_ta_f_bi_t_out = {117437, 117439, 117440, 117441, 117442, 117443};

//Al (8 cm) Filter
std::vector<Int_t> c_ta_f_al8_t_out = {117454, 117455, 117456, 117457, 117458, 117459, 117460}; //117449, 117450, 117451, 117453, 

//Al (5 cm) Filter
std::vector<Int_t> c_ta_f_al5_t_out = {117462, 117463, 117464, 117465, 117466, 117467, 117476, 117477, 117479, 117480, 117481, 117482, 117483, 117484}; //117461, 117478, 117483, 117484
std::vector<Int_t> c_out_f_al5_t_out = {117470, 117471, 117472, 117473, 117474, 117475};
std::vector<Int_t> c_pb_f_al5_t_out = {117497, 117498, 117499, 117500, 117501, 117502};
std::vector<Int_t> c_c_f_al5_t_out = {117503, 117504, 117505, 117506, 117507, 117508, 117509, 117510};
std::vector<Int_t> c_au_f_al5_t_out = {117519, 117520, 117521, 117522, 117530, 117531, 117532};

//Al (3 cm) Filter
std::vector<Int_t> c_ta_f_al3_t_out = {117485, 117486, 117487, 117488, 117489, 117490, 117491, 117492, 117493};
std::vector<Int_t> c_pb_f_al3_t_out = {117496};

//Filter Out runs
std::vector<Int_t> c_au_f_out_t_out = {117350, 117357, 117358, 117359, 117362, 117363, 117364, 117365, 117366, 117367, 117368};
std::vector<Int_t> c_atic_f_out_t_out = {117351, 117355, 117356}; // 117352, 117353,
std::vector<Int_t> c_out_f_out_t_out = {117398, 117405, 117406, 117408, 117409, 117410, 117411, 117412, 117398}; //
std::vector<Int_t> c_pb_f_out_t_out = {117429, 117430, 117431, 117432, 117433, 117434, 117435};
std::vector<Int_t> c_ta_f_out_t_out = {117444, 117445, 117446, 117447, 117448}; //117452
std::vector<Int_t> c_c_f_out_t_out = {117511, 117512, 117513, 117514, 117515, 117516, 117517, 117518};

///////////// Transmisssion setup in
//Filter Out runs
// c_ta_f_out_t_out_ts: 117543, 117544, 117545, 117546, 117547
// c_out_f_out_t_out_ts: 117548, 117549, 117550, 117551, 117552
std::vector<Int_t> f_out_t_out_ts = {117543, 117544, 117545, 117546, 117547, 117548, 117549, 117550, 117551, 117552};

//Carbon Target
// c_ta_f_out_t_c_ts = 117559, 117560, 117561, 117562, 117563, 117564, 117565, 117566
// c_ta_f_al5_t_c_ts = 117567, 117568, 117569, 117570, 117571, 117572, 117573, 117574
// c_out_f_out_t_c_ts = 117575, 117576
std::vector<Int_t> f_out_t_c_ts = {117559, 117560, 117561, 117562, 117563, 117564, 117565, 117566, 117575, 117576};
std::vector<Int_t> f_al5_t_c_ts = {117567, 117568, 117569, 117570, 117571, 117572, 117573, 117574};

void fillCuts(){
    t[0][0] = 780.0;
    a[0][0] = 7700.0;

    t[0][1] = 1090.0;
    a[0][1] = 5451.0;

    t[1][0] = 1090.0;
    a[1][0] = 5451.0;

    t[1][1] = 2605.0;
    a[1][1] = 7167.0;

    t[2][0] = 2605.0;
    a[2][0] = 7960.0;

    t[2][1] = 2856.0;
    a[2][1] = 7960.0;

    t[3][0] = 2856.0;
    a[3][0] = 7432.0;

    t[3][1] = 15290.0;
    a[3][1] = 7432.0;

    t[4][0] = 15290.0;
    a[4][0] = 7432.0;

    t[4][1] = 18708.0;
    a[4][1] = 4000.0; //2416

    t[5][0] = 18708.0;
    a[5][0] = 4000.0; //2416

    t[5][1] = 1e8;
    a[5][1] = 4000.0; //3076
}

void fillRuns(){
    if (!mode.compare("test"))
    {
        filter_in_runs.push_back(117389);
        filter_out_runs.push_back(117367);
        return;
    }
    // Bi 1 cm filter
    if(!filter_name.compare("bi1")){
        cout << "Setting up Bi (1 cm) run list" << endl;
        // filter_in_runs.reserve( c_au_f_bi_t_out.size() + c_out_f_bi_t_out.size() + c_pb_f_bi_t_out.size() + c_ta_f_bi_t_out.size() );
        filter_in_runs.insert( filter_in_runs.end(), c_au_f_bi_t_out.begin(), c_au_f_bi_t_out.end() );
        filter_in_runs.insert( filter_in_runs.end(), c_out_f_bi_t_out.begin(), c_out_f_bi_t_out.end() );
        filter_in_runs.insert( filter_in_runs.end(), c_pb_f_bi_t_out.begin(), c_pb_f_bi_t_out.end() );
        filter_in_runs.insert( filter_in_runs.end(), c_ta_f_bi_t_out.begin(), c_ta_f_bi_t_out.end() );
    }

    // Al 8 cm filter
    if(!filter_name.compare("al8")){
        cout << "Setting up Al (8 cm) run list" << endl;
        // filter_in_runs.reserve( c_ta_f_al8_t_out.size() );
        filter_in_runs.insert( filter_in_runs.end(), c_ta_f_al8_t_out.begin(), c_ta_f_al8_t_out.end() );
    }

    // Al 5 cm filter
    if(!filter_name.compare("al5")){
        cout << "Setting up Al (5 cm) run list" << endl;
        // filter_in_runs.reserve( c_ta_f_al5_t_out.size() + c_out_f_al5_t_out.size() + c_pb_f_al5_t_out.size() + c_c_f_al5_t_out.size() + c_au_f_al5_t_out.size() ); //
        filter_in_runs.insert( filter_in_runs.end(), c_ta_f_al5_t_out.begin(), c_ta_f_al5_t_out.end() );
        filter_in_runs.insert( filter_in_runs.end(), c_out_f_al5_t_out.begin(), c_out_f_al5_t_out.end() );
        filter_in_runs.insert( filter_in_runs.end(), c_pb_f_al5_t_out.begin(), c_pb_f_al5_t_out.end() );
        filter_in_runs.insert( filter_in_runs.end(), c_c_f_al5_t_out.begin(), c_c_f_al5_t_out.end() );
        filter_in_runs.insert( filter_in_runs.end(), c_au_f_al5_t_out.begin(), c_au_f_al5_t_out.end() );
    }

    // Al 3 cm filter
    if(!filter_name.compare("al3")){
        cout << "Setting up Al (3 cm) run list" << endl;
        // filter_in_runs.reserve( c_ta_f_al3_t_out.size() + c_pb_f_al3_t_out.size() );
        filter_in_runs.insert( filter_in_runs.end(), c_ta_f_al3_t_out.begin(), c_ta_f_al3_t_out.end() );
        filter_in_runs.insert( filter_in_runs.end(), c_pb_f_al3_t_out.begin(), c_pb_f_al3_t_out.end() );
    }

    cout << "Setting up filter out run list" << endl;
    // filter_out_runs.reserve( c_au_f_out_t_out.size() + c_atic_f_out_t_out.size() + c_out_f_out_t_out.size() + c_pb_f_out_t_out.size() + c_ta_f_out_t_out.size() + c_c_f_out_t_out.size() );
    filter_out_runs.insert( filter_out_runs.end(), c_au_f_out_t_out.begin(), c_au_f_out_t_out.end() );
    filter_out_runs.insert( filter_out_runs.end(), c_atic_f_out_t_out.begin(), c_atic_f_out_t_out.end() );
    filter_out_runs.insert( filter_out_runs.end(), c_out_f_out_t_out.begin(), c_out_f_out_t_out.end() );
    filter_out_runs.insert( filter_out_runs.end(), c_pb_f_out_t_out.begin(), c_pb_f_out_t_out.end() );
    filter_out_runs.insert( filter_out_runs.end(), c_ta_f_out_t_out.begin(), c_ta_f_out_t_out.end() );
    filter_out_runs.insert( filter_out_runs.end(), c_c_f_out_t_out.begin(), c_c_f_out_t_out.end() );
}

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

// void FillHistograms(Double_t tof, Float_t amp, TH1D* tof_hist, TH1D* energy_hist){
//     //Filling the histograms
//     for (int k = 0; k < 6; k++)
//     {
//         if (tof >= t[k][0] && tof < t[k][1])
//         {
//             if ( (Double_t) amp > yOnTheCutLine(t[k][0], a[k][0], t[k][1], a[k][1], tof) )
//             {
//                 tof_hist->Fill(tof);
//                 energy_hist->Fill( TOFToEnergy(tof * 1e-9) );
//                 break;
//             }
//         }
//     }
// }

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
    if (num == 1)
    {
        return decadePower;
    }
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