/**
 * @file parameters.h
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2024-02-22
 */

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <string>

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

Double_t flight_path_length_PTB = 182.65 - 0.41; //m
Double_t flight_path_length_PTBC = 182.65 - 0.41; //m
Double_t flight_path_length_FIMG = 183.5 - 0.41; //m
Double_t neutron_mass = 939.56542052; //in MeV
Double_t speed_of_light = 299792458.0; //in m/s
Double_t ar_bottle_pressure = 197.385 * 1e5; // in Pa (SI unit)
Double_t ar_bottle_temp = 293.0; // in Kelvin
Double_t delT_pkup_ptbc = 660.0; //in ns
Double_t delT_pkup_fimg = 630.0; //in ns
Double_t min_tof_PTBC = 800.0; //in ns

Double_t n_Bi_1cm = (1.0 /*cm*/) * (9.78 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (208.9804 /*g/mole*/);
Double_t n_Bi_1p2cm = (1.2 /*cm*/) * (9.78 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (208.9804 /*g/mole*/);
Double_t n_Al_3cm = (3.0 /*cm*/) * (2.70 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (26.9815 /*g/mole*/);
Double_t n_Al_5cm = (5.0 /*cm*/) * (2.70 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (26.9815 /*g/mole*/);
Double_t n_Al_8cm = (8.0 /*cm*/) * (2.70 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (26.9815 /*g/mole*/);
Double_t n_C_1p2cm = 0.105;// (1.2 /*cm*/) * (2.267 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (12.011 /*g/mole*/);
Double_t n_CFib_1cm = (1.0 /*cm*/) * (1.8 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (12.011 /*g/mole*/);
Double_t n_Ar_bottle = (11.0 /*cm*/) * ((ar_bottle_pressure)/(8.31446261815324 * ar_bottle_temp * 1e6)) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/);
Double_t n_Al_0p2cm = (0.2 /*cm*/) * (2.70 /*g/cm3*/) * (6.02214076e23 /*atoms/mole*/) * (1e-24 /*cm2/barn*/) / (26.9815 /*g/mole*/);

Double_t t_gamma_PTB = (flight_path_length_PTB / speed_of_light) * 1e9; //converting into ns
Double_t t_gamma_PTBC = (flight_path_length_PTB / speed_of_light) * 1e9; //converting into ns
Double_t t_gamma_FIMG = (flight_path_length_FIMG / speed_of_light) * 1e9; //converting into ns

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

//Filter In runs
//Bi (1 cm) Filter
//fimg starts from 117386
std::vector<Int_t> c_au_f_bi_t_out = {117386, 117387, 117388, 117389, 117390, 117436}; //117369, 117370, 117371, 117372, 117373, 117374, 117375, 117376, 117377, 117378, 117379, 117380, 117381, 117382, 117383, 117384, 117385, 
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
// std::vector<Int_t> c_au_f_out_t_out = {117350, 117357, 117358, 117359, 117362, 117363, 117364, 117365, 117366, 117367, 117368}; //
// std::vector<Int_t> c_atic_f_out_t_out = {117351, 117355, 117356}; // 117352, 117353,
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
std::vector<Int_t> f_out_t_c_ts = {117559, 117560, 117561, 117562, 117563, 117564, 117565, 117566, 117575, 117576, 117577, 117578, 117579, 117580, 117581, 117583};
std::vector<Int_t> f_al5_t_c_ts = {117567, 117568, 117569, 117570, 117571, 117572, 117573, 117574};

// Al5 Target
std::vector<Int_t> f_out_t_al5_ts = {117584, 117586, 117587, 117588, 117589, 117590, 117591};

//Bi 1.2 cm target
std::vector<Int_t> f_out_t_bi_ts = {117593, 117594, 117595, 117596};

//Empty Bottle
std::vector<Int_t> f_out_t_emptyBot_ts = {117607, 117608, 117609, 117612, 117613, 117614, 117615, 117616, 117617, 117618, 117619, 117620};

//Empty Bottle Rotated
std::vector<Int_t> f_out_t_emptyBotRot_ts = {117623, 117624, 117627, 117628, 117629, 117630, 117631, 117632, 117633, 117634, 117635, 117636, 117637, 117638, 117639, 117640, 117641};

//Empty Bottle Rotated
std::vector<Int_t> f_out_t_emptyBotRotBack_ts = {117642, 117645, 117646, 117647, 117648, 117649};

//Argon Bottle
std::vector<Int_t> f_out_t_arBot_ts = {117665, 117680, 117684, 117685, 117688, 117692, 117693, 117710, 117711, 117712, 117713, 117716, 117717, 117718, 117719, 117724, 117725, 117726, 117727, 117728, 117729, 117730, 117731, 117732, 117733, 117734, 117735, 117736, 117737, 117738}; //117714, 117686, 
//117687, 

//Argon Bottle (Argon Out)
std::vector<Int_t> f_out_t_arBotEmpty_ts = {117746, 117747, 117748, 117749, 117750, 117751, 117752, 117753, 117754, 117755, 117756, 117757, 117761, 117762, 117763, 117764, 117765, 117766}; //

////////////////// Detector cuts

//PTBC cuts
Double_t t_det2[4][2];
Double_t a_det2[4][2];
Double_t t_det3to7[5][2][2];
Double_t a_det3to7[5][2][2];
Double_t t_para[4][2];
Double_t a_para[4][2];

Double_t t_det5_early_runs[4][2];
Double_t a_det5_early_runs[4][2];

void fillCutsPTBC(){
    //Det 2
    //End points of cut line 1
    t_det2[0][0] = 800.0;
    a_det2[0][0] = 8000.0;

    t_det2[0][1] = 2605.0;
    a_det2[0][1] = 8000.0;

    //End points of cut line 2
    t_det2[1][0] = 2600.0;
    a_det2[1][0] = 10000.0;

    t_det2[1][1] = 2800.0;
    a_det2[1][1] = 10000.0;

    //End points of cut line 3
    t_det2[2][0] = 2800.0;
    a_det2[2][0] = 8000.0;

    t_det2[2][1] = 5000.0;
    a_det2[2][1] = 8000.0;

    //End points of cut line 4
    t_det2[3][0] = 5000.0;
    a_det2[3][0] = 4000.0; 

    t_det2[3][1] = 1e8;
    a_det2[3][1] = 4000.0;

    //Det 3
    t_det3to7[0][0][0] = 800.0;
    a_det3to7[0][0][0] = 5000.0;

    t_det3to7[0][0][1] = 3000.0;
    a_det3to7[0][0][1] = 5000.0;

    t_det3to7[0][1][0] = 3000.0;
    a_det3to7[0][1][0] = 3500.0;

    t_det3to7[0][1][1] = 1e8;
    a_det3to7[0][1][1] = 3500.0;

    //Det 4
    t_det3to7[1][0][0] = 800.0;
    a_det3to7[1][0][0] = 6000.0;

    t_det3to7[1][0][1] = 2000.0;
    a_det3to7[1][0][1] = 6000.0;

    t_det3to7[1][1][0] = 2000.0;
    a_det3to7[1][1][0] = 3500.0;

    t_det3to7[1][1][1] = 1e8;
    a_det3to7[1][1][1] = 3500.0;

    //Det 5
    t_det3to7[2][0][0] = 800.0;
    a_det3to7[2][0][0] = 7000.0;

    t_det3to7[2][0][1] = 7000.0;
    a_det3to7[2][0][1] = 7000.0;

    t_det3to7[2][1][0] = 7000.0;
    a_det3to7[2][1][0] = 3500.0;

    t_det3to7[2][1][1] = 1e8;
    a_det3to7[2][1][1] = 3500.0;

    //Det 5 - from Run 117386 to 117390
    //End points of cut line 1
    t_det5_early_runs[0][0] = 800.0;
    a_det5_early_runs[0][0] = 7000.0;

    t_det5_early_runs[0][1] = 2300.0;
    a_det5_early_runs[0][1] = 7000.0;

    //End points of cut line 2
    t_det5_early_runs[1][0] = 2300.0;
    a_det5_early_runs[1][0] = 8500.0;

    t_det5_early_runs[1][1] = 4000.0;
    a_det5_early_runs[1][1] = 8500.0;

    //End points of cut line 3
    t_det5_early_runs[2][0] = 4000.0;
    a_det5_early_runs[2][0] = 8500.0;

    t_det5_early_runs[2][1] = 7000.0;
    a_det5_early_runs[2][1] = 3500.0;

    //End points of cut line 3
    t_det5_early_runs[3][0] = 7000.0;
    a_det5_early_runs[3][0] = 3500.0;

    t_det5_early_runs[3][1] = 1e8;
    a_det5_early_runs[3][1] = 3500.0;

    //Det 6
    t_det3to7[3][0][0] = 800.0;
    a_det3to7[3][0][0] = 6000.0;

    t_det3to7[3][0][1] = 6000.0;
    a_det3to7[3][0][1] = 6000.0;

    t_det3to7[3][1][0] = 6000.0;
    a_det3to7[3][1][0] = 4000.0;

    t_det3to7[3][1][1] = 1e8;
    a_det3to7[3][1][1] = 4000.0;

    //Det 7
    t_det3to7[4][0][0] = 800.0;
    a_det3to7[4][0][0] = 4000.0;

    t_det3to7[4][0][1] = 4000.0;
    a_det3to7[4][0][1] = 4000.0;

    t_det3to7[4][1][0] = 4000.0;
    a_det3to7[4][1][0] = 3000.0;

    t_det3to7[4][1][1] = 1e8;
    a_det3to7[4][1][1] = 3000.0;
    
    //Parasitic Cuts
    t_para[0][0] = 800.0;
    a_para[0][0] = 5000.0;

    t_para[0][1] = 3000.0;
    a_para[0][1] = 5000.0;

    t_para[1][0] = 3000.0;
    a_para[1][0] = 4000.0;

    t_para[1][1] = 1e8;
    a_para[1][1] = 4000.0;
}

//My FIMG Cuts
//det number, cut line number, end points of the line
Double_t tof_cut_FIMG[2][5][2];
Double_t amp_cut_FIMG[2][5][2];

void fillCutsFIMG(){
    //Det 1
    tof_cut_FIMG[0][0][0] = 7e3;
    amp_cut_FIMG[0][0][0] = 4600.;
    tof_cut_FIMG[0][0][1] = 1e4;
    amp_cut_FIMG[0][0][1] = 3000.;

    tof_cut_FIMG[0][1][0] = 1e4;
    amp_cut_FIMG[0][1][0] = 3000.;
    tof_cut_FIMG[0][1][1] = 2e4;
    amp_cut_FIMG[0][1][1] = 1525.;

    tof_cut_FIMG[0][2][0] = 2e4;
    amp_cut_FIMG[0][2][0] = 1525.;
    tof_cut_FIMG[0][2][1] = 45090.;
    amp_cut_FIMG[0][2][1] = 1525.;

    tof_cut_FIMG[0][3][0] = 45090.;
    amp_cut_FIMG[0][3][0] = 1525.;
    tof_cut_FIMG[0][3][1] = 1e5;
    amp_cut_FIMG[0][3][1] = 500.;

    tof_cut_FIMG[0][4][0] = 1e5;
    amp_cut_FIMG[0][4][0] = 500.;
    tof_cut_FIMG[0][4][1] = 1e8;
    amp_cut_FIMG[0][4][1] = 500.;

    //Det 2
    tof_cut_FIMG[1][0][0] = 7e3;
    amp_cut_FIMG[1][0][0] = 4000.;
    tof_cut_FIMG[1][0][1] = 8e3;
    amp_cut_FIMG[1][0][1] = 3200.;

    tof_cut_FIMG[1][1][0] = 8e3;
    amp_cut_FIMG[1][1][0] = 3200.;
    tof_cut_FIMG[1][1][1] = 2e4;
    amp_cut_FIMG[1][1][1] = 2250.;

    tof_cut_FIMG[1][2][0] = 2e4;
    amp_cut_FIMG[1][2][0] = 2250.;
    tof_cut_FIMG[1][2][1] = 31800.;
    amp_cut_FIMG[1][2][1] = 2250.;

    tof_cut_FIMG[1][3][0] = 31800.;
    amp_cut_FIMG[1][3][0] = 2250.;
    tof_cut_FIMG[1][3][1] = 7e4;
    amp_cut_FIMG[1][3][1] = 500.;

    tof_cut_FIMG[1][4][0] = 7e4;
    amp_cut_FIMG[1][4][0] = 500.;
    tof_cut_FIMG[1][4][1] = 1e8;
    amp_cut_FIMG[1][4][1] = 500.;
}

/////////////////////////////////////////////////////

//nTOF - PTBC cuts
//det_num, cut num, x or y
Double_t t_det2to4_nTOF[3][2][2];
Double_t a_det2to4_nTOF[3][2][2];
Double_t t_det5to7_nTOF[3][3][2];
Double_t a_det5to7_nTOF[3][3][2];

void fillCutsPTBC_nTOF(){
    //Det 2
    t_det2to4_nTOF[0][0][0] = 700.0;
    a_det2to4_nTOF[0][0][0] = 8000.0;
    t_det2to4_nTOF[0][0][1] = 20000.0;
    a_det2to4_nTOF[0][0][1] = 8000.0;

    t_det2to4_nTOF[0][1][0] = 20000.0;
    a_det2to4_nTOF[0][1][0] = 5000.0;
    t_det2to4_nTOF[0][1][1] = 1e8;
    a_det2to4_nTOF[0][1][1] = 5000.0;

    //Det 3
    t_det2to4_nTOF[1][0][0] = 700.0;
    a_det2to4_nTOF[1][0][0] = 7000.0;
    t_det2to4_nTOF[1][0][1] = 2000.0;
    a_det2to4_nTOF[1][0][1] = 7000.0;

    t_det2to4_nTOF[1][1][0] = 2000.0;
    a_det2to4_nTOF[1][1][0] = 4500.0;
    t_det2to4_nTOF[1][1][1] = 1e8;
    a_det2to4_nTOF[1][1][1] = 4500.0;

    //Det 4
    t_det2to4_nTOF[2][0][0] = 700.0;
    a_det2to4_nTOF[2][0][0] = 10000.0;
    t_det2to4_nTOF[2][0][1] = 2000.0;
    a_det2to4_nTOF[2][0][1] = 10000.0;

    t_det2to4_nTOF[2][1][0] = 2000.0;
    a_det2to4_nTOF[2][1][0] = 5000.0;
    t_det2to4_nTOF[2][1][1] = 1e8;
    a_det2to4_nTOF[2][1][1] = 5000.0;

    //Det 5
    t_det5to7_nTOF[0][0][0] = 700.0;
    a_det5to7_nTOF[0][0][0] = 10000.0;
    t_det5to7_nTOF[0][0][1] = 2000.0;
    a_det5to7_nTOF[0][0][1] = 10000.0;

    t_det5to7_nTOF[0][1][0] = 2000.0;
    a_det5to7_nTOF[0][1][0] = 8000.0;
    t_det5to7_nTOF[0][1][1] = 20000.0;
    a_det5to7_nTOF[0][1][1] = 8000.0;

    t_det5to7_nTOF[0][2][0] = 20000.0;
    a_det5to7_nTOF[0][2][0] = 4500.0;
    t_det5to7_nTOF[0][2][1] = 1e8;
    a_det5to7_nTOF[0][2][1] = 4500.0;

    //Det 6
    t_det5to7_nTOF[1][0][0] = 700.0;
    a_det5to7_nTOF[1][0][0] = 9000.0;
    t_det5to7_nTOF[1][0][1] = 2000.0;
    a_det5to7_nTOF[1][0][1] = 9000.0;

    t_det5to7_nTOF[1][1][0] = 2000.0;
    a_det5to7_nTOF[1][1][0] = 6000.0;
    t_det5to7_nTOF[1][1][1] = 20000.0;
    a_det5to7_nTOF[1][1][1] = 6000.0;

    t_det5to7_nTOF[1][2][0] = 20000.0;
    a_det5to7_nTOF[1][2][0] = 4500.0;
    t_det5to7_nTOF[1][2][1] = 1e8;
    a_det5to7_nTOF[1][2][1] = 4500.0;

    //Det 7
    t_det5to7_nTOF[2][0][0] = 700.0;
    a_det5to7_nTOF[2][0][0] = 8000.0;
    t_det5to7_nTOF[2][0][1] = 2000.0;
    a_det5to7_nTOF[2][0][1] = 8000.0;

    t_det5to7_nTOF[2][1][0] = 2000.0;
    a_det5to7_nTOF[2][1][0] = 4000.0;
    t_det5to7_nTOF[2][1][1] = 20000.0;
    a_det5to7_nTOF[2][1][1] = 4000.0;

    t_det5to7_nTOF[2][2][0] = 20000.0;
    a_det5to7_nTOF[2][2][0] = 3500.0;
    t_det5to7_nTOF[2][2][1] = 1e8;
    a_det5to7_nTOF[2][2][1] = 3500.0;
}

//nTOF - PTBC cuts
Double_t a_para_nTOF[6];

void fillCutsPTBC_para_nTOF(){
    a_para_nTOF[0] = 5000.0;
    a_para_nTOF[1] = 4500.0;
    a_para_nTOF[2] = 5000.0;
    a_para_nTOF[3] = 4500.0;
    a_para_nTOF[4] = 4500.0;
    a_para_nTOF[5] = 3500.0;
}

////////////////////////////////////////////////////// old stuff

// Double_t FIMG_tof_cut_low_det1 = 8102.0; //in ns
// Double_t FIMG_tof_cut_up_det1 = 100000.0; //in ns
// Double_t FIMG_tof_cut_low_det2 = 8174.0; //in ns
// Double_t FIMG_tof_cut_up_det2 = 89408.0; //in ns

// Double_t FIMG_min_amp_cut_det1 = 600.0; //a.u.
// Double_t FIMG_min_amp_cut_det2 = 400.0; //a.u.

// Double_t fimgCutFunction(Double_t x, Int_t det_num){
//     // // A * (log(x))^2 + B + log(x) + C
//     // // Values obtained from cutoffFitter.C
//     // Double_t A = 120556;
//     // Double_t B = -2.75589e6;
//     // Double_t C = 1.57536e7;
//     // return (A * TMath::Log(x) * TMath::Log(x) + B * TMath::Log(x) + C);

//     if(det_num == 1) {
//         // A / ( TMath::Log(x) + B )
//         // Values obtained from cutoffFitter.C
//         Double_t A = 4378.11;
//         Double_t B = -7.695;

//         return (A / (TMath::Log(x) + B));
//     }

//     if(det_num == 2) {
//         // A / ( TMath::Log(x) + B )
//         // Values obtained from cutoffFitter.C
//         Double_t A = 4145.31;
//         Double_t B = -7.67962;

//         return (A / (TMath::Log(x) + B));
//     }
// }

//////////// Old Cuts

// void fillCutsPTBC(){
//     t[0][0] = 780.0;
//     a[0][0] = 7700.0;

//     t[0][1] = 1090.0;
//     a[0][1] = 5451.0;

//     t[1][0] = 1090.0;
//     a[1][0] = 5451.0;

//     t[1][1] = 2605.0;
//     a[1][1] = 7167.0;

//     t[2][0] = 2605.0;
//     a[2][0] = 7960.0;

//     t[2][1] = 2856.0;
//     a[2][1] = 7960.0;

//     t[3][0] = 2856.0;
//     a[3][0] = 7432.0;

//     t[3][1] = 15290.0;
//     a[3][1] = 7432.0;

//     t[4][0] = 15290.0;
//     a[4][0] = 7432.0;

//     t[4][1] = 18708.0;
//     a[4][1] = 4000.0; //2416

//     t[5][0] = 18708.0;
//     a[5][0] = 4000.0; //2416

//     t[5][1] = 1e8;
//     a[5][1] = 4000.0; //3076
// }

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