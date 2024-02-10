/**
 * @file cutoffAnalysis_FIMG.h
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2024-02-08
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
const std::string target_name("no_filter"); //bi1, al3, al5, al8, c1p2_ts, al5_ts, bi1p2_ts, cf_bottle, cf_bottle_rot, cf_bottle_rotBack, ar_bottle_full
const std::string target_out_name("none"); //none, none_ts, ar_bottle
const std::string target_name_title("No target");
//Bi (1 cm), Target Bi (1.2 cm), Al (3 cm), Al (5 cm), Target Al (5 cm), Al (8 cm), Target C (1.2 cm), Empty Bottle, Empty Bottle Rotated
//Argon Tank

// const std::string mode("run");
// bi1, al3, al5, al8 - none

//Root file
TFile *outputRootFile = 0;

// TH1D* trans_loose_cut_det_1 = 0;
// TH1D* trans_mid_cut_det_1 = 0;
// TH1D* trans_tight_cut_det_1 = 0;
// TH1D* trans_loose_cut_det_2 = 0;
// TH1D* trans_mid_cut_det_2 = 0;
// TH1D* trans_tight_cut_det_2 = 0;

TH2D* FIMG_tof_amp_total = 0;
TH2D* FIMG_tof_amp_dedi_det1 = 0;
TH2D* FIMG_tof_amp_dedi_det2 = 0;
TH2D* FIMG_tof_amp_para_det1 = 0;
TH2D* FIMG_tof_amp_para_det2 = 0;

TH2D* FIMG_tof_amp_det1_afterCuts = 0;
TH2D* FIMG_tof_amp_det2_afterCuts = 0;

//cut plots
TCutG* FIMG_tof_amp_cut_dedi_det1;
TCutG* FIMG_tof_amp_cut_dedi_det2;
TCutG* FIMG_tof_amp_cut_para_det1;
TCutG* FIMG_tof_amp_cut_para_det2;

Double_t tof_min = 1e2;
Double_t tof_max = 1e8;
Double_t amp_min = 0.;
Double_t amp_max = 50000.;
Int_t num_bins_amp = 5000;

Int_t bins_per_decade = 1000;
Double_t flight_path_length_PTB = 182.65 - 0.41; //m
Double_t flight_path_length_FIMG = 183.5 - 0.41; //m
Double_t neutron_mass = 939.56542052; //in MeV
Double_t speed_of_light = 299792458.0; //in m/s
Double_t ar_bottle_pressure = 197.385 * 1e5; // in Pa (SI unit)
Double_t ar_bottle_temp = 293.0; // in Kelvin
Double_t delT_pkup_ptbc = 660.0; //in ns
Double_t delT_pkup_fimg = 630.0; //in ns

Double_t FIMG_tof_cut_low_det1 = 8102.0; //in ns
Double_t FIMG_tof_cut_up_det1 = 100000.0; //in ns
Double_t FIMG_tof_cut_low_det2 = 8174.0; //in ns
Double_t FIMG_tof_cut_up_det2 = 89408.0; //in ns

Double_t FIMG_min_amp_cut_det1 = 600.0; //a.u.
Double_t FIMG_min_amp_cut_det2 = 400.0; //a.u.

Double_t t_gamma_PTB = (flight_path_length_PTB / speed_of_light) * 1e9; //converting into ns
Double_t t_gamma_FIMG = (flight_path_length_FIMG / speed_of_light) * 1e9; //converting into ns

Double_t t[6][2];
Double_t a[6][2];

//// vectors to store run numbers
std::vector<Int_t> f_al5_t_out_runs;
std::vector<Int_t> f_out_t_out_runs;

std::vector<Int_t> filter_in_runs;
std::vector<Int_t> filter_out_runs;

std::vector<Int_t> ts_target_in_runs;
std::vector<Int_t> ts_target_out_runs;

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
std::vector<Int_t> f_out_t_arBot_ts = {117665, 117680, 117684, 117685, 117686, 117687, 117688, 117692, 117693, 117710, 117711, 117712, 117713, 117714, 117716, 117717, 117718, 117719, 117724, 117725, 117726, 117727, 117728, 117729, 117730, 117731, 117732, 117733, 117734, 117735, 117736, 117737, 117738};

//Argon Bottle (Argon Out)
std::vector<Int_t> f_out_t_arBotEmpty_ts = {117746, 117747, 117748, 117749, 117750, 117751, 117752, 117753, 117754, 117755, 117756, 117757, 117761, 117762, 117763, 117764, 117765, 117766}; //

void fillRuns(){
    // if (!mode.compare("test"))
    // {
    //     filter_in_runs.push_back(117389);
    //     filter_out_runs.push_back(117367);
    //     return;
    // }

    // Bi 1 cm filter
    if(!target_name.compare("bi1")){
        cout << "Setting up Bi (1 cm) run list" << endl;
        // filter_in_runs.reserve( c_au_f_bi_t_out.size() + c_out_f_bi_t_out.size() + c_pb_f_bi_t_out.size() + c_ta_f_bi_t_out.size() );
        ts_target_in_runs.insert( ts_target_in_runs.end(), c_au_f_bi_t_out.begin(), c_au_f_bi_t_out.end() );
        ts_target_in_runs.insert( ts_target_in_runs.end(), c_out_f_bi_t_out.begin(), c_out_f_bi_t_out.end() );
        ts_target_in_runs.insert( ts_target_in_runs.end(), c_pb_f_bi_t_out.begin(), c_pb_f_bi_t_out.end() );
        ts_target_in_runs.insert( ts_target_in_runs.end(), c_ta_f_bi_t_out.begin(), c_ta_f_bi_t_out.end() );
    }

    // Al 8 cm filter
    if(!target_name.compare("al8")){
        cout << "Setting up Al (8 cm) run list" << endl;
        // ts_target_in_runs.reserve( c_ta_f_al8_t_out.size() );
        ts_target_in_runs.insert( ts_target_in_runs.end(), c_ta_f_al8_t_out.begin(), c_ta_f_al8_t_out.end() );
    }

    // Al 5 cm filter
    if(!target_name.compare("al5")){
        cout << "Setting up Al (5 cm) run list" << endl;
        // ts_target_in_runs.reserve( c_ta_f_al5_t_out.size() + c_out_f_al5_t_out.size() + c_pb_f_al5_t_out.size() + c_c_f_al5_t_out.size() + c_au_f_al5_t_out.size() ); //
        ts_target_in_runs.insert( ts_target_in_runs.end(), c_ta_f_al5_t_out.begin(), c_ta_f_al5_t_out.end() );
        ts_target_in_runs.insert( ts_target_in_runs.end(), c_out_f_al5_t_out.begin(), c_out_f_al5_t_out.end() );
        ts_target_in_runs.insert( ts_target_in_runs.end(), c_pb_f_al5_t_out.begin(), c_pb_f_al5_t_out.end() );
        ts_target_in_runs.insert( ts_target_in_runs.end(), c_c_f_al5_t_out.begin(), c_c_f_al5_t_out.end() );
        ts_target_in_runs.insert( ts_target_in_runs.end(), c_au_f_al5_t_out.begin(), c_au_f_al5_t_out.end() );
    }

    // Al 3 cm filter
    if(!target_name.compare("al3")){
        cout << "Setting up Al (3 cm) run list" << endl;
        // ts_target_in_runs.reserve( c_ta_f_al3_t_out.size() + c_pb_f_al3_t_out.size() );
        ts_target_in_runs.insert( ts_target_in_runs.end(), c_ta_f_al3_t_out.begin(), c_ta_f_al3_t_out.end() );
        ts_target_in_runs.insert( ts_target_in_runs.end(), c_pb_f_al3_t_out.begin(), c_pb_f_al3_t_out.end() );
    }

    //--------------------------------------------------

    if(!target_name.compare("bi1p2_ts")){
        cout << "Setting up Bi (1.2 cm) target run list" << endl;
        ts_target_in_runs.insert( ts_target_in_runs.end(), f_out_t_bi_ts.begin(), f_out_t_bi_ts.end() );
    }

    if(!target_name.compare("c1p2_ts")){
        cout << "Setting up C (1.2 cm) target run list" << endl;
        ts_target_in_runs.insert( ts_target_in_runs.end(), f_out_t_c_ts.begin(), f_out_t_c_ts.end() );
    }

    if(!target_name.compare("al5_ts")){
        cout << "Setting up Al (5 cm) target run list" << endl;
        ts_target_in_runs.insert( ts_target_in_runs.end(), f_out_t_al5_ts.begin(), f_out_t_al5_ts.end() );
    }

    if(!target_name.compare("cf_bottle")){
        cout << "Setting up Empty Bottle target run list" << endl;
        ts_target_in_runs.insert( ts_target_in_runs.end(), f_out_t_emptyBot_ts.begin(), f_out_t_emptyBot_ts.end() );
    }

    if(!target_name.compare("cf_bottle_rot")){
        cout << "Setting up Empty Bottle Rotated target run list" << endl;
        ts_target_in_runs.insert( ts_target_in_runs.end(), f_out_t_emptyBotRot_ts.begin(), f_out_t_emptyBotRot_ts.end() );
    }

    if(!target_name.compare("cf_bottle_rotBack")){
        cout << "Setting up Empty Bottle Rotated Back target run list" << endl;
        ts_target_in_runs.insert( ts_target_in_runs.end(), f_out_t_emptyBotRotBack_ts.begin(), f_out_t_emptyBotRotBack_ts.end() );
    }

    if(!target_name.compare("ar_bottle_full")){
        cout << "Setting up Argon Bottle target run list" << endl;
        ts_target_in_runs.insert( ts_target_in_runs.end(), f_out_t_arBot_ts.begin(), f_out_t_arBot_ts.end() );

        if(!target_out_name.compare("empty_bottle")){
            ts_target_out_runs.insert( ts_target_out_runs.end(), f_out_t_emptyBot_ts.begin(), f_out_t_emptyBot_ts.end() );
            ts_target_out_runs.insert( ts_target_out_runs.end(), f_out_t_emptyBotRot_ts.begin(), f_out_t_emptyBotRot_ts.end() );
            ts_target_out_runs.insert( ts_target_out_runs.end(), f_out_t_emptyBotRotBack_ts.begin(), f_out_t_emptyBotRotBack_ts.end() );
        }

        if(!target_out_name.compare("ar_bottle")){
            ts_target_out_runs.insert( ts_target_out_runs.end(), f_out_t_arBotEmpty_ts.begin(), f_out_t_arBotEmpty_ts.end() );
        }
    }

    // // Al 5 cm filter
    // cout << "Setting up Al (5 cm) filter run list" << endl;
    // // filter_in_runs.reserve( c_ta_f_al5_t_out.size() + c_out_f_al5_t_out.size() + c_pb_f_al5_t_out.size() + c_c_f_al5_t_out.size() + c_au_f_al5_t_out.size() ); //
    // f_al5_t_out_runs.insert( f_al5_t_out_runs.end(), c_ta_f_al5_t_out.begin(), c_ta_f_al5_t_out.end() );
    // f_al5_t_out_runs.insert( f_al5_t_out_runs.end(), c_out_f_al5_t_out.begin(), c_out_f_al5_t_out.end() );
    // f_al5_t_out_runs.insert( f_al5_t_out_runs.end(), c_pb_f_al5_t_out.begin(), c_pb_f_al5_t_out.end() );
    // f_al5_t_out_runs.insert( f_al5_t_out_runs.end(), c_c_f_al5_t_out.begin(), c_c_f_al5_t_out.end() );
    // f_al5_t_out_runs.insert( f_al5_t_out_runs.end(), c_au_f_al5_t_out.begin(), c_au_f_al5_t_out.end() );

    cout << "Setting up filter out run list" << endl;
    // filter_out_runs.reserve( c_au_f_out_t_out.size() + c_atic_f_out_t_out.size() + c_out_f_out_t_out.size() + c_pb_f_out_t_out.size() + c_ta_f_out_t_out.size() + c_c_f_out_t_out.size() );
    // f_out_t_out_runs.insert( f_out_t_out_runs.end(), c_au_f_out_t_out.begin(), c_au_f_out_t_out.end() );
    // f_out_t_out_runs.insert( f_out_t_out_runs.end(), c_atic_f_out_t_out.begin(), c_atic_f_out_t_out.end() );
    f_out_t_out_runs.insert( f_out_t_out_runs.end(), c_out_f_out_t_out.begin(), c_out_f_out_t_out.end() );
    f_out_t_out_runs.insert( f_out_t_out_runs.end(), c_pb_f_out_t_out.begin(), c_pb_f_out_t_out.end() );
    f_out_t_out_runs.insert( f_out_t_out_runs.end(), c_ta_f_out_t_out.begin(), c_ta_f_out_t_out.end() );
    f_out_t_out_runs.insert( f_out_t_out_runs.end(), c_c_f_out_t_out.begin(), c_c_f_out_t_out.end() );
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

// Double_t fimgCutFunction(Double_t x, Int_t det_num, std::string cut_type){

//     Double_t A = 0;
//     Double_t B = 0;
//     Double_t C = 0;

//     if(det_num == 1) {
//         // A / ( TMath::Log(x) + B )
//         if (!cut_type.compare("loose"))
//         {
//             A = 4000;
//             B = -TMath::Log(7750);
//             C = 600;
//         }
//         if (!cut_type.compare("mid"))
//         {
//             A = 4500;
//             B = -TMath::Log(7750);
//             C = 600;
//         }
//         if (!cut_type.compare("tight"))
//         {
//             A = 7000;
//             B = -TMath::Log(7750);
//             C = 600;
//         }
        
//         return (A / (TMath::Log(x) + B));
//     }
//     if(det_num == 2) {
//         // A / ( TMath::Log(x) + B )
//         if (!cut_type.compare("loose"))
//         {
//             A = 3000;
//             B = -TMath::Log(7080);
//             C = 400;
//         }
//         if (!cut_type.compare("mid"))
//         {
//             A = 5000;
//             B = -TMath::Log(7080);
//             C = 400;
//         }
//         if (!cut_type.compare("tight"))
//         {
//             A = 6500;
//             B = -TMath::Log(7080);
//             C = 400;
//         }
        
//         return ((A / (TMath::Log(x) + B)) + C);
//     }
// }

Double_t fimgCutFunction(Double_t tof, Int_t det_num, std::string cut_type, Float_t PulseIntensity){

    if (!cut_type.compare("ntof"))
    {
        if (det_num == 1)
        {
            if (tof >= 1e5)
            {
                return 500;
            }
            if (tof >= 1e4 && tof < 1e5)
            {
                return ((500 - 2350)*(tof - 1e4)/(1e5 - 1e4) + 2350);
            }
        }
        if (det_num == 2)
        {
            if (tof >= 1e5)
            {
                return 500;
            }
            if (tof >= 1e4 && tof < 1e5)
            {
                return ((500 - 3000)*(tof - 1e4)/(1e5 - 1e4) + 3000);
            }
        }
    }
}

void fill_nTOF_cuts(){

    FIMG_tof_amp_cut_dedi_det1 = new TCutG("FIMG_tof_amp_cut_dedi_det1",3);
    FIMG_tof_amp_cut_dedi_det1->SetLineColor(2);
    FIMG_tof_amp_cut_dedi_det1->SetLineWidth(2);
    FIMG_tof_amp_cut_dedi_det1->SetVarX("x");
    FIMG_tof_amp_cut_dedi_det1->SetVarY("y");
    FIMG_tof_amp_cut_dedi_det1->SetPoint(0, 1e4, 70000.);
    FIMG_tof_amp_cut_dedi_det1->SetPoint(1, 1e4, 2350.);
    FIMG_tof_amp_cut_dedi_det1->SetPoint(2, 1e5, 500.);
    FIMG_tof_amp_cut_dedi_det1->SetPoint(3, 1e8, 500.);

    FIMG_tof_amp_cut_para_det1 = new TCutG("FIMG_tof_amp_cut_para_det1",3);
    FIMG_tof_amp_cut_para_det1->SetLineColor(2);
    FIMG_tof_amp_cut_para_det1->SetLineWidth(2);
    FIMG_tof_amp_cut_para_det1->SetVarX("x");
    FIMG_tof_amp_cut_para_det1->SetVarY("y");
    FIMG_tof_amp_cut_para_det1->SetPoint(0, 1e4, 70000.);
    FIMG_tof_amp_cut_para_det1->SetPoint(1, 1e4, 2350.);
    FIMG_tof_amp_cut_para_det1->SetPoint(2, 1e5, 500.);
    FIMG_tof_amp_cut_para_det1->SetPoint(3, 1e8, 500.);
    
    FIMG_tof_amp_cut_dedi_det2 = new TCutG("FIMG_tof_amp_cut_dedi_det2",3);
    FIMG_tof_amp_cut_dedi_det2->SetLineColor(2);
    FIMG_tof_amp_cut_dedi_det2->SetLineWidth(2);
    FIMG_tof_amp_cut_dedi_det2->SetVarX("x");
    FIMG_tof_amp_cut_dedi_det2->SetVarY("y");
    FIMG_tof_amp_cut_dedi_det2->SetPoint(0, 1e4, 70000.);
    FIMG_tof_amp_cut_dedi_det2->SetPoint(1, 1e4, 3000.);
    FIMG_tof_amp_cut_dedi_det2->SetPoint(2, 1e5, 500.);
    FIMG_tof_amp_cut_dedi_det2->SetPoint(3, 1e8, 500.);

    FIMG_tof_amp_cut_para_det2 = new TCutG("FIMG_tof_amp_cut_para_det2",3);
    FIMG_tof_amp_cut_para_det2->SetLineColor(2);
    FIMG_tof_amp_cut_para_det2->SetLineWidth(2);
    FIMG_tof_amp_cut_para_det2->SetVarX("x");
    FIMG_tof_amp_cut_para_det2->SetVarY("y");
    FIMG_tof_amp_cut_para_det2->SetPoint(0, 1e4, 70000.);
    FIMG_tof_amp_cut_para_det2->SetPoint(1, 1e4, 3000.);
    FIMG_tof_amp_cut_para_det2->SetPoint(2, 1e5, 500.);
    FIMG_tof_amp_cut_para_det2->SetPoint(3, 1e8, 500.);
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