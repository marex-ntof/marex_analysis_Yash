/**
 * @file cutoffAnalysis_FIMG.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-11-28
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
const std::string target_name("c1p2_ts"); //bi1, al3, al5, al8, c1p2_ts, al5_ts, bi1p2_ts, cf_bottle, cf_bottle_rot, cf_bottle_rotBack, ar_bottle_full
const std::string target_out_name("none_ts"); //none, none_ts, ar_bottle
// const std::string mode("run");
// bi1, al3, al5, al8 - none

//Root file
TFile *outputRootFile = 0;

TH1D* trans_loose_cut_det_1 = 0;
TH1D* trans_mid_cut_det_1 = 0;
TH1D* trans_tight_cut_det_1 = 0;
TH1D* trans_loose_cut_det_2 = 0;
TH1D* trans_mid_cut_det_2 = 0;
TH1D* trans_tight_cut_det_2 = 0;

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

Double_t fimgCutFunction(Double_t x, Int_t det_num, std::string cut_type){

    Double_t A = 0;
    Double_t B = 0;
    Double_t C = 0;

    if(det_num == 1) {
        // A / ( TMath::Log(x) + B )
        if (!cut_type.compare("loose"))
        {
            A = 4000;
            B = -TMath::Log(7750);
            C = 600;
        }
        if (!cut_type.compare("mid"))
        {
            A = 4500;
            B = -TMath::Log(7750);
            C = 600;
        }
        if (!cut_type.compare("tight"))
        {
            A = 7000;
            B = -TMath::Log(7750);
            C = 600;
        }
        
        return (A / (TMath::Log(x) + B));
    }
    if(det_num == 2) {
        // A / ( TMath::Log(x) + B )
        if (!cut_type.compare("loose"))
        {
            A = 3000;
            B = -TMath::Log(7080);
            C = 400;
        }
        if (!cut_type.compare("mid"))
        {
            A = 5000;
            B = -TMath::Log(7080);
            C = 400;
        }
        if (!cut_type.compare("tight"))
        {
            A = 6500;
            B = -TMath::Log(7080);
            C = 400;
        }
        
        return ((A / (TMath::Log(x) + B)) + C);
    }
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

Double_t GetNormFactor(std::vector<Int_t> run_list){
    
    Double_t NormFactor = 0; //Integral of the pulse intensity

    for (int i = 0; i < run_list.size(); i++)
    {
        TFile *file_ntof = TFile::Open(Form("/eos/experiment/ntof/processing/official/done/run%d.root", run_list.at(i)),"read");
        cout << "Run Number = " << run_list.at(i) << endl;

        //FIMG ---------------------------------------------
        TTree* FIMG;
        Int_t BunchNumber_FIMG = 0;
        Float_t PulseIntensity = 0;

        file_ntof->GetObject("FIMG", FIMG);
        FIMG->SetBranchAddress("BunchNumber", &BunchNumber_FIMG);
        FIMG->SetBranchAddress("PulseIntensity", &PulseIntensity);

        int CurrentBunchNum = 0;

        Long64_t Events_FIMG = FIMG->GetEntriesFast();

        for (int j = 0; j < Events_FIMG; j++)
        {
            FIMG->GetEntry(j);

            if (CurrentBunchNum != BunchNumber_FIMG)
            {
                CurrentBunchNum = BunchNumber_FIMG;
                NormFactor += (Double_t) PulseIntensity;
            }
        }

        file_ntof->Close();
    }

    return NormFactor;
}

void fillEnergyHist(std::vector<Int_t> run_list, std::string cut_type, TH1D* energy_hist_det1, TH1D* energy_hist_det2){

    for (int i = 0; i < run_list.size(); i++)
    {
        TFile *file_ntof = TFile::Open(Form("/eos/experiment/ntof/processing/official/done/run%d.root", run_list.at(i)),"read");
        cout << "Run Number = " << run_list.at(i) << endl;

        //PKUP ---------------------------------------------
        TTree* PKUP;
        Int_t BunchNumber_PKUP = 0;
        Double_t tpkup = 0;

        file_ntof->GetObject("PKUP", PKUP);
        PKUP->SetBranchAddress("BunchNumber", &BunchNumber_PKUP);
        PKUP->SetBranchAddress("tflash", &tpkup);

        std::map<Int_t, Double_t> BNum_tpkup_map;
        Long64_t Events_PKUP = PKUP->GetEntriesFast();

        for (int j = 0; j < Events_PKUP; j++)
        {
            PKUP->GetEntry(j);
            BNum_tpkup_map.emplace(BunchNumber_PKUP, tpkup);
        }

        //FIMG ---------------------------------------------
        TTree* FIMG;
        Double_t tof_FIMG = 0; //tof is in ns
        Float_t amp = 0;
        Float_t area_0 = 0;
        Int_t det_num = 0;
        Int_t BunchNumber_FIMG = 0;

        file_ntof->GetObject("FIMG", FIMG);
        FIMG->SetBranchAddress("BunchNumber", &BunchNumber_FIMG);
        FIMG->SetBranchAddress("tof", &tof_FIMG);
        FIMG->SetBranchAddress("amp", &amp);
        FIMG->SetBranchAddress("area_0", &area_0);
        FIMG->SetBranchAddress("detn", &det_num);

        Long64_t Events_FIMG = FIMG->GetEntriesFast();
        std::cout << "Number of entries - FIMG = " << Events_FIMG << std::endl;

        for (int j = 0; j < Events_FIMG; j++)
        {
            FIMG->GetEntry(j);

            Double_t t_pkup = BNum_tpkup_map[BunchNumber_FIMG];
            Double_t corrected_tof = tof_FIMG - t_pkup + delT_pkup_fimg + t_gamma_FIMG;

            //Filling the histograms after cuts
            if (det_num == 1)
            {
                if ((Double_t) amp >= fimgCutFunction(corrected_tof, det_num, cut_type))
                {
                    energy_hist_det1->Fill( TOFToEnergy(corrected_tof * 1e-9) );
                }
            }

            if (det_num == 2)
            {
                if ((Double_t) amp >= fimgCutFunction(corrected_tof, det_num, cut_type))
                {
                    energy_hist_det2->Fill( TOFToEnergy(corrected_tof * 1e-9) );
                }
            }
        }

        file_ntof->Close();
    }
}

void calcTransmission(TH1D* e_hist_in, TH1D* e_hist_out, TH1D* trans_hist, Double_t Qout_Qin){
    Int_t num_bins_e = e_hist_in->GetNbinsX();
    for (int i = 0; i < num_bins_e; i++)
    {
        Double_t bin_content_in = e_hist_in->GetBinContent(i+1);
        Double_t bin_content_out = e_hist_out->GetBinContent(i+1);
        if (bin_content_in == 0. || bin_content_out == 0.)
        {
            trans_hist->SetBinContent(i+1, 0.);
            trans_hist->SetBinError(i+1, 0.);
        } else {
            Double_t transmission = (bin_content_in * Qout_Qin)/bin_content_out;
            Double_t bin_unc = transmission * std::sqrt( (1./bin_content_in) + (1./bin_content_out) );
            trans_hist->SetBinContent(i+1, transmission);
            trans_hist->SetBinError(i+1, bin_unc);
        }
    }
}

void fillRunHists(Int_t num_bins_e, Double_t bin_edges_e[]){

    TH1D* e_loose_cut_det_1_in = 0;
    TH1D* e_mid_cut_det_1_in = 0;
    TH1D* e_tight_cut_det_1_in = 0;
    TH1D* e_loose_cut_det_2_in = 0;
    TH1D* e_mid_cut_det_2_in = 0;
    TH1D* e_tight_cut_det_2_in = 0;

    TH1D* e_loose_cut_det_1_out = 0;
    TH1D* e_mid_cut_det_1_out = 0;
    TH1D* e_tight_cut_det_1_out = 0;
    TH1D* e_loose_cut_det_2_out = 0;
    TH1D* e_mid_cut_det_2_out = 0;
    TH1D* e_tight_cut_det_2_out = 0;

    //Target In Hists
    e_loose_cut_det_1_in = new TH1D("e_loose_cut_det_1_in","Energy Hist",num_bins_e,bin_edges_e);
    e_mid_cut_det_1_in = new TH1D("e_mid_cut_det_1_in","Energy Hist",num_bins_e,bin_edges_e);
    e_tight_cut_det_1_in = new TH1D("e_tight_cut_det_1_in","Energy Hist",num_bins_e,bin_edges_e);
    e_loose_cut_det_2_in = new TH1D("e_loose_cut_det_2_in","Energy Hist",num_bins_e,bin_edges_e);
    e_mid_cut_det_2_in = new TH1D("e_mid_cut_det_2_in","Energy Hist",num_bins_e,bin_edges_e);
    e_tight_cut_det_2_in = new TH1D("e_tight_cut_det_2_in","Energy Hist",num_bins_e,bin_edges_e);
    //Target Out Hists
    e_loose_cut_det_1_out = new TH1D("e_loose_cut_det_1_out","Energy Hist",num_bins_e,bin_edges_e);
    e_mid_cut_det_1_out = new TH1D("e_mid_cut_det_1_out","Energy Hist",num_bins_e,bin_edges_e);
    e_tight_cut_det_1_out = new TH1D("e_tight_cut_det_1_out","Energy Hist",num_bins_e,bin_edges_e);
    e_loose_cut_det_2_out = new TH1D("e_loose_cut_det_2_out","Energy Hist",num_bins_e,bin_edges_e);
    e_mid_cut_det_2_out = new TH1D("e_mid_cut_det_2_out","Energy Hist",num_bins_e,bin_edges_e);
    e_tight_cut_det_2_out = new TH1D("e_tight_cut_det_2_out","Energy Hist",num_bins_e,bin_edges_e);
    
    //Norm Factors
    Double_t norm_factor_target_in = GetNormFactor(ts_target_in_runs);
    Double_t norm_factor_target_out = 0;
    if(!target_out_name.compare("ar_bottle")) {
        norm_factor_target_out = GetNormFactor(ts_target_out_runs);
    } else if(!target_out_name.compare("none_ts")) {
        norm_factor_target_out = GetNormFactor(f_out_t_out_ts);
    } else if(!target_out_name.compare("none")) {
        norm_factor_target_out = GetNormFactor(f_out_t_out_runs);
    }

    Double_t Qout_Qin = norm_factor_target_out / norm_factor_target_in;
    
    //Filling the Histograms - Target In
    fillEnergyHist(ts_target_in_runs, "loose", e_loose_cut_det_1_in, e_loose_cut_det_2_in);
    fillEnergyHist(ts_target_in_runs, "mid", e_mid_cut_det_1_in, e_mid_cut_det_2_in);
    fillEnergyHist(ts_target_in_runs, "tight", e_tight_cut_det_1_in, e_tight_cut_det_2_in);

    //Filling the Histograms - Target Out
    if(!target_out_name.compare("ar_bottle")) {
        fillEnergyHist(ts_target_out_runs, "loose", e_loose_cut_det_1_out, e_loose_cut_det_2_out);
        fillEnergyHist(ts_target_out_runs, "mid", e_mid_cut_det_1_out, e_mid_cut_det_2_out);
        fillEnergyHist(ts_target_out_runs, "tight", e_tight_cut_det_1_out, e_tight_cut_det_2_out);
    } else if(!target_out_name.compare("none_ts")) {
        fillEnergyHist(f_out_t_out_ts, "loose", e_loose_cut_det_1_out, e_loose_cut_det_2_out);
        fillEnergyHist(f_out_t_out_ts, "mid", e_mid_cut_det_1_out, e_mid_cut_det_2_out);
        fillEnergyHist(f_out_t_out_ts, "tight", e_tight_cut_det_1_out, e_tight_cut_det_2_out);
    } else if(!target_out_name.compare("none")) {
        fillEnergyHist(f_out_t_out_runs, "loose", e_loose_cut_det_1_out, e_loose_cut_det_2_out);
        fillEnergyHist(f_out_t_out_runs, "mid", e_mid_cut_det_1_out, e_mid_cut_det_2_out);
        fillEnergyHist(f_out_t_out_runs, "tight", e_tight_cut_det_1_out, e_tight_cut_det_2_out);
    }

    //Transmission
    for (int i = 0; i < num_bins_e; i++)
    {
        //Loose
        calcTransmission(e_loose_cut_det_1_in, e_loose_cut_det_1_out, trans_loose_cut_det_1, Qout_Qin);
        calcTransmission(e_loose_cut_det_2_in, e_loose_cut_det_2_out, trans_loose_cut_det_2, Qout_Qin);

        //Mid
        calcTransmission(e_mid_cut_det_1_in, e_mid_cut_det_1_out, trans_mid_cut_det_1, Qout_Qin);
        calcTransmission(e_mid_cut_det_2_in, e_mid_cut_det_2_out, trans_mid_cut_det_2, Qout_Qin);

        //Tight
        calcTransmission(e_tight_cut_det_1_in, e_tight_cut_det_1_out, trans_tight_cut_det_1, Qout_Qin);
        calcTransmission(e_tight_cut_det_2_in, e_tight_cut_det_2_out, trans_tight_cut_det_2, Qout_Qin);
    }

    cout << Form("Total Protons %s = ", target_name.c_str()) << norm_factor_target_in << endl;
    cout << "Total Protons Target Out = " << norm_factor_target_out << endl;
}

void cutoffAnalysis_FIMG() {
    
    fillRuns();

    //Calculating Energy bin edges
    Double_t tof_min = 1e3; //ns
    Double_t tof_max = 1e8; //ns
    Double_t e_min = TOFToEnergy(tof_max * 1e-9); //converting into seconds
    Double_t e_max = TOFToEnergy(tof_min * 1e-9); //converting into seconds
    int min_power = FindDecadePower(e_min);
    int max_power = FindDecadePower(e_max);
    int num_decades_e = max_power - min_power;
    int num_bins_e = bins_per_decade * num_decades_e;
    Double_t bin_edges_e[num_bins_e+1];
    Double_t step_e = ((Double_t) 1.0/(Double_t) bins_per_decade);
    for(int i = 0; i < num_bins_e+1; i++)
    {
        Double_t base = 10.;
        Double_t exponent = (step_e * (Double_t) i) + (Double_t) min_power;
        bin_edges_e[i] = (Double_t) std::pow(base, exponent);
    }

    // cout << "Number of tof bins = " << num_bins_tof << endl;
    cout << "Number of e bins = " << num_bins_e << endl;

    //transmission histogram
    trans_loose_cut_det_1 = new TH1D("trans_loose_cut_det_1","Transmission Hist",num_bins_e,bin_edges_e);
    trans_mid_cut_det_1 = new TH1D("trans_mid_cut_det_1","Transmission Hist",num_bins_e,bin_edges_e);
    trans_tight_cut_det_1 = new TH1D("trans_tight_cut_det_1","Transmission Hist",num_bins_e,bin_edges_e);
    trans_loose_cut_det_2 = new TH1D("trans_loose_cut_det_2","Transmission Hist",num_bins_e,bin_edges_e);
    trans_mid_cut_det_2 = new TH1D("trans_mid_cut_det_2","Transmission Hist",num_bins_e,bin_edges_e);
    trans_tight_cut_det_2 = new TH1D("trans_tight_cut_det_2","Transmission Hist",num_bins_e,bin_edges_e);
    
    fillRunHists(num_bins_e,bin_edges_e);

    //Writing to the output file
    outputRootFile = new TFile(Form("../rootFiles/cutoffAnalysis_FIMG_%s.root", target_name.c_str()),"recreate");
    trans_loose_cut_det_1->Write();
    trans_mid_cut_det_1->Write();
    trans_tight_cut_det_1->Write();
    trans_loose_cut_det_2->Write();
    trans_mid_cut_det_2->Write();
    trans_tight_cut_det_2->Write();

    outputRootFile->Close();
    
}