/**
 * @file cutoffSelector.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-11-13
 */

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <functional>

#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TFrame.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TLine.h"
#include "TAxis.h"
#include "TColor.h"
#include "TAttMarker.h"
#include "TRandom3.h"

#include <vector>

const std::string target_name("al5"); ////bi1, al3, al5, al8, c1p2_ts, al5_ts, bi1p2_ts, cf_bottle, cf_bottle_rot, cf_bottle_rotBack, ar_bottle_full
const std::string target_out_name("none"); //none, none_ts, ar_bottle
const std::string target_name_title("Al (5 cm)");
//Bi (1 cm), Target Bi (1.2 cm), Al (3 cm), Al (5 cm), Target Al (5 cm), Al (8 cm), Target C (1.2 cm), Empty Bottle, Empty Bottle Rotated
//Argon Bottle

TH2D* PTBC_tof_amp_fOut_total = 0;
TH2D* PTBC_tof_amp_fOut_det2 = 0;
TH2D* PTBC_tof_amp_fOut_det3 = 0;
TH2D* PTBC_tof_amp_fOut_det4 = 0;
TH2D* PTBC_tof_amp_fOut_det5 = 0;
TH2D* PTBC_tof_amp_fOut_det6 = 0;
TH2D* PTBC_tof_amp_fOut_det7 = 0;
// std::reference_wrapper<TH2D*> PTBC_tof_amp_fOut_dets[] = {PTBC_tof_amp_fOut_det2, PTBC_tof_amp_fOut_det3, PTBC_tof_amp_fOut_det4, PTBC_tof_amp_fOut_det5, PTBC_tof_amp_fOut_det6, PTBC_tof_amp_fOut_det7};

TH2D* PTBC_tof_amp_fIn_total = 0;
TH2D* PTBC_tof_amp_fIn_det2 = 0;
TH2D* PTBC_tof_amp_fIn_det3 = 0;
TH2D* PTBC_tof_amp_fIn_det4 = 0;
TH2D* PTBC_tof_amp_fIn_det5 = 0;
TH2D* PTBC_tof_amp_fIn_det6 = 0;
TH2D* PTBC_tof_amp_fIn_det7 = 0;
// std::reference_wrapper<TH2D*> PTBC_tof_amp_fIn_dets[] = {PTBC_tof_amp_fIn_det2, PTBC_tof_amp_fIn_det3, PTBC_tof_amp_fIn_det4, PTBC_tof_amp_fIn_det5, PTBC_tof_amp_fIn_det6, PTBC_tof_amp_fIn_det7};

// TH2D* PTBC_tof_amp_fOut_det1_cutoff = 0;
// TH2D* PTBC_tof_amp_fOut_det2_cutoff = 0;
// TH2D* PTBC_tof_amp_fIn_det1_cutoff = 0;
// TH2D* PTBC_tof_amp_fIn_det2_cutoff = 0;

//dedicated and parasitic pulses plots
TH2D* PTBC_tof_amp_fIn_dedicated = 0;
TH2D* PTBC_tof_amp_fIn_parasitic = 0;

TH2D* PTBC_tof_amp_fOut_dedicated = 0;
TH2D* PTBC_tof_amp_fOut_parasitic = 0;

//cut plots
TCutG* PTBC_tof_amp_cut;

// Run Vectors
std::vector<Int_t> filter_in_runs;
std::vector<Int_t> filter_out_runs;

//beam off
// TH2D* tof_amp_beam_off_PTBC = 0;
// TH2D* tof_amp_beam_off_FIMG = 0;

//Cut off
// auto cutoff_amp = new TGraph();
// auto cutoff_tof = new TGraph();

Double_t flight_path_length_PTBC = 182.65 - 0.41; //m
Double_t flight_path_length_FIMG = 183.5 - 0.41; //m
Double_t neutron_mass = 939.56542052; //in MeV
Double_t speed_of_light = 299792458.0; //in m/s
Double_t delT_pkup_ptbc = 660.0; //in ns
Double_t delT_pkup_fimg = 630.0; //in ns

Double_t tof_min = 1e2;
Double_t tof_max = 1e8;
Double_t amp_min = 0.;
Double_t amp_max = 70000.;

Double_t t_gamma_PTBC = (flight_path_length_PTBC / speed_of_light) * 1e9; //converting into ns
Double_t t_gamma_FIMG = (flight_path_length_FIMG / speed_of_light) * 1e9; //converting into ns

//PTBC cuts
Double_t t[6][2];
Double_t a[6][2];

void fillCutsPTBC(){
    t[0][0] = 780.0;
    a[0][0] = 7700.0;

    t[0][1] = 1090.0;
    a[0][1] = 5451.0;

    t[1][0] = 1090.0;
    a[1][0] = 5451.0;

    t[1][1] = 2605.0;
    a[1][1] = 7167.0;

    t[2][0] = 2605.0;
    a[2][0] = 9500.0;

    t[2][1] = 2856.0;
    a[2][1] = 9500.0;

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

void fillCutPlot(){
    
    PTBC_tof_amp_cut->SetVarX("x");
    PTBC_tof_amp_cut->SetVarY("y");
    PTBC_tof_amp_cut->SetPoint(0, t[0][0], amp_max);
    PTBC_tof_amp_cut->SetPoint(1, t[0][0], a[0][0]);
    PTBC_tof_amp_cut->SetPoint(2, t[0][1], a[0][1]);
    PTBC_tof_amp_cut->SetPoint(3, t[1][1], a[1][1]);
    PTBC_tof_amp_cut->SetPoint(4, t[2][0], a[2][0]);
    PTBC_tof_amp_cut->SetPoint(5, t[2][1], a[2][1]);
    PTBC_tof_amp_cut->SetPoint(6, t[3][0], a[3][0]);
    PTBC_tof_amp_cut->SetPoint(7, t[3][1], a[3][1]);
    PTBC_tof_amp_cut->SetPoint(8, t[4][1], a[4][1]);
    PTBC_tof_amp_cut->SetPoint(9, t[5][1], a[5][1]);
}

Double_t yOnTheCutLine(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t x3){
    return ((y2 - y1)*(x3 - x1)/(x2 - x1) + y1);
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
std::vector<Int_t> f_out_t_arBot_ts = {117665, 117680, 117684, 117685, 117686, 117687, 117688, 117692, 117693, 117710, 117711, 117712, 117713, 117714, 117716, 117717, 117718, 117719, 117724, 117725, 117726, 117727, 117728, 117729, 117730, 117731, 117732, 117733, 117734, 117735, 117736, 117737, 117738};

//Argon Bottle (Argon Out)
std::vector<Int_t> f_out_t_arBotEmpty_ts = {117746, 117747, 117748, 117749, 117750, 117751, 117752, 117753, 117754, 117755, 117756, 117757, 117761, 117762, 117763, 117764, 117765, 117766}; //

void fillRuns(){
    // Bi 1 cm filter
    if(!target_name.compare("bi1")){
        cout << "Setting up Bi (1 cm) run list" << endl;
        // filter_in_runs.reserve( c_au_f_bi_t_out.size() + c_out_f_bi_t_out.size() + c_pb_f_bi_t_out.size() + c_ta_f_bi_t_out.size() );
        filter_in_runs.insert( filter_in_runs.end(), c_au_f_bi_t_out.begin(), c_au_f_bi_t_out.end() );
        filter_in_runs.insert( filter_in_runs.end(), c_out_f_bi_t_out.begin(), c_out_f_bi_t_out.end() );
        filter_in_runs.insert( filter_in_runs.end(), c_pb_f_bi_t_out.begin(), c_pb_f_bi_t_out.end() );
        filter_in_runs.insert( filter_in_runs.end(), c_ta_f_bi_t_out.begin(), c_ta_f_bi_t_out.end() );
    }

    // Al 8 cm filter
    if(!target_name.compare("al8")){
        cout << "Setting up Al (8 cm) run list" << endl;
        // filter_in_runs.reserve( c_ta_f_al8_t_out.size() );
        filter_in_runs.insert( filter_in_runs.end(), c_ta_f_al8_t_out.begin(), c_ta_f_al8_t_out.end() );
    }

    // Al 5 cm filter
    if(!target_name.compare("al5")){
        cout << "Setting up Al (5 cm) run list" << endl;
        // filter_in_runs.reserve( c_ta_f_al5_t_out.size() + c_out_f_al5_t_out.size() + c_pb_f_al5_t_out.size() + c_c_f_al5_t_out.size() + c_au_f_al5_t_out.size() ); //
        filter_in_runs.insert( filter_in_runs.end(), c_ta_f_al5_t_out.begin(), c_ta_f_al5_t_out.end() );
        filter_in_runs.insert( filter_in_runs.end(), c_out_f_al5_t_out.begin(), c_out_f_al5_t_out.end() );
        filter_in_runs.insert( filter_in_runs.end(), c_pb_f_al5_t_out.begin(), c_pb_f_al5_t_out.end() );
        filter_in_runs.insert( filter_in_runs.end(), c_c_f_al5_t_out.begin(), c_c_f_al5_t_out.end() );
        filter_in_runs.insert( filter_in_runs.end(), c_au_f_al5_t_out.begin(), c_au_f_al5_t_out.end() );
    }

    // Al 3 cm filter
    if(!target_name.compare("al3")){
        cout << "Setting up Al (3 cm) run list" << endl;
        // filter_in_runs.reserve( c_ta_f_al3_t_out.size() + c_pb_f_al3_t_out.size() );
        filter_in_runs.insert( filter_in_runs.end(), c_ta_f_al3_t_out.begin(), c_ta_f_al3_t_out.end() );
        filter_in_runs.insert( filter_in_runs.end(), c_pb_f_al3_t_out.begin(), c_pb_f_al3_t_out.end() );
    }

    //--------------------------------------------------

    if(!target_name.compare("bi1p2_ts")){
        cout << "Setting up Bi (1.2 cm) target run list" << endl;
        filter_in_runs.insert( filter_in_runs.end(), f_out_t_bi_ts.begin(), f_out_t_bi_ts.end() );
    }

    if(!target_name.compare("c1p2_ts")){
        cout << "Setting up C (1.2 cm) target run list" << endl;
        filter_in_runs.insert( filter_in_runs.end(), f_out_t_c_ts.begin(), f_out_t_c_ts.end() );
    }

    if(!target_name.compare("al5_ts")){
        cout << "Setting up Al (5 cm) target run list" << endl;
        filter_in_runs.insert( filter_in_runs.end(), f_out_t_al5_ts.begin(), f_out_t_al5_ts.end() );
    }

    if(!target_name.compare("cf_bottle")){
        cout << "Setting up Empty Bottle target run list" << endl;
        filter_in_runs.insert( filter_in_runs.end(), f_out_t_emptyBot_ts.begin(), f_out_t_emptyBot_ts.end() );
    }

    if(!target_name.compare("cf_bottle_rot")){
        cout << "Setting up Empty Bottle Rotated target run list" << endl;
        filter_in_runs.insert( filter_in_runs.end(), f_out_t_emptyBotRot_ts.begin(), f_out_t_emptyBotRot_ts.end() );
    }

    if(!target_name.compare("cf_bottle_rotBack")){
        cout << "Setting up Empty Bottle Rotated Back target run list" << endl;
        filter_in_runs.insert( filter_in_runs.end(), f_out_t_emptyBotRotBack_ts.begin(), f_out_t_emptyBotRotBack_ts.end() );
    }

    if(!target_name.compare("ar_bottle_full")){
        cout << "Setting up Argon Bottle target run list" << endl;
        filter_in_runs.insert( filter_in_runs.end(), f_out_t_arBot_ts.begin(), f_out_t_arBot_ts.end() );
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
    // filter_out_runs.insert( filter_out_runs.end(), c_au_f_out_t_out.begin(), c_au_f_out_t_out.end() );
    // filter_out_runs.insert( filter_out_runs.end(), c_atic_f_out_t_out.begin(), c_atic_f_out_t_out.end() );

    if(!target_out_name.compare("empty_bottle")){
        filter_out_runs.insert( filter_out_runs.end(), f_out_t_emptyBot_ts.begin(), f_out_t_emptyBot_ts.end() );
        filter_out_runs.insert( filter_out_runs.end(), f_out_t_emptyBotRot_ts.begin(), f_out_t_emptyBotRot_ts.end() );
        filter_out_runs.insert( filter_out_runs.end(), f_out_t_emptyBotRotBack_ts.begin(), f_out_t_emptyBotRotBack_ts.end() );
    }

    if(!target_out_name.compare("ar_bottle")){
        filter_out_runs.insert( filter_out_runs.end(), f_out_t_arBotEmpty_ts.begin(), f_out_t_arBotEmpty_ts.end() );
    }

    if(!target_out_name.compare("none")){
        filter_out_runs.insert( filter_out_runs.end(), c_out_f_out_t_out.begin(), c_out_f_out_t_out.end() );
        filter_out_runs.insert( filter_out_runs.end(), c_pb_f_out_t_out.begin(), c_pb_f_out_t_out.end() );
        filter_out_runs.insert( filter_out_runs.end(), c_ta_f_out_t_out.begin(), c_ta_f_out_t_out.end() );
        filter_out_runs.insert( filter_out_runs.end(), c_c_f_out_t_out.begin(), c_c_f_out_t_out.end() );
    }
}

void FilterIn(){

    for (int i = 0; i < filter_in_runs.size(); i++)
    {
        TFile *file_ntof = TFile::Open(Form("/eos/experiment/ntof/processing/official/done/run%d.root", filter_in_runs.at(i)),"read");
        
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

        //PTBC ---------------------------------------------
        TTree* PTBC;
        Double_t tof_PTBC = 0; //tof is in ns
        Float_t amp_PTBC = 0;
        Int_t BunchNumber_PTBC = 0;
        Int_t det_num_PTBC = 0;
        Float_t PulseIntensity_PTBC = 0;

        file_ntof->GetObject("PTBC", PTBC);
        PTBC->SetBranchAddress("BunchNumber", &BunchNumber_PTBC);
        PTBC->SetBranchAddress("PulseIntensity", &PulseIntensity_PTBC);
        PTBC->SetBranchAddress("tof", &tof_PTBC);
        PTBC->SetBranchAddress("amp", &amp_PTBC);
        PTBC->SetBranchAddress("detn", &det_num_PTBC);

        Long64_t Events_PTBC = PTBC->GetEntriesFast();
        std::cout << "Number of entries - PTBC = " << Events_PTBC << std::endl;
        
        int CurrentBunchNum = 0;

        for (int j = 0; j < Events_PTBC; j++)
        {
            PTBC->GetEntry(j);

            Double_t t_pkup = BNum_tpkup_map[BunchNumber_PTBC];
            Double_t corrected_tof = tof_PTBC - t_pkup + delT_pkup_ptbc + t_gamma_PTBC;

            if (CurrentBunchNum != BunchNumber_PTBC)
            {
                CurrentBunchNum = BunchNumber_PTBC;
            }

            //Filling the histograms
            // for (int k = 0; k < 6; k++)
            // {
            //     if (corrected_tof >= t[k][0] && corrected_tof < t[k][1])
            //     {
            //         if ( (Double_t) amp_PTBC > yOnTheCutLinePTBC(t[k][0], a[k][0], t[k][1], a[k][1], corrected_tof) )
            //         {
            //             energy_hist_PTBC->Fill( TOFToEnergy(corrected_tof * 1e-9) );
            //             break;    
            //         }
            //     }
            // }

            PTBC_tof_amp_fIn_total->Fill(corrected_tof, (Double_t) amp_PTBC);

            if (PulseIntensity_PTBC > 6e12)
            {
                PTBC_tof_amp_fIn_dedicated->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (PulseIntensity_PTBC <= 6e12)
            {
                PTBC_tof_amp_fIn_parasitic->Fill(corrected_tof, (Double_t) amp_PTBC);
            }

            if (det_num_PTBC == 2) {
                PTBC_tof_amp_fIn_det2->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (det_num_PTBC == 3) {
                PTBC_tof_amp_fIn_det3->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (det_num_PTBC == 4) {
                PTBC_tof_amp_fIn_det4->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (det_num_PTBC == 5) {
                PTBC_tof_amp_fIn_det5->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (det_num_PTBC == 6) {
                PTBC_tof_amp_fIn_det6->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (det_num_PTBC == 7) {
                PTBC_tof_amp_fIn_det7->Fill(corrected_tof, (Double_t) amp_PTBC);
            }

            // if (det_num_PTBC != 1 && det_num_PTBC != 8)
            // {
            //     PTBC_tof_amp_fIn_dets[det_num_PTBC-2]->Fill(corrected_tof, (Double_t) amp_PTBC);
            // }
        }

        file_ntof->Close();
    }
}

void FilterOut(){

    for (int i = 0; i < filter_out_runs.size(); i++)
    {
        TFile *file_ntof = TFile::Open(Form("/eos/experiment/ntof/processing/official/done/run%d.root", filter_out_runs.at(i)),"read");
        
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

        //PTBC ---------------------------------------------
        TTree* PTBC;
        Double_t tof_PTBC = 0; //tof is in ns
        Float_t amp_PTBC = 0;
        Int_t BunchNumber_PTBC = 0;
        Int_t det_num_PTBC = 0;
        Float_t PulseIntensity_PTBC = 0;

        file_ntof->GetObject("PTBC", PTBC);
        PTBC->SetBranchAddress("BunchNumber", &BunchNumber_PTBC);
        PTBC->SetBranchAddress("PulseIntensity", &PulseIntensity_PTBC);
        PTBC->SetBranchAddress("tof", &tof_PTBC);
        PTBC->SetBranchAddress("amp", &amp_PTBC);
        PTBC->SetBranchAddress("detn", &det_num_PTBC);

        Long64_t Events_PTBC = PTBC->GetEntriesFast();
        std::cout << "Number of entries - PTBC = " << Events_PTBC << std::endl;
        
        int CurrentBunchNum = 0;

        for (int j = 0; j < Events_PTBC; j++)
        {
            PTBC->GetEntry(j);

            if (det_num_PTBC == 1 || det_num_PTBC == 8) {
                continue;
            }

            Double_t t_pkup = BNum_tpkup_map[BunchNumber_PTBC];
            Double_t corrected_tof = tof_PTBC - t_pkup + delT_pkup_ptbc + t_gamma_PTBC;

            if (CurrentBunchNum != BunchNumber_PTBC)
            {
                CurrentBunchNum = BunchNumber_PTBC;
            }

            //Filling the histograms
            // for (int k = 0; k < 6; k++)
            // {
            //     if (corrected_tof >= t[k][0] && corrected_tof < t[k][1])
            //     {
            //         if ( (Double_t) amp_PTBC > yOnTheCutLinePTBC(t[k][0], a[k][0], t[k][1], a[k][1], corrected_tof) )
            //         {
            //             energy_hist_PTBC->Fill( TOFToEnergy(corrected_tof * 1e-9) );
            //             break;    
            //         }
            //     }
            // }

            PTBC_tof_amp_fOut_total->Fill(corrected_tof, (Double_t) amp_PTBC);

            if (PulseIntensity_PTBC > 6e12)
            {
                PTBC_tof_amp_fOut_dedicated->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (PulseIntensity_PTBC <= 6e12)
            {
                PTBC_tof_amp_fOut_parasitic->Fill(corrected_tof, (Double_t) amp_PTBC);
            }

            if (det_num_PTBC == 2) {
                PTBC_tof_amp_fOut_det2->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (det_num_PTBC == 3) {
                PTBC_tof_amp_fOut_det3->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (det_num_PTBC == 4) {
                PTBC_tof_amp_fOut_det4->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (det_num_PTBC == 5) {
                PTBC_tof_amp_fOut_det5->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (det_num_PTBC == 6) {
                PTBC_tof_amp_fOut_det6->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (det_num_PTBC == 7) {
                PTBC_tof_amp_fOut_det7->Fill(corrected_tof, (Double_t) amp_PTBC);
            }

            // if (det_num_PTBC != 1 && det_num_PTBC != 8)
            // {
            //     PTBC_tof_amp_fOut_dets[det_num_PTBC-2]->Fill(corrected_tof, (Double_t) amp_PTBC);
            // }
        }

        file_ntof->Close();
    }
}

void StoreHist(){
    
    TFile *f = new TFile(Form("../rootFiles/cutoffAnalysis_PTBC_%s.root", target_name.c_str()),"recreate");

    // PTBC_tof_amp_in->Write();
    // PTBC_tof_amp_out->Write();
    PTBC_tof_amp_fOut_total->Write();
    PTBC_tof_amp_fOut_det2->Write();
    PTBC_tof_amp_fOut_det3->Write();
    PTBC_tof_amp_fOut_det4->Write();
    PTBC_tof_amp_fOut_det5->Write();
    PTBC_tof_amp_fOut_det6->Write();
    PTBC_tof_amp_fOut_det7->Write();

    PTBC_tof_amp_fIn_total->Write();
    PTBC_tof_amp_fIn_det2->Write();
    PTBC_tof_amp_fIn_det3->Write();
    PTBC_tof_amp_fIn_det4->Write();
    PTBC_tof_amp_fIn_det5->Write();
    PTBC_tof_amp_fIn_det6->Write();
    PTBC_tof_amp_fIn_det7->Write();

    // PTBC_tof_amp_fIn_det1_cutoff->Write();
    // PTBC_tof_amp_fIn_det2_cutoff->Write();
    // PTBC_tof_amp_fOut_det1_cutoff->Write();
    // PTBC_tof_amp_fOut_det2_cutoff->Write();

    PTBC_tof_amp_fIn_dedicated->Write();
    PTBC_tof_amp_fIn_parasitic->Write();
    PTBC_tof_amp_fOut_dedicated->Write();
    PTBC_tof_amp_fOut_parasitic->Write();

    PTBC_tof_amp_cut->Write();

    f->Close();
}

void cutoffAnalysis_PTBC(){

    gStyle->SetStatX(0.27);
    gStyle->SetStatY(0.9);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.17);

    fillCutsPTBC();
    fillRuns();

    PTBC_tof_amp_cut = new TCutG("PTBC_tof_amp_cut",9);
    fillCutPlot();

    //Calculating TOF (x) bin edges
    int bins_per_decade = 1000;
    int Num_decades = 6;
    int num_bins_tof = bins_per_decade * Num_decades;
    Double_t bin_edges_tof[num_bins_tof+1];
    Double_t step_tof = ((Double_t) 1.0/(Double_t) bins_per_decade);
    for(int i = 0; i < num_bins_tof+1; i++)
    {
        Double_t base = 10.;
        Double_t exponent = (step_tof * (Double_t) i) + 2.;
        bin_edges_tof[i] = (Double_t) std::pow(base, exponent);
    }

    // Calculating amplitude (y) bin edges
    int num_bins_amp = 500;
    Double_t bin_edges_amp[num_bins_amp+1];
    Double_t step_amp = (Double_t) ((amp_max-amp_min)/num_bins_amp);
    // std::cout << "step_amp = " << step_amp << std::endl;
    for(int i = 0; i < num_bins_amp+1; i++)
    {
        bin_edges_amp[i] = step_amp * (Double_t) i;
    }

    //Filter Out
    PTBC_tof_amp_fOut_total = new TH2D("PTBC_tof_amp_fOut_total","ToF vs Amplitude Hist - PTBC All Det - Filter Out",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fOut_det2 = new TH2D("PTBC_tof_amp_fOut_det2","ToF vs Amplitude Hist - PTBC Det 2 - Filter Out",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fOut_det3 = new TH2D("PTBC_tof_amp_fOut_det3","ToF vs Amplitude Hist - PTBC Det 3 - Filter Out",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fOut_det4 = new TH2D("PTBC_tof_amp_fOut_det4","ToF vs Amplitude Hist - PTBC Det 4 - Filter Out",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fOut_det5 = new TH2D("PTBC_tof_amp_fOut_det5","ToF vs Amplitude Hist - PTBC Det 5 - Filter Out",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fOut_det6 = new TH2D("PTBC_tof_amp_fOut_det6","ToF vs Amplitude Hist - PTBC Det 6 - Filter Out",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fOut_det7 = new TH2D("PTBC_tof_amp_fOut_det7","ToF vs Amplitude Hist - PTBC Det 7 - Filter Out",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    
    //Filter In
    PTBC_tof_amp_fIn_total = new TH2D("PTBC_tof_amp_fIn_total",Form("ToF vs Amplitude Hist - PTBC All Det - %s", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fIn_det2 = new TH2D("PTBC_tof_amp_fIn_det2",Form("ToF vs Amplitude Hist - PTBC Det 2 - %s", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fIn_det3 = new TH2D("PTBC_tof_amp_fIn_det3",Form("ToF vs Amplitude Hist - PTBC Det 3 - %s", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fIn_det4 = new TH2D("PTBC_tof_amp_fIn_det4",Form("ToF vs Amplitude Hist - PTBC Det 4 - %s", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fIn_det5 = new TH2D("PTBC_tof_amp_fIn_det5",Form("ToF vs Amplitude Hist - PTBC Det 5 - %s", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fIn_det6 = new TH2D("PTBC_tof_amp_fIn_det6",Form("ToF vs Amplitude Hist - PTBC Det 6 - %s", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fIn_det7 = new TH2D("PTBC_tof_amp_fIn_det7",Form("ToF vs Amplitude Hist - PTBC Det 7 - %s", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);

    // PTBC_tof_amp_fOut_det1_cutoff = new TH2D("PTBC_tof_amp_fOut_det1_cutoff","ToF vs Amp Hist - PTBC Det 1 - Filter Out Cutoff",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    // PTBC_tof_amp_fOut_det2_cutoff = new TH2D("PTBC_tof_amp_fOut_det2_cutoff","ToF vs Amp Hist - PTBC Det 2 - Filter Out Cutoff",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);

    // PTBC_tof_amp_fIn_det1_cutoff = new TH2D("PTBC_tof_amp_fIn_det1_cutoff","ToF vs Amp Hist - PTBC Det 1 - Al (5cm) Filter Cutoff",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    // PTBC_tof_amp_fIn_det2_cutoff = new TH2D("PTBC_tof_amp_fIn_det2_cutoff","ToF vs Amp Hist - PTBC Det 2 - Al (5cm) Filter Cutoff",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);

    PTBC_tof_amp_fIn_dedicated = new TH2D("PTBC_tof_amp_fIn_dedicated",Form("ToF vs Amp Hist - %s - Dedicated", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fIn_parasitic = new TH2D("PTBC_tof_amp_fIn_parasitic",Form("ToF vs Amp Hist - %s - Parasitic", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fOut_dedicated = new TH2D("PTBC_tof_amp_fOut_dedicated","ToF vs Amp Hist - Filter Out - Dedicated",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fOut_parasitic = new TH2D("PTBC_tof_amp_fOut_parasitic","ToF vs Amp Hist - Filter Out - Parasitic",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);

    
    FilterIn();
    FilterOut();

    StoreHist();
}