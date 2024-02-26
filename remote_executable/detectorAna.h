/**
 * @file detectorAna.h
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

#include "parameters.h"
#include "tools.h"

//////// Run variables
const std::string target_name("al5"); //bi1, al3, al5, al8, c1p2_ts, al5_ts, bi1p2_ts, cf_bottle, cf_bottle_rot, cf_bottle_rotBack, ar_bottle_full
const std::string target_out_name("none"); //none, none_ts, ar_bottle
// const std::string mode("run");
// bi1, al3, al5, al8 - none

//Root file
TFile *outputRootFile = 0;

TH1D* energy_hist_target_in_PTB = 0;
TH1D* energy_hist_target_in_FIMG = 0;
TH1D* energy_hist_target_out_PTB = 0;
TH1D* energy_hist_target_out_FIMG = 0;

TH1D* tof_hist_target_in_PTB = 0;
TH1D* tof_hist_target_in_FIMG = 0;
TH1D* tof_hist_target_out_PTB = 0;
TH1D* tof_hist_target_out_FIMG = 0;

std::vector<Double_t> norm_factors;

Int_t bins_per_decade = 1000;

////////////////////////
//// vectors to store run numbers
std::vector<Int_t> f_al5_t_out_runs;
std::vector<Int_t> f_out_t_out_runs;

std::vector<Int_t> filter_in_runs;
std::vector<Int_t> filter_out_runs;

std::vector<Int_t> ts_target_in_runs;
std::vector<Int_t> ts_target_out_runs;

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