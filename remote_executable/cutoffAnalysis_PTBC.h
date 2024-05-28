/**
 * @file cutoffAnalysis_PTBC.h
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2024-01-24
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
#include "parameters.h"

const std::string target_name("bi1"); ////bi1, al3, al5, al8, c1p2_ts, al5_ts, bi1p2_ts, cf_bottle, cf_bottle_rot, cf_bottle_rotBack, ar_bottle_full, bi1sep17
//no targets - none, none_ts, ar_bottle
const std::string target_name_title("Bi (1 cm)");
//Bi (1 cm), Target Bi (1.2 cm), Al (3 cm), Al (5 cm), Target Al (5 cm), Al (8 cm), Target C (1.2 cm), Empty Bottle, Empty Bottle Rotated
//Argon Tank, Bi (1 cm) - Sep 17, No Target, SCUBA Tank

// std::reference_wrapper<TH2D*> PTBC_tof_amp_fOut_dets[] = {PTBC_tof_amp_fOut_det2, PTBC_tof_amp_fOut_det3, PTBC_tof_amp_fOut_det4, PTBC_tof_amp_fOut_det5, PTBC_tof_amp_fOut_det6, PTBC_tof_amp_fOut_det7};

TH2D* PTBC_tof_amp_total = 0;
//dedicated and parasitic pulses plots
TH2D* PTBC_tof_amp_total_dedicated = 0;
TH2D* PTBC_tof_amp_total_parasitic = 0;

TH2D* PTBC_tof_amp_det2 = 0;
TH2D* PTBC_tof_amp_det3 = 0;
TH2D* PTBC_tof_amp_det4 = 0;
TH2D* PTBC_tof_amp_det5 = 0;
TH2D* PTBC_tof_amp_det6 = 0;
TH2D* PTBC_tof_amp_det7 = 0;

TH2D* PTBC_ringing_det2 = 0;
TH2D* PTBC_ringing_det3 = 0;
TH2D* PTBC_ringing_det4 = 0;
TH2D* PTBC_ringing_det5 = 0;
TH2D* PTBC_ringing_det6 = 0;
TH2D* PTBC_ringing_det7 = 0;

// TH2D* PTBC_tof_amp_det2_afterCuts = 0;
// TH2D* PTBC_tof_amp_det3_afterCuts = 0;
// TH2D* PTBC_tof_amp_det4_afterCuts = 0;
// TH2D* PTBC_tof_amp_det5_afterCuts = 0;
// TH2D* PTBC_tof_amp_det6_afterCuts = 0;
// TH2D* PTBC_tof_amp_det7_afterCuts = 0;
// std::reference_wrapper<TH2D*> PTBC_tof_amp_dets[] = {PTBC_tof_amp_det2, PTBC_tof_amp_det3, PTBC_tof_amp_det4, PTBC_tof_amp_det5, PTBC_tof_amp_det6, PTBC_tof_amp_det7};

// TH2D* PTBC_tof_amp_fOut_det1_cutoff = 0;
// TH2D* PTBC_tof_amp_fOut_det2_cutoff = 0;
// TH2D* PTBC_tof_amp_det1_cutoff = 0;
// TH2D* PTBC_tof_amp_det2_cutoff = 0;

//cuts
TH1D* PTBC_cuts_det2;
TH1D* PTBC_cuts_det3;
TH1D* PTBC_cuts_det4;
TH1D* PTBC_cuts_det5;
TH1D* PTBC_cuts_det6;
TH1D* PTBC_cuts_det7;

// Run Vectors
std::vector<Int_t> filter_run_list;
std::vector<Int_t> det5_runs = {117386, 117387, 117388, 117389, 117390, 117391, 117439, 117440, 117441, 117442, 117443};;

//beam off
// TH2D* tof_amp_beam_off_PTBC = 0;
// TH2D* tof_amp_beam_off_FIMG = 0;

//Cut off
// auto cutoff_amp = new TGraph();
// auto cutoff_tof = new TGraph();

Double_t tof_min = 1e2;
Double_t tof_max = 1e8;
Double_t amp_min = 0.;
Double_t amp_max = 50000.;
Int_t num_bins_amp = 5000;
Int_t num_amp_bins_forCuts = 100;

////////////////////////////////////////////////////

// Double_t yOnTheCutLinePTBC(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t x3){
//     return ((y2 - y1)*(x3 - x1)/(x2 - x1) + y1);
// }

void fillRuns(){
    // Bi 1 cm filter
    if(!target_name.compare("bi1")){
        cout << "Setting up Bi (1 cm) run list" << endl;
        // filter_run_list.reserve( c_au_f_bi_t_out.size() + c_out_f_bi_t_out.size() + c_pb_f_bi_t_out.size() + c_ta_f_bi_t_out.size() );
        filter_run_list.insert( filter_run_list.end(), c_au_f_bi_t_out.begin(), c_au_f_bi_t_out.end() );
        filter_run_list.insert( filter_run_list.end(), c_out_f_bi_t_out.begin(), c_out_f_bi_t_out.end() );
        filter_run_list.insert( filter_run_list.end(), c_pb_f_bi_t_out.begin(), c_pb_f_bi_t_out.end() );
        filter_run_list.insert( filter_run_list.end(), c_ta_f_bi_t_out.begin(), c_ta_f_bi_t_out.end() );
    }

    // Bi 1 cm filter - Sep 17
    if(!target_name.compare("bi1sep17")){
        cout << "Setting up Bi (1 cm) Sep 17 run list" << endl;

        std::vector<Int_t> run_list_Bi_sep17 = {117386, 117387, 117388, 117389, 117390};
        filter_run_list.insert( filter_run_list.end(), run_list_Bi_sep17.begin(), run_list_Bi_sep17.end() );
    }

    // Al 8 cm filter
    if(!target_name.compare("al8")){
        cout << "Setting up Al (8 cm) run list" << endl;
        // filter_run_list.reserve( c_ta_f_al8_t_out.size() );
        filter_run_list.insert( filter_run_list.end(), c_ta_f_al8_t_out.begin(), c_ta_f_al8_t_out.end() );
    }

    // Al 5 cm filter
    if(!target_name.compare("al5")){
        cout << "Setting up Al (5 cm) run list" << endl;
        // filter_run_list.reserve( c_ta_f_al5_t_out.size() + c_out_f_al5_t_out.size() + c_pb_f_al5_t_out.size() + c_c_f_al5_t_out.size() + c_au_f_al5_t_out.size() ); //
        filter_run_list.insert( filter_run_list.end(), c_ta_f_al5_t_out.begin(), c_ta_f_al5_t_out.end() );
        filter_run_list.insert( filter_run_list.end(), c_out_f_al5_t_out.begin(), c_out_f_al5_t_out.end() );
        filter_run_list.insert( filter_run_list.end(), c_pb_f_al5_t_out.begin(), c_pb_f_al5_t_out.end() );
        filter_run_list.insert( filter_run_list.end(), c_c_f_al5_t_out.begin(), c_c_f_al5_t_out.end() );
        filter_run_list.insert( filter_run_list.end(), c_au_f_al5_t_out.begin(), c_au_f_al5_t_out.end() );
    }

    // Al 3 cm filter
    if(!target_name.compare("al3")){
        cout << "Setting up Al (3 cm) run list" << endl;
        // filter_run_list.reserve( c_ta_f_al3_t_out.size() + c_pb_f_al3_t_out.size() );
        filter_run_list.insert( filter_run_list.end(), c_ta_f_al3_t_out.begin(), c_ta_f_al3_t_out.end() );
        filter_run_list.insert( filter_run_list.end(), c_pb_f_al3_t_out.begin(), c_pb_f_al3_t_out.end() );
    }

    //--------------------------------------------------

    if(!target_name.compare("bi1p2_ts")){
        cout << "Setting up Bi (1.2 cm) target run list" << endl;
        filter_run_list.insert( filter_run_list.end(), f_out_t_bi_ts.begin(), f_out_t_bi_ts.end() );
    }

    if(!target_name.compare("c1p2_ts")){
        cout << "Setting up C (1.2 cm) target run list" << endl;
        filter_run_list.insert( filter_run_list.end(), f_out_t_c_ts.begin(), f_out_t_c_ts.end() );
    }

    if(!target_name.compare("al5_ts")){
        cout << "Setting up Al (5 cm) target run list" << endl;
        filter_run_list.insert( filter_run_list.end(), f_out_t_al5_ts.begin(), f_out_t_al5_ts.end() );
    }

    if(!target_name.compare("cf_bottle")){
        cout << "Setting up Empty Bottle target run list" << endl;
        filter_run_list.insert( filter_run_list.end(), f_out_t_emptyBot_ts.begin(), f_out_t_emptyBot_ts.end() );
    }

    if(!target_name.compare("cf_bottle_rot")){
        cout << "Setting up Empty Bottle Rotated target run list" << endl;
        filter_run_list.insert( filter_run_list.end(), f_out_t_emptyBotRot_ts.begin(), f_out_t_emptyBotRot_ts.end() );
    }

    if(!target_name.compare("cf_bottle_rotBack")){
        cout << "Setting up Empty Bottle Rotated Back target run list" << endl;
        filter_run_list.insert( filter_run_list.end(), f_out_t_emptyBotRotBack_ts.begin(), f_out_t_emptyBotRotBack_ts.end() );
    }

    if(!target_name.compare("ar_bottle_full")){
        cout << "Setting up Argon Bottle target run list" << endl;
        filter_run_list.insert( filter_run_list.end(), f_out_t_arBot_ts.begin(), f_out_t_arBot_ts.end() );
    }

    if(!target_name.compare("ar_bottle")){
        cout << "Setting up empty argon tank target run list" << endl;
        filter_run_list.insert( filter_run_list.end(), f_out_t_arBotEmpty_ts.begin(), f_out_t_arBotEmpty_ts.end() );
    }

    if(!target_name.compare("none")){
        cout << "Setting up no target/filter run list" << endl;
        filter_run_list.insert( filter_run_list.end(), c_out_f_out_t_out.begin(), c_out_f_out_t_out.end() );
        filter_run_list.insert( filter_run_list.end(), c_pb_f_out_t_out.begin(), c_pb_f_out_t_out.end() );
        filter_run_list.insert( filter_run_list.end(), c_ta_f_out_t_out.begin(), c_ta_f_out_t_out.end() );
        filter_run_list.insert( filter_run_list.end(), c_c_f_out_t_out.begin(), c_c_f_out_t_out.end() );
    }

    if(!target_name.compare("none_ts")){
        cout << "Setting up no target/filter run (transmission station) list" << endl;
        filter_run_list.insert( filter_run_list.end(), f_out_t_out_ts.begin(), f_out_t_out_ts.end() );
    }
}