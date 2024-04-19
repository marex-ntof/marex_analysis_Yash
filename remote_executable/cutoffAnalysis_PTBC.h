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
const std::string target_out_name("none"); //none, none_ts, ar_bottle
const std::string target_name_title("Bi (1 cm)");
//Bi (1 cm), Target Bi (1.2 cm), Al (3 cm), Al (5 cm), Target Al (5 cm), Al (8 cm), Target C (1.2 cm), Empty Bottle, Empty Bottle Rotated
//Argon Tank, Bi (1 cm) - Sep 17

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

TH2D* PTBC_tof_amp_fIn_det2_afterCuts = 0;
TH2D* PTBC_tof_amp_fIn_det3_afterCuts = 0;
TH2D* PTBC_tof_amp_fIn_det4_afterCuts = 0;
TH2D* PTBC_tof_amp_fIn_det5_afterCuts = 0;
TH2D* PTBC_tof_amp_fIn_det6_afterCuts = 0;
TH2D* PTBC_tof_amp_fIn_det7_afterCuts = 0;
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
TCutG* PTBC_tof_amp_cut_det2;
TCutG* PTBC_tof_amp_cut_det3;
TCutG* PTBC_tof_amp_cut_det4;
TCutG* PTBC_tof_amp_cut_det5;
TCutG* PTBC_tof_amp_cut_det5_earlyruns;
TCutG* PTBC_tof_amp_cut_det6;
TCutG* PTBC_tof_amp_cut_det7;
TCutG* PTBC_tof_amp_cut_para;

// Run Vectors
std::vector<Int_t> filter_in_runs;
std::vector<Int_t> filter_out_runs;

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
Double_t t_gamma_PTBC = (flight_path_length_PTBC / speed_of_light) * 1e9; //converting into ns

void fillCutGraph(){
    PTBC_tof_amp_cut_det2 = new TCutG("PTBC_tof_amp_cut_det2",8);
    PTBC_tof_amp_cut_det2->SetLineColor(2);
    PTBC_tof_amp_cut_det2->SetLineWidth(2);
    PTBC_tof_amp_cut_det2->SetVarX("x");
    PTBC_tof_amp_cut_det2->SetVarY("y");
    PTBC_tof_amp_cut_det2->SetPoint(0, t_det2[0][0], amp_max);
    PTBC_tof_amp_cut_det2->SetPoint(1, t_det2[0][0], a_det2[0][0]);
    PTBC_tof_amp_cut_det2->SetPoint(2, t_det2[0][1], a_det2[0][1]);
    PTBC_tof_amp_cut_det2->SetPoint(3, t_det2[1][0], a_det2[1][0]);
    PTBC_tof_amp_cut_det2->SetPoint(4, t_det2[1][1], a_det2[1][1]);
    PTBC_tof_amp_cut_det2->SetPoint(5, t_det2[2][0], a_det2[2][0]);
    PTBC_tof_amp_cut_det2->SetPoint(6, t_det2[2][1], a_det2[2][1]);
    PTBC_tof_amp_cut_det2->SetPoint(7, t_det2[3][0], a_det2[3][0]);
    PTBC_tof_amp_cut_det2->SetPoint(8, t_det2[3][1], a_det2[3][1]);
    
    PTBC_tof_amp_cut_det3 = new TCutG("PTBC_tof_amp_cut_det3",4);
    PTBC_tof_amp_cut_det3->SetLineColor(2);
    PTBC_tof_amp_cut_det3->SetLineWidth(2);
    PTBC_tof_amp_cut_det3->SetVarX("x");
    PTBC_tof_amp_cut_det3->SetVarY("y");
    PTBC_tof_amp_cut_det3->SetPoint(0, t_det3to7[0][0][0], amp_max);
    PTBC_tof_amp_cut_det3->SetPoint(1, t_det3to7[0][0][0], a_det3to7[0][0][0]);
    PTBC_tof_amp_cut_det3->SetPoint(2, t_det3to7[0][0][1], a_det3to7[0][0][1]);
    PTBC_tof_amp_cut_det3->SetPoint(3, t_det3to7[0][1][0], a_det3to7[0][1][0]);
    PTBC_tof_amp_cut_det3->SetPoint(4, t_det3to7[0][1][1], a_det3to7[0][1][1]);
    
    PTBC_tof_amp_cut_det4 = new TCutG("PTBC_tof_amp_cut_det4",4);
    PTBC_tof_amp_cut_det4->SetLineColor(2);
    PTBC_tof_amp_cut_det4->SetLineWidth(2);
    PTBC_tof_amp_cut_det4->SetVarX("x");
    PTBC_tof_amp_cut_det4->SetVarY("y");
    PTBC_tof_amp_cut_det4->SetPoint(0, t_det3to7[1][0][0], amp_max);
    PTBC_tof_amp_cut_det4->SetPoint(1, t_det3to7[1][0][0], a_det3to7[1][0][0]);
    PTBC_tof_amp_cut_det4->SetPoint(2, t_det3to7[1][0][1], a_det3to7[1][0][1]);
    PTBC_tof_amp_cut_det4->SetPoint(3, t_det3to7[1][1][0], a_det3to7[1][1][0]);
    PTBC_tof_amp_cut_det4->SetPoint(4, t_det3to7[1][1][1], a_det3to7[1][1][1]);
    
    PTBC_tof_amp_cut_det5 = new TCutG("PTBC_tof_amp_cut_det5",4);
    PTBC_tof_amp_cut_det5->SetLineColor(2);
    PTBC_tof_amp_cut_det5->SetLineWidth(2);
    PTBC_tof_amp_cut_det5->SetVarX("x");
    PTBC_tof_amp_cut_det5->SetVarY("y");
    PTBC_tof_amp_cut_det5->SetPoint(0, t_det3to7[2][0][0], amp_max);
    PTBC_tof_amp_cut_det5->SetPoint(1, t_det3to7[2][0][0], a_det3to7[2][0][0]);
    PTBC_tof_amp_cut_det5->SetPoint(2, t_det3to7[2][0][1], a_det3to7[2][0][1]);
    PTBC_tof_amp_cut_det5->SetPoint(3, t_det3to7[2][1][0], a_det3to7[2][1][0]);
    PTBC_tof_amp_cut_det5->SetPoint(4, t_det3to7[2][1][1], a_det3to7[2][1][1]);

    // from Run 117386 to 117390
    PTBC_tof_amp_cut_det5_earlyruns = new TCutG("PTBC_tof_amp_cut_det5_earlyruns",6);
    PTBC_tof_amp_cut_det5_earlyruns->SetLineColor(2);
    PTBC_tof_amp_cut_det5_earlyruns->SetLineWidth(2);
    PTBC_tof_amp_cut_det5_earlyruns->SetVarX("x");
    PTBC_tof_amp_cut_det5_earlyruns->SetVarY("y");
    PTBC_tof_amp_cut_det5_earlyruns->SetPoint(0, t_det5_early_runs[0][0], amp_max);
    PTBC_tof_amp_cut_det5_earlyruns->SetPoint(1, t_det5_early_runs[0][0], a_det5_early_runs[0][0]);
    PTBC_tof_amp_cut_det5_earlyruns->SetPoint(2, t_det5_early_runs[0][1], a_det5_early_runs[0][1]);
    PTBC_tof_amp_cut_det5_earlyruns->SetPoint(3, t_det5_early_runs[1][0], a_det5_early_runs[1][0]);
    PTBC_tof_amp_cut_det5_earlyruns->SetPoint(4, t_det5_early_runs[1][1], a_det5_early_runs[1][1]);
    PTBC_tof_amp_cut_det5_earlyruns->SetPoint(5, t_det5_early_runs[2][1], a_det5_early_runs[2][1]);
    PTBC_tof_amp_cut_det5_earlyruns->SetPoint(6, t_det5_early_runs[3][1], a_det5_early_runs[3][1]);
    
    PTBC_tof_amp_cut_det6 = new TCutG("PTBC_tof_amp_cut_det6",4);
    PTBC_tof_amp_cut_det6->SetLineColor(2);
    PTBC_tof_amp_cut_det6->SetLineWidth(2);
    PTBC_tof_amp_cut_det6->SetVarX("x");
    PTBC_tof_amp_cut_det6->SetVarY("y");
    PTBC_tof_amp_cut_det6->SetPoint(0, t_det3to7[3][0][0], amp_max);
    PTBC_tof_amp_cut_det6->SetPoint(1, t_det3to7[3][0][0], a_det3to7[3][0][0]);
    PTBC_tof_amp_cut_det6->SetPoint(2, t_det3to7[3][0][1], a_det3to7[3][0][1]);
    PTBC_tof_amp_cut_det6->SetPoint(3, t_det3to7[3][1][0], a_det3to7[3][1][0]);
    PTBC_tof_amp_cut_det6->SetPoint(4, t_det3to7[3][1][1], a_det3to7[3][1][1]);
    
    PTBC_tof_amp_cut_det7 = new TCutG("PTBC_tof_amp_cut_det7",4);
    PTBC_tof_amp_cut_det7->SetLineColor(2);
    PTBC_tof_amp_cut_det7->SetLineWidth(2);
    PTBC_tof_amp_cut_det7->SetVarX("x");
    PTBC_tof_amp_cut_det7->SetVarY("y");
    PTBC_tof_amp_cut_det7->SetPoint(0, t_det3to7[4][0][0], amp_max);
    PTBC_tof_amp_cut_det7->SetPoint(1, t_det3to7[4][0][0], a_det3to7[4][0][0]);
    PTBC_tof_amp_cut_det7->SetPoint(2, t_det3to7[4][0][1], a_det3to7[4][0][1]);
    PTBC_tof_amp_cut_det7->SetPoint(3, t_det3to7[4][1][0], a_det3to7[4][1][0]);
    PTBC_tof_amp_cut_det7->SetPoint(4, t_det3to7[4][1][1], a_det3to7[4][1][1]);

    PTBC_tof_amp_cut_para = new TCutG("PTBC_tof_amp_cut_para",4);
    PTBC_tof_amp_cut_para->SetLineColor(2);
    PTBC_tof_amp_cut_para->SetLineWidth(2);
    PTBC_tof_amp_cut_para->SetVarX("x");
    PTBC_tof_amp_cut_para->SetVarY("y");
    PTBC_tof_amp_cut_para->SetPoint(0, 800.0, amp_max);
    PTBC_tof_amp_cut_para->SetPoint(1, 800.0, 5000.0);
    PTBC_tof_amp_cut_para->SetPoint(2, 3000.0, 5000.0);
    PTBC_tof_amp_cut_para->SetPoint(3, 3000.0, 4000.0);
    PTBC_tof_amp_cut_para->SetPoint(4, 1e8, 4000.0);
}

////////////////////////////////////////////////////

Double_t yOnTheCutLinePTBC(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t x3){
    return ((y2 - y1)*(x3 - x1)/(x2 - x1) + y1);
}

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

    // Bi 1 cm filter - Sep 17
    if(!target_name.compare("bi1sep17")){
        cout << "Setting up Bi (1 cm) Sep 17 run list" << endl;

        std::vector<Int_t> run_list_Bi_sep17 = {117386, 117387, 117388, 117389, 117390};
        filter_in_runs.insert( filter_in_runs.end(), run_list_Bi_sep17.begin(), run_list_Bi_sep17.end() );
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