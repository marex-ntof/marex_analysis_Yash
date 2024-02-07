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

const std::string target_name("ar_bottle_full"); ////bi1, al3, al5, al8, c1p2_ts, al5_ts, bi1p2_ts, cf_bottle, cf_bottle_rot, cf_bottle_rotBack, ar_bottle_full
const std::string target_out_name("ar_bottle"); //none, none_ts, ar_bottle
const std::string target_name_title("Argon Tank");
//Bi (1 cm), Target Bi (1.2 cm), Al (3 cm), Al (5 cm), Target Al (5 cm), Al (8 cm), Target C (1.2 cm), Empty Bottle, Empty Bottle Rotated
//Argon Tank

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
Double_t t_det2[4][2];
Double_t a_det2[4][2];
Double_t t_det3to7[5][2][2];
Double_t a_det3to7[5][2][2];
Double_t t_para[4][2];
Double_t a_para[4][2];

void fillCutsPTBC(){
    //Det 2
    t_det2[0][0] = 800.0;
    a_det2[0][0] = 8000.0;

    t_det2[0][1] = 2605.0;
    a_det2[0][1] = 8000.0;

    t_det2[1][0] = 2600.0;
    a_det2[1][0] = 9500.0;

    t_det2[1][1] = 2800.0;
    a_det2[1][1] = 9500.0;

    t_det2[2][0] = 2800.0;
    a_det2[2][0] = 8000.0;

    t_det2[2][1] = 5000.0;
    a_det2[2][1] = 8000.0;

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
    
    t_para[0][0] = 800.0;
    a_para[0][0] = 5000.0;

    t_para[0][1] = 3000.0;
    a_para[0][1] = 5000.0;

    t_para[1][0] = 3000.0;
    a_para[1][0] = 4000.0;

    t_para[1][1] = 1e8;
    a_para[1][1] = 4000.0;
}

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

///////////////////////////// nTOF Cuts

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

////////////////////////////////////////////////////

Double_t yOnTheCutLinePTBC(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t x3){
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