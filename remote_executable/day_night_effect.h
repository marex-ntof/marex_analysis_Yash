/**
 * @file day_night_effect.h
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2024-04-11
 */

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <cmath>

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

//Root file
TFile *outputRootFile = 0;

//PTBC cuts
TH1D* PTBC_tof_amp_cuts[6];

//////////// Counts plots
TH1D* counts_filterOut_PTBC = 0;
TH1D* counts_Bi_PTBC = 0;
TH1D* counts_Al5_PTBC = 0;
TH1D* counts_emptyTS_PTBC = 0;
TH1D* counts_emptyTank_PTBC = 0;
TH1D* counts_Argon_PTBC = 0;
TH1D* counts_EmptyArgon_PTBC = 0;

TH1D* counts_filterOut_FIMG = 0;
TH1D* counts_Bi_FIMG = 0;
TH1D* counts_Al5_FIMG = 0;
TH1D* counts_emptyTS_FIMG = 0;
TH1D* counts_emptyTank_FIMG = 0;
TH1D* counts_Argon_FIMG = 0;
TH1D* counts_EmptyArgon_FIMG = 0;

//////////// Pulse Intensity plots
TH1D* pulseIntensity_filterOut_PTBC = 0;
TH1D* pulseIntensity_Bi_PTBC = 0;
TH1D* pulseIntensity_Al5_PTBC = 0;
TH1D* pulseIntensity_emptyTS_PTBC = 0;
TH1D* pulseIntensity_emptyTank_PTBC = 0;
TH1D* pulseIntensity_Argon_PTBC = 0;
TH1D* pulseIntensity_EmptyArgon_PTBC = 0;

TH1D* pulseIntensity_filterOut_FIMG = 0;
TH1D* pulseIntensity_Bi_FIMG = 0;
TH1D* pulseIntensity_Al5_FIMG = 0;
TH1D* pulseIntensity_emptyTS_FIMG = 0;
TH1D* pulseIntensity_emptyTank_FIMG = 0;
TH1D* pulseIntensity_Argon_FIMG = 0;
TH1D* pulseIntensity_EmptyArgon_FIMG = 0;

//////////// PTBC - day-night plots
TH1D* day_night_filterOut_PTBC = 0;
TH1D* day_night_Bi_PTBC = 0;
TH1D* day_night_Al5_PTBC = 0;
TH1D* day_night_emptyTS_PTBC = 0;
TH1D* day_night_emptyTank_PTBC = 0;
TH1D* day_night_Argon_PTBC = 0;
TH1D* day_night_EmptyArgon_PTBC = 0;

//////////// FIMG - day-night plots
TH1D* day_night_filterOut_FIMG = 0;
TH1D* day_night_Bi_FIMG = 0;
TH1D* day_night_Al5_FIMG = 0;
TH1D* day_night_emptyTS_FIMG = 0;
TH1D* day_night_emptyTank_FIMG = 0;
TH1D* day_night_Argon_FIMG = 0;
TH1D* day_night_EmptyArgon_FIMG = 0;

// Run Lists
std::vector<Int_t> filterOut_run_list;
std::vector<Int_t> bi1_run_list;
std::vector<Int_t> al5_run_list;

// Function to convert hhmmss format to seconds
Int_t timeToSeconds(const std::string& timeStr) {

    if (timeStr.length() == 1 || timeStr.length() == 2) {
        return std::stoi(timeStr);
    }

    if (timeStr.length() == 3) {

        Int_t minutes = std::stoi(timeStr.substr(0, 1));
        Int_t seconds = std::stoi(timeStr.substr(1, 2));

        return minutes * 60 + seconds;
    }

    if (timeStr.length() == 4) {

        Int_t minutes = std::stoi(timeStr.substr(0, 2));
        Int_t seconds = std::stoi(timeStr.substr(2, 2));

        return minutes * 60 + seconds;
    }

    if (timeStr.length() == 5) {

        Int_t hours = std::stoi(timeStr.substr(0, 1));
        Int_t minutes = std::stoi(timeStr.substr(1, 2));
        Int_t seconds = std::stoi(timeStr.substr(3, 2));

        return hours * 3600 + minutes * 60 + seconds;
    }

    if (timeStr.length() == 6) {

        Int_t hours = std::stoi(timeStr.substr(0, 2));
        Int_t minutes = std::stoi(timeStr.substr(2, 2));
        Int_t seconds = std::stoi(timeStr.substr(4, 2));

        return hours * 3600 + minutes * 60 + seconds;
    }

    else {
        // Invalid input format
        std::cerr << "Invalid time format: " << timeStr << std::endl;
        return -1;
    }
}

bool select_hit_PTBC(Double_t tof, Float_t amp, Int_t det_num){

    if (tof < min_tof_PTBC)
    {
        return 0;
    }

    Double_t tof_cut_bin = PTBC_tof_amp_cuts[det_num-2]->GetXaxis()->FindBin(tof);
    Double_t amp_cut = PTBC_tof_amp_cuts[det_num-2]->GetBinContent(tof_cut_bin);

    if ( (Double_t) amp > amp_cut)
    {
        return 1;
    }

    return 0;
}

bool select_hit_FIMG(Double_t tof, Float_t amp, Int_t det_num){
    
    if (tof < 7e3)
    {
        return 0;
    }
    
    for (int k = 0; k < 5; k++)
    {
        if (tof >= tof_cut_FIMG[det_num-1][k][0] && tof < tof_cut_FIMG[det_num-1][k][1])
        {
            if ( (Double_t) amp > yOnTheCutLine(tof_cut_FIMG[det_num-1][k][0], amp_cut_FIMG[det_num-1][k][0], tof_cut_FIMG[det_num-1][k][1], amp_cut_FIMG[det_num-1][k][1], tof) )
            {
                return 1;
            }
        }
    }

    return 0;
}