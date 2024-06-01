/**
 * @file data_stability.h
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-12-04
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

//Root file
TFile *outputRootFile = 0;

TH1D* PTBC_tof_amp_cuts[6];

// TH1D* stability_hist_PTBC = 0;
// TH1D* stability_hist_FIMG_det1 = 0;
// TH1D* stability_hist_FIMG_det2 = 0;

//////////// PTBC - counts plots
// TH1D* norm_counts_empty_sep16_PTBC = 0;
TH1D* norm_counts_empty_sep19_PTBC = 0;
TH1D* norm_counts_empty_oct01_PTBC = 0;

TH1D* norm_counts_Bi_sep17_PTBC = 0;
TH1D* norm_counts_Bi_sep23_PTBC = 0;

TH1D* norm_counts_Al5_sep25_PTBC = 0;
TH1D* norm_counts_Al5_oct02_PTBC = 0;

TH1D* norm_counts_emptyTank_oct12_PTBC = 0;
TH1D* norm_counts_emptyTank_oct16_PTBC = 0;

TH1D* norm_counts_ArgonFull_oct18_PTBC = 0;
TH1D* norm_counts_ArgonFull_oct22_PTBC = 0;

//////////// FIMG - counts plots
TH1D* norm_counts_empty_sep19_FIMG = 0;
TH1D* norm_counts_empty_oct01_FIMG = 0;

TH1D* norm_counts_Bi_sep17_FIMG = 0;
TH1D* norm_counts_Bi_sep23_FIMG = 0;

TH1D* norm_counts_Al5_sep25_FIMG = 0;
TH1D* norm_counts_Al5_oct02_FIMG = 0;

TH1D* norm_counts_emptyTank_oct12_FIMG = 0;
TH1D* norm_counts_emptyTank_oct16_FIMG = 0;

TH1D* norm_counts_ArgonFull_oct18_FIMG = 0;
TH1D* norm_counts_ArgonFull_oct22_FIMG = 0;

//////////// PTBC - day-night plots
// TH1D* day_night_Bi_sep18_PTBC = 0;
// TH1D* day_night_Al5_sep27_PTBC = 0;
// TH1D* day_night_emptyTS_oct05_PTBC = 0;
// TH1D* day_night_emptyTank_oct15_PTBC = 0;
// TH1D* day_night_Argon_oct22_PTBC = 0;

//////////// FIMG - day-night plots
// TH1D* day_night_Bi_sep18_FIMG = 0;
// TH1D* day_night_Al5_sep27_FIMG = 0;
// TH1D* day_night_emptyTS_oct05_FIMG = 0;
// TH1D* day_night_emptyTank_oct15_FIMG = 0;
// TH1D* day_night_Argon_oct22_FIMG = 0;

Int_t bins_per_decade = 100;

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

void applyMyCuts_PTBC(Double_t tof, Float_t amp, Int_t det_num, TH1D* hist_tof){

    if (tof < min_tof_PTBC)
    {
        return;
    }

    Double_t tof_cut_bin = PTBC_tof_amp_cuts[det_num-2]->GetXaxis()->FindBin(tof);
    Double_t amp_cut = PTBC_tof_amp_cuts[det_num-2]->GetBinContent(tof_cut_bin);

    if ( (Double_t) amp > amp_cut)
    {
        hist_tof->Fill(tof); //TOFToEnergy(tof * 1e-9, flight_path_length_PTB)
    }

    return;
}

void applyMyCuts_FIMG(Double_t tof, Float_t amp, Int_t det_num, TH1D* hist){
    //Filling the histograms
    if (tof < 7e3)
    {
        return;
    }
    
    for (int k = 0; k < 5; k++)
    {
        if (tof >= tof_cut_FIMG[det_num-1][k][0] && tof < tof_cut_FIMG[det_num-1][k][1])
        {
            if ( (Double_t) amp > yOnTheCutLine(tof_cut_FIMG[det_num-1][k][0], amp_cut_FIMG[det_num-1][k][0], tof_cut_FIMG[det_num-1][k][1], amp_cut_FIMG[det_num-1][k][1], tof) )
            {
                hist->Fill(tof);
                break;    
            }
        }
    }

    return;
}

bool select_hit_PTBC(Double_t tof, Float_t amp, Float_t pulseIntensity, Int_t det_num){
    if (pulseIntensity <= 6e12)
    {
        for (int k = 0; k < 2; k++)
        {
            if (tof >= t_para[k][0] && tof < t_para[k][1])
            {
                if ( (Double_t) amp > yOnTheCutLine(t_para[k][0], a_para[k][0], t_para[k][1], a_para[k][1], tof) )
                {
                    return 1;
                }
            }
        }
    }

    if (det_num == 2) {
        for (int k = 0; k < 4; k++)
        {
            if (tof >= t_det2[k][0] && tof < t_det2[k][1])
            {
                if ( (Double_t) amp > yOnTheCutLine(t_det2[k][0], a_det2[k][0], t_det2[k][1], a_det2[k][1], tof) )
                {
                    return 1;
                }
            }
        }
    } else {
        for (int k = 0; k < 2; k++)
        {
            if (tof >= t_det3to7[det_num-3][k][0] && tof < t_det3to7[det_num-3][k][1])
            {
                if ( (Double_t) amp > yOnTheCutLine(t_det3to7[det_num-3][k][0], a_det3to7[det_num-3][k][0], t_det3to7[det_num-3][k][1], a_det3to7[det_num-3][k][1], tof) )
                {
                    return 1;
                }
            }
        }
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