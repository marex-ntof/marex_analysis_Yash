/**
 * @file tools.h
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2024-02-22
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

Double_t EnergyToVelocity(Double_t e){ //e is in eV
    Double_t KE_M = (e * 1e-6)/neutron_mass;
    Double_t factor = std::sqrt(1.0 - 1.0/((KE_M + 1.0)*(KE_M + 1.0)));
    return speed_of_light * factor; //velocity in m/s
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

Double_t yOnTheCutLine(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t x3){
    if(y1 == y2){
        return y1;
    }
    return ((y2 - y1)*(x3 - x1)/(x2 - x1) + y1);
}

Double_t FindFWHM(TH1D* projection_hist){
    
    Int_t max_bin_num = projection_hist->GetMaximumBin();
    Double_t max_value = projection_hist->GetBinContent(max_bin_num);
    Int_t tot_bins = projection_hist->GetNbinsX();
    Double_t left_edge = 0;
    Double_t right_edge = 0;
    //Find left edge
    for (Int_t i = max_bin_num-1; i > 0; i--)
    {
        Double_t bin_value = projection_hist->GetBinContent(i);
        if (bin_value > max_value/2.0)
        {
            continue;
        } else {
            left_edge = projection_hist->GetXaxis()->GetBinUpEdge(i);
            break;
        }
    }

    //Find right edge
    for (Int_t i = max_bin_num+1; i < tot_bins+1; i++)
    {
        Double_t bin_value = projection_hist->GetBinContent(i);
        if (bin_value > max_value/2.0)
        {
            continue;
        } else {
            right_edge = projection_hist->GetXaxis()->GetBinLowEdge(i);
            break;
        }
    }
    return (right_edge - left_edge);
}