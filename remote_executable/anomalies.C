/**
 * @file anomalies.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2024-04-26
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

TH2D* PTBC_al5_anomaly[6];
TH2D* FIMG_al5_anomaly[2];
TH2D* PTBC_filterOut_anomaly[6];

TH2D* PTBC_filterOut_anomaly_before[6];
TH2D* PTBC_filterOut_anomaly_after[6];

//My cut plots
TCutG* FIMG_my_tof_amp_cut_dedi_det1;
TCutG* FIMG_my_tof_amp_cut_dedi_det2;
TCutG* FIMG_my_tof_amp_cut_para_det1;
TCutG* FIMG_my_tof_amp_cut_para_det2;

Double_t tof_min = 1e2;
Double_t tof_max = 1e8;
Double_t amp_min = 0.;
Double_t amp_max = 50000.;

// Run Lists
std::vector<Int_t> filterOut_run_list;
std::vector<Int_t> al5_run_list;

void fillCutGraph_FIMG(){

    FIMG_my_tof_amp_cut_dedi_det1 = new TCutG("FIMG_my_tof_amp_cut_dedi_det1",6);
    FIMG_my_tof_amp_cut_dedi_det1->SetLineColor(2);
    FIMG_my_tof_amp_cut_dedi_det1->SetLineWidth(2);
    FIMG_my_tof_amp_cut_dedi_det1->SetVarX("x");
    FIMG_my_tof_amp_cut_dedi_det1->SetVarY("y");
    FIMG_my_tof_amp_cut_dedi_det1->SetPoint(0, tof_cut_FIMG[0][0][0], amp_max);
    FIMG_my_tof_amp_cut_dedi_det1->SetPoint(1, tof_cut_FIMG[0][0][0], amp_cut_FIMG[0][0][0]);
    FIMG_my_tof_amp_cut_dedi_det1->SetPoint(2, tof_cut_FIMG[0][1][0], amp_cut_FIMG[0][1][0]);
    FIMG_my_tof_amp_cut_dedi_det1->SetPoint(3, tof_cut_FIMG[0][2][0], amp_cut_FIMG[0][2][0]);
    FIMG_my_tof_amp_cut_dedi_det1->SetPoint(4, tof_cut_FIMG[0][3][0], amp_cut_FIMG[0][3][0]);
    FIMG_my_tof_amp_cut_dedi_det1->SetPoint(5, tof_cut_FIMG[0][4][0], amp_cut_FIMG[0][4][0]);
    FIMG_my_tof_amp_cut_dedi_det1->SetPoint(6, tof_cut_FIMG[0][4][1], amp_cut_FIMG[0][4][1]);

    FIMG_my_tof_amp_cut_dedi_det2 = new TCutG("FIMG_my_tof_amp_cut_dedi_det2",6);
    FIMG_my_tof_amp_cut_dedi_det2->SetLineColor(2);
    FIMG_my_tof_amp_cut_dedi_det2->SetLineWidth(2);
    FIMG_my_tof_amp_cut_dedi_det2->SetVarX("x");
    FIMG_my_tof_amp_cut_dedi_det2->SetVarY("y");
    FIMG_my_tof_amp_cut_dedi_det2->SetPoint(0, tof_cut_FIMG[1][0][0], amp_max);
    FIMG_my_tof_amp_cut_dedi_det2->SetPoint(1, tof_cut_FIMG[1][0][0], amp_cut_FIMG[1][0][0]);
    FIMG_my_tof_amp_cut_dedi_det2->SetPoint(2, tof_cut_FIMG[1][1][0], amp_cut_FIMG[1][1][0]);
    FIMG_my_tof_amp_cut_dedi_det2->SetPoint(3, tof_cut_FIMG[1][2][0], amp_cut_FIMG[1][2][0]);
    FIMG_my_tof_amp_cut_dedi_det2->SetPoint(4, tof_cut_FIMG[1][3][0], amp_cut_FIMG[1][3][0]);
    FIMG_my_tof_amp_cut_dedi_det2->SetPoint(5, tof_cut_FIMG[1][4][0], amp_cut_FIMG[1][4][0]);
    FIMG_my_tof_amp_cut_dedi_det2->SetPoint(6, tof_cut_FIMG[1][4][1], amp_cut_FIMG[1][4][1]);

    FIMG_my_tof_amp_cut_para_det1 = new TCutG("FIMG_my_tof_amp_cut_para_det1",6);
    FIMG_my_tof_amp_cut_para_det1->SetLineColor(2);
    FIMG_my_tof_amp_cut_para_det1->SetLineWidth(2);
    FIMG_my_tof_amp_cut_para_det1->SetVarX("x");
    FIMG_my_tof_amp_cut_para_det1->SetVarY("y");
    FIMG_my_tof_amp_cut_para_det1->SetPoint(0, tof_cut_FIMG[0][0][0], amp_max);
    FIMG_my_tof_amp_cut_para_det1->SetPoint(1, tof_cut_FIMG[0][0][0], amp_cut_FIMG[0][0][0]);
    FIMG_my_tof_amp_cut_para_det1->SetPoint(2, tof_cut_FIMG[0][1][0], amp_cut_FIMG[0][1][0]);
    FIMG_my_tof_amp_cut_para_det1->SetPoint(3, tof_cut_FIMG[0][2][0], amp_cut_FIMG[0][2][0]);
    FIMG_my_tof_amp_cut_para_det1->SetPoint(4, tof_cut_FIMG[0][3][0], amp_cut_FIMG[0][3][0]);
    FIMG_my_tof_amp_cut_para_det1->SetPoint(5, tof_cut_FIMG[0][4][0], amp_cut_FIMG[0][4][0]);
    FIMG_my_tof_amp_cut_para_det1->SetPoint(6, tof_cut_FIMG[0][4][1], amp_cut_FIMG[0][4][1]);

    FIMG_my_tof_amp_cut_para_det2 = new TCutG("FIMG_my_tof_amp_cut_para_det2",6);
    FIMG_my_tof_amp_cut_para_det2->SetLineColor(2);
    FIMG_my_tof_amp_cut_para_det2->SetLineWidth(2);
    FIMG_my_tof_amp_cut_para_det2->SetVarX("x");
    FIMG_my_tof_amp_cut_para_det2->SetVarY("y");
    FIMG_my_tof_amp_cut_para_det2->SetPoint(0, tof_cut_FIMG[1][0][0], amp_max);
    FIMG_my_tof_amp_cut_para_det2->SetPoint(1, tof_cut_FIMG[1][0][0], amp_cut_FIMG[1][0][0]);
    FIMG_my_tof_amp_cut_para_det2->SetPoint(2, tof_cut_FIMG[1][1][0], amp_cut_FIMG[1][1][0]);
    FIMG_my_tof_amp_cut_para_det2->SetPoint(3, tof_cut_FIMG[1][2][0], amp_cut_FIMG[1][2][0]);
    FIMG_my_tof_amp_cut_para_det2->SetPoint(4, tof_cut_FIMG[1][3][0], amp_cut_FIMG[1][3][0]);
    FIMG_my_tof_amp_cut_para_det2->SetPoint(5, tof_cut_FIMG[1][4][0], amp_cut_FIMG[1][4][0]);
    FIMG_my_tof_amp_cut_para_det2->SetPoint(6, tof_cut_FIMG[1][4][1], amp_cut_FIMG[1][4][1]);
}

void fill_run_list(){

    // Filter out runs
    filterOut_run_list.insert(filterOut_run_list.end(), c_out_f_out_t_out.begin(), c_out_f_out_t_out.end());    
    filterOut_run_list.insert(filterOut_run_list.end(), c_pb_f_out_t_out.begin(), c_pb_f_out_t_out.end());
    filterOut_run_list.insert(filterOut_run_list.end(), c_ta_f_out_t_out.begin(), c_ta_f_out_t_out.end());
    filterOut_run_list.insert(filterOut_run_list.end(), c_c_f_out_t_out.begin(), c_c_f_out_t_out.end());

    // Al (5 cm)
    al5_run_list.insert(al5_run_list.end(), c_ta_f_al5_t_out.begin(), c_ta_f_al5_t_out.end());    
    al5_run_list.insert(al5_run_list.end(), c_out_f_al5_t_out.begin(), c_out_f_al5_t_out.end());
    al5_run_list.insert(al5_run_list.end(), c_pb_f_al5_t_out.begin(), c_pb_f_al5_t_out.end());
    al5_run_list.insert(al5_run_list.end(), c_c_f_al5_t_out.begin(), c_c_f_al5_t_out.end());
    al5_run_list.insert(al5_run_list.end(), c_au_f_al5_t_out.begin(), c_au_f_al5_t_out.end());    
}

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

void Fill_tof_amp_hists(std::vector<Int_t> run_list, Int_t det_type, TH2D* hist_array[], Int_t start_time, Int_t end_time){ // 1 = PTBC, 2 = FIMG

    for (int i = 0; i < run_list.size(); i++)
    {
        TFile *file_ntof = TFile::Open(Form("/eos/experiment/ntof/data/rootfiles/2023/ear1/run%d.root", run_list.at(i)),"read");
        
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

        if (det_type == 1)
        {
            //PTBC ---------------------------------------------
            TTree* PTBC;
            Int_t time_PTBC = 0;
            Double_t tof_PTBC = 0; //tof is in ns
            Float_t amp_PTBC = 0;
            Int_t BunchNumber_PTBC = 0;
            Int_t det_num_PTBC = 0;
            // Float_t PulseIntensity_PTBC = 0;

            file_ntof->GetObject("PTBC", PTBC);
            PTBC->SetBranchAddress("time", &time_PTBC);
            PTBC->SetBranchAddress("BunchNumber", &BunchNumber_PTBC);
            // PTBC->SetBranchAddress("PulseIntensity", &PulseIntensity_PTBC);
            PTBC->SetBranchAddress("tof", &tof_PTBC);
            PTBC->SetBranchAddress("amp", &amp_PTBC);
            PTBC->SetBranchAddress("detn", &det_num_PTBC);

            Long64_t Events_PTBC = PTBC->GetEntriesFast();
            std::cout << Form("Number of entries - Run %i - PTBC = ", run_list.at(i)) << Events_PTBC << std::endl;

            for (int j = 0; j < Events_PTBC; j++)
            {
                PTBC->GetEntry(j);

                if (det_num_PTBC == 1 || det_num_PTBC == 8){
                    continue;
                }

                Double_t t_pkup = BNum_tpkup_map[BunchNumber_PTBC];
                Double_t corrected_tof = tof_PTBC - t_pkup + delT_pkup_ptbc + t_gamma_PTBC;

                Int_t hit_time = timeToSeconds(std::to_string(time_PTBC));

                if (hit_time >= start_time && hit_time < end_time)
                {
                    hist_array[det_num_PTBC-2]->Fill(corrected_tof, (Double_t) amp_PTBC);
                }
            }            
        }

        if (det_type == 2)
        {
            //FIMG ---------------------------------------------
            TTree* FIMG;
            Int_t time_FIMG = 0;
            Double_t tof_FIMG = 0; //tof is in ns
            Float_t amp_FIMG = 0;
            Int_t det_num_FIMG = 0;
            Int_t BunchNumber_FIMG = 0;
            // Float_t PulseIntensity_FIMG = 0;

            file_ntof->GetObject("FIMG", FIMG);
            FIMG->SetBranchAddress("time", &time_FIMG);
            FIMG->SetBranchAddress("BunchNumber", &BunchNumber_FIMG);
            FIMG->SetBranchAddress("tof", &tof_FIMG);
            FIMG->SetBranchAddress("amp", &amp_FIMG);
            FIMG->SetBranchAddress("detn", &det_num_FIMG);
            // FIMG->SetBranchAddress("PulseIntensity", &PulseIntensity_FIMG);

            Long64_t Events_FIMG = FIMG->GetEntriesFast();
            std::cout << Form("Number of entries - Run %i - FIMG = ", run_list.at(i)) << Events_FIMG << std::endl;

            for (int j = 0; j < Events_FIMG; j++)
            {
                FIMG->GetEntry(j);

                Double_t t_pkup = BNum_tpkup_map[BunchNumber_FIMG];
                Double_t corrected_tof = tof_FIMG - t_pkup + delT_pkup_fimg + t_gamma_FIMG;

                Int_t hit_time = timeToSeconds(std::to_string(time_FIMG));

                if (hit_time >= start_time && hit_time < end_time)
                {
                    hist_array[det_num_FIMG-1]->Fill(corrected_tof, (Double_t) amp_FIMG);
                }
            }
        }
        file_ntof->Close();
    }
}

void anomalies(){

    // fillCutGraph_PTBC();
    // fillCutGraph_FIMG();
    fill_run_list();

    //Calculating TOF (x) bin edges
    Int_t bins_per_decade = 1000;
    Int_t Num_decades = 6;
    Int_t num_bins_tof = bins_per_decade * Num_decades;
    Double_t bin_edges_tof[num_bins_tof+1];
    Double_t step_tof = ((Double_t) 1.0/(Double_t) bins_per_decade);
    for(Int_t i = 0; i < num_bins_tof+1; i++)
    {
        Double_t base = 10.;
        Double_t exponent = (step_tof * (Double_t) i) + 2.;
        bin_edges_tof[i] = (Double_t) std::pow(base, exponent);
    }

    // Calculating amplitude (y) bin edges
    Int_t num_bins_amp = 5000;
    Double_t bin_edges_amp[num_bins_amp+1];
    Double_t step_amp = (Double_t) ((amp_max-amp_min)/num_bins_amp);
    // std::cout << "step_amp = " << step_amp << std::endl;
    for(Int_t i = 0; i < num_bins_amp+1; i++)
    {
        bin_edges_amp[i] = step_amp * (Double_t) i;
    }

    for (Int_t i = 0; i < 6; i++){
        PTBC_al5_anomaly[i] = new TH2D(Form("PTBC_al5_anomaly_det%i", i+2), Form("ToF vs Amp Hist - Al (5 cm) PTBC Det %i", i+2), num_bins_tof, bin_edges_tof, num_bins_amp, bin_edges_amp);
        // PTBC_filterOut_anomaly[i] = new TH2D(Form("PTBC_filterOut_anomaly_det%i", i+2), Form("ToF vs Amp Hist - Filter Out PTBC Det %i", i+2), num_bins_tof, bin_edges_tof, num_bins_amp, bin_edges_amp);
        // PTBC_filterOut_anomaly_after[i] = new TH2D(Form("PTBC_filterOut_anomaly_after_det%i", i+2), Form("ToF vs Amp Hist - Filter Out PTBC Det %i", i+2), num_bins_tof, bin_edges_tof, num_bins_amp, bin_edges_amp);
        // PTBC_filterOut_anomaly_before[i] = new TH2D(Form("PTBC_filterOut_anomaly_before_det%i", i+2), Form("ToF vs Amp Hist - Filter Out PTBC Det %i", i+2), num_bins_tof, bin_edges_tof, num_bins_amp, bin_edges_amp);
    }

    // for (Int_t i = 0; i < 2; i++){
    //     FIMG_al5_anomaly[i] = new TH2D(Form("FIMG_al5_anomaly_det%i", i+1), Form("ToF vs Amp Hist - Al (5 cm) FIMG Det %i", i+1), num_bins_tof, bin_edges_tof, num_bins_amp, bin_edges_amp);
    // }

    Int_t start_time_al5 = timeToSeconds(std::to_string(121500));
    Int_t end_time_al5 = timeToSeconds(std::to_string(130000));

    // Int_t start_time_filterOut = timeToSeconds(std::to_string(100000));
    // Int_t end_time_filterOut = timeToSeconds(std::to_string(101500));

    // Int_t start_time_filterOut_after = timeToSeconds(std::to_string(101500));
    // Int_t end_time_filterOut_after = timeToSeconds(std::to_string(103000));

    // Int_t start_time_filterOut_before = timeToSeconds(std::to_string(94500));
    // Int_t end_time_filterOut_before = timeToSeconds(std::to_string(100000));

    Fill_tof_amp_hists(al5_run_list, 1, PTBC_al5_anomaly, start_time_al5, end_time_al5); // 1 = PTBC, 2 = FIMG
    Fill_tof_amp_hists(al5_run_list, 2, FIMG_al5_anomaly, start_time_al5, end_time_al5); // 1 = PTBC, 2 = FIMG
    // Fill_tof_amp_hists(filterOut_run_list, 1, PTBC_filterOut_anomaly, start_time_filterOut, end_time_filterOut); // 1 = PTBC, 2 = FIMG
    // Fill_tof_amp_hists(filterOut_run_list, 1, PTBC_filterOut_anomaly_after, start_time_filterOut_after, end_time_filterOut_after); // 1 = PTBC, 2 = FIMG
    // Fill_tof_amp_hists(filterOut_run_list, 1, PTBC_filterOut_anomaly_before, start_time_filterOut_before, end_time_filterOut_before); // 1 = PTBC, 2 = FIMG

    //Writing to the output file
    outputRootFile = new TFile("../rootFiles/anomalies.root","recreate");

    for (Int_t i = 0; i < 6; i++){
        PTBC_al5_anomaly[i]->Write();
    }

    for (Int_t i = 0; i < 2; i++){
        FIMG_al5_anomaly[i]->Write();
    }

    // for (Int_t i = 0; i < 6; i++){
    //     PTBC_filterOut_anomaly[i]->Write();
    // }

    // for (Int_t i = 0; i < 6; i++){
    //     PTBC_filterOut_anomaly_after[i]->Write();
    // }

    // for (Int_t i = 0; i < 6; i++){
    //     PTBC_filterOut_anomaly_before[i]->Write();
    // }

    // FIMG_my_tof_amp_cut_dedi_det1->Write();
    // FIMG_my_tof_amp_cut_dedi_det2->Write();
    // FIMG_my_tof_amp_cut_para_det1->Write();
    // FIMG_my_tof_amp_cut_para_det2->Write();

    outputRootFile->Close();

    std::cout << "Created output file 'anomalies.root'" << std::endl;
}