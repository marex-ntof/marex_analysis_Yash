/**
 * @file PTBCdet5Ana.C
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
#include "tools.h"

TH2D* tof_amp_hists[11];
std::vector<Int_t> run_list = {117386, 117387, 117388, 117389, 117390, 117391, 117439, 117440, 117441, 117442, 117443};

Double_t t_det5[2][2];
Double_t a_det5[2][2];

// Double_t tof_min = 1e2;
// Double_t tof_max = 1e5;
Double_t amp_min = 0.;
Double_t amp_max = 50000.;
Int_t num_bins_amp = 5000;
Double_t t_gamma_PTBC = (flight_path_length_PTBC / speed_of_light) * 1e9; //converting into ns

TCutG* PTBC_tof_amp_cut_det5;

void fill_cuts(){
    t_det5[0][0] = 800.0;
    a_det5[0][0] = 7000.0;

    t_det5[0][1] = 7000.0;
    a_det5[0][1] = 7000.0;

    t_det5[1][0] = 7000.0;
    a_det5[1][0] = 3500.0;

    t_det5[1][1] = 1e8;
    a_det5[1][1] = 3500.0;

    PTBC_tof_amp_cut_det5 = new TCutG("PTBC_tof_amp_cut_det5",4);
    PTBC_tof_amp_cut_det5->SetLineColor(2);
    PTBC_tof_amp_cut_det5->SetLineWidth(2);
    PTBC_tof_amp_cut_det5->SetVarX("x");
    PTBC_tof_amp_cut_det5->SetVarY("y");
    PTBC_tof_amp_cut_det5->SetPoint(0, t_det5[0][0], amp_max);
    PTBC_tof_amp_cut_det5->SetPoint(1, t_det5[0][0], a_det5[0][0]);
    PTBC_tof_amp_cut_det5->SetPoint(2, t_det5[0][1], a_det5[0][1]);
    PTBC_tof_amp_cut_det5->SetPoint(3, t_det5[1][0], a_det5[1][0]);
    PTBC_tof_amp_cut_det5->SetPoint(4, t_det5[1][1], a_det5[1][1]);
}

void Fill_tof_amp_hist(std::vector<Int_t> run_list, Int_t num_bins_tof, Double_t bin_edges_tof[], Int_t num_bins_amp, Double_t bin_edges_amp[]){

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
        std::cout << Form("Number of entries - Run %i - PTBC = ", run_list.at(i)) << Events_PTBC << std::endl;

        tof_amp_hists[i] = new TH2D(Form("tof_amp_hist_run_%i", run_list.at(i)),Form("ToF vs Amp - PTBC Det 5 Bi (1cm) - Run %i", run_list.at(i)),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);

        for (int j = 0; j < Events_PTBC; j++)
        {
            PTBC->GetEntry(j);

            if(det_num_PTBC != 5){
                continue;
            }

            Double_t t_pkup = BNum_tpkup_map[BunchNumber_PTBC];
            Double_t corrected_tof = tof_PTBC - t_pkup + delT_pkup_ptbc + t_gamma_PTBC;

            tof_amp_hists[i]->Fill(corrected_tof, (Double_t) amp_PTBC);
        }
    }
}

void PTBCdet5Ana(){

    fill_cuts();
    
    //Calculating TOF (x) bin edges
    int bins_per_decade = 1000;
    int Num_decades = 3;
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
    Double_t bin_edges_amp[num_bins_amp+1];
    Double_t step_amp = (Double_t) ((amp_max-amp_min)/num_bins_amp);
    // std::cout << "step_amp = " << step_amp << std::endl;
    for(int i = 0; i < num_bins_amp+1; i++)
    {
        bin_edges_amp[i] = step_amp * (Double_t) i;
    }

    Fill_tof_amp_hist(run_list, num_bins_tof, bin_edges_tof, num_bins_amp, bin_edges_amp);

    TFile *f = new TFile("../rootFiles/PTBCdet5Ana.root","recreate");

    for (int i = 0; i < run_list.size(); i++){

        tof_amp_hists[i]->Write();
    }

    PTBC_tof_amp_cut_det5->Write();

    f->Close();

    std::cout << "Created output file 'PTBCdet5Ana.root'" << std::endl;
}