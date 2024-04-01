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

// TH2D* FIMG_tof_amp_in = 0;
// TH2D* FIMG_tof_amp_out = 0;
// TH2D* tof_area_hist_in = 0;

TH2D* FIMG_tof_amp_fOut_det1 = 0;
TH2D* FIMG_tof_amp_fOut_det2 = 0;
TH2D* tof_area_hist_fOut_det1 = 0;
TH2D* tof_area_hist_fOut_det2 = 0;

TH2D* FIMG_tof_amp_fIn_det1 = 0;
TH2D* FIMG_tof_amp_fIn_det2 = 0;
TH2D* tof_area_hist_fIn_det1 = 0;
TH2D* tof_area_hist_fIn_det2 = 0;

// TH2D* tof_area_hist_in_cutoff = 0;

// TH2D* FIMG_tof_amp_fOut_det1_cutoff = 0;
// TH2D* FIMG_tof_amp_fOut_det2_cutoff = 0;
// TH2D* FIMG_tof_amp_fIn_det1_cutoff = 0;
// TH2D* FIMG_tof_amp_fIn_det2_cutoff = 0;

//dedicated and parasitic pulses plots
TH2D* FIMG_tof_amp_fIn_det1_dedicated = 0;
TH2D* FIMG_tof_amp_fIn_det1_parasitic = 0;
TH2D* FIMG_tof_amp_fIn_det2_dedicated = 0;
TH2D* FIMG_tof_amp_fIn_det2_parasitic = 0;

// //beam off
// TH2D* tof_amp_beam_off_PTBC = 0;
// TH2D* tof_amp_beam_off_FIMG = 0;

// // Cut off
// auto cutoff_amp = new TGraph();
// auto cutoff_tof = new TGraph();

Double_t flight_path_length_PTB = 182.65 - 0.41; //m
Double_t flight_path_length_FIMG = 183.5 - 0.41; //m
Double_t neutron_mass = 939.56542052; //in MeV
Double_t speed_of_light = 299792458.0; //in m/s
Double_t delT_pkup_ptbc = 660.0; //in ns
Double_t delT_pkup_fimg = 630.0; //in ns

Double_t FIMG_tof_cut_low_det1 = 8102.0; //in ns
Double_t FIMG_tof_cut_up_det1 = 100000.0; //in ns
Double_t FIMG_tof_cut_low_det2 = 8174.0; //in ns
Double_t FIMG_tof_cut_up_det2 = 89408.0; //in ns

Double_t FIMG_min_amp_cut_det1 = 600.0; //a.u.
Double_t FIMG_min_amp_cut_det2 = 400.0; //a.u.

Double_t t_gamma_PTB = (flight_path_length_PTB / speed_of_light) * 1e9; //converting into ns
Double_t t_gamma_FIMG = (flight_path_length_FIMG / speed_of_light) * 1e9; //converting into ns

Double_t t11 = 780.0;
Double_t a11 = 7700.0;

Double_t t12 = 1090.0;
Double_t a12 = 5451.0;

Double_t t21 = 1090.0;
Double_t a21 = 5451.0;

Double_t t22 = 2605.0;
Double_t a22 = 7167.0;

Double_t t31 = 2605.0;
Double_t a31 = 7960.0;

Double_t t32 = 2856.0;
Double_t a32 = 7960.0;

Double_t t41 = 2856.0;
Double_t a41 = 7432.0;

Double_t t42 = 15290.0;
Double_t a42 = 7432.0;

Double_t t51 = 15290.0;
Double_t a51 = 7432.0;

Double_t t52 = 18708.0;
Double_t a52 = 3600.0; //2416

Double_t t61 = 18708.0;
Double_t a61 = 3600.0; //2416

Double_t t62 = 1e8;
Double_t a62 = 3600.0; //3076

Double_t yOnTheCutLine(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t x3){
    return ((y2 - y1)*(x3 - x1)/(x2 - x1) + y1);
}

///////// FIMG stats from 117386
//Filter In runs
//Bi (1 cm) Filter
std::vector<Int_t> c_au_f_bi_t_out = {117386, 117387, 117388, 117389, 117390, 117436};
std::vector<Int_t> c_out_f_bi_t_out = {117391, 117392, 117393, 117394, 117395, 117396, 117397};
std::vector<Int_t> c_pb_f_bi_t_out = {117422, 117423, 117424, 117425, 117426, 117427, 117428};
std::vector<Int_t> c_ta_f_bi_t_out = {117437, 117439, 117440, 117441, 117442, 117443};
std::vector<Int_t> filter_in_runs;

//Al (8 cm) Filter
std::vector<Int_t> c_ta_f_al8_t_out = {117449, 117450, 117451, 117453, 117454, 117455, 117456, 117457, 117458, 117459, 117460};

//Al (5 cm) Filter
std::vector<Int_t> c_ta_f_al5_t_out = {117462, 117463, 117464, 117465, 117466, 117467, 117476, 117477, 117479, 117480, 117481, 117482, 117483, 117484}; //117461, 117478, 117483, 117484
std::vector<Int_t> c_out_f_al5_t_out = {117470, 117471, 117472, 117473, 117474, 117475};
std::vector<Int_t> c_pb_f_al5_t_out = {117497, 117498, 117499, 117500, 117501, 117502};
std::vector<Int_t> c_c_f_al5_t_out = {117503, 117504, 117505, 117506, 117507, 117508, 117509, 117510};
std::vector<Int_t> c_au_f_al5_t_out = {117519, 117520, 117521, 117522, 117530, 117531, 117532};

//Filter Out runs
std::vector<Int_t> filter_out_runs;
std::vector<Int_t> c_au_f_out_t_out = {117350, 117357, 117358, 117359, 117362, 117363, 117364, 117365, 117366, 117367, 117368};
std::vector<Int_t> c_atic_f_out_t_out = {117351, 117355, 117356}; // 117352, 117353,
std::vector<Int_t> c_out_f_out_t_out = {117398, 117405, 117406, 117408, 117409, 117410, 117411, 117412, 117398}; //
std::vector<Int_t> c_pb_f_out_t_out = {117429, 117430, 117431, 117432, 117433, 117434, 117435};
std::vector<Int_t> c_ta_f_out_t_out = {117444, 117445, 117446, 117447, 117448}; //117452
std::vector<Int_t> c_c_f_out_t_out = {117511, 117512, 117513, 117514, 117515, 117516, 117517, 117518};

//Beam Off runs
std::vector<Int_t> beam_off = {117740, 117741, 117742, 117743, 117744, 117745};

Double_t fimgCutFunction(Double_t x, Int_t det_num){
    // // A * (log(x))^2 + B + log(x) + C
    // // Values obtained from cutoffFitter.C
    // Double_t A = 120556;
    // Double_t B = -2.75589e6;
    // Double_t C = 1.57536e7;
    // return (A * TMath::Log(x) * TMath::Log(x) + B * TMath::Log(x) + C);

    if(det_num == 1) {
        // A / ( TMath::Log(x) + B )
        // Values obtained from cutoffFitter.C
        Double_t A = 4378.11;
        Double_t B = -7.695;

        return (A / (TMath::Log(x) + B));
    }

    if(det_num == 2) {
        // A / ( TMath::Log(x) + B )
        // Values obtained from cutoffFitter.C
        Double_t A = 4145.31;
        Double_t B = -7.67962;

        return (A / (TMath::Log(x) + B));
    }

    return 0;
}

void fillRuns(){
    // filter_in_runs.reserve( c_au_f_bi_t_out.size() + c_out_f_bi_t_out.size() + c_pb_f_bi_t_out.size() + c_ta_f_bi_t_out.size() );
    // filter_in_runs.insert( filter_in_runs.end(), c_au_f_bi_t_out.begin(), c_au_f_bi_t_out.end() );
    // filter_in_runs.insert( filter_in_runs.end(), c_out_f_bi_t_out.begin(), c_out_f_bi_t_out.end() );
    // filter_in_runs.insert( filter_in_runs.end(), c_pb_f_bi_t_out.begin(), c_pb_f_bi_t_out.end() );
    // filter_in_runs.insert( filter_in_runs.end(), c_ta_f_bi_t_out.begin(), c_ta_f_bi_t_out.end() );

    // filter_in_runs.reserve( c_ta_f_al8_t_out.size() );
    // filter_in_runs.insert( filter_in_runs.end(), c_ta_f_al8_t_out.begin(), c_ta_f_al8_t_out.end() );

    filter_in_runs.insert( filter_in_runs.end(), c_ta_f_al5_t_out.begin(), c_ta_f_al5_t_out.end() );
    filter_in_runs.insert( filter_in_runs.end(), c_out_f_al5_t_out.begin(), c_out_f_al5_t_out.end() );
    filter_in_runs.insert( filter_in_runs.end(), c_pb_f_al5_t_out.begin(), c_pb_f_al5_t_out.end() );
    filter_in_runs.insert( filter_in_runs.end(), c_c_f_al5_t_out.begin(), c_c_f_al5_t_out.end() );
    filter_in_runs.insert( filter_in_runs.end(), c_au_f_al5_t_out.begin(), c_au_f_al5_t_out.end() );

    // filter_out_runs.insert( filter_out_runs.end(), c_au_f_out_t_out.begin(), c_au_f_out_t_out.end() );
    // filter_out_runs.insert( filter_out_runs.end(), c_atic_f_out_t_out.begin(), c_atic_f_out_t_out.end() );
    filter_out_runs.insert( filter_out_runs.end(), c_out_f_out_t_out.begin(), c_out_f_out_t_out.end() );
    filter_out_runs.insert( filter_out_runs.end(), c_pb_f_out_t_out.begin(), c_pb_f_out_t_out.end() );
    filter_out_runs.insert( filter_out_runs.end(), c_ta_f_out_t_out.begin(), c_ta_f_out_t_out.end() );
    filter_out_runs.insert( filter_out_runs.end(), c_c_f_out_t_out.begin(), c_c_f_out_t_out.end() );
}

// void beamOff(){
//     for (int i = 0; i < beam_off.size(); i++)
//     {
//         TFile *file_ntof = TFile::Open(Form("/eos/experiment/ntof/data/rootfiles/2023/ear1/run%d.root", beam_off.at(i)),"read");

//         //PTBC ---------------------------------------------
//         TTree* PTBC;
//         Double_t tof_PTB = 0; //tof is in ns
//         Float_t amp_PTB = 0;
//         Int_t BunchNumber_PTB = 0;
//         Float_t PulseIntensity = 0;

//         file_ntof->GetObject("PTBC", PTBC);
//         PTBC->SetBranchAddress("BunchNumber", &BunchNumber_PTB);
//         PTBC->SetBranchAddress("PulseIntensity", &PulseIntensity);
//         PTBC->SetBranchAddress("tof", &tof_PTB);
//         PTBC->SetBranchAddress("amp", &amp_PTB);

//         Long64_t Events_PTB = PTBC->GetEntriesFast();
//         std::cout << "Number of entries - PTBC = " << Events_PTB << std::endl;
        
//         int CurrentBunchNum = 0;

//         for (int j = 0; j < Events_PTB; j++)
//         {
//             PTBC->GetEntry(j);

//             Double_t t_pkup = BNum_tpkup_map[BunchNumber_PTB];
//             Double_t corrected_tof = tof_PTB - t_pkup + delT_pkup_ptbc + t_gamma_PTB;

//             if (CurrentBunchNum != BunchNumber_PTB)
//             {
//                 CurrentBunchNum = BunchNumber_PTB;
//                 NormFactor += (Double_t) PulseIntensity;
//             }

//             //Filling the histograms
//             for (int k = 0; k < 6; k++)
//             {
//                 if (corrected_tof >= t[k][0] && corrected_tof < t[k][1])
//                 {
//                     if ( (Double_t) amp_PTB > yOnTheCutLinePTBC(t[k][0], a[k][0], t[k][1], a[k][1], corrected_tof) )
//                     {
//                         energy_hist_PTB->Fill( TOFToEnergy(corrected_tof * 1e-9) );
//                         break;    
//                     }
//                 }
//             }
//         }

//         //FIMG ---------------------------------------------
//         TTree* FIMG;
//         Double_t tof_FIMG = 0; //tof is in ns
//         Float_t amp = 0;
//         Float_t area_0 = 0;
//         Int_t det_num = 0;
//         Int_t BunchNumber_FIMG = 0;

//         file_ntof->GetObject("FIMG", FIMG);
//         FIMG->SetBranchAddress("BunchNumber", &BunchNumber_FIMG);
//         FIMG->SetBranchAddress("tof", &tof_FIMG);
//         FIMG->SetBranchAddress("amp", &amp);
//         FIMG->SetBranchAddress("area_0", &area_0);
//         FIMG->SetBranchAddress("detn", &det_num);
//     }
// }

void FilterIn(){

    for (int i = 0; i < filter_in_runs.size(); i++)
    {
        TFile *file_ntof = TFile::Open(Form("/eos/experiment/ntof/data/rootfiles/2023/ear1/run%d.root", filter_in_runs.at(i)),"read");
        
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

        //FIMG ---------------------------------------------
        TTree* FIMG;
        Double_t tof = 0; //tof is in ns
        Double_t tflash = 0; //tflash is in ns
        Float_t amp = 0;
        Float_t area_0 = 0;
        Int_t det_num = 0;
        Int_t BunchNumber_FIMG = 0;
        Float_t PulseIntensity = 0;

        file_ntof->GetObject("FIMG", FIMG);
        FIMG->SetBranchAddress("BunchNumber", &BunchNumber_FIMG);
        FIMG->SetBranchAddress("PulseIntensity", &PulseIntensity);
        FIMG->SetBranchAddress("tof", &tof);
        FIMG->SetBranchAddress("tflash", &tflash);
        FIMG->SetBranchAddress("amp", &amp);
        FIMG->SetBranchAddress("area_0", &area_0);
        FIMG->SetBranchAddress("detn", &det_num);

        Long64_t Events_FIMG = FIMG->GetEntriesFast();
        std::cout << "Number of entries - Filter In = " << Events_FIMG << std::endl;
        
        int CurrentBunchNum = 0;

        for (int j = 0; j < Events_FIMG; j++)
        {
            FIMG->GetEntry(j);

            // Double_t corrected_tof = (tof - tflash + t_gamma_FIMG);
            Double_t t_pkup = BNum_tpkup_map[BunchNumber_FIMG];
            Double_t corrected_tof = tof - t_pkup + 630.0 + t_gamma_FIMG;
            // FIMG_tof_amp_out->Fill( corrected_tof , (Double_t) amp);

            // Double_t tof_cut_low = 0;
            // Double_t tof_cut_up = 0;
            // Double_t min_amp_cut = 0;

            if (det_num == 1){
                // tof_cut_low = FIMG_tof_cut_low_det1;
                // tof_cut_up = FIMG_tof_cut_up_det1;
                // min_amp_cut = FIMG_min_amp_cut_det1;

                if (PulseIntensity > 6e12)
                {
                    FIMG_tof_amp_fIn_det1_dedicated->Fill( corrected_tof , (Double_t) amp);
                } else if (PulseIntensity <= 6e12)
                {
                    FIMG_tof_amp_fIn_det1_parasitic->Fill( corrected_tof , (Double_t) amp);
                }

                tof_area_hist_fIn_det1->Fill( corrected_tof , (Double_t) area_0);
                FIMG_tof_amp_fIn_det1->Fill( corrected_tof , (Double_t) amp);
            }
            else if (det_num == 2){
                // tof_cut_low = FIMG_tof_cut_low_det2;
                // tof_cut_up = FIMG_tof_cut_up_det2;
                // min_amp_cut = FIMG_min_amp_cut_det2;

                if (PulseIntensity > 6e12)
                {
                    FIMG_tof_amp_fIn_det2_dedicated->Fill( corrected_tof , (Double_t) amp);
                } else if (PulseIntensity <= 6e12)
                {
                    FIMG_tof_amp_fIn_det2_parasitic->Fill( corrected_tof , (Double_t) amp);
                }

                tof_area_hist_fIn_det2->Fill( corrected_tof , (Double_t) area_0);
                FIMG_tof_amp_fIn_det2->Fill( corrected_tof , (Double_t) amp);
            }

            // //Filling the histograms after cuts
            // if (corrected_tof < tof_cut_low)
            // {
            //     continue;
            // }

            // if ((Double_t) amp < min_amp_cut)
            // {
            //     continue;
            // }

            // if (corrected_tof > tof_cut_up)
            // {
            //     if (det_num == 1)
            //     {
            //         FIMG_tof_amp_fIn_det1_cutoff->Fill( corrected_tof , (Double_t) amp );
            //     }

            //     if (det_num == 2)
            //     {
            //         FIMG_tof_amp_fIn_det2_cutoff->Fill( corrected_tof , (Double_t) amp );
            //     }
            //     continue;
            // }

            // if ((Double_t) amp >= fimgCutFunction(corrected_tof, det_num))
            // {
            //     if (det_num == 1)
            //     {
            //         FIMG_tof_amp_fIn_det1_cutoff->Fill( corrected_tof , (Double_t) amp );
            //     }

            //     if (det_num == 2)
            //     {
            //         FIMG_tof_amp_fIn_det2_cutoff->Fill( corrected_tof , (Double_t) amp );
            //     }
            // }
        }

        file_ntof->Close();
    }
}

void FilterOut(){

    for (int i = 0; i < filter_out_runs.size(); i++)
    {
        TFile *file_ntof = TFile::Open(Form("/eos/experiment/ntof/data/rootfiles/2023/ear1/run%d.root", filter_out_runs.at(i)),"read");
        
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

        //FIMG ---------------------------------------------
        TTree* FIMG;
        Double_t tof = 0; //tof is in ns
        Double_t tflash = 0; //tflash is in ns
        Float_t amp = 0;
        Float_t area_0 = 0;
        Int_t det_num = 0;
        Int_t BunchNumber_FIMG = 0;
        Float_t PulseIntensity = 0;

        file_ntof->GetObject("FIMG", FIMG);
        FIMG->SetBranchAddress("BunchNumber", &BunchNumber_FIMG);
        FIMG->SetBranchAddress("PulseIntensity", &PulseIntensity);
        FIMG->SetBranchAddress("tof", &tof);
        FIMG->SetBranchAddress("tflash", &tflash);
        FIMG->SetBranchAddress("amp", &amp);
        FIMG->SetBranchAddress("area_0", &area_0);
        FIMG->SetBranchAddress("detn", &det_num);

        Long64_t Events_FIMG = FIMG->GetEntriesFast();
        std::cout << "Number of entries - Filter Out = " << Events_FIMG << std::endl;
        
        int CurrentBunchNum = 0;

        for (int j = 0; j < Events_FIMG; j++)
        {
            FIMG->GetEntry(j);

            // Double_t corrected_tof = (tof - tflash + t_gamma_FIMG);
            Double_t t_pkup = BNum_tpkup_map[BunchNumber_FIMG];
            Double_t corrected_tof = tof - t_pkup + 630.0 + t_gamma_FIMG;
            // FIMG_tof_amp_out->Fill( corrected_tof , (Double_t) amp);

            // Double_t tof_cut_low = 0;
            // Double_t tof_cut_up = 0;
            // Double_t min_amp_cut = 0;

            if (det_num == 1){
                // tof_cut_low = FIMG_tof_cut_low_det1;
                // tof_cut_up = FIMG_tof_cut_up_det1;
                // min_amp_cut = FIMG_min_amp_cut_det1;

                tof_area_hist_fOut_det1->Fill( corrected_tof , (Double_t) area_0);
                FIMG_tof_amp_fOut_det1->Fill( corrected_tof , (Double_t) amp);
            }
            else if (det_num == 2){
                // tof_cut_low = FIMG_tof_cut_low_det2;
                // tof_cut_up = FIMG_tof_cut_up_det2;
                // min_amp_cut = FIMG_min_amp_cut_det2;

                tof_area_hist_fOut_det2->Fill( corrected_tof , (Double_t) area_0);
                FIMG_tof_amp_fOut_det2->Fill( corrected_tof , (Double_t) amp);
            }

            // //Filling the histograms after cuts
            // if (corrected_tof < tof_cut_low)
            // {
            //     continue;
            // }

            // if ((Double_t) amp < min_amp_cut)
            // {
            //     continue;
            // }

            // if (corrected_tof > tof_cut_up)
            // {
            //     if (det_num == 1)
            //     {
            //         FIMG_tof_amp_fOut_det1_cutoff->Fill( corrected_tof , (Double_t) amp );
            //     }

            //     if (det_num == 2)
            //     {
            //         FIMG_tof_amp_fOut_det2_cutoff->Fill( corrected_tof , (Double_t) amp );
            //     }
            //     continue;
            // }

            // if ((Double_t) amp >= fimgCutFunction(corrected_tof, det_num))
            // {
            //     if (det_num == 1)
            //     {
            //         FIMG_tof_amp_fOut_det1_cutoff->Fill( corrected_tof , (Double_t) amp );
            //     }

            //     if (det_num == 2)
            //     {
            //         FIMG_tof_amp_fOut_det2_cutoff->Fill( corrected_tof , (Double_t) amp );
            //     }
            // }
        }

        file_ntof->Close();
    }
}

void StoreHist(){

    TFile *f = new TFile("../rootFiles/cutoffSelector_FIMG.root","recreate");

    // FIMG_tof_amp_in->Write();
    // FIMG_tof_amp_out->Write();
    // tof_area_hist_in->Write();
    // tof_area_hist_in_cutoff->Write();
    tof_area_hist_fOut_det1->Write();
    tof_area_hist_fOut_det2->Write();
    FIMG_tof_amp_fOut_det1->Write();
    FIMG_tof_amp_fOut_det2->Write();

    tof_area_hist_fIn_det1->Write();
    tof_area_hist_fIn_det2->Write();
    FIMG_tof_amp_fIn_det1->Write();
    FIMG_tof_amp_fIn_det2->Write();
    // tof_area_hist_out_cutoff->Write();

    // FIMG_tof_amp_fIn_det1_cutoff->Write();
    // FIMG_tof_amp_fIn_det2_cutoff->Write();
    // FIMG_tof_amp_fOut_det1_cutoff->Write();
    // FIMG_tof_amp_fOut_det2_cutoff->Write();

    FIMG_tof_amp_fIn_det1_dedicated->Write();
    FIMG_tof_amp_fIn_det1_parasitic->Write();
    FIMG_tof_amp_fIn_det2_dedicated->Write();
    FIMG_tof_amp_fIn_det2_parasitic->Write();

    f->Close();
}

void plots(){

    TCanvas *c[2];

    int i = 0;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();
    tof_area_hist_fOut_det1->GetXaxis()->SetTitle("Time of Flight (in ns)");
    tof_area_hist_fOut_det1->GetYaxis()->SetTitle("Area (a.u.)");
    tof_area_hist_fOut_det1->SetTitle("ToF vs Area Hist - FIMG Det 1 - Filter Out");
    tof_area_hist_fOut_det1->Draw("colz");
    // tof_area_hist_fOut_det1->SetMarkerStyle(6);
    // tof_area_hist_fOut_det1->SetMarkerSize(0.5);
    gPad->SetLogx();
    gPad->SetLogz();
    gStyle->SetPalette(57);

    c[i]->Print("../plots/h_tof_area_fOut_FIMG_det1.png");

    i++;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();
    tof_area_hist_fOut_det2->GetXaxis()->SetTitle("Time of Flight (in ns)");
    tof_area_hist_fOut_det2->GetYaxis()->SetTitle("Area (a.u.)");
    tof_area_hist_fOut_det2->SetTitle("ToF vs Area Hist - FIMG Det 2 - Filter Out");
    tof_area_hist_fOut_det2->Draw("colz");
    // tof_area_hist_fOut_det2->SetMarkerStyle(6);
    // tof_area_hist_fOut_det2->SetMarkerSize(0.5);
    gPad->SetLogx();
    gPad->SetLogz();
    gStyle->SetPalette(57);

    c[i]->Print("../plots/h_tof_area_fOut_FIMG_det2.png");

}

void cutoffSelector(){

    gStyle->SetStatX(0.27);
    gStyle->SetStatY(0.9);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.17);

    fillRuns();

    //Calculating TOF (x) bin edges
    int bins_per_decade = 100;
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
    int num_bins_amp = 200;
    Double_t amp_min = 0.;
    Double_t amp_max = 4500.;
    Double_t bin_edges_amp[num_bins_amp+1];
    Double_t step_amp = (Double_t) ((amp_max-amp_min)/num_bins_amp);
    // std::cout << "step_amp = " << step_amp << std::endl;
    for(int i = 0; i < num_bins_amp+1; i++)
    {
        bin_edges_amp[i] = step_amp * (Double_t) i;
    }

    //Calculating area (y) bin edges
    int num_bins_area = 200;
    Double_t area_min = 0.;
    Double_t area_max = 8e5;
    Double_t bin_edges_area[num_bins_area+1];
    Double_t step_area = (Double_t) ((area_max-area_min)/num_bins_area);
    // std::cout << "step_area = " << step_area << std::endl;
    for(int i = 0; i < num_bins_area+1; i++)
    {
        bin_edges_area[i] = step_area * (Double_t) i;
    }

    // FIMG_tof_amp_in = new TH2D("FIMG_tof_amp_in","ToF vs Amplitude Hist - FIMG - Al (5cm) Filter In",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    // tof_area_hist_in = new TH2D("tof_area_hist_in","Area vs Area Hist - FIMG - Al (5cm) Filter In",num_bins_tof,bin_edges_tof,num_bins_area,bin_edges_area);
    // tof_area_hist_in_cutoff = new TH2D("tof_area_hist_in_cutoff","ToF vs Area Hist - FIMG - Al (5cm) After Cuts",num_bins_tof,bin_edges_tof,num_bins_area,bin_edges_area);
    // FIMG_tof_amp_out = new TH2D("FIMG_tof_amp_out","ToF vs Amplitude Hist - FIMG - Filter Out",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);

    //Filter Out
    tof_area_hist_fOut_det1 = new TH2D("tof_area_hist_fOut_det1","ToF vs Area Hist - FIMG Det 1 - Filter Out",num_bins_tof,bin_edges_tof,num_bins_area,bin_edges_area);
    tof_area_hist_fOut_det2 = new TH2D("tof_area_hist_fOut_det2","ToF vs Area Hist - FIMG Det 2 - Filter Out",num_bins_tof,bin_edges_tof,num_bins_area,bin_edges_area);

    FIMG_tof_amp_fOut_det1 = new TH2D("FIMG_tof_amp_fOut_det1","ToF vs Amplitude Hist - FIMG Det 1 - Filter Out",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    FIMG_tof_amp_fOut_det2 = new TH2D("FIMG_tof_amp_fOut_det2","ToF vs Amplitude Hist - FIMG Det 2 - Filter Out",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    
    //Filter In
    tof_area_hist_fIn_det1 = new TH2D("tof_area_hist_fIn_det1","ToF vs Area Hist - FIMG Det 1 - Al (5cm) Filter",num_bins_tof,bin_edges_tof,num_bins_area,bin_edges_area);
    tof_area_hist_fIn_det2 = new TH2D("tof_area_hist_fIn_det2","ToF vs Area Hist - FIMG Det 2 - Al (5cm) Filter",num_bins_tof,bin_edges_tof,num_bins_area,bin_edges_area);

    FIMG_tof_amp_fIn_det1 = new TH2D("FIMG_tof_amp_fIn_det1","ToF vs Amplitude Hist - FIMG Det 1 - Al (5cm) Filter",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    FIMG_tof_amp_fIn_det2 = new TH2D("FIMG_tof_amp_fIn_det2","ToF vs Amplitude Hist - FIMG Det 2 - Al (5cm) Filter",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    
    // tof_area_hist_out_cutoff = new TH2D("tof_area_hist_out_cutoff","ToF vs Area Hist - FIMG - Filter Out After Cuts",num_bins_tof,bin_edges_tof,num_bins_area,bin_edges_area);

    // FIMG_tof_amp_fOut_det1_cutoff = new TH2D("FIMG_tof_amp_fOut_det1_cutoff","ToF vs Amp Hist - FIMG Det 1 - Filter Out Cutoff",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    // FIMG_tof_amp_fOut_det2_cutoff = new TH2D("FIMG_tof_amp_fOut_det2_cutoff","ToF vs Amp Hist - FIMG Det 2 - Filter Out Cutoff",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);

    // FIMG_tof_amp_fIn_det1_cutoff = new TH2D("FIMG_tof_amp_fIn_det1_cutoff","ToF vs Amp Hist - FIMG Det 1 - Al (5cm) Filter Cutoff",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    // FIMG_tof_amp_fIn_det2_cutoff = new TH2D("FIMG_tof_amp_fIn_det2_cutoff","ToF vs Amp Hist - FIMG Det 2 - Al (5cm) Filter Cutoff",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);

    FIMG_tof_amp_fIn_det1_dedicated = new TH2D("FIMG_tof_amp_fIn_det1_dedicated","ToF vs Amp Hist - FIMG Det 1 - Al(5cm) - Dedicated",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    FIMG_tof_amp_fIn_det1_parasitic = new TH2D("FIMG_tof_amp_fIn_det1_parasitic","ToF vs Amp Hist - FIMG Det 1 - Al(5cm) - Parasitic",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    FIMG_tof_amp_fIn_det2_dedicated = new TH2D("FIMG_tof_amp_fIn_det2_dedicated","ToF vs Amp Hist - FIMG Det 2 - Al(5cm) - Dedicated",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    FIMG_tof_amp_fIn_det2_parasitic = new TH2D("FIMG_tof_amp_fIn_det2_parasitic","ToF vs Amp Hist - FIMG Det 2 - Al(5cm) - Parasitic",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);

    
    FilterIn();
    FilterOut();

    StoreHist();
}