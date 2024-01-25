/**
 * @file cutoffAnalysis_PTBC.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2024-01-17
 */

#include "cutoffAnalysis_PTBC.h"

void FilterIn(){

    for (int i = 0; i < filter_in_runs.size(); i++)
    {
        TFile *file_ntof = TFile::Open(Form("/eos/experiment/ntof/processing/official/done/run%d.root", filter_in_runs.at(i)),"read");
        
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
        std::cout << "Number of entries - PTBC = " << Events_PTBC << std::endl;
        
        int CurrentBunchNum = 0;

        for (int j = 0; j < Events_PTBC; j++)
        {
            PTBC->GetEntry(j);

            Double_t t_pkup = BNum_tpkup_map[BunchNumber_PTBC];
            Double_t corrected_tof = tof_PTBC - t_pkup + delT_pkup_ptbc + t_gamma_PTBC;

            if (CurrentBunchNum != BunchNumber_PTBC)
            {
                CurrentBunchNum = BunchNumber_PTBC;
            }

            PTBC_tof_amp_fIn_total->Fill(corrected_tof, (Double_t) amp_PTBC);

            if (PulseIntensity_PTBC > 6e12)
            {
                PTBC_tof_amp_fIn_dedicated->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (PulseIntensity_PTBC <= 6e12)
            {
                PTBC_tof_amp_fIn_parasitic->Fill(corrected_tof, (Double_t) amp_PTBC);
            }

            if (det_num_PTBC == 2) {
                PTBC_tof_amp_fIn_det2->Fill(corrected_tof, (Double_t) amp_PTBC);
                for (int k = 0; k < 4; k++)
                {
                    if (corrected_tof >= t_det2[k][0] && corrected_tof < t_det2[k][1])
                    {
                        if ( (Double_t) amp_PTBC > yOnTheCutLinePTBC(t_det2[k][0], a_det2[k][0], t_det2[k][1], a_det2[k][1], corrected_tof) )
                        {
                            PTBC_tof_amp_fIn_det2_afterCuts->Fill(corrected_tof, (Double_t) amp_PTBC);
                            break;    
                        }
                    }
                }
            } else if (det_num_PTBC == 3) {
                PTBC_tof_amp_fIn_det3->Fill(corrected_tof, (Double_t) amp_PTBC);
                for (int k = 0; k < 2; k++)
                {
                    if (corrected_tof >= t_det3to7[det_num_PTBC-3][k][0] && corrected_tof < t_det3to7[det_num_PTBC-3][k][1])
                    {
                        if ( (Double_t) amp_PTBC > yOnTheCutLinePTBC(t_det3to7[det_num_PTBC-3][k][0], a_det3to7[det_num_PTBC-3][k][0], t_det3to7[det_num_PTBC-3][k][1], a_det3to7[det_num_PTBC-3][k][1], corrected_tof) )
                        {
                            PTBC_tof_amp_fIn_det3_afterCuts->Fill(corrected_tof, (Double_t) amp_PTBC);
                            break;    
                        }
                    }
                }
            } else if (det_num_PTBC == 4) {
                PTBC_tof_amp_fIn_det4->Fill(corrected_tof, (Double_t) amp_PTBC);
                for (int k = 0; k < 2; k++)
                {
                    if (corrected_tof >= t_det3to7[det_num_PTBC-3][k][0] && corrected_tof < t_det3to7[det_num_PTBC-3][k][1])
                    {
                        if ( (Double_t) amp_PTBC > yOnTheCutLinePTBC(t_det3to7[det_num_PTBC-3][k][0], a_det3to7[det_num_PTBC-3][k][0], t_det3to7[det_num_PTBC-3][k][1], a_det3to7[det_num_PTBC-3][k][1], corrected_tof) )
                        {
                            PTBC_tof_amp_fIn_det4_afterCuts->Fill(corrected_tof, (Double_t) amp_PTBC);
                            break;    
                        }
                    }
                }
            } else if (det_num_PTBC == 5) {
                PTBC_tof_amp_fIn_det5->Fill(corrected_tof, (Double_t) amp_PTBC);
                for (int k = 0; k < 2; k++)
                {
                    if (corrected_tof >= t_det3to7[det_num_PTBC-3][k][0] && corrected_tof < t_det3to7[det_num_PTBC-3][k][1])
                    {
                        if ( (Double_t) amp_PTBC > yOnTheCutLinePTBC(t_det3to7[det_num_PTBC-3][k][0], a_det3to7[det_num_PTBC-3][k][0], t_det3to7[det_num_PTBC-3][k][1], a_det3to7[det_num_PTBC-3][k][1], corrected_tof) )
                        {
                            PTBC_tof_amp_fIn_det5_afterCuts->Fill(corrected_tof, (Double_t) amp_PTBC);
                            break;    
                        }
                    }
                }
            } else if (det_num_PTBC == 6) {
                PTBC_tof_amp_fIn_det6->Fill(corrected_tof, (Double_t) amp_PTBC);
                for (int k = 0; k < 2; k++)
                {
                    if (corrected_tof >= t_det3to7[det_num_PTBC-3][k][0] && corrected_tof < t_det3to7[det_num_PTBC-3][k][1])
                    {
                        if ( (Double_t) amp_PTBC > yOnTheCutLinePTBC(t_det3to7[det_num_PTBC-3][k][0], a_det3to7[det_num_PTBC-3][k][0], t_det3to7[det_num_PTBC-3][k][1], a_det3to7[det_num_PTBC-3][k][1], corrected_tof) )
                        {
                            PTBC_tof_amp_fIn_det6_afterCuts->Fill(corrected_tof, (Double_t) amp_PTBC);
                            break;    
                        }
                    }
                }
            } else if (det_num_PTBC == 7) {
                PTBC_tof_amp_fIn_det7->Fill(corrected_tof, (Double_t) amp_PTBC);
                for (int k = 0; k < 2; k++)
                {
                    if (corrected_tof >= t_det3to7[det_num_PTBC-3][k][0] && corrected_tof < t_det3to7[det_num_PTBC-3][k][1])
                    {
                        if ( (Double_t) amp_PTBC > yOnTheCutLinePTBC(t_det3to7[det_num_PTBC-3][k][0], a_det3to7[det_num_PTBC-3][k][0], t_det3to7[det_num_PTBC-3][k][1], a_det3to7[det_num_PTBC-3][k][1], corrected_tof) )
                        {
                            PTBC_tof_amp_fIn_det7_afterCuts->Fill(corrected_tof, (Double_t) amp_PTBC);
                            break;    
                        }
                    }
                }
            }

            // if (det_num_PTBC != 1 && det_num_PTBC != 8)
            // {
            //     PTBC_tof_amp_fIn_dets[det_num_PTBC-2]->Fill(corrected_tof, (Double_t) amp_PTBC);
            // }
        }

        file_ntof->Close();
    }
}

void FilterOut(){

    for (int i = 0; i < filter_out_runs.size(); i++)
    {
        TFile *file_ntof = TFile::Open(Form("/eos/experiment/ntof/processing/official/done/run%d.root", filter_out_runs.at(i)),"read");
        
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
        std::cout << "Number of entries - PTBC = " << Events_PTBC << std::endl;
        
        int CurrentBunchNum = 0;

        for (int j = 0; j < Events_PTBC; j++)
        {
            PTBC->GetEntry(j);

            if (det_num_PTBC == 1 || det_num_PTBC == 8) {
                continue;
            }

            Double_t t_pkup = BNum_tpkup_map[BunchNumber_PTBC];
            Double_t corrected_tof = tof_PTBC - t_pkup + delT_pkup_ptbc + t_gamma_PTBC;

            if (CurrentBunchNum != BunchNumber_PTBC)
            {
                CurrentBunchNum = BunchNumber_PTBC;
            }

            //Filling the histograms
            // for (int k = 0; k < 6; k++)
            // {
            //     if (corrected_tof >= t[k][0] && corrected_tof < t[k][1])
            //     {
            //         if ( (Double_t) amp_PTBC > yOnTheCutLinePTBC(t[k][0], a[k][0], t[k][1], a[k][1], corrected_tof) )
            //         {
            //             energy_hist_PTBC->Fill( TOFToEnergy(corrected_tof * 1e-9) );
            //             break;    
            //         }
            //     }
            // }

            PTBC_tof_amp_fOut_total->Fill(corrected_tof, (Double_t) amp_PTBC);

            if (PulseIntensity_PTBC > 6e12)
            {
                PTBC_tof_amp_fOut_dedicated->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (PulseIntensity_PTBC <= 6e12)
            {
                PTBC_tof_amp_fOut_parasitic->Fill(corrected_tof, (Double_t) amp_PTBC);
            }

            if (det_num_PTBC == 2) {
                PTBC_tof_amp_fOut_det2->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (det_num_PTBC == 3) {
                PTBC_tof_amp_fOut_det3->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (det_num_PTBC == 4) {
                PTBC_tof_amp_fOut_det4->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (det_num_PTBC == 5) {
                PTBC_tof_amp_fOut_det5->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (det_num_PTBC == 6) {
                PTBC_tof_amp_fOut_det6->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (det_num_PTBC == 7) {
                PTBC_tof_amp_fOut_det7->Fill(corrected_tof, (Double_t) amp_PTBC);
            }

            // if (det_num_PTBC != 1 && det_num_PTBC != 8)
            // {
            //     PTBC_tof_amp_fOut_dets[det_num_PTBC-2]->Fill(corrected_tof, (Double_t) amp_PTBC);
            // }
        }

        file_ntof->Close();
    }
}

void StoreHist(){
    
    TFile *f = new TFile(Form("../rootFiles/cutoffAnalysis_PTBC_%s.root", target_name.c_str()),"recreate");

    // PTBC_tof_amp_in->Write();
    // PTBC_tof_amp_out->Write();
    PTBC_tof_amp_fOut_total->Write();
    PTBC_tof_amp_fOut_det2->Write();
    PTBC_tof_amp_fOut_det3->Write();
    PTBC_tof_amp_fOut_det4->Write();
    PTBC_tof_amp_fOut_det5->Write();
    PTBC_tof_amp_fOut_det6->Write();
    PTBC_tof_amp_fOut_det7->Write();

    PTBC_tof_amp_fIn_total->Write();
    PTBC_tof_amp_fIn_det2->Write();
    PTBC_tof_amp_fIn_det3->Write();
    PTBC_tof_amp_fIn_det4->Write();
    PTBC_tof_amp_fIn_det5->Write();
    PTBC_tof_amp_fIn_det6->Write();
    PTBC_tof_amp_fIn_det7->Write();

    PTBC_tof_amp_fIn_det2_afterCuts->Write();
    PTBC_tof_amp_fIn_det3_afterCuts->Write();
    PTBC_tof_amp_fIn_det4_afterCuts->Write();
    PTBC_tof_amp_fIn_det5_afterCuts->Write();
    PTBC_tof_amp_fIn_det6_afterCuts->Write();
    PTBC_tof_amp_fIn_det7_afterCuts->Write();

    // PTBC_tof_amp_fIn_det1_cutoff->Write();
    // PTBC_tof_amp_fIn_det2_cutoff->Write();
    // PTBC_tof_amp_fOut_det1_cutoff->Write();
    // PTBC_tof_amp_fOut_det2_cutoff->Write();

    PTBC_tof_amp_fIn_dedicated->Write();
    PTBC_tof_amp_fIn_parasitic->Write();
    PTBC_tof_amp_fOut_dedicated->Write();
    PTBC_tof_amp_fOut_parasitic->Write();

    PTBC_tof_amp_cut_det2->Write();
    PTBC_tof_amp_cut_det3->Write();
    PTBC_tof_amp_cut_det4->Write();
    PTBC_tof_amp_cut_det5->Write();
    PTBC_tof_amp_cut_det6->Write();
    PTBC_tof_amp_cut_det7->Write();
    PTBC_tof_amp_cut_para->Write();

    f->Close();
}

void cutoffAnalysis_PTBC(){

    gStyle->SetStatX(0.27);
    gStyle->SetStatY(0.9);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.17);

    fillCutsPTBC_det2();
    fillCutsPTBC_det3();
    fillCutsPTBC_det4();
    fillCutsPTBC_det5();
    fillCutsPTBC_det6();
    fillCutsPTBC_det7();
    fillCut_para();
    fillCut_det2();
    fillCut_det3();
    fillCut_det4();
    fillCut_det5();
    fillCut_det6();
    fillCut_det7();
    fillRuns();

    //Calculating TOF (x) bin edges
    int bins_per_decade = 1000;
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
    int num_bins_amp = 500;
    Double_t bin_edges_amp[num_bins_amp+1];
    Double_t step_amp = (Double_t) ((amp_max-amp_min)/num_bins_amp);
    // std::cout << "step_amp = " << step_amp << std::endl;
    for(int i = 0; i < num_bins_amp+1; i++)
    {
        bin_edges_amp[i] = step_amp * (Double_t) i;
    }

    //Filter Out
    PTBC_tof_amp_fOut_total = new TH2D("PTBC_tof_amp_fOut_total","ToF vs Amplitude Hist - PTBC All Det - Filter Out",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fOut_det2 = new TH2D("PTBC_tof_amp_fOut_det2","ToF vs Amplitude Hist - PTBC Det 2 - Filter Out",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fOut_det3 = new TH2D("PTBC_tof_amp_fOut_det3","ToF vs Amplitude Hist - PTBC Det 3 - Filter Out",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fOut_det4 = new TH2D("PTBC_tof_amp_fOut_det4","ToF vs Amplitude Hist - PTBC Det 4 - Filter Out",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fOut_det5 = new TH2D("PTBC_tof_amp_fOut_det5","ToF vs Amplitude Hist - PTBC Det 5 - Filter Out",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fOut_det6 = new TH2D("PTBC_tof_amp_fOut_det6","ToF vs Amplitude Hist - PTBC Det 6 - Filter Out",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fOut_det7 = new TH2D("PTBC_tof_amp_fOut_det7","ToF vs Amplitude Hist - PTBC Det 7 - Filter Out",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    
    //Filter In
    PTBC_tof_amp_fIn_total = new TH2D("PTBC_tof_amp_fIn_total",Form("ToF vs Amplitude Hist - PTBC All Det - %s", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fIn_det2 = new TH2D("PTBC_tof_amp_fIn_det2",Form("ToF vs Amplitude Hist - PTBC Det 2 - %s", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fIn_det3 = new TH2D("PTBC_tof_amp_fIn_det3",Form("ToF vs Amplitude Hist - PTBC Det 3 - %s", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fIn_det4 = new TH2D("PTBC_tof_amp_fIn_det4",Form("ToF vs Amplitude Hist - PTBC Det 4 - %s", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fIn_det5 = new TH2D("PTBC_tof_amp_fIn_det5",Form("ToF vs Amplitude Hist - PTBC Det 5 - %s", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fIn_det6 = new TH2D("PTBC_tof_amp_fIn_det6",Form("ToF vs Amplitude Hist - PTBC Det 6 - %s", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fIn_det7 = new TH2D("PTBC_tof_amp_fIn_det7",Form("ToF vs Amplitude Hist - PTBC Det 7 - %s", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);

    PTBC_tof_amp_fIn_det2_afterCuts = new TH2D("PTBC_tof_amp_fIn_det2_afterCuts",Form("ToF vs Amp - PTBC Det 2 - %s - After Cuts", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);;
    PTBC_tof_amp_fIn_det3_afterCuts = new TH2D("PTBC_tof_amp_fIn_det3_afterCuts",Form("ToF vs Amp - PTBC Det 3 - %s - After Cuts", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);;
    PTBC_tof_amp_fIn_det4_afterCuts = new TH2D("PTBC_tof_amp_fIn_det4_afterCuts",Form("ToF vs Amp - PTBC Det 4 - %s - After Cuts", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);;
    PTBC_tof_amp_fIn_det5_afterCuts = new TH2D("PTBC_tof_amp_fIn_det5_afterCuts",Form("ToF vs Amp - PTBC Det 5 - %s - After Cuts", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);;
    PTBC_tof_amp_fIn_det6_afterCuts = new TH2D("PTBC_tof_amp_fIn_det6_afterCuts",Form("ToF vs Amp - PTBC Det 6 - %s - After Cuts", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);;
    PTBC_tof_amp_fIn_det7_afterCuts = new TH2D("PTBC_tof_amp_fIn_det7_afterCuts",Form("ToF vs Amp - PTBC Det 7 - %s - After Cuts", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);;

    // PTBC_tof_amp_fOut_det1_cutoff = new TH2D("PTBC_tof_amp_fOut_det1_cutoff","ToF vs Amp Hist - PTBC Det 1 - Filter Out Cutoff",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    // PTBC_tof_amp_fOut_det2_cutoff = new TH2D("PTBC_tof_amp_fOut_det2_cutoff","ToF vs Amp Hist - PTBC Det 2 - Filter Out Cutoff",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);

    // PTBC_tof_amp_fIn_det1_cutoff = new TH2D("PTBC_tof_amp_fIn_det1_cutoff","ToF vs Amp Hist - PTBC Det 1 - Al (5cm) Filter Cutoff",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    // PTBC_tof_amp_fIn_det2_cutoff = new TH2D("PTBC_tof_amp_fIn_det2_cutoff","ToF vs Amp Hist - PTBC Det 2 - Al (5cm) Filter Cutoff",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);

    PTBC_tof_amp_fIn_dedicated = new TH2D("PTBC_tof_amp_fIn_dedicated",Form("ToF vs Amp Hist - %s - Dedicated", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fIn_parasitic = new TH2D("PTBC_tof_amp_fIn_parasitic",Form("ToF vs Amp Hist - %s - Parasitic", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fOut_dedicated = new TH2D("PTBC_tof_amp_fOut_dedicated","ToF vs Amp Hist - Filter Out - Dedicated",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_fOut_parasitic = new TH2D("PTBC_tof_amp_fOut_parasitic","ToF vs Amp Hist - Filter Out - Parasitic",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);

    
    FilterIn();
    FilterOut();

    StoreHist();
}