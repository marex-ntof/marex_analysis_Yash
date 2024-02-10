/**
 * @file cutoffAnalysis_FIMG.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-11-28
 */

#include "cutoffAnalysis_FIMG.h"

void applynTOFCuts_FIMG(Double_t tof, Float_t amp, Float_t pulseIntensity, Int_t det_num){
    //Filling the histograms
    if (tof < 1e4)
    {
        return;
    }

    if (tof >= 1e4)
    {
        if (det_num == 1)
        {
            if ( (Double_t) amp >= fimgCutFunction(tof, det_num, "ntof", pulseIntensity) )
            {
                FIMG_tof_amp_det1_afterCuts->Fill(tof, (Double_t) amp);
            }
        }

        if (det_num == 2)
        {
            if ( (Double_t) amp >= fimgCutFunction(tof, det_num, "ntof", pulseIntensity) )
            {
                FIMG_tof_amp_det2_afterCuts->Fill(tof, (Double_t) amp);
            }
        }
    }
    return;
}

Double_t GetNormFactor(std::vector<Int_t> run_list){
    
    Double_t NormFactor = 0; //Integral of the pulse intensity

    for (int i = 0; i < run_list.size(); i++)
    {
        TFile *file_ntof = TFile::Open(Form("/eos/experiment/ntof/processing/official/done/run%d.root", run_list.at(i)),"read");
        cout << "Run Number = " << run_list.at(i) << endl;

        //FIMG ---------------------------------------------
        TTree* FIMG;
        Int_t BunchNumber_FIMG = 0;
        Float_t PulseIntensity = 0;

        file_ntof->GetObject("FIMG", FIMG);
        FIMG->SetBranchAddress("BunchNumber", &BunchNumber_FIMG);
        FIMG->SetBranchAddress("PulseIntensity", &PulseIntensity);

        int CurrentBunchNum = 0;

        Long64_t Events_FIMG = FIMG->GetEntriesFast();

        for (int j = 0; j < Events_FIMG; j++)
        {
            FIMG->GetEntry(j);

            if (CurrentBunchNum != BunchNumber_FIMG)
            {
                CurrentBunchNum = BunchNumber_FIMG;
                NormFactor += (Double_t) PulseIntensity;
            }
        }

        file_ntof->Close();
    }

    return NormFactor;
}

void fillEnergyHist(std::vector<Int_t> run_list){

    for (int i = 0; i < run_list.size(); i++)
    {
        TFile *file_ntof = TFile::Open(Form("/eos/experiment/ntof/processing/official/done/run%d.root", run_list.at(i)),"read");
        cout << "Run Number = " << run_list.at(i) << endl;

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
        Double_t tof_FIMG = 0; //tof is in ns
        Float_t amp = 0;
        Float_t area_0 = 0;
        Int_t det_num = 0;
        Int_t BunchNumber_FIMG = 0;
        Float_t PulseIntensity = 0;

        file_ntof->GetObject("FIMG", FIMG);
        FIMG->SetBranchAddress("BunchNumber", &BunchNumber_FIMG);
        FIMG->SetBranchAddress("tof", &tof_FIMG);
        FIMG->SetBranchAddress("amp", &amp);
        FIMG->SetBranchAddress("area_0", &area_0);
        FIMG->SetBranchAddress("detn", &det_num);
        FIMG->SetBranchAddress("PulseIntensity", &PulseIntensity);

        Long64_t Events_FIMG = FIMG->GetEntriesFast();
        std::cout << "Number of entries - FIMG = " << Events_FIMG << std::endl;

        for (int j = 0; j < Events_FIMG; j++)
        {
            FIMG->GetEntry(j);

            Double_t t_pkup = BNum_tpkup_map[BunchNumber_FIMG];
            Double_t corrected_tof = tof_FIMG - t_pkup + delT_pkup_fimg + t_gamma_FIMG;

            FIMG_tof_amp_total->Fill(corrected_tof, (Double_t) amp);

            if (PulseIntensity > 6e12)
            {
                if (det_num == 1)
                {
                    FIMG_tof_amp_dedi_det1->Fill(corrected_tof, (Double_t) amp);
                }

                if (det_num == 2)
                {
                    FIMG_tof_amp_dedi_det2->Fill(corrected_tof, (Double_t) amp);
                }
            } else if (PulseIntensity <= 6e12)
            {
                if (det_num == 1)
                {
                    FIMG_tof_amp_para_det1->Fill(corrected_tof, (Double_t) amp);
                }

                if (det_num == 2)
                {
                    FIMG_tof_amp_para_det2->Fill(corrected_tof, (Double_t) amp);
                }
            }

            applynTOFCuts_FIMG(corrected_tof, amp, PulseIntensity, det_num);
        }

        file_ntof->Close();
    }
}

void calcTransmission(TH1D* e_hist_in, TH1D* e_hist_out, TH1D* trans_hist, Double_t Qout_Qin){
    Int_t num_bins_e = e_hist_in->GetNbinsX();
    for (int i = 0; i < num_bins_e; i++)
    {
        Double_t bin_content_in = e_hist_in->GetBinContent(i+1);
        Double_t bin_content_out = e_hist_out->GetBinContent(i+1);
        if (bin_content_in == 0. || bin_content_out == 0.)
        {
            trans_hist->SetBinContent(i+1, 0.);
            trans_hist->SetBinError(i+1, 0.);
        } else {
            Double_t transmission = (bin_content_in * Qout_Qin)/bin_content_out;
            Double_t bin_unc = transmission * std::sqrt( (1./bin_content_in) + (1./bin_content_out) );
            trans_hist->SetBinContent(i+1, transmission);
            trans_hist->SetBinError(i+1, bin_unc);
        }
    }
}

// void fillRunHists(Int_t num_bins_e, Double_t bin_edges_e[]){

//     TH1D* e_loose_cut_det_1_in = 0;
//     TH1D* e_mid_cut_det_1_in = 0;
//     TH1D* e_tight_cut_det_1_in = 0;
//     TH1D* e_loose_cut_det_2_in = 0;
//     TH1D* e_mid_cut_det_2_in = 0;
//     TH1D* e_tight_cut_det_2_in = 0;

//     TH1D* e_loose_cut_det_1_out = 0;
//     TH1D* e_mid_cut_det_1_out = 0;
//     TH1D* e_tight_cut_det_1_out = 0;
//     TH1D* e_loose_cut_det_2_out = 0;
//     TH1D* e_mid_cut_det_2_out = 0;
//     TH1D* e_tight_cut_det_2_out = 0;

//     //Target In Hists
//     e_loose_cut_det_1_in = new TH1D("e_loose_cut_det_1_in","Energy Hist",num_bins_e,bin_edges_e);
//     e_mid_cut_det_1_in = new TH1D("e_mid_cut_det_1_in","Energy Hist",num_bins_e,bin_edges_e);
//     e_tight_cut_det_1_in = new TH1D("e_tight_cut_det_1_in","Energy Hist",num_bins_e,bin_edges_e);
//     e_loose_cut_det_2_in = new TH1D("e_loose_cut_det_2_in","Energy Hist",num_bins_e,bin_edges_e);
//     e_mid_cut_det_2_in = new TH1D("e_mid_cut_det_2_in","Energy Hist",num_bins_e,bin_edges_e);
//     e_tight_cut_det_2_in = new TH1D("e_tight_cut_det_2_in","Energy Hist",num_bins_e,bin_edges_e);
//     //Target Out Hists
//     e_loose_cut_det_1_out = new TH1D("e_loose_cut_det_1_out","Energy Hist",num_bins_e,bin_edges_e);
//     e_mid_cut_det_1_out = new TH1D("e_mid_cut_det_1_out","Energy Hist",num_bins_e,bin_edges_e);
//     e_tight_cut_det_1_out = new TH1D("e_tight_cut_det_1_out","Energy Hist",num_bins_e,bin_edges_e);
//     e_loose_cut_det_2_out = new TH1D("e_loose_cut_det_2_out","Energy Hist",num_bins_e,bin_edges_e);
//     e_mid_cut_det_2_out = new TH1D("e_mid_cut_det_2_out","Energy Hist",num_bins_e,bin_edges_e);
//     e_tight_cut_det_2_out = new TH1D("e_tight_cut_det_2_out","Energy Hist",num_bins_e,bin_edges_e);
    
//     //Norm Factors
//     Double_t norm_factor_target_in = GetNormFactor(ts_target_in_runs);
//     Double_t norm_factor_target_out = 0;
//     if(!target_out_name.compare("ar_bottle")) {
//         norm_factor_target_out = GetNormFactor(ts_target_out_runs);
//     } else if(!target_out_name.compare("none_ts")) {
//         norm_factor_target_out = GetNormFactor(f_out_t_out_ts);
//     } else if(!target_out_name.compare("none")) {
//         norm_factor_target_out = GetNormFactor(f_out_t_out_runs);
//     }

//     Double_t Qout_Qin = norm_factor_target_out / norm_factor_target_in;
    
//     //Filling the Histograms - Target In
//     fillEnergyHist(ts_target_in_runs, "loose", e_loose_cut_det_1_in, e_loose_cut_det_2_in);
//     fillEnergyHist(ts_target_in_runs, "mid", e_mid_cut_det_1_in, e_mid_cut_det_2_in);
//     fillEnergyHist(ts_target_in_runs, "tight", e_tight_cut_det_1_in, e_tight_cut_det_2_in);

//     //Filling the Histograms - Target Out
//     if(!target_out_name.compare("ar_bottle")) {
//         fillEnergyHist(ts_target_out_runs, "loose", e_loose_cut_det_1_out, e_loose_cut_det_2_out);
//         fillEnergyHist(ts_target_out_runs, "mid", e_mid_cut_det_1_out, e_mid_cut_det_2_out);
//         fillEnergyHist(ts_target_out_runs, "tight", e_tight_cut_det_1_out, e_tight_cut_det_2_out);
//     } else if(!target_out_name.compare("none_ts")) {
//         fillEnergyHist(f_out_t_out_ts, "loose", e_loose_cut_det_1_out, e_loose_cut_det_2_out);
//         fillEnergyHist(f_out_t_out_ts, "mid", e_mid_cut_det_1_out, e_mid_cut_det_2_out);
//         fillEnergyHist(f_out_t_out_ts, "tight", e_tight_cut_det_1_out, e_tight_cut_det_2_out);
//     } else if(!target_out_name.compare("none")) {
//         fillEnergyHist(f_out_t_out_runs, "loose", e_loose_cut_det_1_out, e_loose_cut_det_2_out);
//         fillEnergyHist(f_out_t_out_runs, "mid", e_mid_cut_det_1_out, e_mid_cut_det_2_out);
//         fillEnergyHist(f_out_t_out_runs, "tight", e_tight_cut_det_1_out, e_tight_cut_det_2_out);
//     }

//     //Transmission
//     for (int i = 0; i < num_bins_e; i++)
//     {
//         //Loose
//         calcTransmission(e_loose_cut_det_1_in, e_loose_cut_det_1_out, trans_loose_cut_det_1, Qout_Qin);
//         calcTransmission(e_loose_cut_det_2_in, e_loose_cut_det_2_out, trans_loose_cut_det_2, Qout_Qin);

//         //Mid
//         calcTransmission(e_mid_cut_det_1_in, e_mid_cut_det_1_out, trans_mid_cut_det_1, Qout_Qin);
//         calcTransmission(e_mid_cut_det_2_in, e_mid_cut_det_2_out, trans_mid_cut_det_2, Qout_Qin);

//         //Tight
//         calcTransmission(e_tight_cut_det_1_in, e_tight_cut_det_1_out, trans_tight_cut_det_1, Qout_Qin);
//         calcTransmission(e_tight_cut_det_2_in, e_tight_cut_det_2_out, trans_tight_cut_det_2, Qout_Qin);
//     }

//     cout << Form("Total Protons %s = ", target_name.c_str()) << norm_factor_target_in << endl;
//     cout << "Total Protons Target Out = " << norm_factor_target_out << endl;
// }

void cutoffAnalysis_FIMG() {
    
    fillRuns();

    //Calculating Energy bin edges
    // Double_t tof_min = 1e3; //ns
    // Double_t tof_max = 1e8; //ns
    // Double_t e_min = TOFToEnergy(tof_max * 1e-9, flight_path_length_FIMG); //converting into seconds
    // Double_t e_max = TOFToEnergy(tof_min * 1e-9, flight_path_length_FIMG); //converting into seconds
    // int min_power = FindDecadePower(e_min);
    // int max_power = FindDecadePower(e_max);
    // int num_decades_e = max_power - min_power;
    // int num_bins_e = bins_per_decade * num_decades_e;
    // Double_t bin_edges_e[num_bins_e+1];
    // Double_t step_e = ((Double_t) 1.0/(Double_t) bins_per_decade);
    // for(int i = 0; i < num_bins_e+1; i++)
    // {
    //     Double_t base = 10.;
    //     Double_t exponent = (step_e * (Double_t) i) + (Double_t) min_power;
    //     bin_edges_e[i] = (Double_t) std::pow(base, exponent);
    // }

    //Calculating TOF (x) bin edges
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
    Double_t bin_edges_amp[num_bins_amp+1];
    Double_t step_amp = (Double_t) ((amp_max-amp_min)/num_bins_amp);
    // std::cout << "step_amp = " << step_amp << std::endl;
    for(int i = 0; i < num_bins_amp+1; i++)
    {
        bin_edges_amp[i] = step_amp * (Double_t) i;
    }

    // cout << "Number of tof bins = " << num_bins_tof << endl;
    // cout << "Number of e bins = " << num_bins_e << endl;

    //Defining histograms
    FIMG_tof_amp_total = new TH2D("FIMG_tof_amp_total",Form("ToF vs Amp - FIMG All Det - %s", target_name_title.c_str()), num_bins_tof, bin_edges_tof, num_bins_amp, bin_edges_amp);
    FIMG_tof_amp_dedi_det1 = new TH2D("FIMG_tof_amp_dedi_det1",Form("ToF vs Amp - FIMG Det 1 - Dedicated - %s", target_name_title.c_str()), num_bins_tof, bin_edges_tof, num_bins_amp, bin_edges_amp);
    FIMG_tof_amp_dedi_det2 = new TH2D("FIMG_tof_amp_dedi_det2",Form("ToF vs Amp - FIMG Det 2 - Dedicated - %s", target_name_title.c_str()), num_bins_tof, bin_edges_tof, num_bins_amp, bin_edges_amp);
    FIMG_tof_amp_para_det1 = new TH2D("FIMG_tof_amp_para_det1",Form("ToF vs Amp - FIMG Det 1 - Parasitic - %s", target_name_title.c_str()), num_bins_tof, bin_edges_tof, num_bins_amp, bin_edges_amp);
    FIMG_tof_amp_para_det2 = new TH2D("FIMG_tof_amp_para_det2",Form("ToF vs Amp - FIMG Det 2 - Parasitic - %s", target_name_title.c_str()), num_bins_tof, bin_edges_tof, num_bins_amp, bin_edges_amp);

    FIMG_tof_amp_det1_afterCuts = new TH2D("FIMG_tof_amp_det1_afterCuts",Form("ToF vs Amp - FIMG Det 1 - After Cuts - %s", target_name_title.c_str()), num_bins_tof, bin_edges_tof, num_bins_amp, bin_edges_amp);
    FIMG_tof_amp_det2_afterCuts = new TH2D("FIMG_tof_amp_det2_afterCuts",Form("ToF vs Amp - FIMG Det 2 - After Cuts - %s", target_name_title.c_str()), num_bins_tof, bin_edges_tof, num_bins_amp, bin_edges_amp);

    //transmission histogram
    // trans_loose_cut_det_1 = new TH1D("trans_loose_cut_det_1","Transmission Hist",num_bins_e,bin_edges_e);
    // trans_mid_cut_det_1 = new TH1D("trans_mid_cut_det_1","Transmission Hist",num_bins_e,bin_edges_e);
    // trans_tight_cut_det_1 = new TH1D("trans_tight_cut_det_1","Transmission Hist",num_bins_e,bin_edges_e);
    // trans_loose_cut_det_2 = new TH1D("trans_loose_cut_det_2","Transmission Hist",num_bins_e,bin_edges_e);
    // trans_mid_cut_det_2 = new TH1D("trans_mid_cut_det_2","Transmission Hist",num_bins_e,bin_edges_e);
    // trans_tight_cut_det_2 = new TH1D("trans_tight_cut_det_2","Transmission Hist",num_bins_e,bin_edges_e);
    
    // fillRunHists(num_bins_e,bin_edges_e);
    if(!target_name.compare("no_filter")){
        fillEnergyHist(f_out_t_out_runs);
    } else {
        fillEnergyHist(ts_target_in_runs);
    }

    fill_nTOF_cuts();

    //Writing to the output file
    outputRootFile = new TFile(Form("../rootFiles/cutoffAnalysis_FIMG_%s_nTOF_cuts.root", target_name.c_str()),"recreate");
    FIMG_tof_amp_total->Write();
    FIMG_tof_amp_dedi_det1->Write();
    FIMG_tof_amp_dedi_det2->Write();
    FIMG_tof_amp_para_det1->Write();
    FIMG_tof_amp_para_det2->Write();
    FIMG_tof_amp_det1_afterCuts->Write();
    FIMG_tof_amp_det2_afterCuts->Write();

    FIMG_tof_amp_cut_dedi_det1->Write();
    FIMG_tof_amp_cut_dedi_det2->Write();
    FIMG_tof_amp_cut_para_det1->Write();
    FIMG_tof_amp_cut_para_det2->Write();

    outputRootFile->Close();

    std::cout << Form("Created output file 'cutoffAnalysis_FIMG_%s_nTOF_cuts.root'", target_name.c_str()) << std::endl;
    
}