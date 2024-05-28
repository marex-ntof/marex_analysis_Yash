/**
 * @file cutoffAnalysis_PTBC.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2024-01-17
 */

#include "cutoffAnalysis_PTBC.h"

void Fill_tof_amp_hists(std::vector<Int_t> run_list){

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

            // PTBC_tof_amp_total->Fill(corrected_tof, (Double_t) amp_PTBC);

            if (PulseIntensity_PTBC > 6e12)
            {
                PTBC_tof_amp_total_dedicated->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (PulseIntensity_PTBC <= 6e12)
            {
                PTBC_tof_amp_total_parasitic->Fill(corrected_tof, (Double_t) amp_PTBC);
            }

            if (det_num_PTBC == 2) {
                PTBC_tof_amp_det2->Fill(corrected_tof, (Double_t) amp_PTBC);
                PTBC_ringing_det2->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (det_num_PTBC == 3) {
                PTBC_tof_amp_det3->Fill(corrected_tof, (Double_t) amp_PTBC);
                PTBC_ringing_det3->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (det_num_PTBC == 4) {
                PTBC_tof_amp_det4->Fill(corrected_tof, (Double_t) amp_PTBC);
                PTBC_ringing_det4->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (det_num_PTBC == 5) {
                PTBC_tof_amp_det5->Fill(corrected_tof, (Double_t) amp_PTBC);
                PTBC_ringing_det5->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (det_num_PTBC == 6) {
                PTBC_tof_amp_det6->Fill(corrected_tof, (Double_t) amp_PTBC);
                PTBC_ringing_det6->Fill(corrected_tof, (Double_t) amp_PTBC);
            } else if (det_num_PTBC == 7) {
                PTBC_tof_amp_det7->Fill(corrected_tof, (Double_t) amp_PTBC);
                PTBC_ringing_det7->Fill(corrected_tof, (Double_t) amp_PTBC);
            }
        }

        file_ntof->Close();
    }
}

void determine_cuts(){

    //Fill cut histograms
    TH2D* PTBC_tof_amp_det2_forCuts = (TH2D*)PTBC_tof_amp_det2->Rebin2D(50, 25, "PTBC_tof_amp_det2_forCuts");
    // TH2D* PTBC_tof_amp_det3_forCuts = (TH2D*)PTBC_tof_amp_det3->Rebin2D(50, 50, "PTBC_tof_amp_det3_forCuts");
    // TH2D* PTBC_tof_amp_det4_forCuts = (TH2D*)PTBC_tof_amp_det4->Rebin2D(50, 50, "PTBC_tof_amp_det4_forCuts");
    // TH2D* PTBC_tof_amp_det5_forCuts = (TH2D*)PTBC_tof_amp_det5->Rebin2D(50, 50, "PTBC_tof_amp_det5_forCuts");
    // TH2D* PTBC_tof_amp_det6_forCuts = (TH2D*)PTBC_tof_amp_det6->Rebin2D(50, 50, "PTBC_tof_amp_det6_forCuts");
    // TH2D* PTBC_tof_amp_det7_forCuts = (TH2D*)PTBC_tof_amp_det7->Rebin2D(50, 50, "PTBC_tof_amp_det7_forCuts");

    Int_t num_tof_bins = PTBC_tof_amp_det2_forCuts->GetNbinsX();

    for (Int_t i = 1; i <= num_tof_bins; i++) 
    {
        if (i < 19) //starting from xbin = 19 (around tof = 800ns)
        {
            PTBC_cuts_det2->SetBinContent(i, 50000.);
            continue;
        }
        
        std::string projection_name_det2 = "profile_det2_bin_" + std::to_string(i);
        TH1D* projection_det2 = (TH1D*)PTBC_tof_amp_det2_forCuts->ProjectionY(projection_name_det2.c_str(),i, i);

        if (i >= 19 && i <= 40) // For tof below 10^4 ns
        {
            std::vector<Int_t> peaks_bin_num;
            std::vector<Int_t> valleys_bin_num;
            std::vector<Int_t> peaks_bin_value;
            std::vector<Int_t> valleys_bin_value;

            for (Int_t j = 3; j < 41; j++) // scanning only till 10k ampluitude to get the min value
            {
                Int_t prevprevBinContent = projection_det2->GetBinContent(j-2);
                Int_t prevBinContent = projection_det2->GetBinContent(j-1);
                Int_t currentBinContent = projection_det2->GetBinContent(j);
                Int_t nextBinContent = projection_det2->GetBinContent(j+1);
                Int_t nextnextBinContent = projection_det2->GetBinContent(j+2);

                // Check for peak
                if (currentBinContent > prevprevBinContent && currentBinContent > prevBinContent && currentBinContent > nextBinContent && currentBinContent > nextnextBinContent) {
                    peaks_bin_num.push_back(j);
                    peaks_bin_value.push_back(currentBinContent);
                }
                
                // Check for valley
                if (currentBinContent < prevprevBinContent && currentBinContent < prevBinContent && currentBinContent < nextBinContent && currentBinContent < nextnextBinContent) {
                    valleys_bin_num.push_back(j);
                    valleys_bin_value.push_back(currentBinContent);
                }
            }

            // if (*max_element(peaks_bin_value.begin(), peaks_bin_value.end()) < ) 
            // {
            //     /* code */
            // }

            Int_t min_bin_num = 0; //Bin after which the cut needs to be placed
            for (Int_t j = 0; j < peaks_bin_value.size(); j++)
            {
                if (peaks_bin_value.at(j) > 500)
                {
                    min_bin_num = peaks_bin_num.at(j);
                }
            }

            Int_t amp_cut_bin = 0; //Bin where the cut is determined
            for (Int_t j = 0; j < valleys_bin_num.size(); j++)
            {
                if (valleys_bin_num.at(j) > min_bin_num && valleys_bin_value.at(j) < 100)
                {
                    amp_cut_bin = valleys_bin_num.at(j);
                    break;
                }
            }
            PTBC_cuts_det2->SetBinContent(i, projection_det2->GetBinCenter(amp_cut_bin));
            continue;
        }
        

        if (i > 40){ // For tof above 10^4 ns
            Int_t low_bin = -1;
            Int_t low_value = std::numeric_limits<int>::max();
            // Int_t num_bins_projection = projection_det2->GetNbinsX();
            for (Int_t j = 3; j < 21; j++) // scanning only till 5k ampluitude to get the min value
            {
                Int_t binContent = projection_det2->GetBinContent(j);

                // Check for low point
                if (binContent < low_value) {
                    low_value = binContent;
                    low_bin = j;
                }
            }
            PTBC_cuts_det2->SetBinContent(i, projection_det2->GetBinCenter(low_bin));
            continue;
        }
        
    }
}

void StoreHist(){
    
    TFile *f = new TFile(Form("../rootFiles/cutoffAnalysis_PTBC_%s.root", target_name.c_str()),"recreate");

    // PTBC_tof_amp_in->Write();
    // PTBC_tof_amp_out->Write();

    // PTBC_tof_amp_total->Write();
    PTBC_tof_amp_det2->Write();
    PTBC_tof_amp_det3->Write();
    PTBC_tof_amp_det4->Write();
    PTBC_tof_amp_det5->Write();
    PTBC_tof_amp_det6->Write();
    PTBC_tof_amp_det7->Write();

    PTBC_ringing_det2->Write();
    PTBC_ringing_det3->Write();
    PTBC_ringing_det4->Write();
    PTBC_ringing_det5->Write();
    PTBC_ringing_det6->Write();
    PTBC_ringing_det7->Write();

    //cuts
    // PTBC_cuts_det2->Write();
    // PTBC_cuts_det3->Write();
    // PTBC_cuts_det4->Write();
    // PTBC_cuts_det5->Write();
    // PTBC_cuts_det6->Write();
    // PTBC_cuts_det7->Write();

    // PTBC_tof_amp_det2_afterCuts->Write();
    // PTBC_tof_amp_det3_afterCuts->Write();
    // PTBC_tof_amp_det4_afterCuts->Write();
    // PTBC_tof_amp_det5_afterCuts->Write();
    // PTBC_tof_amp_det6_afterCuts->Write();
    // PTBC_tof_amp_det7_afterCuts->Write();

    PTBC_tof_amp_total_dedicated->Write();
    PTBC_tof_amp_total_parasitic->Write();

    f->Close();

    std::cout << "Created output file '" << Form("cutoffAnalysis_PTBC_%s.root", target_name.c_str()) << "'" << std::endl;
}

void cutoffAnalysis_PTBC(){

    gStyle->SetStatX(0.27);
    gStyle->SetStatY(0.9);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.17);

    // fillCutsPTBC();
    // fillCutGraph();
    fillRuns();

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

    //Calculating TOF (x) bin edges FOR CUTS
    // Int_t bins_per_decade_cuts = 20;
    // Int_t num_bins_tof_cuts = bins_per_decade_cuts * Num_decades;
    // Double_t bin_edges_tof_cuts[num_bins_tof_cuts+1];
    // Double_t step_tof_cuts = ((Double_t) 1.0/(Double_t) bins_per_decade_cuts);
    // for(Int_t i = 0; i < num_bins_tof_cuts+1; i++)
    // {
    //     Double_t base = 10.;
    //     Double_t exponent = (step_tof_cuts * (Double_t) i) + 2.;
    //     bin_edges_tof_cuts[i] = (Double_t) std::pow(base, exponent);
    // }

    // Calculating amplitude (y) bin edges
    Double_t bin_edges_amp[num_bins_amp+1];
    Double_t step_amp = (Double_t) ((amp_max-amp_min)/num_bins_amp);
    // std::cout << "step_amp = " << step_amp << std::endl;
    for(Int_t i = 0; i < num_bins_amp+1; i++)
    {
        bin_edges_amp[i] = step_amp * (Double_t) i;
    }

    //Initializing histograms
    // PTBC_tof_amp_total = new TH2D("PTBC_tof_amp_total",Form("ToF vs Amplitude Hist - PTBC All Det - %s", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_det2 = new TH2D("PTBC_tof_amp_det2",Form("ToF vs Amplitude Hist - PTBC Det 2 - %s", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_det3 = new TH2D("PTBC_tof_amp_det3",Form("ToF vs Amplitude Hist - PTBC Det 3 - %s", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_det4 = new TH2D("PTBC_tof_amp_det4",Form("ToF vs Amplitude Hist - PTBC Det 4 - %s", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_det5 = new TH2D("PTBC_tof_amp_det5",Form("ToF vs Amplitude Hist - PTBC Det 5 - %s", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_det6 = new TH2D("PTBC_tof_amp_det6",Form("ToF vs Amplitude Hist - PTBC Det 6 - %s", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_det7 = new TH2D("PTBC_tof_amp_det7",Form("ToF vs Amplitude Hist - PTBC Det 7 - %s", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);

    PTBC_ringing_det2 = new TH2D("PTBC_ringing_det2",Form("ToF vs Amp - Ringing - PTBC Det 2 - %s", target_name_title.c_str()), 920, 800., 10000., 2000, 0., 20000.);
    PTBC_ringing_det3 = new TH2D("PTBC_ringing_det3",Form("ToF vs Amp - Ringing - PTBC Det 3 - %s", target_name_title.c_str()), 920, 800., 10000., 2000, 0., 20000.);
    PTBC_ringing_det4 = new TH2D("PTBC_ringing_det4",Form("ToF vs Amp - Ringing - PTBC Det 4 - %s", target_name_title.c_str()), 920, 800., 10000., 2000, 0., 20000.);
    PTBC_ringing_det5 = new TH2D("PTBC_ringing_det5",Form("ToF vs Amp - Ringing - PTBC Det 5 - %s", target_name_title.c_str()), 920, 800., 10000., 2000, 0., 20000.);
    PTBC_ringing_det6 = new TH2D("PTBC_ringing_det6",Form("ToF vs Amp - Ringing - PTBC Det 6 - %s", target_name_title.c_str()), 920, 800., 10000., 2000, 0., 20000.);
    PTBC_ringing_det7 = new TH2D("PTBC_ringing_det7",Form("ToF vs Amp - Ringing - PTBC Det 7 - %s", target_name_title.c_str()), 920, 800., 10000., 2000, 0., 20000.);

    //Cut histograms
    // PTBC_cuts_det2 = new TH1D("PTBC_cuts_det2", Form("ToF-Amp cut Hist - PTBC Det 2 - %s", target_name_title.c_str()), num_bins_tof_cuts, bin_edges_tof_cuts);
    // PTBC_cuts_det3 = new TH1D("PTBC_cuts_det3", Form("ToF-Amp cut Hist - PTBC Det 3 - %s", target_name_title.c_str()), num_bins_tof_cuts, bin_edges_tof_cuts);
    // PTBC_cuts_det4 = new TH1D("PTBC_cuts_det4", Form("ToF-Amp cut Hist - PTBC Det 4 - %s", target_name_title.c_str()), num_bins_tof_cuts, bin_edges_tof_cuts);
    // PTBC_cuts_det5 = new TH1D("PTBC_cuts_det5", Form("ToF-Amp cut Hist - PTBC Det 5 - %s", target_name_title.c_str()), num_bins_tof_cuts, bin_edges_tof_cuts);
    // PTBC_cuts_det6 = new TH1D("PTBC_cuts_det6", Form("ToF-Amp cut Hist - PTBC Det 6 - %s", target_name_title.c_str()), num_bins_tof_cuts, bin_edges_tof_cuts);
    // PTBC_cuts_det7 = new TH1D("PTBC_cuts_det7", Form("ToF-Amp cut Hist - PTBC Det 7 - %s", target_name_title.c_str()), num_bins_tof_cuts, bin_edges_tof_cuts);
    
    
    // PTBC_tof_amp_det2_afterCuts = new TH2D("PTBC_tof_amp_det2_afterCuts",Form("ToF vs Amp - PTBC Det 2 - %s - After Cuts", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);;
    // PTBC_tof_amp_det3_afterCuts = new TH2D("PTBC_tof_amp_det3_afterCuts",Form("ToF vs Amp - PTBC Det 3 - %s - After Cuts", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);;
    // PTBC_tof_amp_det4_afterCuts = new TH2D("PTBC_tof_amp_det4_afterCuts",Form("ToF vs Amp - PTBC Det 4 - %s - After Cuts", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);;
    // PTBC_tof_amp_det5_afterCuts = new TH2D("PTBC_tof_amp_det5_afterCuts",Form("ToF vs Amp - PTBC Det 5 - %s - After Cuts", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);;
    // PTBC_tof_amp_det6_afterCuts = new TH2D("PTBC_tof_amp_det6_afterCuts",Form("ToF vs Amp - PTBC Det 6 - %s - After Cuts", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);;
    // PTBC_tof_amp_det7_afterCuts = new TH2D("PTBC_tof_amp_det7_afterCuts",Form("ToF vs Amp - PTBC Det 7 - %s - After Cuts", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);;

    // PTBC_tof_amp_fOut_det1_cutoff = new TH2D("PTBC_tof_amp_fOut_det1_cutoff","ToF vs Amp Hist - PTBC Det 1 - Filter Out Cutoff",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    // PTBC_tof_amp_fOut_det2_cutoff = new TH2D("PTBC_tof_amp_fOut_det2_cutoff","ToF vs Amp Hist - PTBC Det 2 - Filter Out Cutoff",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);

    // PTBC_tof_amp_det1_cutoff = new TH2D("PTBC_tof_amp_det1_cutoff","ToF vs Amp Hist - PTBC Det 1 - Al (5cm) Filter Cutoff",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    // PTBC_tof_amp_det2_cutoff = new TH2D("PTBC_tof_amp_det2_cutoff","ToF vs Amp Hist - PTBC Det 2 - Al (5cm) Filter Cutoff",num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);

    PTBC_tof_amp_total_dedicated = new TH2D("PTBC_tof_amp_total_dedicated",Form("ToF vs Amp Hist - %s - Total Dedicated", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    PTBC_tof_amp_total_parasitic = new TH2D("PTBC_tof_amp_total_parasitic",Form("ToF vs Amp Hist - %s - Total Parasitic", target_name_title.c_str()),num_bins_tof,bin_edges_tof,num_bins_amp,bin_edges_amp);
    
    Fill_tof_amp_hists(filter_run_list);
    // determine_cuts();
    // Fill_tof_amp_hists_det5(det5_runs);

    StoreHist();
}