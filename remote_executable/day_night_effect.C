/**
 * @file day_night_effect.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2024-04-11
 */

#include "day_night_effect.h"

void fill_dn_hists(std::vector<Int_t> run_list, TH1D* counts_hist_PTBC, TH1D* counts_hist_FIMG, TH1D* PI_hist_PTBC, TH1D* PI_hist_FIMG, TH1D* dn_hist_PTBC, TH1D* dn_hist_FIMG){

    Int_t num_bins_PTBC = dn_hist_PTBC->GetNbinsX();
    Int_t num_bins_FIMG = dn_hist_FIMG->GetNbinsX();

    Float_t PTBC_PulseIntensity[96] = {}; // 96 for 15 min bins
    Float_t FIMG_PulseIntensity[96] = {}; // 96 for 15 min bins

    for (int i = 0; i < run_list.size(); i++)
    {
        TFile *file_ntof = TFile::Open(Form("/eos/experiment/ntof/data/rootfiles/2023/ear1/run%d.root", run_list.at(i)),"read");
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

        //PTBC ---------------------------------------------
        TTree* PTBC;
        Double_t tof_PTBC = 0; //tof is in ns
        Float_t amp_PTBC = 0;
        Int_t time_PTBC = 0;
        Int_t BunchNumber_PTBC = 0;
        Float_t PulseIntensity_PTBC = 0;
        Int_t det_num_PTBC = 0;

        file_ntof->GetObject("PTBC", PTBC);
        PTBC->SetBranchAddress("time", &time_PTBC);
        PTBC->SetBranchAddress("tof", &tof_PTBC);
        PTBC->SetBranchAddress("amp", &amp_PTBC);
        PTBC->SetBranchAddress("BunchNumber", &BunchNumber_PTBC);
        PTBC->SetBranchAddress("PulseIntensity", &PulseIntensity_PTBC);
        PTBC->SetBranchAddress("detn", &det_num_PTBC);

        Long64_t Events_PTBC = PTBC->GetEntriesFast();
        std::cout << "Number of entries - PTBC = " << Events_PTBC << std::endl;
        
        Int_t CurrentBunchNum = 0;
        Int_t hit_time = 0;
        Int_t entry_bin_num = 0;
        
        for (Int_t j = 0; j < Events_PTBC; j++)
        {
            PTBC->GetEntry(j);

            if (det_num_PTBC == 1 || det_num_PTBC == 8) {
                continue;
            }

            Double_t t_pkup = BNum_tpkup_map[BunchNumber_PTBC];
            Double_t corrected_tof = tof_PTBC - t_pkup + delT_pkup_ptbc + t_gamma_PTBC;

            if (!select_hit_PTBC(corrected_tof, amp_PTBC, det_num_PTBC))
            {
                continue;
            }
            
            hit_time = timeToSeconds(std::to_string(time_PTBC));

            if (j == 0)
            {
                CurrentBunchNum = BunchNumber_PTBC;
                dn_hist_PTBC->Fill( hit_time );
                counts_hist_PTBC->Fill( hit_time );
                entry_bin_num = dn_hist_PTBC->FindBin( hit_time );
                PTBC_PulseIntensity[entry_bin_num-1] += PulseIntensity_PTBC;
                continue;
            }   

            dn_hist_PTBC->Fill( hit_time );
            counts_hist_PTBC->Fill( hit_time );

            if (CurrentBunchNum != BunchNumber_PTBC)
            {
                entry_bin_num = dn_hist_PTBC->FindBin( hit_time );
                CurrentBunchNum = BunchNumber_PTBC;
                PTBC_PulseIntensity[entry_bin_num-1] += PulseIntensity_PTBC;
            }
        }

        //FIMG ---------------------------------------------
        TTree* FIMG;
        Double_t tof_FIMG = 0; //tof is in ns
        Float_t amp_FIMG = 0;
        Int_t det_num_FIMG = 0;
        Int_t time_FIMG = 0;
        Int_t BunchNumber_FIMG = 0;
        Float_t PulseIntensity_FIMG = 0;

        file_ntof->GetObject("FIMG", FIMG);
        FIMG->SetBranchAddress("time", &time_FIMG);
        FIMG->SetBranchAddress("tof", &tof_FIMG);
        FIMG->SetBranchAddress("amp", &amp_FIMG);
        FIMG->SetBranchAddress("BunchNumber", &BunchNumber_FIMG);
        FIMG->SetBranchAddress("detn", &det_num_FIMG);
        FIMG->SetBranchAddress("PulseIntensity", &PulseIntensity_FIMG);

        Long64_t Events_FIMG = FIMG->GetEntriesFast();
        std::cout << "Number of entries - FIMG = " << Events_FIMG << std::endl;

        //reset variables
        CurrentBunchNum = 0;
        hit_time = 0;
        entry_bin_num = 0;

        for (int j = 0; j < Events_FIMG; j++)
        {
            FIMG->GetEntry(j);

            Double_t t_pkup = BNum_tpkup_map[BunchNumber_FIMG];
            Double_t corrected_tof = tof_FIMG - t_pkup + delT_pkup_fimg + t_gamma_FIMG;

            if (!select_hit_FIMG(corrected_tof, amp_FIMG, det_num_FIMG))
            {
                continue;
            }

            hit_time = timeToSeconds(std::to_string(time_FIMG));

            if (j == 0)
            {
                CurrentBunchNum = BunchNumber_FIMG;
                dn_hist_FIMG->Fill( hit_time );
                counts_hist_FIMG->Fill( hit_time );
                entry_bin_num = dn_hist_FIMG->FindBin( hit_time );
                FIMG_PulseIntensity[entry_bin_num-1] += PulseIntensity_FIMG;
                continue;
            }   

            dn_hist_FIMG->Fill( hit_time );
            counts_hist_FIMG->Fill( hit_time );
            // std::cout << "FIMG - Filled a bin! (bin_num_FIMG = " << bin_num_FIMG << ", bin_content = " << bin_content << ")" << std::endl;

            if (CurrentBunchNum != BunchNumber_FIMG)
            {
                entry_bin_num = dn_hist_FIMG->FindBin( hit_time );
                CurrentBunchNum = BunchNumber_FIMG;
                FIMG_PulseIntensity[entry_bin_num-1] += PulseIntensity_FIMG;
            }
        }

        file_ntof->Close();
    }

    for (Int_t i = 0; i < num_bins_PTBC; i++)
    {
        Int_t bin_content = dn_hist_PTBC->GetBinContent(i+1);
        Double_t bin_error = sqrt(bin_content);
        if (PTBC_PulseIntensity[i] != 0){
            PI_hist_PTBC->SetBinContent(i+1, (Double_t) PTBC_PulseIntensity[i] );
            counts_hist_PTBC->SetBinError(i+1, bin_error);
            dn_hist_PTBC->SetBinContent(i+1, (Double_t) bin_content/(Double_t) PTBC_PulseIntensity[i] );
            dn_hist_PTBC->SetBinError(i+1, bin_error/(Double_t) PTBC_PulseIntensity[i] );
        }
    }

    for (Int_t i = 0; i < num_bins_FIMG; i++)
    {
        Int_t bin_content = dn_hist_FIMG->GetBinContent(i+1);
        Double_t bin_error = sqrt(bin_content);
        if (FIMG_PulseIntensity[i] != 0){
            PI_hist_FIMG->SetBinContent(i+1, (Double_t) FIMG_PulseIntensity[i] );
            counts_hist_FIMG->SetBinError(i+1, bin_error);
            dn_hist_FIMG->SetBinContent(i+1, (Double_t) bin_content/(Double_t) FIMG_PulseIntensity[i] );
            dn_hist_FIMG->SetBinError(i+1, bin_error/(Double_t) FIMG_PulseIntensity[i] );
        }
    }
    
}

void fill_run_list(){

    // Filter out runs
    filterOut_run_list.insert(filterOut_run_list.end(), c_out_f_out_t_out.begin(), c_out_f_out_t_out.end());    
    filterOut_run_list.insert(filterOut_run_list.end(), c_pb_f_out_t_out.begin(), c_pb_f_out_t_out.end());
    filterOut_run_list.insert(filterOut_run_list.end(), c_ta_f_out_t_out.begin(), c_ta_f_out_t_out.end());
    filterOut_run_list.insert(filterOut_run_list.end(), c_c_f_out_t_out.begin(), c_c_f_out_t_out.end());

    // Bi (1 cm)
    bi1_run_list.insert(bi1_run_list.end(), c_au_f_bi_t_out.begin(), c_au_f_bi_t_out.end());    
    bi1_run_list.insert(bi1_run_list.end(), c_out_f_bi_t_out.begin(), c_out_f_bi_t_out.end());
    bi1_run_list.insert(bi1_run_list.end(), c_pb_f_bi_t_out.begin(), c_pb_f_bi_t_out.end());
    bi1_run_list.insert(bi1_run_list.end(), c_ta_f_bi_t_out.begin(), c_ta_f_bi_t_out.end());

    // Al (5 cm)
    al5_run_list.insert(al5_run_list.end(), c_ta_f_al5_t_out.begin(), c_ta_f_al5_t_out.end());    
    al5_run_list.insert(al5_run_list.end(), c_out_f_al5_t_out.begin(), c_out_f_al5_t_out.end());
    al5_run_list.insert(al5_run_list.end(), c_pb_f_al5_t_out.begin(), c_pb_f_al5_t_out.end());
    al5_run_list.insert(al5_run_list.end(), c_c_f_al5_t_out.begin(), c_c_f_al5_t_out.end());
    al5_run_list.insert(al5_run_list.end(), c_au_f_al5_t_out.begin(), c_au_f_al5_t_out.end());    
}

void day_night_effect(){

    fillCutsFIMG();
    fillNumDensityMap();
    fill_run_list();

    //Get PTBC Cuts
    for (Int_t i = 0; i < 6; i++)
    {
        PTBC_tof_amp_cuts[i] = GetHist1D("../inputFiles/PTBC_cuts.root", Form("PTBC_cuts_det%i", i+2));
    }

    // Calculating day-night effect bins
    Int_t num_bins_dn = 96; // 96 bins for 15 min bins // 24 hours
    Double_t x_min_dn = 0.; //in seconds
    Double_t x_max_dn = 86400.; //in seconds
    Double_t bin_edges_dn[num_bins_dn+1];
    Double_t bin_step_dn = (x_max_dn - x_min_dn)/ (Double_t) num_bins_dn;
    for(Int_t i = 0; i < num_bins_dn+1; i++)
    {
        bin_edges_dn[i] = x_min_dn + (Double_t) i * bin_step_dn;
    }

    // Initializing counts plots
    counts_filterOut_PTBC = new TH1D("counts_filterOut_PTBC", "Counts vs Time - No Filter - PTBC", num_bins_dn, bin_edges_dn);
    counts_Bi_PTBC = new TH1D("counts_Bi_PTBC", "Counts vs Time - Bi (1 cm) - PTBC", num_bins_dn, bin_edges_dn);
    counts_Al5_PTBC = new TH1D("counts_Al5_PTBC", "Counts vs Time - Al (5 cm) - PTBC", num_bins_dn, bin_edges_dn);
    counts_emptyTS_PTBC = new TH1D("counts_emptyTS_PTBC", "Counts vs Time - Empty (TS) - PTBC", num_bins_dn, bin_edges_dn);
    counts_emptyTank_PTBC = new TH1D("counts_emptyTank_PTBC", "Counts vs Time - Empty Tank - PTBC", num_bins_dn, bin_edges_dn);
    counts_Argon_PTBC = new TH1D("counts_Argon_PTBC", "Counts vs Time - Argon Tank - PTBC", num_bins_dn, bin_edges_dn);
    counts_EmptyArgon_PTBC = new TH1D("counts_EmptyArgon_PTBC", "Counts vs Time - Empty Ar Tank - PTBC", num_bins_dn, bin_edges_dn);

    counts_filterOut_FIMG = new TH1D("counts_filterOut_FIMG", "Counts vs Time - No Filter - FIMG", num_bins_dn, bin_edges_dn);
    counts_Bi_FIMG = new TH1D("counts_Bi_FIMG", "Counts vs Time - Bi (1 cm) - FIMG", num_bins_dn, bin_edges_dn);
    counts_Al5_FIMG = new TH1D("counts_Al5_FIMG", "Counts vs Time - Al (5 cm) - FIMG", num_bins_dn, bin_edges_dn);
    counts_emptyTS_FIMG = new TH1D("counts_emptyTS_FIMG", "Counts vs Time - Empty (TS) - FIMG", num_bins_dn, bin_edges_dn);
    counts_emptyTank_FIMG = new TH1D("counts_emptyTank_FIMG", "Counts vs Time - Empty Tank - FIMG", num_bins_dn, bin_edges_dn);
    counts_Argon_FIMG = new TH1D("counts_Argon_FIMG", "Counts vs Time - Argon Tank - FIMG", num_bins_dn, bin_edges_dn);
    counts_EmptyArgon_FIMG = new TH1D("counts_EmptyArgon_FIMG", "Counts vs Time - Empty Ar Tank - FIMG", num_bins_dn, bin_edges_dn);

    // Initializing Pulse Intensity plots
    pulseIntensity_filterOut_PTBC = new TH1D("pulseIntensity_filterOut_PTBC", "Pulse Intensity vs Time - No Filter - PTBC", num_bins_dn, bin_edges_dn);
    pulseIntensity_Bi_PTBC = new TH1D("pulseIntensity_Bi_PTBC", "Pulse Intensity vs Time - Bi (1 cm) - PTBC", num_bins_dn, bin_edges_dn);
    pulseIntensity_Al5_PTBC = new TH1D("pulseIntensity_Al5_PTBC", "Pulse Intensity vs Time - Al (5 cm) - PTBC", num_bins_dn, bin_edges_dn);
    pulseIntensity_emptyTS_PTBC = new TH1D("pulseIntensity_emptyTS_PTBC", "Pulse Intensity vs Time - Empty (TS) - PTBC", num_bins_dn, bin_edges_dn);
    pulseIntensity_emptyTank_PTBC = new TH1D("pulseIntensity_emptyTank_PTBC", "Pulse Intensity vs Time - Empty Tank - PTBC", num_bins_dn, bin_edges_dn);
    pulseIntensity_Argon_PTBC = new TH1D("pulseIntensity_Argon_PTBC", "Pulse Intensity vs Time - Argon Tank - PTBC", num_bins_dn, bin_edges_dn);
    pulseIntensity_EmptyArgon_PTBC = new TH1D("pulseIntensity_EmptyArgon_PTBC", "Pulse Intensity vs Time - Empty Ar Tank - PTBC", num_bins_dn, bin_edges_dn);

    pulseIntensity_filterOut_FIMG = new TH1D("pulseIntensity_filterOut_FIMG", "Pulse Intensity vs Time - No Filter - FIMG", num_bins_dn, bin_edges_dn);
    pulseIntensity_Bi_FIMG = new TH1D("pulseIntensity_Bi_FIMG", "Pulse Intensity vs Time - Bi (1 cm) - FIMG", num_bins_dn, bin_edges_dn);
    pulseIntensity_Al5_FIMG = new TH1D("pulseIntensity_Al5_FIMG", "Pulse Intensity vs Time - Al (5 cm) - FIMG", num_bins_dn, bin_edges_dn);
    pulseIntensity_emptyTS_FIMG = new TH1D("pulseIntensity_emptyTS_FIMG", "Pulse Intensity vs Time - Empty (TS) - FIMG", num_bins_dn, bin_edges_dn);
    pulseIntensity_emptyTank_FIMG = new TH1D("pulseIntensity_emptyTank_FIMG", "Pulse Intensity vs Time - Empty Tank - FIMG", num_bins_dn, bin_edges_dn);
    pulseIntensity_Argon_FIMG = new TH1D("pulseIntensity_Argon_FIMG", "Pulse Intensity vs Time - Argon Tank - FIMG", num_bins_dn, bin_edges_dn);
    pulseIntensity_EmptyArgon_FIMG = new TH1D("pulseIntensity_EmptyArgon_FIMG", "Pulse Intensity vs Time - Empty Ar Tank - FIMG", num_bins_dn, bin_edges_dn);

    // Initializing day night plots
    // PTBC
    day_night_filterOut_PTBC = new TH1D("day_night_filterOut_PTBC", "Day-Night Effect - No Filter - PTBC", num_bins_dn, bin_edges_dn);
    day_night_Bi_PTBC = new TH1D("day_night_Bi_PTBC", "Day-Night Effect - Bi (1 cm) - PTBC", num_bins_dn, bin_edges_dn);
    day_night_Al5_PTBC = new TH1D("day_night_Al5_PTBC", "Day-Night Effect - Al (5 cm) - PTBC", num_bins_dn, bin_edges_dn);
    day_night_emptyTS_PTBC = new TH1D("day_night_emptyTS_PTBC", "Day-Night Effect - Empty (TS) - PTBC", num_bins_dn, bin_edges_dn);
    day_night_emptyTank_PTBC = new TH1D("day_night_emptyTank_PTBC", "Day-Night Effect - Empty Tank - PTBC", num_bins_dn, bin_edges_dn);
    day_night_Argon_PTBC = new TH1D("day_night_Argon_PTBC", "Day-Night Effect - Argon Tank - PTBC", num_bins_dn, bin_edges_dn);
    day_night_EmptyArgon_PTBC = new TH1D("day_night_EmptyArgon_PTBC", "Day-Night Effect - Empty Ar Tank - PTBC", num_bins_dn, bin_edges_dn);

    // FIMG
    day_night_filterOut_FIMG = new TH1D("day_night_filterOut_FIMG", "Day-Night Effect - No Filter - FIMG", num_bins_dn, bin_edges_dn);
    day_night_Bi_FIMG = new TH1D("day_night_Bi_FIMG", "Day-Night Effect - Bi (1 cm) - FIMG", num_bins_dn, bin_edges_dn);
    day_night_Al5_FIMG = new TH1D("day_night_Al5_FIMG", "Day-Night Effect - Al (5 cm) - FIMG", num_bins_dn, bin_edges_dn);
    day_night_emptyTS_FIMG = new TH1D("day_night_emptyTS_FIMG", "Day-Night Effect - Empty (TS) - FIMG", num_bins_dn, bin_edges_dn);
    day_night_emptyTank_FIMG = new TH1D("day_night_emptyTank_FIMG", "Day-Night Effect - Empty Tank - FIMG", num_bins_dn, bin_edges_dn);
    day_night_Argon_FIMG = new TH1D("day_night_Argon_FIMG", "Day-Night Effect - Argon Tank - FIMG", num_bins_dn, bin_edges_dn);
    day_night_EmptyArgon_FIMG = new TH1D("day_night_EmptyArgon_FIMG", "Day-Night Effect - Empty Ar Tank - FIMG", num_bins_dn, bin_edges_dn);

    // Run Lists (FIMG starts from 117386)
    // Day-Night effect
    // first and last runs also include previous and the next day runs. This has to be removed
    // std::vector<Int_t> run_list_Bi_sep18 = {117389, 117390, 117391, 117392, 117393, 117394};
    // std::vector<Int_t> run_list_Al5_sep27 = {117472, 117473, 117474, 117475, 117476, 117477, 117479, 117480, 117481};
    // std::vector<Int_t> run_list_emptyTS_oct05 = {117543, 117544, 117545, 117546, 117547, 117548, 117549, 117550};
    // std::vector<Int_t> run_list_emptyTank_oct15 = {117632, 117633, 117634, 117635, 117636, 117637, 117638, 117639};
    // std::vector<Int_t> run_list_Argon_oct22 = {117729, 117730, 117731, 117732, 117733, 117734, 117735, 117736};

    // Filling the histograms
    // Day Night hists
    fill_dn_hists(filterOut_run_list, counts_filterOut_PTBC, counts_filterOut_FIMG, pulseIntensity_filterOut_PTBC, pulseIntensity_filterOut_FIMG, day_night_filterOut_PTBC, day_night_filterOut_FIMG);
    fill_dn_hists(bi1_run_list, counts_Bi_PTBC, counts_Bi_FIMG, pulseIntensity_Bi_PTBC, pulseIntensity_Bi_FIMG, day_night_Bi_PTBC, day_night_Bi_FIMG);
    fill_dn_hists(al5_run_list, counts_Al5_PTBC, counts_Al5_FIMG, pulseIntensity_Al5_PTBC, pulseIntensity_Al5_FIMG, day_night_Al5_PTBC, day_night_Al5_FIMG);
    fill_dn_hists(f_out_t_out_ts, counts_emptyTS_PTBC, counts_emptyTS_FIMG, pulseIntensity_emptyTS_PTBC, pulseIntensity_emptyTS_FIMG, day_night_emptyTS_PTBC, day_night_emptyTS_FIMG);
    fill_dn_hists(f_out_t_emptyBotRot_ts, counts_emptyTank_PTBC, counts_emptyTank_FIMG, pulseIntensity_emptyTank_PTBC, pulseIntensity_emptyTank_FIMG, day_night_emptyTank_PTBC, day_night_emptyTank_FIMG);
    fill_dn_hists(f_out_t_arBot_ts, counts_Argon_PTBC, counts_Argon_FIMG, pulseIntensity_Argon_PTBC, pulseIntensity_Argon_FIMG, day_night_Argon_PTBC, day_night_Argon_FIMG);
    fill_dn_hists(f_out_t_arBotEmpty_ts, counts_EmptyArgon_PTBC, counts_EmptyArgon_FIMG, pulseIntensity_EmptyArgon_PTBC, pulseIntensity_EmptyArgon_FIMG, day_night_EmptyArgon_PTBC, day_night_EmptyArgon_FIMG);

    //Writing to the output file
    outputRootFile = new TFile("../rootFiles/day_night_15mins_bins.root","recreate");

    // Count plots
    counts_filterOut_PTBC->Write();
    counts_Bi_PTBC->Write();
    counts_Al5_PTBC->Write();
    counts_emptyTS_PTBC->Write();
    counts_emptyTank_PTBC->Write();
    counts_Argon_PTBC->Write();
    counts_EmptyArgon_PTBC->Write();

    counts_filterOut_FIMG->Write();
    counts_Bi_FIMG->Write();
    counts_Al5_FIMG->Write();
    counts_emptyTS_FIMG->Write();
    counts_emptyTank_FIMG->Write();
    counts_Argon_FIMG->Write();
    counts_EmptyArgon_FIMG->Write();

    // Pulse Intensity plots
    pulseIntensity_filterOut_PTBC->Write();
    pulseIntensity_Bi_PTBC->Write();
    pulseIntensity_Al5_PTBC->Write();
    pulseIntensity_emptyTS_PTBC->Write();
    pulseIntensity_emptyTank_PTBC->Write();
    pulseIntensity_Argon_PTBC->Write();
    pulseIntensity_EmptyArgon_PTBC->Write();

    pulseIntensity_filterOut_FIMG->Write();
    pulseIntensity_Bi_FIMG->Write();
    pulseIntensity_Al5_FIMG->Write();
    pulseIntensity_emptyTS_FIMG->Write();
    pulseIntensity_emptyTank_FIMG->Write();
    pulseIntensity_Argon_FIMG->Write();
    pulseIntensity_EmptyArgon_FIMG->Write();

    // Day-night plots
    day_night_filterOut_PTBC->Write();
    day_night_Bi_PTBC->Write();
    day_night_Al5_PTBC->Write();
    day_night_emptyTS_PTBC->Write();
    day_night_emptyTank_PTBC->Write();
    day_night_Argon_PTBC->Write();
    day_night_EmptyArgon_PTBC->Write();

    day_night_filterOut_FIMG->Write();
    day_night_Bi_FIMG->Write();
    day_night_Al5_FIMG->Write();
    day_night_emptyTS_FIMG->Write();
    day_night_emptyTank_FIMG->Write();
    day_night_Argon_FIMG->Write();
    day_night_EmptyArgon_FIMG->Write();

    outputRootFile->Close();

    std::cout << "Created output file 'day_night_15mins_bins.root'" << std::endl;
}