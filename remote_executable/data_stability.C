/**
 * @file data_stability.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-12-04
 */

#include "data_stability.h"

void fillTOFHist(std::vector<Int_t> run_list, TH1D* energy_hist_PTBC, TH1D* energy_hist_FIMG, Double_t bin_edges_tof_PTBC[], Double_t bin_edges_tof_FIMG[]){

    Double_t NormFactor = 0; //Integral of the pulse intensity

    Int_t num_bins_PTBC = energy_hist_PTBC->GetNbinsX();
    Int_t num_bins_FIMG = energy_hist_FIMG->GetNbinsX();

    TH1D* tof_hist_PTBC = new TH1D("tof_hist_PTBC","ToF Hist PTBC", num_bins_PTBC, bin_edges_tof_PTBC);
    TH1D* tof_hist_FIMG = new TH1D("tof_hist_FIMG","ToF Hist FIMG", num_bins_FIMG, bin_edges_tof_FIMG);

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
        Int_t det_num_PTBC = 0;
        Int_t BunchNumber_PTBC = 0;
        Float_t PulseIntensity = 0;

        file_ntof->GetObject("PTBC", PTBC);
        PTBC->SetBranchAddress("BunchNumber", &BunchNumber_PTBC);
        PTBC->SetBranchAddress("PulseIntensity", &PulseIntensity);
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
                NormFactor += (Double_t) PulseIntensity;
            }

            applyMyCuts_PTBC(corrected_tof, amp_PTBC, PulseIntensity, det_num_PTBC, tof_hist_PTBC, run_list.at(i));
            // applynTOFCuts_PTBC(corrected_tof, amp_PTBC, PulseIntensity, det_num_PTBC, tof_hist_PTBC);
        }

        //FIMG ---------------------------------------------
        TTree* FIMG;
        Double_t tof_FIMG = 0; //tof is in ns
        Float_t amp_FIMG = 0;
        Int_t det_num_FIMG = 0;
        Int_t BunchNumber_FIMG = 0;

        file_ntof->GetObject("FIMG", FIMG);
        FIMG->SetBranchAddress("BunchNumber", &BunchNumber_FIMG);
        FIMG->SetBranchAddress("tof", &tof_FIMG);
        FIMG->SetBranchAddress("amp", &amp_FIMG);
        FIMG->SetBranchAddress("detn", &det_num_FIMG);

        Long64_t Events_FIMG = FIMG->GetEntriesFast();
        std::cout << "Number of entries - FIMG = " << Events_FIMG << std::endl;

        for (int j = 0; j < Events_FIMG; j++)
        {
            FIMG->GetEntry(j);

            Double_t t_pkup = BNum_tpkup_map[BunchNumber_FIMG];
            Double_t corrected_tof = tof_FIMG - t_pkup + delT_pkup_fimg + t_gamma_FIMG;

            applyMyCuts_FIMG(corrected_tof, amp_FIMG, det_num_FIMG, tof_hist_FIMG);
        }

        file_ntof->Close();
    }

    // Int_t num_bins_PTBC = energy_hist_PTBC->GetNbinsX();
    for(Int_t i = 0; i < num_bins_PTBC; i++){
        energy_hist_PTBC->SetBinContent(i+1, tof_hist_PTBC->GetBinContent(num_bins_PTBC-i));
    }

    // Int_t num_bins_FIMG = energy_hist_FIMG->GetNbinsX();
    for(Int_t i = 0; i < num_bins_FIMG; i++){
        energy_hist_FIMG->SetBinContent(i+1, tof_hist_FIMG->GetBinContent(num_bins_FIMG-i));
    }

    delete tof_hist_PTBC;
    delete tof_hist_FIMG;

    energy_hist_PTBC->Scale(1.0/NormFactor);
    energy_hist_FIMG->Scale(1.0/NormFactor);
}

void fill_dn_hists(std::vector<Int_t> run_list, TH1D* dn_hist_PTBC, TH1D* dn_hist_FIMG){

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

            if (!select_hit_PTBC(corrected_tof, amp_PTBC, PulseIntensity_PTBC, det_num_PTBC))
            {
                continue;
            }
            
            hit_time = timeToSeconds(std::to_string(time_PTBC));

            if (j == 0)
            {
                CurrentBunchNum = BunchNumber_PTBC;
                dn_hist_PTBC->Fill( hit_time );
                entry_bin_num = dn_hist_PTBC->FindBin( hit_time );
                PTBC_PulseIntensity[entry_bin_num-1] += PulseIntensity_PTBC;
                continue;
            }   

            dn_hist_PTBC->Fill( hit_time );

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
                entry_bin_num = dn_hist_FIMG->FindBin( hit_time );
                FIMG_PulseIntensity[entry_bin_num-1] += PulseIntensity_FIMG;
                continue;
            }   

            dn_hist_FIMG->Fill( hit_time );
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
        Int_t bin_error = dn_hist_PTBC->GetBinError(i+1);
        if (PTBC_PulseIntensity[i] != 0){
            dn_hist_PTBC->SetBinContent(i+1, (Double_t) bin_content/(Double_t) PTBC_PulseIntensity[i] );
            dn_hist_PTBC->SetBinError(i+1, (Double_t) bin_error/(Double_t) PTBC_PulseIntensity[i] );
        }
    }

    for (Int_t i = 0; i < num_bins_FIMG; i++)
    {
        Int_t bin_content = dn_hist_FIMG->GetBinContent(i+1);
        Int_t bin_error = dn_hist_FIMG->GetBinError(i+1);
        if (FIMG_PulseIntensity[i] != 0){
            dn_hist_FIMG->SetBinContent(i+1, bin_content/FIMG_PulseIntensity[i] );
            dn_hist_FIMG->SetBinError(i+1, bin_error/FIMG_PulseIntensity[i] );
        }
    }
    
}

void data_stability(){

    fillCutsPTBC();
    fillCutsFIMG();
    fillNumDensityMap();

    // Int_t num_bins = 240; // 4 hours
    // Double_t x_min = 0.;
    // Double_t y_min = 14400.;

    // stability_hist_PTBC = new TH1D("stability_hist_PTBC","Stability Hist - PTBC",num_bins,x_min,y_min);
    // stability_hist_FIMG_det1 = new TH1D("stability_hist_FIMG_det1","Stability Hist - FIMG - Det 1",num_bins,x_min,y_min);
    // stability_hist_FIMG_det2 = new TH1D("stability_hist_FIMG_det2","Stability Hist - FIMG - Det 1",num_bins,x_min,y_min);

    // Calculating TOF (x) bin edges
    Int_t num_decades_tof = 6;
    Int_t num_bins_tof = bins_per_decade * num_decades_tof;
    Double_t bin_edges_tof[num_bins_tof+1];
    Double_t step_tof = ((Double_t) 1.0/(Double_t) bins_per_decade);
    for(Int_t i = 0; i < num_bins_tof+1; i++)
    {
        Double_t base = 10.;
        Double_t exponent = (step_tof * (Double_t) i) + 2.;
        bin_edges_tof[i] = (Double_t) std::pow(base, exponent);
    }

    //Calculating Energy bin edges
    Double_t tof_min = 1e3; //ns
    Double_t tof_max = 1e8; //ns
    Double_t e_min = TOFToEnergy(tof_max * 1e-9, flight_path_length_FIMG); //converting into seconds
    Double_t e_max = TOFToEnergy(tof_min * 1e-9, flight_path_length_FIMG); //converting into seconds
    int min_power = FindDecadePower(e_min);
    int max_power = FindDecadePower(e_max);
    int num_decades_e = max_power - min_power;
    int num_bins_e = bins_per_decade * num_decades_e;
    Double_t bin_edges_e[num_bins_e+1];
    Double_t step_e = ((Double_t) 1.0/(Double_t) bins_per_decade);
    for(int i = 0; i < num_bins_e+1; i++)
    {
        Double_t base = 10.;
        Double_t exponent = (step_e * (Double_t) i) + (Double_t) min_power;
        bin_edges_e[i] = (Double_t) std::pow(base, exponent);
    }

    //Getting tof edges
    Double_t bin_edges_tof_PTBC[num_bins_e+1];
    for(Int_t i = 0; i < num_bins_e+1; i++){
        bin_edges_tof_PTBC[i] = EnergyToTOF(bin_edges_e[num_bins_e-i], flight_path_length_PTB) * 1e9;
    }

    Double_t bin_edges_tof_FIMG[num_bins_e+1];
    for(Int_t i = 0; i < num_bins_e+1; i++){
        bin_edges_tof_FIMG[i] = EnergyToTOF(bin_edges_e[num_bins_e-i], flight_path_length_FIMG) * 1e9;
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

    // Initializing Stability plots
    // PTBC
    // norm_counts_empty_sep16_PTBC = new TH1D("norm_counts_empty_sep16_PTBC","Norm Counts - Empty Sep 16 - PTBC", num_bins_e, bin_edges_e);
    norm_counts_empty_sep19_PTBC = new TH1D("norm_counts_empty_sep19_PTBC","Norm Counts - Empty Sep 19 - PTBC", num_bins_e, bin_edges_e);
    norm_counts_empty_oct01_PTBC = new TH1D("norm_counts_empty_oct01_PTBC","Norm Counts - Empty Oct 01 - PTBC", num_bins_e, bin_edges_e);

    norm_counts_Bi_sep17_PTBC = new TH1D("norm_counts_Bi_sep17_PTBC","Norm Counts - Bi Sep 17 - PTBC", num_bins_e, bin_edges_e);
    norm_counts_Bi_sep23_PTBC = new TH1D("norm_counts_Bi_sep23_PTBC","Norm Counts - Bi Sep 23 - PTBC", num_bins_e, bin_edges_e);

    norm_counts_Al5_sep25_PTBC = new TH1D("norm_counts_Al5_sep25_PTBC","Norm Counts - Al (5cm) Sep 25 - PTBC", num_bins_e, bin_edges_e);
    norm_counts_Al5_oct02_PTBC = new TH1D("norm_counts_Al5_oct02_PTBC","Norm Counts - Al (5cm) Oct 02 - PTBC", num_bins_e, bin_edges_e);

    norm_counts_emptyTank_oct12_PTBC = new TH1D("norm_counts_emptyTank_oct12_PTBC","Norm Counts - Empty Tank Oct 12 - PTBC", num_bins_e, bin_edges_e);
    norm_counts_emptyTank_oct16_PTBC = new TH1D("norm_counts_emptyTank_oct16_PTBC","Norm Counts - Empty Tank Oct 16 - PTBC", num_bins_e, bin_edges_e);

    norm_counts_ArgonFull_oct18_PTBC = new TH1D("norm_counts_ArgonFull_oct18_PTBC","Norm Counts - Argon Tank Oct 18 - PTBC", num_bins_e, bin_edges_e);
    norm_counts_ArgonFull_oct22_PTBC = new TH1D("norm_counts_ArgonFull_oct22_PTBC","Norm Counts - Argon Tank Oct 22 - PTBC", num_bins_e, bin_edges_e);

    // FIMG
    norm_counts_empty_sep19_FIMG = new TH1D("norm_counts_empty_sep19_FIMG","Norm Counts - Empty Sep 19 - FIMG", num_bins_e, bin_edges_e);
    norm_counts_empty_oct01_FIMG = new TH1D("norm_counts_empty_oct01_FIMG","Norm Counts - Empty Oct 01 - FIMG", num_bins_e, bin_edges_e);

    norm_counts_Bi_sep17_FIMG = new TH1D("norm_counts_Bi_sep17_FIMG","Norm Counts - Bi Sep 17 - FIMG", num_bins_e, bin_edges_e);
    norm_counts_Bi_sep23_FIMG = new TH1D("norm_counts_Bi_sep23_FIMG","Norm Counts - Bi Sep 23 - FIMG", num_bins_e, bin_edges_e);

    norm_counts_Al5_sep25_FIMG = new TH1D("norm_counts_Al5_sep25_FIMG","Norm Counts - Al (5cm) Sep 25 - FIMG", num_bins_e, bin_edges_e);
    norm_counts_Al5_oct02_FIMG = new TH1D("norm_counts_Al5_oct02_FIMG","Norm Counts - Al (5cm) Oct 02 - FIMG", num_bins_e, bin_edges_e);

    norm_counts_emptyTank_oct12_FIMG = new TH1D("norm_counts_emptyTank_oct12_FIMG","Norm Counts - Empty Tank Oct 12 - FIMG", num_bins_e, bin_edges_e);
    norm_counts_emptyTank_oct16_FIMG = new TH1D("norm_counts_emptyTank_oct16_FIMG","Norm Counts - Empty Tank Oct 16 - FIMG", num_bins_e, bin_edges_e);

    norm_counts_ArgonFull_oct18_FIMG = new TH1D("norm_counts_ArgonFull_oct18_FIMG","Norm Counts - Argon Tank Oct 18 - FIMG", num_bins_e, bin_edges_e);
    norm_counts_ArgonFull_oct22_FIMG = new TH1D("norm_counts_ArgonFull_oct22_FIMG","Norm Counts - Argon Tank Oct 22 - FIMG", num_bins_e, bin_edges_e);

    // Initializing day night plots
    // PTBC
    // day_night_Bi_sep18_PTBC = new TH1D("day_night_Bi_sep18_PTBC","Day-Night Effect - Bi Sep 18 - PTBC", num_bins_dn, bin_edges_dn);
    // day_night_Al5_sep27_PTBC = new TH1D("day_night_Al5_sep27_PTBC","Day-Night Effect - Al (5 cm) Sep 27 - PTBC", num_bins_dn, bin_edges_dn);
    // day_night_emptyTS_oct05_PTBC = new TH1D("day_night_emptyTS_oct05_PTBC","Day-Night Effect - Empty (TS) Oct 05 - PTBC", num_bins_dn, bin_edges_dn);
    // day_night_emptyTank_oct15_PTBC = new TH1D("day_night_emptyTank_oct15_PTBC","Day-Night Effect - Empty Tank Oct 15 - PTBC", num_bins_dn, bin_edges_dn);
    // day_night_Argon_oct22_PTBC = new TH1D("day_night_Argon_oct22_PTBC","Day-Night Effect - Argon Oct 22 - PTBC", num_bins_dn, bin_edges_dn);

    // FIMG
    // day_night_Bi_sep18_FIMG = new TH1D("day_night_Bi_sep18_FIMG","Day-Night Effect - Bi Sep 18 - FIMG", num_bins_dn, bin_edges_dn);
    // day_night_Al5_sep27_FIMG = new TH1D("day_night_Al5_sep27_FIMG","Day-Night Effect - Al (5 cm) Sep 27 - FIMG", num_bins_dn, bin_edges_dn);
    // day_night_emptyTS_oct05_FIMG = new TH1D("day_night_emptyTS_oct05_FIMG","Day-Night Effect - Empty (TS) Oct 05 - FIMG", num_bins_dn, bin_edges_dn);
    // day_night_emptyTank_oct15_FIMG = new TH1D("day_night_emptyTank_oct15_FIMG","Day-Night Effect - Empty Tank Oct 15 - FIMG", num_bins_dn, bin_edges_dn);
    // day_night_Argon_oct22_FIMG = new TH1D("day_night_Argon_oct22_FIMG","Day-Night Effect - Argon Oct 22 - FIMG", num_bins_dn, bin_edges_dn);

    // Run Lists (FIMG starts from 117386)
    // Beam Stability
    // std::vector<Int_t> run_list_empty_sep16 = {117357, 117358, 117359, 117362, 117363}; // No FIMG in these runs
    std::vector<Int_t> run_list_empty_sep19 = {117405, 117406, 117408, 117409, 117410, 117411}; //Sep 19 and 20th
    std::vector<Int_t> run_list_empty_oct01 = {117511, 117512, 117513, 117514, 117515};
    std::vector<Int_t> run_list_Bi_sep17 = {117386, 117387, 117388, 117389, 117390};
    std::vector<Int_t> run_list_Bi_sep23 = {117439, 117440, 117441, 117442, 117443};
    std::vector<Int_t> run_list_Al5_sep25 = {117462, 117463, 117464, 117465};
    std::vector<Int_t> run_list_Al5_oct02 = {117519, 117520, 117521, 117522};
    std::vector<Int_t> run_list_emptyTank_oct12 = {117612, 117613, 117614, 117615, 117616};
    std::vector<Int_t> run_list_emptyTank_oct16 = {117645, 117646, 117647, 117648, 117649};
    std::vector<Int_t> run_list_ArgonFull_oct18 = {117680, 117684, 117685, 117686, 117687};
    std::vector<Int_t> run_list_ArgonFull_oct22 = {117730, 117731, 117732, 117733, 117734, 117735};

    // Day-Night effect
    // first and last runs also include previous and the next day runs. This has to be removed
    // std::vector<Int_t> run_list_Bi_sep18 = {117389, 117390, 117391, 117392, 117393, 117394};
    // std::vector<Int_t> run_list_Al5_sep27 = {117472, 117473, 117474, 117475, 117476, 117477, 117479, 117480, 117481};
    // std::vector<Int_t> run_list_emptyTS_oct05 = {117543, 117544, 117545, 117546, 117547, 117548, 117549, 117550};
    // std::vector<Int_t> run_list_emptyTank_oct15 = {117632, 117633, 117634, 117635, 117636, 117637, 117638, 117639};
    // std::vector<Int_t> run_list_Argon_oct22 = {117729, 117730, 117731, 117732, 117733, 117734, 117735, 117736};

    // Filling the histograms
    // fillTOFHist(run_list_empty_sep16, norm_counts_empty_sep16_PTBC);

    // Stability hists
    fillTOFHist(run_list_empty_sep19, norm_counts_empty_sep19_PTBC, norm_counts_empty_sep19_FIMG, bin_edges_tof_PTBC, bin_edges_tof_FIMG);
    fillTOFHist(run_list_empty_oct01, norm_counts_empty_oct01_PTBC, norm_counts_empty_oct01_FIMG, bin_edges_tof_PTBC, bin_edges_tof_FIMG);
    fillTOFHist(run_list_Bi_sep17, norm_counts_Bi_sep17_PTBC, norm_counts_Bi_sep17_FIMG, bin_edges_tof_PTBC, bin_edges_tof_FIMG);
    fillTOFHist(run_list_Bi_sep23, norm_counts_Bi_sep23_PTBC, norm_counts_Bi_sep23_FIMG, bin_edges_tof_PTBC, bin_edges_tof_FIMG);
    fillTOFHist(run_list_Al5_sep25, norm_counts_Al5_sep25_PTBC, norm_counts_Al5_sep25_FIMG, bin_edges_tof_PTBC, bin_edges_tof_FIMG);
    fillTOFHist(run_list_Al5_oct02, norm_counts_Al5_oct02_PTBC, norm_counts_Al5_oct02_FIMG, bin_edges_tof_PTBC, bin_edges_tof_FIMG);
    fillTOFHist(run_list_emptyTank_oct12, norm_counts_emptyTank_oct12_PTBC, norm_counts_emptyTank_oct12_FIMG, bin_edges_tof_PTBC, bin_edges_tof_FIMG);
    fillTOFHist(run_list_emptyTank_oct16, norm_counts_emptyTank_oct16_PTBC, norm_counts_emptyTank_oct16_FIMG, bin_edges_tof_PTBC, bin_edges_tof_FIMG);
    fillTOFHist(run_list_ArgonFull_oct18, norm_counts_ArgonFull_oct18_PTBC, norm_counts_ArgonFull_oct18_FIMG, bin_edges_tof_PTBC, bin_edges_tof_FIMG);
    fillTOFHist(run_list_ArgonFull_oct22, norm_counts_ArgonFull_oct22_PTBC, norm_counts_ArgonFull_oct22_FIMG, bin_edges_tof_PTBC, bin_edges_tof_FIMG);

    // Day Night hists
    // fill_dn_hists(run_list_Bi_sep18, day_night_Bi_sep18_PTBC, day_night_Bi_sep18_FIMG);
    // fill_dn_hists(run_list_Al5_sep27, day_night_Al5_sep27_PTBC, day_night_Al5_sep27_FIMG);
    // fill_dn_hists(run_list_emptyTS_oct05, day_night_emptyTS_oct05_PTBC, day_night_emptyTS_oct05_FIMG);
    // fill_dn_hists(run_list_emptyTank_oct15, day_night_emptyTank_oct15_PTBC, day_night_emptyTank_oct15_FIMG);
    // fill_dn_hists(run_list_Argon_oct22, day_night_Argon_oct22_PTBC, day_night_Argon_oct22_FIMG);

    // Draw atributes
    norm_counts_empty_oct01_PTBC->SetLineColor(2);
    norm_counts_Bi_sep23_PTBC->SetLineColor(2);
    norm_counts_Al5_oct02_PTBC->SetLineColor(2);
    norm_counts_emptyTank_oct16_PTBC->SetLineColor(2);
    norm_counts_ArgonFull_oct22_PTBC->SetLineColor(2);

    norm_counts_empty_oct01_FIMG->SetLineColor(2);
    norm_counts_Bi_sep23_FIMG->SetLineColor(2);
    norm_counts_Al5_oct02_FIMG->SetLineColor(2);
    norm_counts_emptyTank_oct16_FIMG->SetLineColor(2);
    norm_counts_ArgonFull_oct22_FIMG->SetLineColor(2);

    //Writing to the output file
    outputRootFile = new TFile("../rootFiles/data_stability.root","recreate");

    // norm_counts_empty_sep16_PTBC->Write();

    norm_counts_empty_sep19_PTBC->Write();
    norm_counts_empty_oct01_PTBC->Write();
    norm_counts_Bi_sep17_PTBC->Write();
    norm_counts_Bi_sep23_PTBC->Write();
    norm_counts_Al5_sep25_PTBC->Write();
    norm_counts_Al5_oct02_PTBC->Write();
    norm_counts_emptyTank_oct12_PTBC->Write();
    norm_counts_emptyTank_oct16_PTBC->Write();
    norm_counts_ArgonFull_oct18_PTBC->Write();
    norm_counts_ArgonFull_oct22_PTBC->Write();

    norm_counts_empty_sep19_FIMG->Write();
    norm_counts_empty_oct01_FIMG->Write();
    norm_counts_Bi_sep17_FIMG->Write();
    norm_counts_Bi_sep23_FIMG->Write();
    norm_counts_Al5_sep25_FIMG->Write();
    norm_counts_Al5_oct02_FIMG->Write();
    norm_counts_emptyTank_oct12_FIMG->Write();
    norm_counts_emptyTank_oct16_FIMG->Write();
    norm_counts_ArgonFull_oct18_FIMG->Write();
    norm_counts_ArgonFull_oct22_FIMG->Write();

    // day_night_Bi_sep18_PTBC->Write();
    // day_night_Al5_sep27_PTBC->Write();
    // day_night_emptyTS_oct05_PTBC->Write();
    // day_night_emptyTank_oct15_PTBC->Write();
    // day_night_Argon_oct22_PTBC->Write();

    // day_night_Bi_sep18_FIMG->Write();
    // day_night_Al5_sep27_FIMG->Write();
    // day_night_emptyTS_oct05_FIMG->Write();
    // day_night_emptyTank_oct15_FIMG->Write();
    // day_night_Argon_oct22_FIMG->Write();

    outputRootFile->Close();

    std::cout << "Created output file 'data_stability.root'" << std::endl;
}