/**
 * @file PTBC_cuts_ana.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2024-09-06
 */

 #include "PTBC_cuts_ana.h"

Double_t fillEnergyHist(std::vector<Int_t> run_list, TH1D* energy_hist[]){

    Double_t NormFactor = 0; //Integral of the pulse intensity

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
            Double_t corrected_tof = tof_PTBC - t_pkup + delT_pkup_ptbc + t_gamma_PTB;

            if (CurrentBunchNum != BunchNumber_PTBC)
            {
                CurrentBunchNum = BunchNumber_PTBC;
                NormFactor += (Double_t) PulseIntensity;
            }

            if (corrected_tof < min_tof_PTBC)
            {
                continue;
            }

            Double_t tof_cut_bin = PTBC_tof_amp_cuts[det_num_PTBC-2]->GetXaxis()->FindBin(corrected_tof);
            Double_t amp_cut = PTBC_tof_amp_cuts[det_num_PTBC-2]->GetBinContent(tof_cut_bin);
            Int_t tf_bin_num = transfer_function_PTBC->GetXaxis()->FindBin(corrected_tof);
            Double_t neutron_e = transfer_function_PTBC->GetBinContent(tf_bin_num);

            if ( (Double_t) amp_PTBC > amp_cut)
            {
                energy_hist[0]->Fill(neutron_e);
            }

            if ( (Double_t) amp_PTBC > (amp_cut + amp_cut*0.2))
            {
                energy_hist[1]->Fill(neutron_e);
            }

            if ( (Double_t) amp_PTBC > (amp_cut + amp_cut*0.5))
            {
                energy_hist[2]->Fill(neutron_e);
            }

            if ( (Double_t) amp_PTBC > (amp_cut - amp_cut*0.2))
            {
                energy_hist[3]->Fill(neutron_e);
            }

            if ( (Double_t) amp_PTBC > (amp_cut - amp_cut*0.5))
            {
                energy_hist[4]->Fill(neutron_e);
            }
        }
        file_ntof->Close();
    }

    return NormFactor;
}

void fillRunHists(){

    Double_t norm_factor_target_in = 0;
    Double_t norm_factor_target_out = 0;
    
    //Filling the Histograms
    norm_factor_target_in = fillEnergyHist(ts_target_in_runs, energy_hist_target_in);
    norm_factor_target_out = 0;
    if(!target_out_name.compare("ar_bottle")) {
        norm_factor_target_out = fillEnergyHist(ts_target_out_runs, energy_hist_target_out);
    } else if(!target_out_name.compare("none_ts")) {
        norm_factor_target_out = fillEnergyHist(f_out_t_out_ts, energy_hist_target_out);
    } else if(!target_out_name.compare("none")) {
        norm_factor_target_out = fillEnergyHist(f_out_t_out_runs, energy_hist_target_out);
    }

    norm_factors.push_back(norm_factor_target_in);
    norm_factors.push_back(norm_factor_target_out);
    // Double_t Qout_Qin = norm_factor_target_out / norm_factor_target_in;

    cout << Form("Total Protons %s = ", target_name.c_str()) << norm_factor_target_in << endl;
    cout << "Total Protons Target Out = " << norm_factor_target_out << endl;
}

void PTBC_cuts_ana(){

    fillRuns();

    //Get Transfer Functions
    transfer_function_PTBC = GetHist1D("../inputFiles/transfer_function.root", "transfer_func_mean_PTBC_5itr");

    //Get PTBC Cuts
    for (Int_t i = 0; i < 6; i++)
    {
        PTBC_tof_amp_cuts[i] = GetHist1D("../inputFiles/PTBC_cuts.root", Form("PTBC_cuts_det%i", i+2));
    }

    //Calculating Energy bin edges
    Double_t tof_min = 1e3; //ns
    Double_t tof_max = 1e8; //ns
    Double_t e_min = TOFToEnergy(tof_max * 1e-9, flight_path_length_PTBC); //converting into seconds
    Double_t e_max = TOFToEnergy(tof_min * 1e-9, flight_path_length_PTBC); //converting into seconds
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

    cout << "Number of e bins = " << num_bins_e << endl;

    energy_hist_target_in[0] = new TH1D("energy_hist_target_in_thecut","Energy Hist Target In - The Cut", num_bins_e, bin_edges_e);
    energy_hist_target_in[1] = new TH1D("energy_hist_target_in_20more","Energy Hist Target In - 20%% More", num_bins_e, bin_edges_e);
    energy_hist_target_in[2] = new TH1D("energy_hist_target_in_50more","Energy Hist Target In - 50%% More", num_bins_e, bin_edges_e);
    energy_hist_target_in[3] = new TH1D("energy_hist_target_in_20less","Energy Hist Target In - 20%% Less", num_bins_e, bin_edges_e);
    energy_hist_target_in[4] = new TH1D("energy_hist_target_in_50less","Energy Hist Target In - 50%% Less", num_bins_e, bin_edges_e);

    energy_hist_target_out[0] = new TH1D("energy_hist_target_out_thecut","Energy Hist Target Out - The Cut", num_bins_e, bin_edges_e);
    energy_hist_target_out[1] = new TH1D("energy_hist_target_out_20more","Energy Hist Target Out - 20%% More", num_bins_e, bin_edges_e);
    energy_hist_target_out[2] = new TH1D("energy_hist_target_out_50more","Energy Hist Target Out - 50%% More", num_bins_e, bin_edges_e);
    energy_hist_target_out[3] = new TH1D("energy_hist_target_out_20less","Energy Hist Target Out - 20%% Less", num_bins_e, bin_edges_e);
    energy_hist_target_out[4] = new TH1D("energy_hist_target_out_50less","Energy Hist Target Out - 50%% Less", num_bins_e, bin_edges_e);

    fillRunHists();

    //Writing to the output file
    outputRootFile = new TFile(Form("../rootFiles/PTBC_cuts_ana_%s_%iBPD.root", target_name.c_str(), bins_per_decade),"recreate");

    for (Int_t i = 0; i < 5; i++)
    {
        energy_hist_target_in[i]->Write();
    }
    for (Int_t i = 0; i < 5; i++)
    {
        energy_hist_target_out[i]->Write();
    }

    outputRootFile->WriteObject(&norm_factors, "norm_factors");
    outputRootFile->Close();

    std::cout << Form("Created output file 'PTBC_cuts_ana_%s_%iBPD.root'", target_name.c_str(), bins_per_decade) << std::endl;
}