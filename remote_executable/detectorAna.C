/**
 * @file detectorAna.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-11-13
 */

#include "detectorAna.h"

Double_t fillEnergyHist(std::vector<Int_t> run_list, TH1D* energy_hist_PTB, TH1D* energy_hist_FIMG){

    Double_t NormFactor = 0; //Integral of the pulse intensity

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

        //PTBC ---------------------------------------------
        TTree* PTBC;
        Double_t tof_PTB = 0; //tof is in ns
        Float_t amp_PTB = 0;
        Int_t BunchNumber_PTB = 0;
        Float_t PulseIntensity = 0;

        file_ntof->GetObject("PTBC", PTBC);
        PTBC->SetBranchAddress("BunchNumber", &BunchNumber_PTB);
        PTBC->SetBranchAddress("PulseIntensity", &PulseIntensity);
        PTBC->SetBranchAddress("tof", &tof_PTB);
        PTBC->SetBranchAddress("amp", &amp_PTB);

        Long64_t Events_PTB = PTBC->GetEntriesFast();
        std::cout << "Number of entries - PTBC = " << Events_PTB << std::endl;
        
        int CurrentBunchNum = 0;

        for (int j = 0; j < Events_PTB; j++)
        {
            PTBC->GetEntry(j);

            Double_t t_pkup = BNum_tpkup_map[BunchNumber_PTB];
            Double_t corrected_tof = tof_PTB - t_pkup + delT_pkup_ptbc + t_gamma_PTB;

            if (CurrentBunchNum != BunchNumber_PTB)
            {
                CurrentBunchNum = BunchNumber_PTB;
                NormFactor += (Double_t) PulseIntensity;
            }

            //Filling the histograms
            for (int k = 0; k < 6; k++)
            {
                if (corrected_tof >= t[k][0] && corrected_tof < t[k][1])
                {
                    if ( (Double_t) amp_PTB > yOnTheCutLinePTBC(t[k][0], a[k][0], t[k][1], a[k][1], corrected_tof) )
                    {
                        energy_hist_PTB->Fill( TOFToEnergy(corrected_tof * 1e-9) );
                        break;    
                    }
                }
            }
        }

        //FIMG ---------------------------------------------
        TTree* FIMG;
        Double_t tof_FIMG = 0; //tof is in ns
        Float_t amp = 0;
        Float_t area_0 = 0;
        Int_t det_num = 0;
        Int_t BunchNumber_FIMG = 0;

        file_ntof->GetObject("FIMG", FIMG);
        FIMG->SetBranchAddress("BunchNumber", &BunchNumber_FIMG);
        FIMG->SetBranchAddress("tof", &tof_FIMG);
        FIMG->SetBranchAddress("amp", &amp);
        FIMG->SetBranchAddress("area_0", &area_0);
        FIMG->SetBranchAddress("detn", &det_num);

        Long64_t Events_FIMG = FIMG->GetEntriesFast();
        std::cout << "Number of entries - FIMG = " << Events_FIMG << std::endl;

        for (int j = 0; j < Events_FIMG; j++)
        {
            FIMG->GetEntry(j);

            Double_t t_pkup = BNum_tpkup_map[BunchNumber_FIMG];
            Double_t corrected_tof = tof_FIMG - t_pkup + delT_pkup_fimg + t_gamma_FIMG;

            Double_t tof_cut_low = 0;
            Double_t tof_cut_up = 0;
            Double_t min_amp_cut = 0;

            if (det_num == 1){
                tof_cut_low = FIMG_tof_cut_low_det1;
                tof_cut_up = FIMG_tof_cut_up_det1;
                min_amp_cut = FIMG_min_amp_cut_det1;
            }
            else if (det_num == 2){
                tof_cut_low = FIMG_tof_cut_low_det2;
                tof_cut_up = FIMG_tof_cut_up_det2;
                min_amp_cut = FIMG_min_amp_cut_det2;
            }

            //Filling the histograms after cuts
            if (corrected_tof < tof_cut_low)
            {
                continue;
            }

            if ((Double_t) amp < min_amp_cut)
            {
                continue;
            }

            if (corrected_tof > tof_cut_up)
            {
                energy_hist_FIMG->Fill( TOFToEnergy(corrected_tof * 1e-9) );
                continue;
            }

            if ((Double_t) amp >= fimgCutFunction(corrected_tof, det_num))
            {
                energy_hist_FIMG->Fill( TOFToEnergy(corrected_tof * 1e-9) );
            }
        }

        file_ntof->Close();
    }

    return NormFactor;

    // tof_hist_filter_out->Scale(1.0/NormFactor);
    // energy_hist_filter_out->Scale(1.0/NormFactor);
}

void fillRunHists(Int_t num_bins_e, Double_t bin_edges_e[]){

    TH1D* energy_hist_target_in_PTB = 0;
    TH1D* energy_hist_target_in_FIMG = 0;
    TH1D* energy_hist_target_out_PTB = 0;
    TH1D* energy_hist_target_out_FIMG = 0;

    //Target In Hists
    energy_hist_target_in_PTB = new TH1D("energy_hist_target_in_PTB","Energy Hist Target In - PTBC",num_bins_e,bin_edges_e);
    energy_hist_target_in_FIMG = new TH1D("energy_hist_target_in_FIMG","Energy Hist Target In - FIMG",num_bins_e,bin_edges_e);
    //Target Out Hists
    energy_hist_target_out_PTB = new TH1D("energy_hist_target_out_PTB","Energy Hist Target Out - PTBC",num_bins_e,bin_edges_e);
    energy_hist_target_out_FIMG = new TH1D("energy_hist_target_out_FIMG","Energy Hist Target Out - FIMG",num_bins_e,bin_edges_e);
    
    //Filling the Histograms
    Double_t norm_factor_target_in = fillEnergyHist(ts_target_in_runs, energy_hist_target_in_PTB, energy_hist_target_in_FIMG);
    Double_t norm_factor_target_out = 0;
    if(!target_out_name.compare("ar_bottle")) {
        norm_factor_target_out = fillEnergyHist(ts_target_out_runs, energy_hist_target_out_PTB, energy_hist_target_out_FIMG);
    } else if(!target_out_name.compare("none_ts")) {
        norm_factor_target_out = fillEnergyHist(f_out_t_out_ts, energy_hist_target_out_PTB, energy_hist_target_out_FIMG);
    } else if(!target_out_name.compare("none")) {
        norm_factor_target_out = fillEnergyHist(f_out_t_out_runs, energy_hist_target_out_PTB, energy_hist_target_out_FIMG);
    }
    Double_t Qout_Qin = norm_factor_target_out / norm_factor_target_in;

    //Transmission
    for (int i = 0; i < num_bins_e; i++)
    {
        //PTBC
        Double_t bin_content_in_PTB = energy_hist_target_in_PTB->GetBinContent(i+1);
        Double_t bin_content_out_PTB = energy_hist_target_out_PTB->GetBinContent(i+1);
        if (bin_content_in_PTB == 0. || bin_content_out_PTB == 0.)
        {
            transmission_hist_e_PTB->SetBinContent(i+1, 0.);
            transmission_hist_e_PTB->SetBinError(i+1, 0.);
        } else {
            Double_t transmission_PTB = (bin_content_in_PTB * Qout_Qin)/bin_content_out_PTB;
            Double_t bin_unc_PTB = transmission_PTB * std::sqrt( (1./bin_content_in_PTB) + (1./bin_content_out_PTB) );
            transmission_hist_e_PTB->SetBinContent(i+1, transmission_PTB);
            transmission_hist_e_PTB->SetBinError(i+1, bin_unc_PTB);
        }

        //FIMG
        Double_t bin_content_in_FIMG = energy_hist_target_in_FIMG->GetBinContent(i+1);
        Double_t bin_content_out_FIMG = energy_hist_target_out_FIMG->GetBinContent(i+1);
        if (bin_content_in_FIMG == 0. || bin_content_out_FIMG == 0.)
        {
            transmission_hist_e_FIMG->SetBinContent(i+1, 0.);
            transmission_hist_e_FIMG->SetBinError(i+1, 0.);
        } else {
            Double_t transmission_FIMG = (bin_content_in_FIMG * Qout_Qin)/bin_content_out_FIMG;
            Double_t bin_unc_FIMG = transmission_FIMG * std::sqrt( (1./bin_content_in_FIMG) + (1./bin_content_out_FIMG) );
            transmission_hist_e_FIMG->SetBinContent(i+1, transmission_FIMG);
            transmission_hist_e_FIMG->SetBinError(i+1, bin_unc_FIMG);
        }
    }

    //Cross Section
    Double_t num_density = num_density_map[target_name];
    cout << Form("Number density of %s = ", target_name.c_str()) << num_density << endl;
    Double_t n_inverse = ( (Double_t) 1.0/num_density);
    cout << Form("1/n of %s = ", target_name.c_str()) << n_inverse << endl;

    for (int i = 0; i < num_bins_e; i++)
    {
        //PTBC
        Double_t trans_bin_content_PTB = transmission_hist_e_PTB->GetBinContent(i+1);
        Double_t trans_bin_error_PTB = transmission_hist_e_PTB->GetBinError(i+1);
        if (trans_bin_content_PTB == 0)
        {
            cross_section_hist_e_PTB->SetBinContent(i+1, 0);
            cross_section_hist_e_PTB->SetBinError(i+1, 0);
        } else {
            Double_t cross_section_PTB = - n_inverse * std::log(trans_bin_content_PTB);
            Double_t bin_unc_PTB = n_inverse * (1./trans_bin_content_PTB) * trans_bin_error_PTB;
            cross_section_hist_e_PTB->SetBinContent(i+1, cross_section_PTB);
            cross_section_hist_e_PTB->SetBinError(i+1, bin_unc_PTB);
        }

        //FIMG
        Double_t trans_bin_content_FIMG = transmission_hist_e_FIMG->GetBinContent(i+1);
        Double_t trans_bin_error_FIMG = transmission_hist_e_FIMG->GetBinError(i+1);
        if (trans_bin_content_FIMG == 0)
        {
            cross_section_hist_e_FIMG->SetBinContent(i+1, 0);
            cross_section_hist_e_FIMG->SetBinError(i+1, 0);
        } else {
            Double_t cross_section_FIMG = - n_inverse * std::log(trans_bin_content_FIMG);
            Double_t bin_unc_FIMG = n_inverse * (1./trans_bin_content_FIMG) * trans_bin_error_FIMG;
            cross_section_hist_e_FIMG->SetBinContent(i+1, cross_section_FIMG);
            cross_section_hist_e_FIMG->SetBinError(i+1, bin_unc_FIMG);
        }
    }

    cout << Form("Total Protons %s = ", target_name.c_str()) << norm_factor_target_in << endl;
    cout << "Total Protons Target Out = " << norm_factor_target_out << endl;
}

void detectorAna(){

    fillCutsPTBC();
    fillRuns();
    fillNumDensityMap();

    //Calculating TOF bin edges
    int num_decades = 5;
    int num_bins_tof = bins_per_decade * num_decades;
    Double_t bin_edges_tof[num_bins_tof+1];
    Double_t step_tof = ((Double_t) 1.0/(Double_t) bins_per_decade);
    for(int i = 0; i < num_bins_tof+1; i++)
    {
        Double_t base = 10.;
        Double_t exponent = (step_tof * (Double_t) i) + 3.;
        bin_edges_tof[i] = (Double_t) std::pow(base, exponent);
    }

    //Calculating Energy bin edges
    Double_t tof_min = 1e3; //ns
    Double_t tof_max = 1e8; //ns
    Double_t e_min = TOFToEnergy(tof_max * 1e-9); //converting into seconds
    Double_t e_max = TOFToEnergy(tof_min * 1e-9); //converting into seconds
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

    cout << "Number of tof bins = " << num_bins_tof << endl;
    cout << "Number of e bins = " << num_bins_e << endl;

    //transmission histogram
    transmission_hist_e_PTB = new TH1D("transmission_hist_e_PTB","Transmission Hist - PTBC",num_bins_e,bin_edges_e);
    transmission_hist_e_FIMG = new TH1D("transmission_hist_e_FIMG","Transmission Hist - FIMG",num_bins_e,bin_edges_e);
    //cross section histogram
    cross_section_hist_e_PTB = new TH1D("cross_section_hist_e_PTB","Cross Section Hist - PTBC",num_bins_e,bin_edges_e);
    cross_section_hist_e_FIMG = new TH1D("cross_section_hist_e_FIMG","Cross Section Hist - FIMG",num_bins_e,bin_edges_e);
    
    // fillTSRunsHists(num_bins_e,bin_edges_e);

    fillRunHists(num_bins_e,bin_edges_e);

    //Writing to the output file
    outputRootFile = new TFile(Form("../rootFiles/crossSectionAna_%s.root", target_name.c_str()),"recreate");
    // trans_hist_fOut->Write();
    // trans_hist_fIn->Write();
    // trans_hist_fOut_endf->Write();
    // trans_hist_fIn_endf->Write();
    transmission_hist_e_PTB->Write();
    transmission_hist_e_FIMG->Write();
    cross_section_hist_e_PTB->Write();
    cross_section_hist_e_FIMG->Write();

    outputRootFile->Close();
}