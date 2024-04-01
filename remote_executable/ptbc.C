/**
 * @file ptbc.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-11-13
 */

#include "ptbc.h"
/*
While changing the filter target
- check root file name
- check filter in runs
- check txt file name in endf() (not needed if onl changing thickness)
- check max line count in endf() (not needed if onl changing thickness)
- check names and titles of the plots
- check n_inverse
- check n passed to endf()
*/

// void StoreHist(){

//     TFile *f = new TFile(Form("rootFiles/crossSectionAna_%s.root", filter_name.c_str()),"recreate");

//     tof_hist_filter_in->Write();
//     energy_hist_filter_in->Write();
//     tof_hist_filter_out->Write();
//     energy_hist_filter_out->Write();
//     transmission_hist_e->Write();
//     cross_section_hist_e->Write();

//     f->Close();
// }

Double_t fillTofEnergyHists(std::vector<Int_t> run_list, TH1D* tof_hist, TH1D* energy_hist){

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
        Double_t tof = 0; //tof is in ns
        Double_t tflash = 0; //tflash is in ns
        Float_t amp = 0;
        Int_t BunchNumber_PTB = 0;
        Float_t PulseIntensity = 0;

        file_ntof->GetObject("PTBC", PTBC);
        PTBC->SetBranchAddress("BunchNumber", &BunchNumber_PTB);
        PTBC->SetBranchAddress("PulseIntensity", &PulseIntensity);
        PTBC->SetBranchAddress("tof", &tof);
        PTBC->SetBranchAddress("tflash", &tflash);
        PTBC->SetBranchAddress("amp", &amp);

        Long64_t Events_PTB = PTBC->GetEntriesFast();
        std::cout << "Number of entries = " << Events_PTB << std::endl;
        
        int CurrentBunchNum = 0;

        for (int j = 0; j < Events_PTB; j++)
        {
            PTBC->GetEntry(j);

            Double_t t_pkup = BNum_tpkup_map[BunchNumber_PTB];
            Double_t corrected_tof = tof - t_pkup + 660.0 + t_gamma_PTB;

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
                    if ( (Double_t) amp > yOnTheCutLine(t[k][0], a[k][0], t[k][1], a[k][1], corrected_tof) )
                    {
                        tof_hist->Fill(corrected_tof);
                        energy_hist->Fill( TOFToEnergy(corrected_tof * 1e-9) );
                        break;    
                    }
                }
            }
        }
        file_ntof->Close();
    }

    return NormFactor;

    // tof_hist_filter_out->Scale(1.0/NormFactor);
    // energy_hist_filter_out->Scale(1.0/NormFactor);
}

void ptbc(){

    fillCuts();
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
    transmission_hist_e = new TH1D("transmission_hist_e","Transmission Hist",num_bins_e,bin_edges_e);
    //cross section histogram
    cross_section_hist_e = new TH1D("cross_section_hist_e","Cross Section Hist",num_bins_e,bin_edges_e);
    //Filter In Hists
    tof_hist_filter_in = new TH1D("tof_hist_filter_in","Time of Flight Hist Filter In",num_bins_tof,bin_edges_tof);
    energy_hist_filter_in = new TH1D("energy_hist_filter_in","Energy Hist Filter In",num_bins_e,bin_edges_e);
    //Filter Out Hists
    tof_hist_filter_out = new TH1D("tof_hist_filter_out","Time of Flight Hist Filter Out",num_bins_tof,bin_edges_tof);
    energy_hist_filter_out = new TH1D("energy_hist_filter_out","Energy Hist Filter Out",num_bins_e,bin_edges_e);
    
    //Filling the Histograms
    Double_t norm_factor_filter_in = fillTofEnergyHists(filter_in_runs, tof_hist_filter_in, energy_hist_filter_in); 
    Double_t norm_factor_filter_out = fillTofEnergyHists(filter_out_runs, tof_hist_filter_out, energy_hist_filter_out);
    Double_t Qout_Qin = norm_factor_filter_out / norm_factor_filter_in;

    //Transmission
    for (int i = 0; i < num_bins_e; i++)
    {
        Double_t bin_content_in = energy_hist_filter_in->GetBinContent(i+1);
        Double_t bin_content_out = energy_hist_filter_out->GetBinContent(i+1);
        if (bin_content_in == 0 || bin_content_out == 0)
        {
            transmission_hist_e->SetBinContent(i+1, 0.);
            transmission_hist_e->SetBinError(i+1, 0.);
        } else {
            Double_t transmission = (bin_content_in * Qout_Qin)/bin_content_out;
            Double_t bin_unc = transmission * std::sqrt( (1./bin_content_in) + (1./bin_content_out) );
            transmission_hist_e->SetBinContent(i+1, transmission);
            transmission_hist_e->SetBinError(i+1, bin_unc);
        }
    }

    //Normalizing the Histograms
    tof_hist_filter_in->Scale(1.0/norm_factor_filter_in);
    energy_hist_filter_in->Scale(1.0/norm_factor_filter_in);
    tof_hist_filter_out->Scale(1.0/norm_factor_filter_out);
    energy_hist_filter_out->Scale(1.0/norm_factor_filter_out);

    //Cross Section
    Double_t num_density = num_density_map[filter_name];
    cout << Form("Number density of %s = ", filter_name.c_str()) << num_density << endl;
    Double_t n_inverse = ( (Double_t) 1.0/num_density);
    cout << Form("1/n of %s = ", filter_name.c_str()) << n_inverse << endl;

    for (int i = 0; i < num_bins_e; i++)
    {
        Double_t trans_bin_content = transmission_hist_e->GetBinContent(i+1);
        Double_t trans_bin_error = transmission_hist_e->GetBinError(i+1);
        if (trans_bin_content == 0)
        {
            cross_section_hist_e->SetBinContent(i+1, 0);
            cross_section_hist_e->SetBinError(i+1, 0);
        } else {
            Double_t cross_section = - n_inverse * std::log(trans_bin_content);
            Double_t bin_unc = n_inverse * (1./trans_bin_content) * trans_bin_error;
            cross_section_hist_e->SetBinContent(i+1, cross_section);
            cross_section_hist_e->SetBinError(i+1, bin_unc);
        }
    }

    cout << Form("Total Protons %s = ", filter_name.c_str()) << norm_factor_filter_in << endl;
    cout << "Total Protons Filter Out = " << norm_factor_filter_out << endl;

    //Writing to the output file
    outputRootFile = new TFile(Form("../rootFiles/crossSectionAna_%s.root", filter_name.c_str()),"recreate");
    tof_hist_filter_in->Write();
    energy_hist_filter_in->Write();
    tof_hist_filter_out->Write();
    energy_hist_filter_out->Write();
    transmission_hist_e->Write();
    cross_section_hist_e->Write();

    outputRootFile->Close();
}