/**
 * @file ptbc_ts.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-11-13
 */

#include "ptbc_ts.h"
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

//     TFile *f = new TFile(Form("rootFiles/crossSectionAna_%s.root", target_name.c_str()),"recreate");

//     tof_hist_filter_in->Write();
//     energy_hist_filter_in->Write();
//     tof_hist_filter_out->Write();
//     energy_hist_filter_out->Write();
//     transmission_hist_e->Write();
//     cross_section_hist_e->Write();

//     f->Close();
// }

Double_t fillEnergyHist(std::vector<Int_t> run_list, TH1D* energy_hist){

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

void fillENDFHist(Double_t energy_bin_edges[]){
    //Extracting ENDF Cross Section and transmission

    std::ifstream inputFile("../evalData/Al_tot_xsec.txt"); // endf_file_name.c_str()

    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the file.\n";
    }

    //Extracting the data from the text file
    Double_t val1, val2; //val1 -> Energy (eV); val2 -> Cross Section (barns)
    std::string line;
    Double_t line_count = 0;
    std::vector<Double_t> energy;
    std::vector<Double_t> xsec;

    while (std::getline(inputFile, line)) {

        line_count++;

        if (line_count == 1)
        {
            continue;
        }

        if (line_count == 9745) //Bi - 27019; Al - 9745
        {
            break;
        }

        std::istringstream iss(line);

        if (iss >> val1 >> val2) {
            energy.push_back(val1);
            xsec.push_back(val2);
        } else {
            std::cerr << "Invalid data format.\n";
        }
    }

    //Filling the histograms
    Double_t xsec_sum = 0;
    Int_t sum_counter = 0;
    Int_t bin_counter = 1;

    for (int i = 0; i < energy.size(); i++)
    {
        if (energy[i] > energy_bin_edges[bin_counter])
        {
            if (sum_counter == 0)
            {
                endf_xsec_hist->SetBinContent(bin_counter, 0);
            } else {
                endf_xsec_hist->SetBinContent(bin_counter, xsec_sum/sum_counter);
            }
            
            xsec_sum = xsec[i];
            sum_counter = 1;
            bin_counter++;
            i--;
            continue;
        }

        if (energy[i] < energy_bin_edges[bin_counter])
        {
            xsec_sum += xsec[i];
            sum_counter++;
        }
    }
}

void fillTSRunsHists(Int_t num_bins_e, Double_t bin_edges_e[]){    

    //Transmission - Method 1 //////////////////////////////////////////////////////////////
    TH1D* n_fOut_tIn_ts = 0;
    TH1D* n_fOut_tOut_ts = 0;

    n_fOut_tIn_ts = new TH1D("n_fOut_tIn_ts","Energy Hist Filter Out and Target In - with transmission setup",num_bins_e,bin_edges_e);
    Double_t norm_factor_fOut_tIn = fillEnergyHist(f_out_t_c_ts, n_fOut_tIn_ts);

    n_fOut_tOut_ts = new TH1D("n_fOut_tOut_ts","Energy Hist Filter Out and Target Out - with transmission setup",num_bins_e,bin_edges_e);
    Double_t norm_factor_fOut_tOut_ts = fillEnergyHist(f_out_t_out_ts, n_fOut_tOut_ts);

    Double_t norm_factor_method1 = (norm_factor_fOut_tOut_ts) / (norm_factor_fOut_tIn);

    //Transmission
    trans_hist_fOut = new TH1D("trans_hist_fOut","Transmission Hist - Filter Out method",num_bins_e,bin_edges_e);

    for (int i = 0; i < num_bins_e; i++)
    {
        Double_t bin_content_fOut_tIn = n_fOut_tIn_ts->GetBinContent(i+1);
        Double_t bin_content_fOut_tOut_ts = n_fOut_tOut_ts->GetBinContent(i+1);
        if (bin_content_fOut_tIn == 0. || bin_content_fOut_tOut_ts == 0.)
        {
            trans_hist_fOut->SetBinContent(i+1, 0.);
            trans_hist_fOut->SetBinError(i+1, 0.);
        } else {
            Double_t transmission = (bin_content_fOut_tIn * norm_factor_method1) / (bin_content_fOut_tOut_ts);
            Double_t bin_unc = transmission * std::sqrt( (1./bin_content_fOut_tIn) + (1./bin_content_fOut_tOut_ts) );
            trans_hist_fOut->SetBinContent(i+1, transmission);
            trans_hist_fOut->SetBinError(i+1, bin_unc);
        }
    }

    //Transmission - Method 2 //////////////////////////////////////////////////////////////
    TH1D* n_fIn_tIn_ts = 0;
    TH1D* n_fIn_tOut_ts_out = 0;
    TH1D* n_fOut_tOut_ts_out = 0;

    n_fIn_tIn_ts = new TH1D("n_fIn_tIn_ts","Energy Hist Filter In and Target In - with transmission setup",num_bins_e,bin_edges_e);
    Double_t norm_factor_fIn_tIn = fillEnergyHist(f_al5_t_c_ts, n_fIn_tIn_ts);

    n_fIn_tOut_ts_out = new TH1D("n_fIn_tOut_ts_out","Energy Hist Filter In and Target Out- without transmission setup",num_bins_e,bin_edges_e);
    Double_t norm_factor_fIn_tOut = fillEnergyHist(f_al5_t_out_runs, n_fIn_tOut_ts_out);

    n_fOut_tOut_ts_out = new TH1D("n_fOut_tOut_ts_out","Energy Hist Filter Out and Target Out - without transmission setup",num_bins_e,bin_edges_e);
    Double_t norm_factor_ts_out = fillEnergyHist(f_out_t_out_runs, n_fOut_tOut_ts_out);

    Double_t norm_factor_method2 = (norm_factor_fIn_tOut * norm_factor_fOut_tOut_ts) / (norm_factor_fIn_tIn * norm_factor_ts_out);

    //Transmission
    trans_hist_fIn = new TH1D("trans_hist_fIn","Transmission Hist - Filter In method",num_bins_e,bin_edges_e);
    for (int i = 0; i < num_bins_e; i++)
    {
        Double_t bin_content_fIn_tIn_ts = n_fIn_tIn_ts->GetBinContent(i+1);
        Double_t bin_content_fIn_tOut_ts_out = n_fIn_tOut_ts_out->GetBinContent(i+1);
        Double_t bin_content_fOut_tOut_ts_out = n_fOut_tOut_ts_out->GetBinContent(i+1);
        Double_t bin_content_fOut_tOut_ts = n_fOut_tOut_ts->GetBinContent(i+1);
        if (bin_content_fIn_tIn_ts == 0. || bin_content_fIn_tOut_ts_out == 0. || bin_content_fOut_tOut_ts_out == 0. || bin_content_fOut_tOut_ts == 0.)
        {
            trans_hist_fIn->SetBinContent(i+1, 0.);
            trans_hist_fIn->SetBinError(i+1, 0.);
        } else {
            Double_t transmission = (bin_content_fIn_tIn_ts * bin_content_fOut_tOut_ts_out * norm_factor_method2) / (bin_content_fIn_tOut_ts_out * bin_content_fOut_tOut_ts);
            Double_t bin_unc = transmission * std::sqrt( (1./bin_content_fIn_tIn_ts) + (1./bin_content_fOut_tOut_ts_out) + (1./bin_content_fIn_tOut_ts_out) + (1./bin_content_fOut_tOut_ts) );
            trans_hist_fIn->SetBinContent(i+1, transmission);
            trans_hist_fIn->SetBinError(i+1, bin_unc);
        }
    }

    //Transmission - Method 3 //////////////////////////////////////////////////////////////
    endf_xsec_hist = new TH1D("endf_xsec_hist","ENDF Cross Section Hist",num_bins_e,bin_edges_e);
    fillENDFHist(bin_edges_e);

    Double_t norm_factor_method3 = (norm_factor_ts_out) / (norm_factor_fOut_tIn);

    trans_hist_fOut_endf = new TH1D("trans_hist_fOut_endf","Transmission Hist - Filter Out ENDF method",num_bins_e,bin_edges_e);

    for (int i = 0; i < num_bins_e; i++)
    {
        Double_t bin_content_fOut_tIn = n_fOut_tIn_ts->GetBinContent(i+1);
        Double_t bin_content_fOut_tOut_ts_out = n_fOut_tOut_ts_out->GetBinContent(i+1);
        if (bin_content_fOut_tIn == 0. || bin_content_fOut_tOut_ts_out == 0.)
        {
            trans_hist_fOut_endf->SetBinContent(i+1, 0.);
            trans_hist_fOut_endf->SetBinError(i+1, 0.);
        } else {
            Double_t endf_factor = std::exp(- n_Al_0p2cm * endf_xsec_hist->GetBinContent(i+1));
            Double_t transmission = (bin_content_fOut_tIn * norm_factor_method3) / (bin_content_fOut_tOut_ts_out * endf_factor);
            Double_t bin_unc = transmission * std::sqrt( (1./bin_content_fOut_tIn) + (1./bin_content_fOut_tOut_ts_out) );
            trans_hist_fOut_endf->SetBinContent(i+1, transmission);
            trans_hist_fOut_endf->SetBinError(i+1, bin_unc);
        }
    }

    //Transmission - Method 4 //////////////////////////////////////////////////////////////
    fillENDFHist(bin_edges_e);

    Double_t norm_factor_method4 = (norm_factor_fIn_tOut) / (norm_factor_fIn_tIn);

    trans_hist_fIn_endf = new TH1D("trans_hist_fIn_endf","Transmission Hist - Filter In ENDF method",num_bins_e,bin_edges_e);

    for (int i = 0; i < num_bins_e; i++)
    {
        Double_t bin_content_fIn_tIn = n_fIn_tIn_ts->GetBinContent(i+1);
        Double_t bin_content_fIn_tOut_ts_out = n_fIn_tOut_ts_out->GetBinContent(i+1);
        if (bin_content_fIn_tIn == 0. || bin_content_fIn_tOut_ts_out == 0.)
        {
            trans_hist_fIn_endf->SetBinContent(i+1, 0.);
            trans_hist_fIn_endf->SetBinError(i+1, 0.);
        } else {
            Double_t endf_factor = std::exp(- n_Al_0p2cm * endf_xsec_hist->GetBinContent(i+1));
            Double_t transmission = (bin_content_fIn_tIn * norm_factor_method4) / (bin_content_fIn_tOut_ts_out * endf_factor);
            Double_t bin_unc = transmission * std::sqrt( (1./bin_content_fIn_tIn) + (1./bin_content_fIn_tOut_ts_out) );
            trans_hist_fIn_endf->SetBinContent(i+1, transmission);
            trans_hist_fIn_endf->SetBinError(i+1, bin_unc);
        }
    }

    //Combining the transmissions //////////////////////////////////////////////////////////////
    for (int i = 0; i < num_bins_e; i++)
    {
        Double_t bin_content_fIn = trans_hist_fIn->GetBinContent(i+1);
        Double_t bin_error_fIn = trans_hist_fIn->GetBinError(i+1);
        Double_t bin_content_fOut = trans_hist_fOut->GetBinContent(i+1);
        Double_t bin_error_fOut = trans_hist_fOut->GetBinError(i+1);
        Double_t transmission = (bin_content_fIn + bin_content_fOut) * 0.5;
        Double_t bin_unc = 0.5 * std::sqrt( (bin_error_fIn * bin_error_fIn) + (bin_error_fOut * bin_error_fOut) );
        transmission_hist_e->SetBinContent(i+1, transmission);
        transmission_hist_e->SetBinError(i+1, bin_unc);
    }

    //Cross Section Calculation
    Double_t num_density = num_density_map[target_name];
    cout << Form("Number density of %s = ", target_name.c_str()) << num_density << endl;
    Double_t n_inverse = ( (Double_t) 1.0/num_density);
    cout << Form("1/n of %s = ", target_name.c_str()) << n_inverse << endl;

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
            Double_t bin_unc = n_inverse * (1/trans_bin_content) * trans_bin_error;
            cross_section_hist_e->SetBinContent(i+1, cross_section);
            cross_section_hist_e->SetBinError(i+1, bin_unc);
        }
    }

    cout << "Total Protons F: Out, T: C (1.2 cm) - with transmission setup = " << norm_factor_fOut_tIn << endl;
    cout << "Total Protons F: Out, T: Out - with transmission setup = " << norm_factor_fOut_tOut_ts << endl;
    cout << "Total Protons F: Al (5 cm), T: C (1.2 cm) - with transmission setup = " << norm_factor_fIn_tIn << endl;
    cout << "Total Protons F: Al (5 cm), T: Out - without transmission setup = " << norm_factor_fIn_tOut << endl;
    cout << "Total Protons F: Out, T: Out - without transmission setup = " << norm_factor_ts_out << endl;
}

void fillNoFilterRunHists(Int_t num_bins_e, Double_t bin_edges_e[]){

    TH1D* energy_hist_target_in = 0;
    TH1D* energy_hist_target_out = 0;

    //Target In Hists
    energy_hist_target_in = new TH1D("energy_hist_target_in","Energy Hist Target In",num_bins_e,bin_edges_e);
    //Target Out Hists
    energy_hist_target_out = new TH1D("energy_hist_target_out","Energy Hist Target Out",num_bins_e,bin_edges_e);
    
    //Filling the Histograms
    Double_t norm_factor_target_in = fillEnergyHist(ts_target_in_runs, energy_hist_target_in);
    Double_t norm_factor_target_out = 0;
    if(!target_name.compare("ar_bottle_full")){
        norm_factor_target_out = fillEnergyHist(ts_target_out_runs, energy_hist_target_out);
    } else {
        norm_factor_target_out = fillEnergyHist(f_out_t_out_ts, energy_hist_target_out);
    }
    Double_t Qout_Qin = norm_factor_target_out / norm_factor_target_in;

    //Transmission
    for (int i = 0; i < num_bins_e; i++)
    {
        Double_t bin_content_in = energy_hist_target_in->GetBinContent(i+1);
        Double_t bin_content_out = energy_hist_target_out->GetBinContent(i+1);
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

    //Cross Section
    Double_t num_density = num_density_map[target_name];
    cout << Form("Number density of %s = ", target_name.c_str()) << num_density << endl;
    Double_t n_inverse = ( (Double_t) 1.0/num_density);
    cout << Form("1/n of %s = ", target_name.c_str()) << n_inverse << endl;

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

    cout << Form("Total Protons %s = ", target_name.c_str()) << norm_factor_target_in << endl;
    cout << "Total Protons Target Out = " << norm_factor_target_out << endl;

}

void ptbc_ts(){

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
    
    // fillTSRunsHists(num_bins_e,bin_edges_e);

    fillNoFilterRunHists(num_bins_e,bin_edges_e);

    //Writing to the output file
    outputRootFile = new TFile(Form("../rootFiles/crossSectionAna_%s.root", target_name.c_str()),"recreate");
    // trans_hist_fOut->Write();
    // trans_hist_fIn->Write();
    // trans_hist_fOut_endf->Write();
    // trans_hist_fIn_endf->Write();
    transmission_hist_e->Write();
    cross_section_hist_e->Write();

    outputRootFile->Close();
}