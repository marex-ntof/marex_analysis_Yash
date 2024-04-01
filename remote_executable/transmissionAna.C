/**
 * @file transmissionAna.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-11-13
 */

#include "transmissionAna.h"

void StoreHist(){

    TFile *f = new TFile("../rootFiles/transmissionAna_al8_al5.root","recreate");

    tof_hist_filter_al5->Write();
    energy_hist_filter_al5->Write();
    tof_hist_filter_al8->Write();
    energy_hist_filter_al8->Write();

    f->Close();
}

void fillRuns(){
    // if (mode == "test")
    // {
    //     filter_in_runs.push_back(117389);
    //     filter_out_runs.push_back(117367);
    // }
    // Bi 1 cm filter
    // filter_in_runs.reserve( c_au_f_bi_t_out.size() + c_out_f_bi_t_out.size() + c_pb_f_bi_t_out.size() + c_ta_f_bi_t_out.size() );
    // filter_in_runs.insert( filter_in_runs.end(), c_au_f_bi_t_out.begin(), c_au_f_bi_t_out.end() );
    // filter_in_runs.insert( filter_in_runs.end(), c_out_f_bi_t_out.begin(), c_out_f_bi_t_out.end() );
    // filter_in_runs.insert( filter_in_runs.end(), c_pb_f_bi_t_out.begin(), c_pb_f_bi_t_out.end() );
    // filter_in_runs.insert( filter_in_runs.end(), c_ta_f_bi_t_out.begin(), c_ta_f_bi_t_out.end() );

    // Al 8 cm filter 
    filter_al8_runs.reserve( c_ta_f_al8_t_out.size() );
    filter_al8_runs.insert( filter_al8_runs.end(), c_ta_f_al8_t_out.begin(), c_ta_f_al8_t_out.end() );

    // Al 5 cm filter
    filter_al5_runs.reserve( c_ta_f_al5_t_out.size() + c_out_f_al5_t_out.size() );
    filter_al5_runs.insert( filter_al5_runs.end(), c_ta_f_al5_t_out.begin(), c_ta_f_al5_t_out.end() );
    filter_al5_runs.insert( filter_al5_runs.end(), c_out_f_al5_t_out.begin(), c_out_f_al5_t_out.end() );

    // filter_out_runs.reserve( c_au_f_out_t_out.size() + c_out_f_out_t_out.size() + c_pb_f_out_t_out.size() + c_ta_f_out_t_out.size() );
    // filter_out_runs.insert( filter_out_runs.end(), c_au_f_out_t_out.begin(), c_au_f_out_t_out.end() );
    // filter_out_runs.insert( filter_out_runs.end(), c_out_f_out_t_out.begin(), c_out_f_out_t_out.end() );
    // filter_out_runs.insert( filter_out_runs.end(), c_pb_f_out_t_out.begin(), c_pb_f_out_t_out.end() );
    // filter_out_runs.insert( filter_out_runs.end(), c_ta_f_out_t_out.begin(), c_ta_f_out_t_out.end() );
}

Double_t FilterAl5(){

    Double_t NormFactor = 0; //Integral of the pulse intensity

    for (int i = 0; i < filter_al5_runs.size(); i++)
    {
        TTree* PTBC;

        Double_t tof = 0; //tof is in ns
        Double_t tflash = 0; //tflash is in ns
        Float_t amp = 0;
        Int_t BunchNumber = 0;
        Float_t PulseIntensity = 0;

        TFile *file_ntof = TFile::Open(Form("/eos/experiment/ntof/data/rootfiles/2023/ear1/run%d.root", filter_al5_runs.at(i)),"read");
        file_ntof->GetObject("PTBC", PTBC);

        PTBC->SetBranchAddress("BunchNumber", &BunchNumber);
        PTBC->SetBranchAddress("PulseIntensity", &PulseIntensity);
        PTBC->SetBranchAddress("tof", &tof);
        PTBC->SetBranchAddress("tflash", &tflash);
        PTBC->SetBranchAddress("amp", &amp); 

        Long64_t Events = PTBC->GetEntriesFast();
        std::cout << "Number of entries - Filter In = " << Events << std::endl;
        
        int CurrentBunchNum = 0;

        for (int i = 0; i < Events; i++)
        {
            PTBC->GetEntry(i);

            Double_t corrected_tof = (tof - tflash + t_gamma_PTB);

            if (CurrentBunchNum != BunchNumber)
            {
                CurrentBunchNum = BunchNumber;
                NormFactor += (Double_t) PulseIntensity;
            }
            
            //Line 1
            if (corrected_tof >= t11 && corrected_tof < t12)
            {
                if ( (Double_t) amp > yOnTheCutLine(t11, (Double_t) a11, t12, (Double_t) a12, corrected_tof) )
                {
                    tof_hist_filter_al5->Fill(corrected_tof);
                    energy_hist_filter_al5->Fill( TOFToEnergy(corrected_tof * 1e-9) );
                }
            }
            
            //Line 2
            if (corrected_tof >= t21 && corrected_tof < t22)
            {
                if ( (Double_t) amp > yOnTheCutLine(t21, (Double_t) a21, t22, (Double_t) a22, corrected_tof) )
                {
                    tof_hist_filter_al5->Fill(corrected_tof);
                    energy_hist_filter_al5->Fill( TOFToEnergy(corrected_tof * 1e-9) );
                }
            }

            //Line 3
            if (corrected_tof >= t31 && corrected_tof < t32)
            {
                if ( (Double_t) amp > yOnTheCutLine(t31, (Double_t) a31, t32, (Double_t) a32, corrected_tof) )
                {
                    tof_hist_filter_al5->Fill(corrected_tof);
                    energy_hist_filter_al5->Fill( TOFToEnergy(corrected_tof * 1e-9) );
                }
            }

            //Line 4
            if (corrected_tof >= t41 && corrected_tof < t42)
            {
                if ( (Double_t) amp > yOnTheCutLine(t41, (Double_t) a41, t42, (Double_t) a42, corrected_tof) )
                {
                    tof_hist_filter_al5->Fill(corrected_tof);
                    energy_hist_filter_al5->Fill( TOFToEnergy(corrected_tof * 1e-9) );
                }
            }

            //Line 5
            if (corrected_tof >= t51 && corrected_tof < t52)
            {
                if ( (Double_t) amp > yOnTheCutLine(t51, (Double_t) a51, t52, (Double_t) a52, corrected_tof) )
                {
                    tof_hist_filter_al5->Fill(corrected_tof);
                    energy_hist_filter_al5->Fill( TOFToEnergy(corrected_tof * 1e-9) );
                }
            }

            //Line 6
            if (corrected_tof >= t61 && corrected_tof < t62)
            {
                if ( (Double_t) amp > yOnTheCutLine(t61, (Double_t) a61, t62, (Double_t) a62, corrected_tof) ) 
                {
                    tof_hist_filter_al5->Fill(corrected_tof);
                    energy_hist_filter_al5->Fill( TOFToEnergy(corrected_tof * 1e-9) );
                }
            }
        }

        std::cout << "Integrated pulse intensity - Filter In = " << NormFactor << std::endl;
    }

    return NormFactor;

    // tof_hist_filter_al5->Scale(1.0/NormFactor);
    // energy_hist_filter_al5->Scale(1.0/NormFactor);
}

Double_t FilterAl8(){

    Double_t NormFactor = 0; //Integral of the pulse intensity

    for (int i = 0; i < filter_al8_runs.size(); i++)
    {
        TTree* PTBC;

        Double_t tof = 0; //tof is in ns
        Double_t tflash = 0; //tflash is in ns
        Float_t amp = 0;
        Int_t BunchNumber = 0;
        Float_t PulseIntensity = 0;

        TFile *file_ntof = TFile::Open(Form("/eos/experiment/ntof/data/rootfiles/2023/ear1/run%d.root", filter_al8_runs.at(i)),"read");

        file_ntof->GetObject("PTBC", PTBC);

        PTBC->SetBranchAddress("BunchNumber", &BunchNumber);
        PTBC->SetBranchAddress("PulseIntensity", &PulseIntensity);
        PTBC->SetBranchAddress("tof", &tof);
        PTBC->SetBranchAddress("tflash", &tflash);
        PTBC->SetBranchAddress("amp", &amp);

        Long64_t Events = PTBC->GetEntriesFast();
        std::cout << "Number of entries - Filter Out = " << Events << std::endl;
        
        int CurrentBunchNum = 0;

        for (int i = 0; i < Events; i++)
        {
            PTBC->GetEntry(i);

            Double_t corrected_tof = (tof - tflash + t_gamma_PTB);

            if (CurrentBunchNum != BunchNumber)
            {
                CurrentBunchNum = BunchNumber;
                NormFactor += (Double_t) PulseIntensity;
            }

            //Line 1
            if (corrected_tof >= t11 && corrected_tof < t12)
            {
                if ( (Double_t) amp > yOnTheCutLine(t11, (Double_t) a11, t12, (Double_t) a12, corrected_tof) )
                {
                    tof_hist_filter_al8->Fill(corrected_tof);
                    energy_hist_filter_al8->Fill( TOFToEnergy(corrected_tof * 1e-9) );
                }
            }
            
            //Line 2
            if (corrected_tof >= t21 && corrected_tof < t22)
            {
                if ( (Double_t) amp > yOnTheCutLine(t21, (Double_t) a21, t22, (Double_t) a22, corrected_tof) )
                {
                    tof_hist_filter_al8->Fill(corrected_tof);
                    energy_hist_filter_al8->Fill( TOFToEnergy(corrected_tof * 1e-9) );
                }
            }

            //Line 3
            if (corrected_tof >= t31 && corrected_tof < t32)
            {
                if ( (Double_t) amp > yOnTheCutLine(t31, (Double_t) a31, t32, (Double_t) a32, corrected_tof) )
                {
                    tof_hist_filter_al8->Fill(corrected_tof);
                    energy_hist_filter_al8->Fill( TOFToEnergy(corrected_tof * 1e-9) );
                }
            }

            //Line 4
            if (corrected_tof >= t41 && corrected_tof < t42)
            {
                if ( (Double_t) amp > yOnTheCutLine(t41, (Double_t) a41, t42, (Double_t) a42, corrected_tof) )
                {
                    tof_hist_filter_al8->Fill(corrected_tof);
                    energy_hist_filter_al8->Fill( TOFToEnergy(corrected_tof * 1e-9) );
                }
            }

            //Line 5
            if (corrected_tof >= t51 && corrected_tof < t52)
            {
                if ( (Double_t) amp > yOnTheCutLine(t51, (Double_t) a51, t52, (Double_t) a52, corrected_tof) )
                {
                    tof_hist_filter_al8->Fill(corrected_tof);
                    energy_hist_filter_al8->Fill( TOFToEnergy(corrected_tof * 1e-9) );
                }
            }

            //Line 6
            if (corrected_tof >= t61 && corrected_tof < t62)
            {
                if ( (Double_t) amp > yOnTheCutLine(t61, (Double_t) a61, t62, (Double_t) a62, corrected_tof) ) 
                {
                    tof_hist_filter_al8->Fill(corrected_tof);
                    energy_hist_filter_al8->Fill( TOFToEnergy(corrected_tof * 1e-9) );
                }
            }
        }

        std::cout << "Integrated pulse intensity - Filter In = " << NormFactor << std::endl;
    }

    return NormFactor;

    // tof_hist_filter_al8->Scale(1.0/NormFactor);
    // energy_hist_filter_al8->Scale(1.0/NormFactor);
}

void endf(Double_t n, Double_t energy_bin_edges[], bool fillENDF, bool fillENDFSmeared){
    //Extracting ENDF Bi Cross Section and transmission
    //number density - Bismuth Filter

    std::ifstream inputFile("Al_tot_xsec.txt"); //Bi_tot_xsec.txt

    TFile rfFile("RF.root");
    TH2D *rf_hist = (TH2D*)rfFile.Get("histfluka");

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
            line_count++;
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
    Double_t trans_sum = 0;
    Int_t sum_counter = 0;
    Int_t bin_counter = 1;
    
    Double_t xsec_sum_rf = 0;
    Double_t trans_sum_rf = 0;
    Int_t sum_counter_rf = 0;
    Int_t bin_counter_rf = 1;

    for (int i = 0; i < energy.size(); i++)
    {
        Double_t new_e = 0;
        if (fillENDFSmeared == true)
        {
            //Convoluting endf with rf
            Int_t e_bin_num = rf_hist->GetXaxis()->FindBin(energy[i]);
            std::string projection_name = "profile_" + std::to_string(energy[i]);
            TH1D* projection = (TH1D*)rf_hist->ProjectionY(
                projection_name.c_str(),
                e_bin_num, e_bin_num
            );
            // Double_t fwhm = FindFWHM(projection); //in cm
            Double_t rf_length = projection->GetMean(1) * 0.01; //projection->GetBinCenter( projection->GetMaximumBin() ) * 0.01; //converting to m
            Double_t e_tof = EnergyToTOF(energy[i]);
            new_e = TOFToEnergy(e_tof, rf_length); //in eV
        }
        
        if (fillENDF == true && energy[i] > energy_bin_edges[bin_counter])
        {
            if (sum_counter == 0)
            {
                endf_xsec_hist->SetBinContent(bin_counter, 0);
                endf_trans_hist->SetBinContent(bin_counter, 0);
            } else {
                endf_xsec_hist->SetBinContent(bin_counter, xsec_sum/sum_counter);
                endf_trans_hist->SetBinContent(bin_counter, trans_sum/sum_counter);
            }
            
            xsec_sum = xsec[i];
            trans_sum = std::exp(- n * xsec[i]);
            sum_counter = 1;
            bin_counter++;
            i--;
            continue;
        }

        if (fillENDFSmeared == true && new_e > energy_bin_edges[bin_counter_rf])
        {
            if (sum_counter_rf == 0)
            {
                endf_rf_xsec_hist->SetBinContent(bin_counter_rf, 0);
                endf_rf_trans_hist->SetBinContent(bin_counter_rf, 0);
            } else {
                endf_rf_xsec_hist->SetBinContent(bin_counter_rf, xsec_sum_rf/sum_counter_rf);
                endf_rf_trans_hist->SetBinContent(bin_counter_rf, trans_sum_rf/sum_counter_rf);
            }
            
            xsec_sum_rf = xsec[i];
            trans_sum_rf = std::exp(- n * xsec[i]);
            sum_counter_rf = 1;
            bin_counter_rf++;
            i--;
            continue;
        }

        if (fillENDF == true && energy[i] < energy_bin_edges[bin_counter])
        {
            xsec_sum += xsec[i];
            trans_sum += std::exp(- n * xsec[i]);
            sum_counter++;
        }

        if (fillENDFSmeared == true && new_e < energy_bin_edges[bin_counter])
        {
            xsec_sum_rf += xsec[i];
            trans_sum_rf += std::exp(- n * xsec[i]);
            sum_counter_rf++;
        }
  
        // //Filling the endf
        // endf_trans->AddPoint(energy[i], std::exp(- n * val2));
        // endf_xsec->AddPoint(energy[i], val2);

        // endf_rf_trans->AddPoint(new_e, std::exp(- n * val2));
        // endf_rf_xsec->AddPoint(new_e, val2);
    }
}

void plots(){
    
    gStyle->SetStatX(0.27);
    gStyle->SetStatY(0.9);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.17);

    TCanvas *c[4];
    TLegend *l[4];

    int i = 0;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();

    l[i] = new TLegend(0.72,0.8,0.90,0.9);
    l[i]->AddEntry(energy_hist_filter_al8,"Al (8 cm)","l");
    l[i]->AddEntry(energy_hist_filter_al5,"Al (5 cm)","l");
    // l[i]->AddEntry(endf_rf_trans_hist,"ENDF smeared","l");

    energy_hist_filter_al8->GetXaxis()->SetTitle("Energy (in eV)");
    energy_hist_filter_al8->GetYaxis()->SetTitle("Counts / total protons");
    energy_hist_filter_al8->SetTitle("Normalized Events - Al (8 cm) and Al (5 cm)");
    energy_hist_filter_al8->GetXaxis()->SetRange(2e4,3e5);
    energy_hist_filter_al8->Draw(); //"HISTE"
    energy_hist_filter_al8->SetStats(0);
    // energy_hist_filter_al8->SetMarkerStyle(6);
    // energy_hist_filter_al8->SetMarkerSize(0.5);
    gPad->SetGrid();
    gPad->SetLogx();
    // gStyle->SetPalette(57);

    energy_hist_filter_al5->SetLineColor(2);
    // energy_hist_filter_al5->SetLineWidth(4);
    energy_hist_filter_al5->Draw("SAME");

    // endf_rf_trans_hist->SetLineColor(8);
    // // endf_rf_trans_hist->SetLineWidth(4);
    // endf_rf_trans_hist->Draw("SAME");
    l[i]->Draw();

    c[i]->Print("/afs/cern.ch/user/y/yabezawa/nTOF_data_analysis/plots/norm_events_e_al8_al5.png");

    i++;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();

    l[i] = new TLegend(0.72,0.8,0.90,0.9);
    l[i]->AddEntry(tof_hist_filter_al8,"Al (8 cm)","l");
    l[i]->AddEntry(tof_hist_filter_al5,"Al (5 cm)","l");
    // l[i]->AddEntry(endf_rf_trans_hist,"ENDF smeared","l");

    tof_hist_filter_al8->GetXaxis()->SetTitle("Time of Flight (in ns)");
    tof_hist_filter_al8->GetYaxis()->SetTitle("Counts / total protons");
    tof_hist_filter_al8->SetTitle("Normalized Events - Al (8 cm) and Al (5 cm)");
    tof_hist_filter_al8->GetXaxis()->SetRange(24061,93168);
    tof_hist_filter_al8->Draw(); //"HISTE"
    tof_hist_filter_al8->SetStats(0);
    // tof_hist_filter_al8->SetMarkerStyle(6);
    // tof_hist_filter_al8->SetMarkerSize(0.5);
    gPad->SetGrid();
    gPad->SetLogx();
    // gStyle->SetPalette(57);

    tof_hist_filter_al5->SetLineColor(2);
    // tof_hist_filter_al5->SetLineWidth(4);
    tof_hist_filter_al5->Draw("SAME");

    // endf_rf_trans_hist->SetLineColor(8);
    // // endf_rf_trans_hist->SetLineWidth(4);
    // endf_rf_trans_hist->Draw("SAME");
    l[i]->Draw();

    c[i]->Print("/afs/cern.ch/user/y/yabezawa/nTOF_data_analysis/plots/norm_events_tof_al8_al5.png");

}

void transmissionAna(){

    fillRuns();

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

    //Filter In Hists
    tof_hist_filter_al5 = new TH1D("tof_hist_filter_al5","Time of Flight Hist Al 5 cm",num_bins_tof,bin_edges_tof);
    energy_hist_filter_al5 = new TH1D("energy_hist_filter_al5","Energy Hist Al 5 cm",num_bins_e,bin_edges_e);
    //Filter Out Hists
    tof_hist_filter_al8 = new TH1D("tof_hist_filter_al8","Time of Flight Hist Al 8 cm",num_bins_tof,bin_edges_tof);
    energy_hist_filter_al8 = new TH1D("energy_hist_filter_al8","Energy Hist Al 8 cm",num_bins_e,bin_edges_e);
    //ENDF Hists
    // endf_trans_hist = new TH1D("endf_trans_hist","ENDF Transmission Hist",num_bins_e,bin_edges_e);
    // endf_xsec_hist = new TH1D("endf_xsec_hist","ENDF Cross Section Hist",num_bins_e,bin_edges_e);
    // endf_rf_trans_hist = new TH1D("endf_rf_trans_hist","ENDF Transmission Hist - RF Convoluted",num_bins_e,bin_edges_e);
    // endf_rf_xsec_hist = new TH1D("endf_rf_xsec_hist","ENDF Cross Section Hist - RF Convoluted",num_bins_e,bin_edges_e);

    //Filling the Histograms
    Double_t norm_factor_filter_al5 = FilterAl5();
    Double_t norm_factor_filter_al8 = FilterAl8();

    //Normalizing the Histograms
    tof_hist_filter_al5->Scale(1.0/norm_factor_filter_al5);
    energy_hist_filter_al5->Scale(1.0/norm_factor_filter_al5);
    tof_hist_filter_al8->Scale(1.0/norm_factor_filter_al8);
    energy_hist_filter_al8->Scale(1.0/norm_factor_filter_al8);

    cout << "Total Protons Al 5 cm = " << norm_factor_filter_al5 << endl;
    cout << "Total Protons Al 8 cm = " << norm_factor_filter_al8 << endl;

    //endf(n_Al_8cm, bin_edges_e, true, false);
    StoreHist();
    plots();
}