/**
 * @file crossSectionPlots.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-11-13
 */

#include "crossSectionPlots.h"

void endf(Double_t n, Double_t energy_bin_edges[], bool fillENDF, bool fillENDFSmeared){
    //Extracting ENDF Cross Section and transmission

    fillMaxLineCount();

    std::ifstream inputFile(Form("../evalData/%s", eval_file_name_map[filter_name].c_str())); // endf_file_name.c_str()

    if (fillENDFSmeared == true){
        TFile *rfFile = TFile::Open("../inputFiles/RF.root", "READ");
        rf_hist = (TH2D*)rfFile->Get("histfluka");
    }

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

        if (line_count == evalFile_max_line_num) //Bi - 27019; Al - 9745
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
            Double_t e_tof = EnergyToTOF(energy[i], flight_path_length_PTB);
            new_e = TOFToEnergy(e_tof, flight_path_length_PTB, rf_length); //in eV
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
    }
}

void jendl(Double_t n, Double_t energy_bin_edges[]){
    //Extracting ENDF Cross Section and transmission

    std::ifstream inputFile("../evalData/JENDL_Ar_tot_xsec.txt");

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

        if (line_count == 51631)
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

    for (int i = 0; i < energy.size(); i++)
    {
        
        if (energy[i] > energy_bin_edges[bin_counter])
        {
            if (sum_counter == 0)
            {
                jendl_xsec_hist->SetBinContent(bin_counter, 0);
                jendl_trans_hist->SetBinContent(bin_counter, 0);
            } else {
                jendl_xsec_hist->SetBinContent(bin_counter, xsec_sum/sum_counter);
                jendl_trans_hist->SetBinContent(bin_counter, trans_sum/sum_counter);
            }
            
            xsec_sum = xsec[i];
            trans_sum = std::exp(- n * xsec[i]);
            sum_counter = 1;
            bin_counter++;
            i--;
            continue;
        }

        if (energy[i] < energy_bin_edges[bin_counter])
        {
            xsec_sum += xsec[i];
            trans_sum += std::exp(- n * xsec[i]);
            sum_counter++;
        }
    }
}

TH1D* retriveHistogramsChangeBPD(const char *fname, const char *hist_name, Int_t bpd_old, Int_t bpd_new){

    // if(bpd_old < bpd_new) {

    // }

    cout << "Opening File " << fname << endl;
    TFile *hist_file = TFile::Open(fname, "READ");
    // if (!hist_file || hist_file->IsZombie()) {
    //     cout << "Unable to open " << fname << " for reading..." << endl;
    //     return;
    // }

    // for (auto&& keyAsObj : *hist_file->GetListOfKeys()){
    //     auto key = (TKey*) keyAsObj;
    //     cout << key->GetName() << " " << key->GetClassName() << endl;
    // }

    cout << "Extracting Histogram " << hist_name << endl;
    TH1D* hist_old = (TH1D*)hist_file->Get(hist_name);
    Int_t num_bins_old = hist_old->GetNbinsX();
    // cout << "Number of old bins = " << num_bins_old << endl;
    Int_t max_sum_bin_count = (bpd_old/bpd_new); // Number of old bins that need to be summed to make a new bin
    Int_t num_bins_new = num_bins_old / max_sum_bin_count;
    // cout << "Number of new bins = " << num_bins_new << endl;

    Double_t bin_edges_new[num_bins_new + 1];
    Double_t bin_content[num_bins_new];
    Double_t bin_error[num_bins_new];
    Int_t sum_bin_count = 0;
    Int_t new_bin_counter = 0;
    for (int i = 1; i < num_bins_old + 1; i++)
    {
        if (i == 1)
        {
            bin_edges_new[new_bin_counter] = hist_old->GetXaxis()->GetBinLowEdge(i);
            bin_content[new_bin_counter] = 0;
            bin_error[new_bin_counter] = 0;
        }
        bin_content[new_bin_counter] += hist_old->GetBinContent(i);
        bin_error[new_bin_counter] += hist_old->GetBinError(i) * hist_old->GetBinError(i);
        sum_bin_count++;

        if (i == num_bins_old)
        {
            // sum_bin_count = 0;
            bin_error[new_bin_counter] = sqrt(bin_error[new_bin_counter])/max_sum_bin_count;
            bin_content[new_bin_counter] = bin_content[new_bin_counter]/max_sum_bin_count;
            bin_edges_new[new_bin_counter+1] = hist_old->GetXaxis()->GetBinUpEdge(i);
        }
        else if (sum_bin_count == max_sum_bin_count)
        {
            sum_bin_count = 0;
            bin_error[new_bin_counter] = sqrt(bin_error[new_bin_counter])/max_sum_bin_count;
            bin_content[new_bin_counter] = bin_content[new_bin_counter]/max_sum_bin_count;
            new_bin_counter++;
            bin_edges_new[new_bin_counter] = hist_old->GetXaxis()->GetBinUpEdge(i);
            bin_content[new_bin_counter] = 0;
            bin_error[new_bin_counter] = 0;
        } 
    }

    auto hist_new = new TH1D(hist_name,hist_old->GetTitle(),num_bins_new,bin_edges_new);
    for(int i = 0; i < num_bins_new; i++){
        hist_new->SetBinContent(i+1, bin_content[i]);
        hist_new->SetBinError(i+1, bin_error[i]);
    }

    return hist_new;
}

void retriveHistograms(const char *fname){
    
    TFile* hist_file = TFile::Open(fname, "READ");
    if (!hist_file || hist_file->IsZombie()) {
        cout << "Unable to open " << fname << " for reading..." <<endl;
        return;
    }

    // tof_hist_filter_in = (TH1D*)hist_file->Get("tof_hist_filter_in");
    // energy_hist_filter_in = (TH1D*)hist_file->Get("energy_hist_filter_in");
    // tof_hist_filter_out = (TH1D*)hist_file->Get("tof_hist_filter_out");
    // energy_hist_filter_out = (TH1D*)hist_file->Get("energy_hist_filter_out");

    // trans_hist_fOut = (TH1D*)hist_file->Get("trans_hist_fOut");
    // trans_hist_fIn = (TH1D*)hist_file->Get("trans_hist_fIn");
    // trans_hist_fOut_endf = (TH1D*)hist_file->Get("trans_hist_fOut_endf");
    // trans_hist_fIn_endf = (TH1D*)hist_file->Get("trans_hist_fIn_endf");
    transmission_hist_e_PTB = (TH1D*)hist_file->Get("transmission_hist_e_PTB");
    transmission_hist_e_FIMG = (TH1D*)hist_file->Get("transmission_hist_e_FIMG");
    cross_section_hist_e_PTB = (TH1D*)hist_file->Get("cross_section_hist_e_PTB");
    cross_section_hist_e_FIMG = (TH1D*)hist_file->Get("cross_section_hist_e_FIMG");

    // hist_file->Close();
}

void crossSectionPlots(){

    fillNumDensityMap();
    fillEValFileNameMap();

    // retriveHistograms(Form("../rootFiles/crossSectionAna_%s.root", filter_name.c_str())); //root_file_name.c_str()
    transmission_hist_e_PTB = retriveHistogramsChangeBPD(Form("../rootFiles/crossSectionAna_%s.root", filter_name.c_str()),"transmission_hist_e_PTB", 1000, bins_per_decade);
    transmission_hist_e_FIMG = retriveHistogramsChangeBPD(Form("../rootFiles/crossSectionAna_%s.root", filter_name.c_str()),"transmission_hist_e_FIMG", 1000, bins_per_decade);
    cross_section_hist_e_PTB = retriveHistogramsChangeBPD(Form("../rootFiles/crossSectionAna_%s.root", filter_name.c_str()),"cross_section_hist_e_PTB", 1000, bins_per_decade);
    cross_section_hist_e_FIMG = retriveHistogramsChangeBPD(Form("../rootFiles/crossSectionAna_%s.root", filter_name.c_str()),"cross_section_hist_e_FIMG", 1000, bins_per_decade);

    bool fillENDF = true;
    bool fillENDFSmeared = false;
    bool fillJENDL = true;

    //Getting energy bin edges
    Int_t num_bins_e = transmission_hist_e_PTB->GetNbinsX();
    Double_t bin_edges_e[num_bins_e + 1];
    for (int i = 0; i < num_bins_e; i++)
    {
        bin_edges_e[i] = transmission_hist_e_PTB->GetXaxis()->GetBinLowEdge(i+1);

        if (i == num_bins_e - 1)
        {
            bin_edges_e[i+1] = transmission_hist_e_PTB->GetXaxis()->GetBinUpEdge(i+1);
        }
    }

    //ENDF Hists
    if (fillENDF){
        endf_trans_hist = new TH1D("endf_trans_hist","ENDF Transmission Hist",num_bins_e,bin_edges_e);
        endf_xsec_hist = new TH1D("endf_xsec_hist","ENDF Cross Section Hist",num_bins_e,bin_edges_e);  
    }

    if (fillENDFSmeared){
        endf_rf_trans_hist = new TH1D("endf_rf_trans_hist","ENDF Transmission Hist - RF Convoluted",num_bins_e,bin_edges_e);
        endf_rf_xsec_hist = new TH1D("endf_rf_xsec_hist","ENDF Cross Section Hist - RF Convoluted",num_bins_e,bin_edges_e);
    }

    if(fillJENDL){
        jendl_trans_hist = new TH1D("jendl_trans_hist","JENDL Transmission Hist",num_bins_e,bin_edges_e);
        jendl_xsec_hist = new TH1D("jendl_xsec_hist","JENDL Cross Section Hist",num_bins_e,bin_edges_e);  
    }

    //Filling ENDF hist
    Double_t num_density = num_density_map[filter_name];
    endf(num_density, bin_edges_e, fillENDF, fillENDFSmeared);

    if(fillJENDL)
    {
        jendl(num_density, bin_edges_e);
    }
    
    
    //Plotting
    SetMArEXStyle();

    gStyle->SetStatX(0.27);
    gStyle->SetStatY(0.9);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.17);

    TCanvas *c[4];
    TLegend *l[4];

    int i = 0;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();

    l[i] = new TLegend(0.77,0.7,0.86,0.8); //0.68,0.7,0.86,0.8       ;         0.72,0.8,0.90,0.9
    l[i]->AddEntry(transmission_hist_e_PTB,"PTBC","l");
    
    transmission_hist_e_PTB->GetXaxis()->SetTitle("Energy (in eV)");
    transmission_hist_e_PTB->GetYaxis()->SetTitle("Transmission");
    transmission_hist_e_PTB->SetTitle(Form("Transmission Histogram - %s", filter_name_title.c_str()));
    transmission_hist_e_PTB->SetLineWidth(2);
    transmission_hist_e_PTB->Draw(); //"HISTE"
    transmission_hist_e_PTB->SetStats(0);
    // transmission_hist_e_PTB->SetMarkerStyle(6);
    // transmission_hist_e_PTB->SetMarkerSize(0.5);
    // gPad->SetGrid();
    gPad->SetLogx();
    // gStyle->SetPalette(57);

    l[i]->AddEntry(transmission_hist_e_FIMG,"FIMG","l");
    transmission_hist_e_FIMG->SetLineColor(3);
    transmission_hist_e_FIMG->SetLineWidth(2);
    transmission_hist_e_FIMG->GetXaxis()->SetRangeUser(1e-2,1e3);
    transmission_hist_e_FIMG->Draw("SAME");

    // l[i]->AddEntry(trans_hist_fOut,"No Al5","l");
    // trans_hist_fOut->SetLineColor(2);
    // trans_hist_fOut->Draw("SAME");

    // l[i]->AddEntry(trans_hist_fIn,"With Al5","l");
    // trans_hist_fIn->SetLineColor(1);
    // trans_hist_fIn->Draw("SAME");

    // l[i]->AddEntry(trans_hist_fOut_endf,"No Al5 with ENDF","l");
    // trans_hist_fOut_endf->SetLineColor(6);
    // trans_hist_fOut_endf->Draw("SAME");

    // l[i]->AddEntry(trans_hist_fIn_endf,"With Al5 and ENDF","l");
    // trans_hist_fIn_endf->SetLineColor(7);
    // trans_hist_fIn_endf->Draw("SAME");

    if (fillENDF){
        l[i]->AddEntry(endf_trans_hist,"ENDF","l");
        endf_trans_hist->SetLineColor(2);
        endf_trans_hist->SetLineWidth(2);
        // endf_trans_hist->GetXaxis()->SetRange(1e-2,2e7);
        endf_trans_hist->Draw("SAME");
    }

    if (fillENDFSmeared){
        l[i]->AddEntry(endf_rf_trans_hist,"ENDF smeared","l");
        endf_rf_trans_hist->SetLineColor(7);
        endf_rf_trans_hist->SetLineWidth(2);
        endf_rf_trans_hist->Draw("SAME");
    }

    if (fillJENDL)
    {
        l[i]->AddEntry(jendl_trans_hist,"JENDL-5","l");
        jendl_trans_hist->SetLineColor(1);
        jendl_trans_hist->SetLineWidth(2);
        // jendl_trans_hist->GetXaxis()->SetRange(1e-2,2e7);
        jendl_trans_hist->Draw("SAME");
    }
    

    l[i]->SetMargin(0.4);
    l[i]->Draw();
    // c[i]->Print(Form("../plots/h_trans_e_%s_%iBPD.png", filter_name.c_str(), bins_per_decade));

    i++;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();

    l[i] = new TLegend(0.80,0.7,0.86,0.8);
    l[i]->AddEntry(cross_section_hist_e_PTB,"PTBC","l");

    cross_section_hist_e_PTB->GetXaxis()->SetTitle("Energy (in eV)");
    cross_section_hist_e_PTB->GetYaxis()->SetTitle("Cross Section (in barns)");
    cross_section_hist_e_PTB->SetTitle(Form("Cross Section Histogram - %s", filter_name_title.c_str()));
    cross_section_hist_e_PTB->SetLineWidth(2);
    cross_section_hist_e_PTB->Draw(); //"HISTE"
    cross_section_hist_e_PTB->SetStats(0);
    // cross_section_hist_e_PTB->SetMarkerStyle(6);
    // cross_section_hist_e_PTB->SetMarkerSize(0.5);
    // gPad->SetGrid();
    gPad->SetLogx();
    // gStyle->SetPalette(57);

    l[i]->AddEntry(cross_section_hist_e_FIMG,"FIMG","l");
    cross_section_hist_e_FIMG->SetLineColor(3);
    cross_section_hist_e_FIMG->SetLineWidth(2);
    cross_section_hist_e_FIMG->Draw("SAME");

    if (fillENDF){
        l[i]->AddEntry(endf_xsec_hist,"ENDF","l");
        endf_xsec_hist->SetLineColor(2);
        endf_xsec_hist->SetLineWidth(2);
        // endf_xsec_hist->GetXaxis()->SetRange(1e-2,2e7);
        endf_xsec_hist->Draw("SAME");
    }
    
    if (fillENDFSmeared){
        l[i]->AddEntry(endf_rf_xsec_hist,"ENDF smeared","l");
        endf_rf_xsec_hist->SetLineColor(7);
        endf_rf_xsec_hist->SetLineWidth(2);
        endf_rf_xsec_hist->Draw("SAME");
    }

    l[i]->SetMargin(0.4);
    l[i]->Draw();
    c[i]->Print(Form("../plots/h_xsec_e_%s_%iBPD.png", filter_name.c_str(), bins_per_decade));

    // i++;

    // c[i] = new TCanvas(Form("c%d", i)," ");
    // c[i]->cd();

    // l[i] = new TLegend(0.72,0.8,0.90,0.9);
    // l[i]->AddEntry(tof_hist_filter_in,"Filter In","l");
    // l[i]->AddEntry(tof_hist_filter_out,"Filter Out","l");

    // tof_hist_filter_in->GetXaxis()->SetTitle("Time of Flight (in ns)");
    // tof_hist_filter_in->GetYaxis()->SetTitle("Counts (per proton)");
    // tof_hist_filter_in->SetTitle("Time of Flight Histograms");
    // tof_hist_filter_in->Draw(); //"HISTE"
    // gPad->SetGrid();
    // gPad->SetLogx();
    // gPad->SetLogy();
    // tof_hist_filter_out->SetLineColor(2);
    // tof_hist_filter_out->Draw("HISTESAME");
    // l[i]->Draw();

    // c[i]->Print("../plots/tof_hist.png");

    // i++;

    // c[i] = new TCanvas(Form("c%d", i)," ");
    // c[i]->cd();

    // l[i] = new TLegend(0.72,0.8,0.90,0.9);
    // l[i]->AddEntry(energy_hist_filter_in,"Filter In","l");
    // l[i]->AddEntry(energy_hist_filter_out,"Filter Out","l");

    // energy_hist_filter_in->GetXaxis()->SetTitle("Energy (in eV)");
    // energy_hist_filter_in->GetYaxis()->SetTitle("Counts (per proton)");
    // energy_hist_filter_in->SetTitle("Energy Histograms");
    // energy_hist_filter_in->Draw("HISTE");
    // gPad->SetGrid();
    // gPad->SetLogx();
    // gPad->SetLogy();
    // energy_hist_filter_out->SetLineColor(2);
    // energy_hist_filter_out->Draw("HISTESAME");
    // l[i]->Draw();

    // c[i]->Print("../plots/energy_hist.png");
}