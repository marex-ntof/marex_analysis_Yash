/**
 * @file pkup.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-11-13
 */

#include "pkup.h"
/*
While changing the filter target
- change root file name
- change filter in runs
- change txt file name in endf() (not needed if onl changing thickness)
- change max line count in endf() (not needed if onl changing thickness)
- change names and titles of the plots
- change n_inverse
- change n passed to endf()
*/

void StoreHist(){
    
    TFile *f = new TFile("../rootFiles/pkup.root","recreate");

    tpkup_beam_intensity_hist->Write();
    delT_PTBC_beam_intensity_hist->Write();
    delT_FIMG_beam_intensity_hist->Write();

    f->Close();
}

void plots(){
    
    TCanvas *c[3];

    int i = 0;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();
    tpkup_beam_intensity_hist->GetXaxis()->SetTitle("Pulse Intensity (in Num Protons)");
    tpkup_beam_intensity_hist->GetYaxis()->SetTitle("T PKUP (in ns)");
    tpkup_beam_intensity_hist->SetTitle("Pulse Intensity vs T pkup Hist");
    tpkup_beam_intensity_hist->Draw("colz");
    // tpkup_beam_intensity_hist->SetMarkerStyle(6);
    // tpkup_beam_intensity_hist->SetMarkerSize(0.5);
    // gPad->SetLogx();
    gPad->SetLogz();
    gStyle->SetPalette(57);

    c[i]->Print("../plots/tpkup_pulse_intensity_hist.png");

    i++;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();
    delT_PTBC_beam_intensity_hist->GetXaxis()->SetTitle("Pulse Intensity (in Num Protons)");
    delT_PTBC_beam_intensity_hist->GetYaxis()->SetTitle("del T (PKUP - PTBC) (in ns)");
    delT_PTBC_beam_intensity_hist->SetTitle("Pulse Intensity vs del T (PKUP - PTBC) Hist");
    delT_PTBC_beam_intensity_hist->Draw("colz");
    // delT_PTBC_beam_intensity_hist->SetMarkerStyle(6);
    // delT_PTBC_beam_intensity_hist->SetMarkerSize(0.5);
    // gPad->SetLogx();
    gPad->SetLogz();
    gStyle->SetPalette(57);

    c[i]->Print("../plots/delT_PTBC_pulse_intensity_hist.png");

    i++;

    c[i] = new TCanvas(Form("c%d", i)," ");
    c[i]->cd();
    delT_FIMG_beam_intensity_hist->GetXaxis()->SetTitle("Pulse Intensity (in Num Protons)");
    delT_FIMG_beam_intensity_hist->GetYaxis()->SetTitle("del T (PKUP - FIMG) (in ns)");
    delT_FIMG_beam_intensity_hist->SetTitle("Pulse Intensity vs del T (PKUP - FIMG) Hist");
    delT_FIMG_beam_intensity_hist->Draw("colz");
    // delT_FIMG_beam_intensity_hist->SetMarkerStyle(6);
    // delT_FIMG_beam_intensity_hist->SetMarkerSize(0.5);
    // gPad->SetLogx();
    gPad->SetLogz();
    gStyle->SetPalette(57);

    c[i]->Print("../plots/delT_FIMG_pulse_intensity_hist.png");
}

void pkup(){

    fillRuns();

    tpkup_beam_intensity_hist = new TH2D("tpkup_beam_intensity_hist", "Pulse Intensity vs T pkup Hist",50,1e12,9e12,50,11200,12200);
    delT_PTBC_beam_intensity_hist = new TH2D("delT_PTBC_beam_intensity_hist", "Pulse Intensity vs del T (PKUP - PTBC) Hist",50,1e12,9e12,100,200,1000);
    delT_FIMG_beam_intensity_hist = new TH2D("delT_FIMG_beam_intensity_hist", "Pulse Intensity vs del T (PKUP - FIMG) Hist",50,1e12,9e12,100,200,1000);

    for (int i = 0; i < list_of_runs.size(); i++)
    {
        TFile *file_ntof = TFile::Open(Form("/eos/experiment/ntof/processing/official/done/run%d.root", list_of_runs.at(i)),"read");
        cout << "Run Number = " << list_of_runs.at(i) << endl;

        //PKUP ---------------------------------------------
        TTree* PKUP;
        Int_t BunchNumber_PKUP = 0;
        Double_t tpkup = 0;
        Float_t PulseIntensity_PKUP = 0;

        file_ntof->GetObject("PKUP", PKUP);
        PKUP->SetBranchAddress("PulseIntensity", &PulseIntensity_PKUP);
        PKUP->SetBranchAddress("BunchNumber", &BunchNumber_PKUP);
        PKUP->SetBranchAddress("tflash", &tpkup);

        std::map<Int_t, Double_t> BNum_tpkup_map;
        Long64_t Events_PKUP = PKUP->GetEntriesFast();

        for (int j = 0; j < Events_PKUP; j++)
        {
            PKUP->GetEntry(j);

            BNum_tpkup_map.emplace(BunchNumber_PKUP, tpkup);
            tpkup_beam_intensity_hist->Fill( (Double_t) PulseIntensity_PKUP, tpkup);
        }

        //PTBC ---------------------------------------------
        TTree* PTBC;
        Double_t tflash_PTBC = 0; //tflash_PTBC is in ns
        Int_t BunchNumber_PTBC = 0;
        Float_t PulseIntensity_PTBC = 0;

        file_ntof->GetObject("PTBC", PTBC);
        PTBC->SetBranchAddress("BunchNumber", &BunchNumber_PTBC);
        PTBC->SetBranchAddress("PulseIntensity", &PulseIntensity_PTBC);
        PTBC->SetBranchAddress("tflash", &tflash_PTBC);

        Long64_t Events_PTBC = PTBC->GetEntriesFast();
        std::cout << "Number of entries PTBC = " << Events_PTBC << std::endl;

        int CurrentBunchNum_PTBC = 0;

        for (int j = 0; j < Events_PTBC; j++)
        {
            PTBC->GetEntry(j);

            if (CurrentBunchNum_PTBC != BunchNumber_PTBC)
            {
                CurrentBunchNum_PTBC = BunchNumber_PTBC;
                // Double_t t_pkup = BNum_tpkup_map[BunchNumber_PTBC];
                Double_t del_tflash_PTBC = BNum_tpkup_map[BunchNumber_PTBC] - tflash_PTBC;
                delT_PTBC_beam_intensity_hist->Fill( (Double_t) PulseIntensity_PTBC, del_tflash_PTBC);
            }
        }

        //FIMG ---------------------------------------------
        TTree* FIMG;
        Double_t tflash_FIMG = 0; //tflash_FIMG is in ns
        Int_t BunchNumber_FIMG = 0;
        Float_t PulseIntensity_FIMG = 0;

        file_ntof->GetObject("FIMG", FIMG);
        FIMG->SetBranchAddress("BunchNumber", &BunchNumber_FIMG);
        FIMG->SetBranchAddress("PulseIntensity", &PulseIntensity_FIMG);
        FIMG->SetBranchAddress("tflash", &tflash_FIMG);

        Long64_t Events_FIMG = FIMG->GetEntriesFast();
        std::cout << "Number of entries FIMG = " << Events_FIMG << std::endl;

        int CurrentBunchNum_FIMG = 0;

        for (int j = 0; j < Events_FIMG; j++)
        {
            FIMG->GetEntry(j);

            if (CurrentBunchNum_FIMG != BunchNumber_FIMG)
            {
                CurrentBunchNum_FIMG = BunchNumber_FIMG;
                // Double_t t_pkup = BNum_tpkup_map[BunchNumber_FIMG];
                Double_t del_tflash_FIMG = BNum_tpkup_map[BunchNumber_FIMG] - tflash_FIMG;
                delT_FIMG_beam_intensity_hist->Fill( (Double_t) PulseIntensity_FIMG, del_tflash_FIMG);
            }
        }
    }
    // endf(n_Al_8cm, bin_edges_e, true, false);
    StoreHist();
    plots();
}