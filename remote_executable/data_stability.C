/**
 * @file data_stability.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2023-12-04
 */

#include "data_stability.h"

void data_stability(){

    Int_t num_bins = 240; // 4 hours
    Double_t x_min = 0.;
    Double_t y_min = 240.;

    stability_hist_PTBC = new TH1D("stability_hist_PTBC","Stability Hist - PTBC",num_bins,x_min,y_min);
    stability_hist_FIMG_det1 = new TH1D("stability_hist_FIMG_det1","Stability Hist - FIMG - Det 1",num_bins,x_min,y_min);
    stability_hist_FIMG_det2 = new TH1D("stability_hist_FIMG_det2","Stability Hist - FIMG - Det 1",num_bins,x_min,y_min);

    std::vector<Int_t> run_list = {117471};

    for (int i = 0; i < run_list.size(); i++)
    {
        TFile *file_ntof = TFile::Open(Form("/eos/experiment/ntof/processing/official/done/run%d.root", run_list.at(i)),"read");
        cout << "Run Number = " << run_list.at(i) << endl;

        //PTBC ---------------------------------------------
        TTree* PTBC;
        Int_t time_PTBC = 0;
        Int_t BunchNumber_PTBC = 0;
        Float_t PulseIntensity_PTBC = 0;

        file_ntof->GetObject("PTBC", PTBC);
        PTBC->SetBranchAddress("time", &time_PTBC);
        PTBC->SetBranchAddress("BunchNumber", &BunchNumber_PTBC);
        PTBC->SetBranchAddress("PulseIntensity", &PulseIntensity_PTBC);

        Long64_t Events_PTB = PTBC->GetEntriesFast();
        std::cout << "Number of entries - PTBC = " << Events_PTB << std::endl;
        
        Int_t CurrentBunchNum = 0;
        Int_t CurrentTime = 0;
        Float_t sum_pulse_intensity = 0;
        Int_t sum_num_events = 0;
        Int_t bin_num = 1;

        for (int j = 0; j < Events_PTB; j++)
        {
            PTBC->GetEntry(j);

            if (j == 0)
            {
                CurrentTime = time_PTBC;
                CurrentBunchNum = BunchNumber_PTBC;
            }        

            if (time_PTBC - CurrentTime > 60)
            {
                Double_t bin_content = sum_num_events/sum_pulse_intensity;
                stability_hist_PTBC->SetBinContent(bin_num, bin_content);
                bin_num++;
                sum_num_events = 0;
                sum_pulse_intensity = 0;
                CurrentTime = time_PTBC;
            }

            if (CurrentBunchNum != BunchNumber_PTBC)
            {
                CurrentBunchNum = BunchNumber_PTBC;
                CurrentTime = time_PTBC;
                sum_pulse_intensity += PulseIntensity_PTBC;
            }

            if (CurrentBunchNum == BunchNumber_PTBC)
            {
                sum_num_events++;
            }
        }

        //FIMG ---------------------------------------------
        TTree* FIMG;
        Int_t time_FIMG = 0;
        Int_t BunchNumber_FIMG = 0;
        Int_t det_num = 0;
        Float_t PulseIntensity_FIMG = 0;

        file_ntof->GetObject("FIMG", FIMG);
        FIMG->SetBranchAddress("time", &time_FIMG);
        FIMG->SetBranchAddress("BunchNumber", &BunchNumber_FIMG);
        FIMG->SetBranchAddress("detn", &det_num);
        FIMG->SetBranchAddress("PulseIntensity", &PulseIntensity_FIMG);

        Long64_t Events_FIMG = FIMG->GetEntriesFast();
        std::cout << "Number of entries - FIMG = " << Events_FIMG << std::endl;

        //reset variables
        CurrentBunchNum = 0;
        CurrentTime = 0;
        bin_num = 1;
        sum_pulse_intensity = 0;
        Int_t sum_num_events_det1 = 0;
        Int_t sum_num_events_det2 = 0;

        for (int j = 0; j < Events_FIMG; j++)
        {
            FIMG->GetEntry(j);

            if (j == 0)
            {
                CurrentTime = time_FIMG;
                CurrentBunchNum = BunchNumber_FIMG;
            }        

            if (time_FIMG - CurrentTime > 60)
            {
                Double_t bin_content_det1 = sum_num_events_det1/sum_pulse_intensity;
                stability_hist_FIMG_det1->SetBinContent(bin_num, bin_content_det1);
                stability_hist_FIMG_det1->SetBinError(bin_num, TMath::Sqrt(sum_num_events_det1)/sum_pulse_intensity);
                Double_t bin_content_det2 = sum_num_events_det2/sum_pulse_intensity;
                stability_hist_FIMG_det2->SetBinContent(bin_num, bin_content_det2);
                stability_hist_FIMG_det2->SetBinError(bin_num, TMath::Sqrt(sum_num_events_det2)/sum_pulse_intensity);
                bin_num++;
                sum_num_events_det1 = 0;
                sum_num_events_det2 = 0;
                sum_pulse_intensity = 0;
                CurrentTime = time_FIMG;
            }

            if (CurrentBunchNum != BunchNumber_FIMG)
            {
                CurrentBunchNum = BunchNumber_FIMG;
                CurrentTime = time_FIMG;
                sum_pulse_intensity += PulseIntensity_FIMG;
            }

            if (CurrentBunchNum == BunchNumber_FIMG)
            {
                if (det_num == 1)
                {
                    sum_num_events_det1++;
                }
                if (det_num == 2)
                {
                    sum_num_events_det2++;
                }
            }
        }

        file_ntof->Close();
    }

    //Writing to the output file
    outputRootFile = new TFile("../rootFiles/data_stability.root","recreate");
    stability_hist_PTBC->Write();
    stability_hist_FIMG_det1->Write();
    stability_hist_FIMG_det2->Write();

    outputRootFile->Close();
}