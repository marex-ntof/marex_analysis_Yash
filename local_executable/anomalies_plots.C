/**
 * @file anomalies_plots.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2024-05-14
 */

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TFrame.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TLine.h"
#include "TAxis.h"
#include "TColor.h"
#include "TLegend.h"
#include "TAttMarker.h"
#include "TRandom3.h"

#include "parameters.h"
#include "tools.h"

TH2D* PTBC_al5_anomaly_det2 = 0;
TH2D* PTBC_al5_anomaly_det3 = 0;
TH2D* PTBC_al5_anomaly_det4 = 0;
TH2D* PTBC_al5_anomaly_det5 = 0;
TH2D* PTBC_al5_anomaly_det6 = 0;
TH2D* PTBC_al5_anomaly_det7 = 0;

TH2D* FIMG_al5_anomaly_det1 = 0;
TH2D* FIMG_al5_anomaly_det2 = 0;

TH2D* PTBC_filterOut_anomaly_det2 = 0;
TH2D* PTBC_filterOut_anomaly_det3 = 0;
TH2D* PTBC_filterOut_anomaly_det4 = 0;
TH2D* PTBC_filterOut_anomaly_det5 = 0;
TH2D* PTBC_filterOut_anomaly_det6 = 0;
TH2D* PTBC_filterOut_anomaly_det7 = 0;

TCutG* PTBC_tof_amp_cut_det2 = 0;
TCutG* PTBC_tof_amp_cut_det3 = 0;
TCutG* PTBC_tof_amp_cut_det4 = 0;
TCutG* PTBC_tof_amp_cut_det5 = 0;
TCutG* PTBC_tof_amp_cut_det6 = 0;
TCutG* PTBC_tof_amp_cut_det7 = 0;
TCutG* PTBC_tof_amp_cut_para = 0;
TCutG* FIMG_my_tof_amp_cut_dedi_det1 = 0;
TCutG* FIMG_my_tof_amp_cut_dedi_det2 = 0;
TCutG* FIMG_my_tof_amp_cut_para_det1 = 0;
TCutG* FIMG_my_tof_amp_cut_para_det2 = 0;

TH1D* retriveHistograms(const char *file_name, const char *hist_name){
    
    TFile* hist_file = TFile::Open(file_name, "READ");
    TH1D* hist_new = (TH1D*)hist_file->Get(hist_name);
    // hist_file->Close();

    return hist_new;
}

