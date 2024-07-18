import numpy as np
import ROOT
import math
import matplotlib.pyplot as plt

def calculate_skewness(hist, mean_val, std_val):
    num_bins = hist.GetNbinsX()
    x = 0
    sum = 0
    np = 0
    for i in range(1,num_bins+1):
        x = hist.GetBinCenter(i)
        freq = hist.GetBinContent(i)
        np += freq
        sum += freq*(x-mean_val)*(x-mean_val)*(x-mean_val)
    sum = sum / (np*std_val*std_val*std_val)
    return sum

def calculate_kurtosis(hist, mean_val, std_val):
    num_bins = hist.GetNbinsX()
    x = 0
    sum = 0
    np = 0
    for i in range(1,num_bins+1):
        x = hist.GetBinCenter(i)
        freq = hist.GetBinContent(i)
        np += freq
        sum += freq*(x-mean_val)*(x-mean_val)*(x-mean_val)*(x-mean_val)
    sum = sum / (np*std_val*std_val*std_val*std_val)
    return sum

def determine_skew_kurt(det_num, tof_amp_hist):
    PTBC_tof_amp_hist_forCuts = tof_amp_hist.Rebin2D(int(1000 / 5), 5, f"tof_amp_forCuts_det{det_num}")
    bin_center = []
    skewness = []
    kurtosis = []
    for i in range(5, 11):
        projection_name = f"profile_forCuts_bin_{i}"
        proj_hist = PTBC_tof_amp_hist_forCuts.ProjectionY(projection_name, i, i)
        proj_hist.GetXaxis().SetRangeUser(0., 10000.)
        gaus_fit_total = ROOT.TF1("gaus_fit_total", "gaus", 0, 10000)
        proj_hist.Fit(gaus_fit_total, "0R")

        mean_val = gaus_fit_total.GetParameter(1)
        std_dev = gaus_fit_total.GetParameter(2)

        x_val = PTBC_tof_amp_hist_forCuts.GetXaxis().GetBinCenter(i)
        y_val_skew = calculate_skewness(proj_hist, mean_val, std_dev)
        y_val_kurt = calculate_kurtosis(proj_hist, mean_val, std_dev)

        bin_center.append(x_val)
        skewness.append(y_val_skew)
        kurtosis.append(y_val_kurt)
    return bin_center, skewness, kurtosis

def main():
    rootFile = ROOT.TFile.Open("../rootFiles/cutoffAnalysis_PTBC_none.root")
    for i in range(0,6):
        hist = rootFile.Get(f"PTBC_tof_amp_det{i+2}")
        bin_c, skew, kurt = determine_skew_kurt(i+2, hist)

        fig, ax1 = plt.subplots()

        skew_plot = ax1.scatter(bin_c, skew, color='b', label='Skewness', marker="o", s=50, alpha=0.5)
        ax1.set_xlabel('ToF (in ns)')
        ax1.set_xscale('log')
        ax1.set_ylabel('Skewness', color='b')
        ax1.set_yscale('log')
        ax1.tick_params(axis='y', labelcolor='b')

        ax2 = ax1.twinx()
        kurt_plot = ax2.scatter(bin_c, kurt, color='r', label='Kurtosis', marker="o", s=50, alpha=0.5)
        ax2.set_ylabel('Kurtosis', color='r')
        ax2.set_yscale('log')
        ax2.tick_params(axis='y', labelcolor='r')

        ax1.set_title(f"Skewness and Kurtosis - Det {i+2} - No Target")

        fig.show()

if __name__ == "__main__":
    main()
    input("Press Enter to exit...")