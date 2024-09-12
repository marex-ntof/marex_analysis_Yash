import numpy as np
import ROOT
import math
import array
import matplotlib.pyplot as plt

c_data_path = '../evalData/C_tot_xsec.txt'
h_data_path = '../evalData/JENDL_H_tot_xsec.txt'
ar_t = 10 #cm
cf_t = 0.04 #cm
bins_per_decade = 50
ar_vary_t = 1 #cm
cf_vary_t = 0.01 #cm

ar_bottle_pressure = 197.385 * 1e5
ar_bottle_temp = 293.0

n_value_ar = ar_t * ((ar_bottle_pressure)/(8.31446261815324 * ar_bottle_temp * 1e6)) * (6.02214076e23) * (1e-24)
n_value_ar_more = (ar_t+1) * ((ar_bottle_pressure)/(8.31446261815324 * ar_bottle_temp * 1e6)) * (6.02214076e23) * (1e-24)
n_value_ar_less = (ar_t-1) * ((ar_bottle_pressure)/(8.31446261815324 * ar_bottle_temp * 1e6)) * (6.02214076e23) * (1e-24)

cf_density = 1.93 #g/cm3
#2.267 - c density
n_value_c = cf_t * cf_density * 0.9 * (6.02214076e23) * (1e-24) / (12.011)
n_value_c_more = (cf_t+cf_vary_t) * cf_density * 0.9 * (6.02214076e23) * (1e-24) / (12.011)
n_value_c_less = (cf_t-cf_vary_t) * cf_density * 0.9 * (6.02214076e23) * (1e-24) / (12.011)

n_value_h = cf_t * cf_density * 0.1 * (6.02214076e23) * (1e-24) / (1.008)
n_value_h_more = (cf_t+cf_vary_t) * cf_density * 0.1 * (6.02214076e23) * (1e-24) / (1.008)
n_value_h_less = (cf_t-cf_vary_t) * cf_density * 0.1 * (6.02214076e23) * (1e-24) / (1.008)

def extract_xsec(file_path):
    e_arr = []
    xsec_arr = []
    with open(file_path, 'r') as text_file:
        for line in text_file:
            values = line.split()
            if values[0] == '#E,eV':
                continue
            if values[0] == '#END':
                break
            e_arr.append(float(values[0]))
            xsec_arr.append(float(values[1]))
    return e_arr, xsec_arr

def fill_xsec_hist(e_arr, xsec_arr, e_bin_edges, xsec_hist):

    bin_content_dict = {}

    for i in range(1,len(e_bin_edges)):
        bin_content_dict[i] = []

    for i in range(0,len(e_arr)):
        if e_arr[i] < e_bin_edges[0]:
            continue
        if e_arr[i] >= e_bin_edges[-1]:
            break
        bin_num = xsec_hist.GetXaxis().FindBin(e_arr[i])
        if bin_num > (len(e_bin_edges)-1):
            break
        bin_content_dict[bin_num].append(xsec_arr[i])

    bin_contents = []

    for bin_n, xsecs in bin_content_dict.items():
        if not xsecs:
            xsec_hist.SetBinContent(bin_n, 0)
            bin_contents.append(0)
        else:
            bin_content = sum(xsecs)/len(xsecs)
            xsec_hist.SetBinContent(bin_n, bin_content)
            bin_contents.append(bin_content)

    bin_contents = np.array(bin_contents)
    num_bins_e = len(e_bin_edges) - 1
    
    # Find indices of zero bins
    zero_indices = np.where(bin_contents == 0)[0]

    for idx in zero_indices:
        left_idx = idx - 1
        right_idx = idx + 1
        
        # Find the nearest non-zero bins to the left
        while left_idx >= 0 and bin_contents[left_idx] == 0:
            left_idx -= 1
        
        # Find the nearest non-zero bins to the right
        while right_idx < num_bins_e and bin_contents[right_idx] == 0:
            right_idx += 1
        
        if left_idx >= 0 and right_idx < num_bins_e:
            left_value = bin_contents[left_idx]
            left_e_val = xsec_hist.GetBinCenter(int(left_idx + 1))
            right_value = bin_contents[right_idx]
            right_e_val = xsec_hist.GetBinCenter(int(right_idx + 1))
            #linear interpolation
            slope = (right_value-left_value)/(right_e_val-left_e_val)
            num_bins_to_update = right_idx - left_idx - 1
            for j in range(1,num_bins_to_update+1):
                bin_center = xsec_hist.GetBinCenter(int(left_idx + j + 1))
                bin_contents[left_idx + j] = slope * (bin_center - left_e_val) + left_value

    # Update histogram with interpolated values
    for i in range(1, num_bins_e + 1):
        xsec_hist.SetBinContent(i, bin_contents[i - 1])

def fill_trans_array(xsec_hist_array, n_array, trans_array):
    num_bins = xsec_hist_array[0].GetXaxis().GetNbins()
    for i in range(0,num_bins):
        bin_content = 1
        for hist, n_val in zip(xsec_hist_array, n_array):
            bin_content = bin_content * math.exp(- n_val * hist.GetBinContent(i+1))
        trans_array.append(bin_content)

def calc_trans(xsec_hist, n_value, trans_array):
    num_bins = xsec_hist.GetXaxis().GetNbins()
    for i in range(0,num_bins):
        bin_content = math.exp(- n_value * xsec_hist.GetBinContent(i+1))
        trans_array.append(bin_content)

def varrying_ar_c_thickness():
    arFile = ROOT.TFile.Open("../inputFiles/natAr_xsec_hists.root")
    ar_xsec_hist = arFile.Get("natAr_xsec_hist_50bpd_jendl")

    c_e_arr, c_xsec_arr = extract_xsec(c_data_path)
    h_e_arr, h_xsec_arr = extract_xsec(h_data_path)

    #calculating bin edges
    num_decades = 10
    num_bins = bins_per_decade * num_decades
    bin_edges = []
    bin_step = 1./bins_per_decade
    for i in range(0,num_bins+1):
        base = 10.
        exponent = (bin_step * i) - 2.
        bin_edges.append(base**exponent)

    bin_edges_array = array.array('d', bin_edges)

    #cross section hists
    c_xsec_hist = ROOT.TH1D("c_xsec_hist", "Carbon Cross Section Hist", num_bins, bin_edges_array)
    h_xsec_hist = ROOT.TH1D("h_xsec_hist", "Hydrogen Cross Section Hist", num_bins, bin_edges_array)

    fill_xsec_hist(c_e_arr, c_xsec_arr, bin_edges, c_xsec_hist)
    fill_xsec_hist(h_e_arr, h_xsec_arr, bin_edges, h_xsec_hist)

    trans_actual = []
    trans_more_ar = []
    trans_more_cf = []
    trans_less_ar = []
    trans_less_cf = []

    fill_trans_array([ar_xsec_hist, c_xsec_hist, h_xsec_hist], [n_value_ar, n_value_c, n_value_h], trans_actual)
    fill_trans_array([ar_xsec_hist, c_xsec_hist, h_xsec_hist], [n_value_ar_more, n_value_c, n_value_h], trans_more_ar)
    fill_trans_array([ar_xsec_hist, c_xsec_hist, h_xsec_hist], [n_value_ar_less, n_value_c, n_value_h], trans_less_ar)
    fill_trans_array([ar_xsec_hist, c_xsec_hist, h_xsec_hist], [n_value_ar, n_value_c_more, n_value_h_more], trans_more_cf)
    fill_trans_array([ar_xsec_hist, c_xsec_hist, h_xsec_hist], [n_value_ar, n_value_c_less, n_value_h_less], trans_less_cf)
    
    # fill_trans_array(ar_xsec_hist, n_value_ar, c_xsec_hist, n_value_c, trans_actual)
    # fill_trans_array(ar_xsec_hist, n_value_ar_more, c_xsec_hist, n_value_c, trans_more_ar)
    # fill_trans_array(ar_xsec_hist, n_value_ar_less, c_xsec_hist, n_value_c, trans_less_ar)
    # fill_trans_array(ar_xsec_hist, n_value_ar, c_xsec_hist, n_value_c_more, trans_more_cf)
    # fill_trans_array(ar_xsec_hist, n_value_ar, c_xsec_hist, n_value_c_less, trans_less_cf)

    energy_bin_centers = []
    for i in range(0,num_bins):
        energy_bin_centers.append(ar_xsec_hist.GetBinCenter(i+1))

    plt.figure(figsize=(8, 6))
    plt.plot(energy_bin_centers[2:], trans_actual[2:], linestyle = '-', color='black', label=f'Transmission - {ar_t}cm Ar, {cf_t}cm CF')
    plt.fill_between(energy_bin_centers[2:], trans_less_ar[2:], trans_more_ar[2:], alpha=0.4, color='blue', label=f'Varrying Ar {ar_vary_t}cm')
    plt.fill_between(energy_bin_centers[2:], trans_less_cf[2:], trans_more_cf[2:], alpha=0.2, color='red', label=f'Varrying CF {cf_vary_t}cm')
    # plt.xlim(1e5,1e7)
    plt.xscale('log')
    plt.xlabel('Incident Energy (keV)')
    plt.ylabel('Transmission')
    plt.legend()
    plt.show()

def different_trans_ar_c():
    arFile = ROOT.TFile.Open("../inputFiles/natAr_xsec_hists.root")
    ar_xsec_hist = arFile.Get("natAr_xsec_hist_50bpd_jendl")

    c_e_arr, c_xsec_arr = extract_xsec(c_data_path)
    h_e_arr, h_xsec_arr = extract_xsec(h_data_path)

    #calculating bin edges
    num_decades = 10
    num_bins = bins_per_decade * num_decades
    bin_edges = []
    bin_step = 1./bins_per_decade
    for i in range(0,num_bins+1):
        base = 10.
        exponent = (bin_step * i) - 2.
        bin_edges.append(base**exponent)

    bin_edges_array = array.array('d', bin_edges)

    #cross section hists
    c_xsec_hist = ROOT.TH1D("c_xsec_hist", "Carbon Cross Section Hist", num_bins, bin_edges_array)
    h_xsec_hist = ROOT.TH1D("h_xsec_hist", "Hydrogen Cross Section Hist", num_bins, bin_edges_array)

    fill_xsec_hist(c_e_arr, c_xsec_arr, bin_edges, c_xsec_hist)
    fill_xsec_hist(h_e_arr, h_xsec_arr, bin_edges, h_xsec_hist)

    trans_actual = []
    trans_more_ar = []
    trans_more_cf = []
    trans_less_ar = []
    trans_less_cf = []

    calc_trans(ar_xsec_hist, n_value_ar, trans_actual)
    calc_trans(ar_xsec_hist, n_value_ar_more, trans_more_ar)
    calc_trans(ar_xsec_hist, n_value_ar_less, trans_less_ar)

    # num_bins = xsec_hist_1.GetXaxis().GetNbins()
    for i in range(0,num_bins):
        numa = math.exp(- n_value_ar * ar_xsec_hist.GetBinContent(i+1)) * math.exp(- n_value_c * c_xsec_hist.GetBinContent(i+1)) * math.exp(- n_value_h * h_xsec_hist.GetBinContent(i+1))
        denom_1 = math.exp(- n_value_c_more * c_xsec_hist.GetBinContent(i+1)) * math.exp(- n_value_h_more * h_xsec_hist.GetBinContent(i+1))
        denom_2 = math.exp(- n_value_c_less * c_xsec_hist.GetBinContent(i+1)) * math.exp(- n_value_h_less * h_xsec_hist.GetBinContent(i+1))
        trans_more_cf.append(numa/denom_1)
        trans_less_cf.append(numa/denom_2)

    energy_bin_centers = []
    for i in range(0,num_bins):
        energy_bin_centers.append(ar_xsec_hist.GetBinCenter(i+1))

    plt.figure(figsize=(8, 6))
    plt.plot(energy_bin_centers[2:], trans_actual[2:], linestyle = '-', color='black', label=f'Transmission - {ar_t}cm Ar, {cf_t}cm CF')
    plt.fill_between(energy_bin_centers[2:], trans_less_ar[2:], trans_more_ar[2:], alpha=0.4, color='blue', label=f'Varrying Ar {ar_vary_t}cm')
    plt.fill_between(energy_bin_centers[2:], trans_less_cf[2:], trans_more_cf[2:], alpha=0.2, color='red', label=f'Varrying CF {cf_vary_t}cm')
    # plt.xlim(1e4,1e7)
    # plt.ylim(0.8,1.1)
    plt.xscale('log')
    plt.xlabel('Incident Energy (keV)')
    plt.ylabel('Transmission')
    plt.legend()
    plt.show()

def cf_trans():
    c_e_arr, c_xsec_arr = extract_xsec(c_data_path)
    h_e_arr, h_xsec_arr = extract_xsec(h_data_path)

    #calculating bin edges
    num_decades = 10
    num_bins = bins_per_decade * num_decades
    bin_edges = []
    bin_step = 1./bins_per_decade
    for i in range(0,num_bins+1):
        base = 10.
        exponent = (bin_step * i) - 2.
        bin_edges.append(base**exponent)

    bin_edges_array = array.array('d', bin_edges)

    #cross section hists
    c_xsec_hist = ROOT.TH1D("c_xsec_hist", "Carbon Cross Section Hist", num_bins, bin_edges_array)
    h_xsec_hist = ROOT.TH1D("h_xsec_hist", "Hydrogen Cross Section Hist", num_bins, bin_edges_array)

    fill_xsec_hist(c_e_arr, c_xsec_arr, bin_edges, c_xsec_hist)
    fill_xsec_hist(h_e_arr, h_xsec_arr, bin_edges, h_xsec_hist)

    trans_cf = []
    trans_cf_more = []
    trans_cf_less = []

    fill_trans_array([c_xsec_hist, h_xsec_hist], [n_value_c, n_value_h], trans_cf)
    fill_trans_array([c_xsec_hist, h_xsec_hist], [n_value_c_more, n_value_h_more], trans_cf_more)
    fill_trans_array([c_xsec_hist, h_xsec_hist], [n_value_c_less, n_value_h_less], trans_cf_less)

    energy_bin_centers = []
    for i in range(0,num_bins):
        energy_bin_centers.append(c_xsec_hist.GetBinCenter(i+1))

    plt.figure(figsize=(8, 6))
    plt.plot(energy_bin_centers[2:], trans_cf[2:], linestyle = '-', color='black', label=f'Transmission - {cf_t}cm CF')
    plt.fill_between(energy_bin_centers[2:], trans_cf_less[2:], trans_cf_more[2:], alpha=0.2, color='red', label=f'Varrying CF {cf_vary_t}cm')
    plt.xlim(1e-1,1e2)
    plt.ylim(0.2,0.4)
    plt.xscale('log')
    plt.xlabel('Incident Energy (keV)')
    plt.ylabel('Transmission')
    plt.legend()
    plt.show()

def main():
    varrying_ar_c_thickness()
    # different_trans_ar_c()
    # cf_trans()

if __name__ == "__main__":
    main()
    input("Press Enter to exit...")