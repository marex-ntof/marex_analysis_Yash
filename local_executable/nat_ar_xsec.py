import numpy as np
import ROOT
import math
import array

ar36_data_path_endf = '../evalData/Ar36_tot_xsec.txt'
ar38_data_path_endf = '../evalData/Ar38_tot_xsec.txt'
ar40_data_path_endf = '../evalData/Ar40_tot_xsec.txt'

ar36_data_path_jendl = '../evalData/JENDL_Ar36_tot_xsec.txt'
ar38_data_path_jendl = '../evalData/JENDL_Ar38_tot_xsec.txt'
ar40_data_path_jendl = '../evalData/JENDL_Ar40_tot_xsec.txt'

# min_e = 1e-2 #eV
# max_e = 2e7 #eV

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

def calculate_ar_xsec_hist(bins_per_decade, eval_data_base):
    if eval_data_base == 'endf':
        ar36_e_arr, ar36_xsec_arr = extract_xsec(ar36_data_path_endf)
        ar38_e_arr, ar38_xsec_arr = extract_xsec(ar38_data_path_endf)
        ar40_e_arr, ar40_xsec_arr = extract_xsec(ar40_data_path_endf)
        num_decades = 9
    
    if eval_data_base == 'jendl':
        ar36_e_arr, ar36_xsec_arr = extract_xsec(ar36_data_path_jendl)
        ar38_e_arr, ar38_xsec_arr = extract_xsec(ar38_data_path_jendl)
        ar40_e_arr, ar40_xsec_arr = extract_xsec(ar40_data_path_jendl)
        num_decades = 10

    #calculating bin edges
    num_bins = bins_per_decade * num_decades
    bin_edges = []
    bin_step = 1./bins_per_decade
    for i in range(0,num_bins+1):
        base = 10.
        exponent = (bin_step * i) - 2.
        bin_edges.append(base**exponent)

    bin_edges_array = array.array('d', bin_edges)

    #cross section hists
    ar36_xsec_hist = ROOT.TH1D("ar36_xsec_hist", "Ar 36 Cross Section Hist", num_bins, bin_edges_array)
    ar38_xsec_hist = ROOT.TH1D("ar38_xsec_hist", "Ar 38 Cross Section Hist", num_bins, bin_edges_array)
    ar40_xsec_hist = ROOT.TH1D("ar40_xsec_hist", "Ar 40 Cross Section Hist", num_bins, bin_edges_array)

    # natAr_xsec_hist = ROOT.TH1D("natAr_xsec_hist", "Natural Ar Cross Section Hist", num_bins, bin_edges_array)

    fill_xsec_hist(ar36_e_arr, ar36_xsec_arr, bin_edges, ar36_xsec_hist)
    fill_xsec_hist(ar38_e_arr, ar38_xsec_arr, bin_edges, ar38_xsec_hist)
    fill_xsec_hist(ar40_e_arr, ar40_xsec_arr, bin_edges, ar40_xsec_hist)

    natAr_bin_content = []

    for i in range(1,num_bins+1):
        new_bin_content = 0.996 * ar40_xsec_hist.GetBinContent(i) + 0.00063 * ar38_xsec_hist.GetBinContent(i) + 0.00334 * ar36_xsec_hist.GetBinContent(i)
        natAr_bin_content.append(new_bin_content)

    return bin_edges_array, natAr_bin_content

#ENDF hists
natAr_binEdges_20bpd_endf, natAr_binContent_20bpd_endf = calculate_ar_xsec_hist(20, 'endf')
natAr_binEdges_50bpd_endf, natAr_binContent_50bpd_endf = calculate_ar_xsec_hist(50, 'endf')
natAr_binEdges_100bpd_endf, natAr_binContent_100bpd_endf = calculate_ar_xsec_hist(100, 'endf')

natAr_xsec_hist_20bpd_endf = ROOT.TH1D("natAr_xsec_hist_20bpd_endf", "Natural Ar Cross Section Hist", len(natAr_binContent_20bpd_endf), natAr_binEdges_20bpd_endf)
natAr_xsec_hist_50bpd_endf = ROOT.TH1D("natAr_xsec_hist_50bpd_endf", "Natural Ar Cross Section Hist", len(natAr_binContent_50bpd_endf), natAr_binEdges_50bpd_endf)
natAr_xsec_hist_100bpd_endf = ROOT.TH1D("natAr_xsec_hist_100bpd_endf", "Natural Ar Cross Section Hist", len(natAr_binContent_100bpd_endf), natAr_binEdges_100bpd_endf)

for i in range(1,len(natAr_binContent_20bpd_endf)+1):
    natAr_xsec_hist_20bpd_endf.SetBinContent(i, natAr_binContent_20bpd_endf[i-1])

for i in range(1,len(natAr_binContent_50bpd_endf)+1):
    natAr_xsec_hist_50bpd_endf.SetBinContent(i, natAr_binContent_50bpd_endf[i-1])

for i in range(1,len(natAr_binContent_100bpd_endf)+1):
    natAr_xsec_hist_100bpd_endf.SetBinContent(i, natAr_binContent_100bpd_endf[i-1])

#jendl hists
natAr_binEdges_20bpd_jendl, natAr_binContent_20bpd_jendl = calculate_ar_xsec_hist(20, 'jendl')
natAr_binEdges_50bpd_jendl, natAr_binContent_50bpd_jendl = calculate_ar_xsec_hist(50, 'jendl')
natAr_binEdges_100bpd_jendl, natAr_binContent_100bpd_jendl = calculate_ar_xsec_hist(100, 'jendl')

natAr_xsec_hist_20bpd_jendl = ROOT.TH1D("natAr_xsec_hist_20bpd_jendl", "Natural Ar Cross Section Hist", len(natAr_binContent_20bpd_jendl), natAr_binEdges_20bpd_jendl)
natAr_xsec_hist_50bpd_jendl = ROOT.TH1D("natAr_xsec_hist_50bpd_jendl", "Natural Ar Cross Section Hist", len(natAr_binContent_50bpd_jendl), natAr_binEdges_50bpd_jendl)
natAr_xsec_hist_100bpd_jendl = ROOT.TH1D("natAr_xsec_hist_100bpd_jendl", "Natural Ar Cross Section Hist", len(natAr_binContent_100bpd_jendl), natAr_binEdges_100bpd_jendl)

for i in range(1,len(natAr_binContent_20bpd_jendl)+1):
    natAr_xsec_hist_20bpd_jendl.SetBinContent(i, natAr_binContent_20bpd_jendl[i-1])

for i in range(1,len(natAr_binContent_50bpd_jendl)+1):
    natAr_xsec_hist_50bpd_jendl.SetBinContent(i, natAr_binContent_50bpd_jendl[i-1])

for i in range(1,len(natAr_binContent_100bpd_jendl)+1):
    natAr_xsec_hist_100bpd_jendl.SetBinContent(i, natAr_binContent_100bpd_jendl[i-1])

# c1 = ROOT.TCanvas("", "", 800, 500)
# c1.SetGridx(1)
# c1.SetGridy(1)

# l1 = ROOT.TLegend(0.2,0.25,0.35,0.35)
# l1.AddEntry(natAr_xsec_hist_20bpd_endf, "Nat Ar 20 BPD", "l")
# l1.AddEntry(natAr_xsec_hist_50bpd_endf, "Nat Ar 50 BPD", "l")
# l1.AddEntry(natAr_xsec_hist_100bpd_endf, "Nat Ar 100 BPD", "l")

# natAr_xsec_hist_20bpd_endf.SetStats(0)
# natAr_xsec_hist_20bpd_endf.SetLineColor(1)
# natAr_xsec_hist_20bpd_endf.Draw() #"P"
# natAr_xsec_hist_50bpd_endf.SetLineColor(2)
# natAr_xsec_hist_50bpd_endf.Draw("SAME")
# natAr_xsec_hist_100bpd_endf.SetLineColor(4)
# natAr_xsec_hist_100bpd_endf.Draw("SAME")

# ar36_xsec_hist.SetStats(0)
# ar36_xsec_hist.GetYaxis().SetRangeUser(1e-4,100)
# ar36_xsec_hist.SetLineColor(1)
# # ar36_xsec_hist.SetMarkerStyle(2)
# # ar36_xsec_hist.SetMarkerSize(1)
# ar36_xsec_hist.Draw() #"P"
# ar38_xsec_hist.SetLineColor(2)
# ar38_xsec_hist.Draw("SAME")
# ar40_xsec_hist.SetLineColor(4)
# ar40_xsec_hist.Draw("SAME")

output_file = ROOT.TFile("../inputFiles/natAr_xsec_hists.root", "RECREATE")
natAr_xsec_hist_20bpd_endf.Write()
natAr_xsec_hist_50bpd_endf.Write()
natAr_xsec_hist_100bpd_endf.Write()
natAr_xsec_hist_20bpd_jendl.Write()
natAr_xsec_hist_50bpd_jendl.Write()
natAr_xsec_hist_100bpd_jendl.Write()
output_file.Close()

# l1.Draw()
# c1.SetLogx()
# c1.SetLogy()
# c1.Print("../plots/natAr_xsec_hists.png")

input("Press Enter to exit...")