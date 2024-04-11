import numpy as np
import ROOT
import math

ca_data_path = '../evalData/Ca_tot_xsec.txt'
si_data_path = '../evalData/INDEN_Si_tot_xsec.txt'
o_data_path = '../evalData/INDEN_O_tot_xsec.txt'
al_data_path = '../evalData/Al_tot_xsec.txt'
fe_data_path = '../evalData/INDEN_Fe_tot_xsec.txt'
h_data_path = '../evalData/JENDL_H_tot_xsec.txt'
s_endf_data_path = '../evalData/S_tot_xsec.txt'
s_jeff_data_path = '../evalData/JEFF_S_tot_xsec.txt'

min_energy = 1e6 #eV
max_energy = 1.5e8 #eV

def extract_xsec(file_path, is_sulfur):
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
    
    if is_sulfur == True:
        with open(s_jeff_data_path, 'r') as text_file:
            for line in text_file:
                values = line.split()
                if values[0] == '#E,eV':
                    continue
                if values[0] == '#END':
                    break
                if float(values[0]) > 2e7:
                    e_arr.append(float(values[0]))
                    xsec_arr.append(float(values[1]))

    return e_arr, xsec_arr

def fill_xsec_hist(e_arr, xsec_arr, e_bin_edges, xsec_hist):
    xsec_sum = 0
    sum_counter = 0
    bin_counter = 1

    for i in range(0,len(e_arr)):
        if e_arr[i] < min_energy:
            continue
        if e_arr[i] > max_energy:
            break
        if e_arr[i] < e_bin_edges[bin_counter]:
            xsec_sum+=xsec_arr[i]
            sum_counter+=1
        if e_arr[i] >= e_bin_edges[bin_counter]:
            if sum_counter == 0:
                xsec_hist.SetBinContent(bin_counter, 0)
            else:
                xsec_hist.SetBinContent(bin_counter, xsec_sum/sum_counter)
            xsec_sum = 0
            sum_counter = 0
            bin_counter+=1
            i-=1

def fill_trans_hist(xsec_hist, n_value, trans_hist):
    num_bins = xsec_hist.GetXaxis().GetNbins()
    for i in range(0,num_bins):
        bin_content = math.exp(- n_value * xsec_hist.GetBinContent(i+1))
        trans_hist.SetBinContent(i+1, bin_content)

ca_e_arr, ca_xsec_arr = extract_xsec(ca_data_path, False)
si_e_arr, si_xsec_arr = extract_xsec(si_data_path, False)
o_e_arr, o_xsec_arr = extract_xsec(o_data_path, False)
al_e_arr, al_xsec_arr = extract_xsec(al_data_path, False)
fe_e_arr, fe_xsec_arr = extract_xsec(fe_data_path, False)
h_e_arr, h_xsec_arr = extract_xsec(h_data_path, False)
s_e_arr, s_xsec_arr = extract_xsec(s_endf_data_path, True)

wall_thickness = 200 #cm
concrete_density = 2.3 #g/cm3

n_ca = wall_thickness * concrete_density * 0.4565 * 0.602214076 / 40.078
n_si = wall_thickness * concrete_density * 0.08724 * 0.602214076 / 28.085
n_o = wall_thickness * concrete_density * 0.35783 * 0.602214076 / 15.999
n_al = wall_thickness * concrete_density * 0.0447 * 0.602214076 / 26.982
n_fe = wall_thickness * concrete_density * 0.04626 * 0.602214076 / 55.845
n_h = wall_thickness * concrete_density * 0.00083 * 0.602214076 / 1.008
n_s = wall_thickness * concrete_density * 0.00664 * 0.602214076 / 32.06

num_bins = 30
bin_step = (max_energy-min_energy)/num_bins
bin_edges = []
for i in range(0,num_bins+1):
    bin_edges.append(min_energy + i*bin_step)
    if i == num_bins:
        bin_edges.append(max_energy)

#cross section hists
ca_xsec_hist = ROOT.TH1D("ca_xsec_hist", "Ca Cross Section Hist", num_bins, min_energy, max_energy)
si_xsec_hist = ROOT.TH1D("si_xsec_hist", "Si Cross Section Hist", num_bins, min_energy, max_energy)
o_xsec_hist = ROOT.TH1D("o_xsec_hist", "O Cross Section Hist", num_bins, min_energy, max_energy)
al_xsec_hist = ROOT.TH1D("al_xsec_hist", "Al Cross Section Hist", num_bins, min_energy, max_energy)
fe_xsec_hist = ROOT.TH1D("fe_xsec_hist", "Fe Cross Section Hist", num_bins, min_energy, max_energy)
h_xsec_hist = ROOT.TH1D("h_xsec_hist", "H Cross Section Hist", num_bins, min_energy, max_energy)
s_xsec_hist = ROOT.TH1D("s_xsec_hist", "S Cross Section Hist", num_bins, min_energy, max_energy)

#filling cross section hists
fill_xsec_hist(ca_e_arr, ca_xsec_arr, bin_edges, ca_xsec_hist)
fill_xsec_hist(si_e_arr, si_xsec_arr, bin_edges, si_xsec_hist)
fill_xsec_hist(o_e_arr, o_xsec_arr, bin_edges, o_xsec_hist)
fill_xsec_hist(al_e_arr, al_xsec_arr, bin_edges, al_xsec_hist)
fill_xsec_hist(fe_e_arr, fe_xsec_arr, bin_edges, fe_xsec_hist)
fill_xsec_hist(h_e_arr, h_xsec_arr, bin_edges, h_xsec_hist)
fill_xsec_hist(s_e_arr, s_xsec_arr, bin_edges, s_xsec_hist)

#transmission hists
ca_trans_hist = ROOT.TH1D("ca_trans_hist", "Ca Transmission Hist", num_bins, min_energy, max_energy)
si_trans_hist = ROOT.TH1D("si_trans_hist", "Si Transmission Hist", num_bins, min_energy, max_energy)
o_trans_hist = ROOT.TH1D("o_trans_hist", "O Transmission Hist", num_bins, min_energy, max_energy)
al_trans_hist = ROOT.TH1D("al_trans_hist", "Al Transmission Hist", num_bins, min_energy, max_energy)
fe_trans_hist = ROOT.TH1D("fe_trans_hist", "Fe Transmission Hist", num_bins, min_energy, max_energy)
h_trans_hist = ROOT.TH1D("h_trans_hist", "H Transmission Hist", num_bins, min_energy, max_energy)
s_trans_hist = ROOT.TH1D("s_trans_hist", "S Transmission Hist", num_bins, min_energy, max_energy)

#filling transmission hist
fill_trans_hist(ca_xsec_hist, n_ca, ca_trans_hist)
fill_trans_hist(si_xsec_hist, n_si, si_trans_hist)
fill_trans_hist(o_xsec_hist, n_o, o_trans_hist)
fill_trans_hist(al_xsec_hist, n_al, al_trans_hist)
fill_trans_hist(fe_xsec_hist, n_fe, fe_trans_hist)
fill_trans_hist(h_xsec_hist, n_h, h_trans_hist)
fill_trans_hist(s_xsec_hist, n_s, s_trans_hist)

concrete_trans_hist = ROOT.TH1D("concrete_trans_hist", "Concrete Transmission Histogram", num_bins, min_energy, max_energy)

for i in range(0, num_bins):
    bin_content = ca_trans_hist.GetBinContent(i+1) * si_trans_hist.GetBinContent(i+1) * o_trans_hist.GetBinContent(i+1) * al_trans_hist.GetBinContent(i+1) * fe_trans_hist.GetBinContent(i+1) * h_trans_hist.GetBinContent(i+1) * s_trans_hist.GetBinContent(i+1)
    concrete_trans_hist.SetBinContent(i+1, bin_content)

concrete_trans_hist.GetXaxis().SetTitle("Energy (in eV)")
concrete_trans_hist.GetYaxis().SetTitle("Transmission")
concrete_trans_hist.SetStats(0)
concrete_trans_hist.SetLineWidth(2)
concrete_trans_hist.GetYaxis().SetRangeUser(1e-13, 1.)

ca_trans_hist.SetLineColor(1)   
si_trans_hist.SetLineColor(2)   
o_trans_hist.SetLineColor(3)    
al_trans_hist.SetLineColor(6)   
fe_trans_hist.SetLineColor(7)   
h_trans_hist.SetLineColor(8)    
s_trans_hist.SetLineColor(9)    

c1 = ROOT.TCanvas("", "", 800, 500)
c1.SetGridx(1)
c1.SetGridy(1)

l1 = ROOT.TLegend(0.7,0.15,0.85,0.55)
l1.AddEntry(concrete_trans_hist, "Concrete", "l")
l1.AddEntry(ca_trans_hist, "Ca (45.650%)", "l")
l1.AddEntry(si_trans_hist, "Si (8.724%)", "l")
l1.AddEntry(o_trans_hist, "O (35.783%)", "l")
l1.AddEntry(al_trans_hist, "Al (4.470%)", "l")
l1.AddEntry(fe_trans_hist, "Fe (4.626%)", "l")
l1.AddEntry(h_trans_hist, "H (0.664%)", "l")
l1.AddEntry(s_trans_hist, "S (0.083%)", "l")

concrete_trans_hist.Draw()
ca_trans_hist.Draw("SAME")
si_trans_hist.Draw("SAME")
o_trans_hist.Draw("SAME")
al_trans_hist.Draw("SAME")
fe_trans_hist.Draw("SAME")
h_trans_hist.Draw("SAME")
s_trans_hist.Draw("SAME")
l1.Draw()
c1.SetLogy()
c1.Print("../plots/stability_plots/concrete_trans.png")