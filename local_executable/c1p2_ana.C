/**
 * @file c1p2_ana.C
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2024-09-24
 */

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <map>
#include <vector>

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
#include "TKey.h"
#include "TList.h"

#include "MArEXStyle.C"

Int_t plot_index=0;
Double_t n_C_1p2cm = 0.105;

TH1D* trans_hist_c1p2_ts_ptbc = 0;
TH1D* trans_hist_c1p2_ts_fimg = 0;
TH1D* xsec_hist_c1p2_ts_ptbc = 0;
TH1D* xsec_hist_c1p2_ts_fimg = 0;
TH1D* endf_trans_hist = 0;
TH1D* endf_xsec_hist = 0;

///////////////////////////////////
//Systematics histogram
//backgrounds
TH1D* trans_bkgd_hist_ptbc = 0;
TH1D* trans_bkgd_hist_fimg = 0;
TH1D* xsec_bkgd_hist_ptbc = 0;
TH1D* xsec_bkgd_hist_fimg = 0;

//////////////////////////////// PTBC
Double_t bin_edges_ptbc[] = {1.00000000e-01, 1.12201845e-01, 1.25892541e-01, 1.41253754e-01,
    1.58489319e-01, 1.77827941e-01, 1.99526231e-01, 2.23872114e-01,
    2.51188643e-01, 2.81838293e-01, 3.16227766e-01, 3.54813389e-01,
    3.98107171e-01, 4.46683592e-01, 5.01187234e-01, 5.62341325e-01,
    6.30957344e-01, 7.07945784e-01, 7.94328235e-01, 8.91250938e-01,
    1.00000000e+00, 1.12201845e+00, 1.25892541e+00, 1.41253754e+00,
    1.58489319e+00, 1.77827941e+00, 1.99526231e+00, 2.23872114e+00,
    2.51188643e+00, 2.81838293e+00, 3.16227766e+00, 3.54813389e+00,
    3.98107171e+00, 4.46683592e+00, 5.01187234e+00, 5.62341325e+00,
    6.30957344e+00, 7.07945784e+00, 7.94328235e+00, 8.91250938e+00,
    1.00000000e+01, 1.12201845e+01, 1.25892541e+01, 1.41253754e+01,
    1.58489319e+01, 1.77827941e+01, 1.99526231e+01, 2.23872114e+01,
    2.51188643e+01, 2.81838293e+01, 3.16227766e+01, 3.54813389e+01,
    3.98107171e+01, 4.46683592e+01, 5.01187234e+01, 5.62341325e+01,
    6.30957344e+01, 7.07945784e+01, 7.94328235e+01, 8.91250938e+01,
    1.00000000e+02, 1.12201845e+02, 1.25892541e+02, 1.41253754e+02,
    1.58489319e+02, 1.77827941e+02, 1.99526231e+02, 2.23872114e+02,
    2.51188643e+02, 2.81838293e+02, 3.16227766e+02, 3.54813389e+02,
    3.98107171e+02, 4.46683592e+02, 5.01187234e+02, 5.62341325e+02,
    6.30957344e+02, 7.07945784e+02, 7.94328235e+02, 8.91250938e+02,
    1.00000000e+03, 1.12201845e+03, 1.25892541e+03, 1.41253754e+03,
    1.58489319e+03, 1.77827941e+03, 1.99526231e+03, 2.23872114e+03,
    2.51188643e+03, 2.81838293e+03, 3.16227766e+03, 3.54813389e+03,
    3.98107171e+03, 4.46683592e+03, 5.01187234e+03, 5.62341325e+03,
    6.30957344e+03, 7.07945784e+03, 7.94328235e+03, 8.91250938e+03,
    1.00000000e+04, 1.12201845e+04, 1.25892541e+04, 1.41253754e+04,
    1.58489319e+04, 1.77827941e+04, 1.99526231e+04, 2.23872114e+04,
    2.51188643e+04, 2.81838293e+04, 3.16227766e+04, 3.54813389e+04,
    3.98107171e+04, 4.46683592e+04, 5.01187234e+04, 5.62341325e+04,
    6.30957344e+04, 7.07945784e+04, 7.94328235e+04, 8.91250938e+04,
    1.00000000e+05, 1.12201845e+05, 1.25892541e+05, 1.41253754e+05,
    1.58489319e+05, 1.77827941e+05, 1.99526231e+05, 2.23872114e+05,
    2.51188643e+05, 2.81838293e+05, 3.16227766e+05, 3.54813389e+05,
    3.98107171e+05, 4.46683592e+05, 5.01187234e+05, 5.62341325e+05,
    6.30957344e+05, 7.07945784e+05, 7.94328235e+05, 8.91250938e+05,
    1.00000000e+06, 1.12201845e+06, 1.25892541e+06, 1.41253754e+06,
    1.58489319e+06, 1.77827941e+06, 1.99526231e+06, 2.23872114e+06,
    2.51188643e+06, 2.81838293e+06, 3.16227766e+06, 3.54813389e+06,
    3.98107171e+06, 4.46683592e+06, 5.01187234e+06, 5.62341325e+06,
    6.30957344e+06, 7.07945784e+06, 7.94328235e+06, 8.91250938e+06,
    1.00000000e+07, 1.12201845e+07, 1.25892541e+07, 1.41253754e+07,
    1.58489319e+07, 1.77827941e+07, 1.99526231e+07, 2.23872114e+07,
    2.51188643e+07, 2.81838293e+07, 3.16227766e+07, 3.54813389e+07,
    3.98107171e+07, 4.46683592e+07, 5.01187234e+07, 5.62341325e+07,
    6.30957344e+07, 7.07945784e+07, 7.94328235e+07, 8.91250938e+07,
    1.00000000e+08};

//    , 1.12201845e+08, 1.25892541e+08, 1.41253754e+08,
//    1.58489319e+08, 1.77827941e+08, 1.99526231e+08, 2.23872114e+08,
//    2.51188643e+08, 2.81838293e+08, 3.16227766e+08, 3.54813389e+08,
//    3.98107171e+08, 4.46683592e+08, 5.01187234e+08, 5.62341325e+08,
//    6.30957344e+08, 7.07945784e+08, 7.94328235e+08, 8.91250938e+08,
//    1.00000000e+09

Double_t bin_err_up_ptbc[] = {1.74975385e-04, 1.77789899e-03, 9.71639382e-03, 5.56957747e-03,
    2.69431208e-03, 4.12744606e-03, 1.00495401e-03, 2.45248720e-03,
    2.81277881e-03, 1.34403154e-03, 2.50337706e-03, 1.48388191e-03,
    2.36810348e-03, 1.80907618e-03, 4.28997504e-03, 2.18613267e-03,
    4.72117255e-03, 1.58047925e-03, 3.86693302e-03, 2.68686275e-03,
    1.32161731e-03, 2.46367309e-03, 8.21935907e-03, 1.23130993e-03,
    1.05014331e-02, 1.70658057e-03, 7.27017805e-03, 7.21772294e-03,
    6.69200856e-04, 1.55958526e-03, 1.20836149e-03, 1.65827630e-03,
    1.14569325e-01, 2.90015383e-03, 3.29929218e-03, 3.59862445e-04,
    9.87573958e-04, 1.00189834e-03, 1.95575975e-04, 2.02291816e-04,
    9.63478151e-04, 2.34843618e-04, 9.58569894e-04, 2.86313183e-04,
    9.22101350e-04, 2.06553465e-04, 3.74874447e-04, 5.16159652e-05,
    1.65406810e-04, 9.54722483e-04, 8.91320235e-05, 3.08045127e-04,
    2.17570088e-04, 2.31220078e-04, 1.71652459e-04, 3.07543684e-04,
    4.35111551e-04, 3.60897765e-04, 2.12452386e-04, 5.66901038e-04,
    3.64017238e-04, 1.84380388e-04, 1.45841756e-04, 2.12486645e-04,
    2.31705287e-04, 4.94828979e-04, 7.08063187e-05, 1.93933111e-04,
    1.97142964e-04, 8.29116471e-05, 2.01391949e-04, 6.42524690e-05,
    1.11514450e-04, 2.01948063e-04, 9.98550605e-05, 1.72170732e-04,
    1.31828851e-04, 1.88894054e-04, 1.08855910e-04, 1.07843913e-04,
    2.26549752e-04, 3.33293011e-04, 1.61659234e-04, 2.16744763e-04,
    9.97418349e-05, 6.09398578e-05, 1.25834089e-04, 2.90760530e-04,
    2.70687684e-04, 1.35615236e-04, 9.33381869e-05, 7.29358156e-05,
    2.12653672e-04, 1.02215951e-05, 3.06298044e-04, 5.21523311e-04,
    9.87779370e-05, 4.28191845e-07, 2.50528742e-04, 3.40154235e-04,
    1.41756145e-04, 1.23459159e-04, 4.31182615e-06, 1.85610971e-04,
    1.69432923e-04, 1.82373695e-04, 3.24564091e-05, 1.29538178e-05,
    7.74220279e-06, 1.58485569e-04, 1.90026553e-04, 3.85466059e-04,
    1.22308758e-04, 1.84279145e-04, 5.56658375e-05, 2.06122899e-05,
    3.49506471e-05, 6.64841479e-04, 7.49939882e-04, 1.50512069e-04,
    1.09866898e-04, 1.19606539e-04, 4.40337759e-05, 7.54972331e-05,
    1.79830442e-04, 7.02792755e-05, 6.42393174e-05, 3.12599532e-05,
    5.34998318e-05, 1.90533985e-05, 2.33736622e-05, 1.88211916e-05,
    2.49424727e-05, 7.20867117e-06, 3.96395096e-06, 1.26810714e-05,
    5.64645717e-06, 7.35283219e-06, 1.09353243e-05, 6.19190233e-06,
    8.71933360e-06, 5.72417175e-06, 9.55355899e-06, 7.27569405e-06,
    4.92171027e-06, 7.87184964e-06, 7.44361760e-06, 8.29569112e-06,
    9.16604449e-06, 1.09247109e-05, 1.20869299e-05, 1.44878017e-05,
    1.35292453e-05, 1.28972390e-05, 1.77776500e-05, 1.69761568e-05,
    1.25994188e-05, 1.30911494e-05, 8.97797850e-06, 4.02501867e-06,
    9.35157876e-06, 7.38719985e-06, 6.56560022e-06, 4.13202023e-06,
    1.01763062e-06, 4.43522854e-06, 5.75938723e-07, 6.11262235e-06,
    3.52198404e-06, 4.69817423e-06, 4.69913888e-06, 2.56850110e-06,
    2.12042572e-06, 2.17356540e-06, 2.97497399e-06, 2.17695811e-06,
    3.29274020e-06, 9.21747967e-07, 7.70255985e-07, 1.60965758e-06,
    1.11168331e-06};

//    , 9.53959687e-07, 8.81867841e-07, 2.99813557e-08,
//    4.14260140e-07, 1.35251259e-07, 8.19599534e-07, 1.03032699e-06,
//    5.68604086e-07, 2.89348316e-07, 1.58732444e-06, 1.06449694e-06,
//    7.68250228e-07, 8.05180113e-07, 8.58082643e-07, 7.19911744e-07,
//    6.33704798e-07, 1.20902356e-07, 1.07916729e-06, 6.76301102e-08

Double_t bin_err_low_ptbc[] = {1.20182835e-03, 1.08698661e-02, 5.18117960e-02, 2.74876590e-02,
    1.26554658e-02, 1.83645266e-02, 4.21287600e-03, 9.40441549e-03,
    1.00497835e-02, 4.64644591e-03, 8.67181713e-03, 5.29992799e-03,
    8.72923857e-03, 6.81444422e-03, 1.63956043e-02, 8.49694951e-03,
    1.81442889e-02, 6.04656025e-03, 1.43503067e-02, 9.46519853e-03,
    4.22415558e-03, 8.09713328e-03, 4.10339929e-02, 7.98108564e-03,
    6.74228015e-02, 1.02140496e-02, 3.45399305e-02, 5.90529740e-02,
    5.79941265e-03, 5.27717182e-03, 4.10156228e-03, 5.51545920e-03,
    7.38866274e-01, 1.79333599e-02, 1.39088395e-02, 1.13492572e-03,
    2.99525086e-03, 3.83448891e-03, 5.53181265e-04, 5.85989578e-04,
    3.10987163e-03, 6.71934333e-04, 2.79832956e-03, 8.56647380e-04,
    2.75792274e-03, 5.82680392e-04, 1.09181972e-03, 1.47802969e-04,
    4.74638372e-04, 2.95163067e-03, 2.50677911e-04, 8.86361597e-04,
    6.24257782e-04, 6.61545122e-04, 4.83407044e-04, 8.72657634e-04,
    1.25470223e-03, 1.03134863e-03, 6.04954654e-04, 1.62164536e-03,
    1.05032493e-03, 5.24449784e-04, 4.14219266e-04, 6.06003198e-04,
    6.57166557e-04, 1.41170027e-03, 2.00627829e-04, 5.48968026e-04,
    5.56069911e-04, 2.37502370e-04, 5.74288957e-04, 1.84086222e-04,
    3.16720543e-04, 5.71792682e-04, 2.82181569e-04, 4.86577942e-04,
    3.73688705e-04, 5.34381494e-04, 3.09685392e-04, 3.06406573e-04,
    6.43520016e-04, 9.42353051e-04, 4.57864755e-04, 6.15174379e-04,
    2.82434543e-04, 1.72630920e-04, 3.57488839e-04, 8.24607198e-04,
    7.67728994e-04, 3.84282202e-04, 2.64226767e-04, 2.06530572e-04,
    6.00752223e-04, 2.89240879e-05, 8.65187084e-04, 1.47535079e-03,
    2.79308736e-04, 1.21130919e-06, 7.07265110e-04, 9.59612555e-04,
    4.00126307e-04, 3.48200484e-04, 1.21441838e-05, 5.22547386e-04,
    4.76910542e-04, 5.12222447e-04, 9.12638357e-05, 3.63752288e-05,
    2.17191260e-05, 4.44145312e-04, 5.35644569e-04, 1.08634359e-03,
    3.42732960e-04, 5.15849346e-04, 1.55742577e-04, 5.76303572e-05,
    9.76811793e-05, 1.85674796e-03, 2.09922577e-03, 4.20905656e-04,
    3.06586298e-04, 3.33693558e-04, 1.22758673e-04, 2.10657917e-04,
    5.01277210e-04, 1.95808186e-04, 1.78981457e-04, 8.70558478e-05,
    1.48964214e-04, 5.30513548e-05, 6.50614177e-05, 5.23866175e-05,
    6.94370051e-05, 2.00604505e-05, 1.10296353e-05, 3.52829167e-05,
    1.57093540e-05, 2.04560894e-05, 3.04220717e-05, 1.72255767e-05,
    2.42556340e-05, 1.59229465e-05, 2.65747374e-05, 2.02379294e-05,
    1.36901727e-05, 2.18958336e-05, 2.07049261e-05, 2.30749702e-05,
    2.54968108e-05, 3.03898128e-05, 3.36250053e-05, 4.03068003e-05,
    3.76405217e-05, 3.58831215e-05, 4.94622866e-05, 4.72324561e-05,
    3.50518061e-05, 3.64183139e-05, 2.49758519e-05, 1.11973811e-05,
    2.60153912e-05, 2.05508896e-05, 1.82645431e-05, 1.14943341e-05,
    2.83085002e-06, 1.23373978e-05, 1.60205724e-06, 1.70025563e-05,
    9.79657610e-06, 1.30679261e-05, 1.30705150e-05, 7.14413771e-06,
    5.89779252e-06, 6.04561278e-06, 8.27457286e-06, 6.05497964e-06,
    9.15831739e-06, 2.56370772e-06, 2.14235068e-06, 4.47698096e-06,
    3.09197540e-06};

//    , 2.65326521e-06, 2.45273749e-06, 8.33875616e-08,
//    1.15218442e-06, 3.76174938e-07, 2.27954448e-06, 2.86563006e-06,
//    1.58145199e-06, 8.04761240e-07, 4.41478266e-06, 2.96064260e-06,
//    2.13670267e-06, 2.23941077e-06, 2.38653693e-06, 2.00224687e-06,
//    1.76248027e-06, 3.36257369e-07, 3.00141680e-06, 1.88096616e-07

//////////////////////////////// FIMG
Double_t bin_edges_fimg[] = {1.00000000e-01, 1.12201845e-01, 1.25892541e-01, 1.41253754e-01,
    1.58489319e-01, 1.77827941e-01, 1.99526231e-01, 2.23872114e-01,
    2.51188643e-01, 2.81838293e-01, 3.16227766e-01, 3.54813389e-01,
    3.98107171e-01, 4.46683592e-01, 5.01187234e-01, 5.62341325e-01,
    6.30957344e-01, 7.07945784e-01, 7.94328235e-01, 8.91250938e-01,
    1.00000000e+00, 1.12201845e+00, 1.25892541e+00, 1.41253754e+00,
    1.58489319e+00, 1.77827941e+00, 1.99526231e+00, 2.23872114e+00,
    2.51188643e+00, 2.81838293e+00, 3.16227766e+00, 3.54813389e+00,
    3.98107171e+00, 4.46683592e+00, 5.01187234e+00, 5.62341325e+00,
    6.30957344e+00, 7.07945784e+00, 7.94328235e+00, 8.91250938e+00,
    1.00000000e+01, 1.12201845e+01, 1.25892541e+01, 1.41253754e+01,
    1.58489319e+01, 1.77827941e+01, 1.99526231e+01, 2.23872114e+01,
    2.51188643e+01, 2.81838293e+01, 3.16227766e+01, 3.54813389e+01,
    3.98107171e+01, 4.46683592e+01, 5.01187234e+01, 5.62341325e+01,
    6.30957344e+01, 7.07945784e+01, 7.94328235e+01, 8.91250938e+01,
    1.00000000e+02, 1.12201845e+02, 1.25892541e+02, 1.41253754e+02,
    1.58489319e+02, 1.77827941e+02, 1.99526231e+02, 2.23872114e+02,
    2.51188643e+02, 2.81838293e+02, 3.16227766e+02, 3.54813389e+02,
    3.98107171e+02, 4.46683592e+02, 5.01187234e+02, 5.62341325e+02,
    6.30957344e+02, 7.07945784e+02, 7.94328235e+02, 8.91250938e+02,
    1.00000000e+03, 1.12201845e+03, 1.25892541e+03, 1.41253754e+03,
    1.58489319e+03, 1.77827941e+03, 1.99526231e+03, 2.23872114e+03,
    2.51188643e+03, 2.81838293e+03, 3.16227766e+03, 3.54813389e+03,
    3.98107171e+03, 4.46683592e+03, 5.01187234e+03, 5.62341325e+03,
    6.30957344e+03, 7.07945784e+03, 7.94328235e+03, 8.91250938e+03,
    1.00000000e+04, 1.12201845e+04, 1.25892541e+04, 1.41253754e+04,
    1.58489319e+04, 1.77827941e+04, 1.99526231e+04, 2.23872114e+04,
    2.51188643e+04, 2.81838293e+04, 3.16227766e+04, 3.54813389e+04,
    3.98107171e+04, 4.46683592e+04, 5.01187234e+04, 5.62341325e+04,
    6.30957344e+04, 7.07945784e+04, 7.94328235e+04, 8.91250938e+04,
    1.00000000e+05, 1.12201845e+05, 1.25892541e+05, 1.41253754e+05,
    1.58489319e+05, 1.77827941e+05, 1.99526231e+05, 2.23872114e+05,
    2.51188643e+05, 2.81838293e+05, 3.16227766e+05, 3.54813389e+05,
    3.98107171e+05, 4.46683592e+05, 5.01187234e+05, 5.62341325e+05,
    6.30957344e+05, 7.07945784e+05, 7.94328235e+05, 8.91250938e+05,
    1.00000000e+06};

Double_t bin_err_low_fimg[] = {7.44330791e-03, 5.62878727e-03, 5.26444842e-03, 6.75438194e-03,
    3.12791401e-03, 2.90086201e-03, 4.48315669e-03, 3.91170255e-03,
    2.14444987e-03, 2.30664008e-03, 1.13053623e-03, 1.64090160e-03,
    2.05720995e-03, 2.37283602e-03, 2.18010796e-03, 1.93648164e-03,
    1.38560311e-03, 1.35092112e-03, 1.17317283e-03, 1.54036094e-03,
    1.04081453e-03, 1.20651728e-03, 6.39178853e-04, 1.89387109e-05,
    1.09842063e-03, 1.47494875e-03, 8.74458268e-04, 1.23845742e-03,
    5.47558594e-04, 3.84481131e-04, 1.34714309e-03, 8.30113399e-04,
    1.13366202e-03, 9.35199597e-04, 9.66706323e-04, 5.71441755e-04,
    9.54208129e-04, 5.68995270e-04, 8.47901032e-04, 6.16907593e-04,
    1.23996441e-03, 8.67421216e-04, 1.03594350e-03, 2.73848297e-04,
    6.29219978e-04, 1.24356944e-03, 7.70711124e-04, 6.47578149e-04,
    1.23030498e-03, 5.73690005e-04, 6.05492436e-04, 9.83491158e-04,
    1.07551893e-03, 8.16745635e-04, 9.64682134e-04, 9.76865872e-04,
    9.73238967e-04, 3.76589005e-04, 3.29487648e-04, 7.31992660e-04,
    8.61220272e-04, 2.87126833e-04, 8.23604322e-04, 1.36343406e-04,
    8.80332073e-04, 6.43982274e-04, 1.00180520e-03, 4.70811949e-04,
    8.31391457e-04, 3.12675253e-04, 1.72144642e-04, 8.69067034e-04,
    1.33347248e-03, 4.17241134e-05, 5.42569418e-05, 7.48816671e-04,
    1.56404135e-04, 5.04689944e-04, 2.43135019e-05, 1.14578072e-03,
    6.92957834e-05, 5.01273617e-04, 8.90840463e-04, 2.20257346e-04,
    1.73387700e-04, 2.79023653e-04, 5.27487403e-04, 1.05466487e-04,
    8.39149235e-04, 2.10129028e-04, 1.06399849e-03, 6.63144122e-04,
    9.09787300e-04, 2.50709881e-04, 8.60593813e-04, 1.74378632e-04,
    1.93325981e-03, 1.32451436e-03, 7.30178226e-04, 5.04132632e-05,
    1.60487341e-05, 1.02149804e-03, 3.42765165e-04, 9.28974046e-04,
    4.04032425e-04, 2.68252863e-04, 4.93231507e-04, 2.79072213e-04,
    1.02956254e-04, 8.22532754e-04, 1.07982202e-03, 8.80689818e-04,
    5.19607262e-04, 1.31751188e-03, 1.86395227e-04, 8.33640188e-04,
    3.25817010e-04, 3.48080933e-03, 4.47408307e-03, 1.73881958e-03,
    6.04278650e-04, 3.67749595e-04, 1.37849805e-04, 5.90002759e-05,
    6.90638312e-04, 5.08368916e-04, 1.22614701e-04, 2.35397438e-04,
    1.88284630e-04, 3.07398958e-04, 3.13007855e-04, 4.55336069e-05,
    7.13025457e-05, 1.45768354e-04, 1.38653151e-04, 2.26687761e-04,
    1.52132690e-04, 3.35091105e-04, 7.77289033e-06, 3.90340452e-04,
    3.61172261e-04, 5.25950513e-04, 2.33374148e-04, 2.08962209e-04};


TH1D* Get_1D_hist(const char *fname, const char *hist_name){
    
    TFile* hist_file = TFile::Open(fname, "READ");
    TH1D* hist_new = (TH1D*)hist_file->Get(hist_name);

    return hist_new;
}

TH1D* subtractHists_theoryEval(TH1D* h1, TH1D* h2){
    TH1D *h_sub = (TH1D*)h1->Clone();
    Int_t num_bins = h1->GetNbinsX();
    for (int i = 0; i < num_bins; i++){
        Double_t h1_bin_content = h1->GetBinContent(i+1);
        Double_t h2_bin_content = h2->GetBinContent(i+1);
        if (h1_bin_content == 0 || h2_bin_content == 0)
        {
            h_sub->SetBinContent(i+1, 0);
            h_sub->SetBinError(i+1, 0);
            continue;
        }
        Double_t h1_err = h1->GetBinError(i+1);
        h_sub->SetBinContent(i+1, h1_bin_content - h2_bin_content);
        h_sub->SetBinError(i+1, h1_err);
    }
    return h_sub;
}


void fill_bkgd_sys(){

    Int_t num_bins_ptbc = static_cast<Int_t>(sizeof(bin_edges_ptbc)/sizeof(*bin_edges_ptbc) - 1);
    trans_bkgd_hist_ptbc = new TH1D("trans_bkgd_hist_ptbc", "Background errors transmission - ptbc", num_bins_ptbc, bin_edges_ptbc);
    xsec_bkgd_hist_ptbc = new TH1D("xsec_bkgd_hist_ptbc", "Background errors cross section - ptbc", num_bins_ptbc, bin_edges_ptbc);

    for (Int_t i = 1; i < num_bins_ptbc+1; i++)
    {
        Double_t trans_bin_lim_up = trans_hist_c1p2_ts_ptbc->GetBinContent(20+i) + bin_err_up_ptbc[i-1];
        Double_t trans_bin_lim_low = trans_hist_c1p2_ts_ptbc->GetBinContent(20+i) - bin_err_low_ptbc[i-1];
        Double_t trans_bin_content = (trans_bin_lim_low + trans_bin_lim_up)/2;
        Double_t trans_bin_error = (trans_bin_lim_up - trans_bin_lim_low)/2;
        trans_bkgd_hist_ptbc->SetBinContent(i, trans_bin_content);
        trans_bkgd_hist_ptbc->SetBinError(i, trans_bin_error);

        Double_t xsec_bin_lim_low = 0;
        Double_t xsec_bin_lim_up = 0;
        if (trans_bin_lim_up > 0)
        {
            xsec_bin_lim_low = - std::log(trans_bin_lim_up) / n_C_1p2cm;
        }
        if (trans_bin_lim_low > 0)
        {
            xsec_bin_lim_up = - std::log(trans_bin_lim_low) / n_C_1p2cm;
        }
        // cout << "( " << trans_bin_lim_low << ", " << xsec_bin_lim_up << " )" << endl;
        Double_t xsec_bin_content = (xsec_bin_lim_up + xsec_bin_lim_low)/2;
        Double_t xsec_bin_err = (xsec_bin_lim_up - xsec_bin_lim_low)/2;
        xsec_bkgd_hist_ptbc->SetBinContent(i, xsec_bin_content);
        xsec_bkgd_hist_ptbc->SetBinError(i, xsec_bin_err);
    }

    Int_t num_bins_fimg = static_cast<Int_t>(sizeof(bin_edges_fimg)/sizeof(*bin_edges_fimg) - 1);
    trans_bkgd_hist_fimg = new TH1D("trans_bkgd_hist_fimg", "Background errors transmission - fimg", num_bins_fimg, bin_edges_fimg);
    xsec_bkgd_hist_fimg = new TH1D("xsec_bkgd_hist_fimg", "Background errors cross section - fimg", num_bins_fimg, bin_edges_fimg);

    for (Int_t i = 1; i < num_bins_fimg+1; i++)
    {
        Double_t trans_bin_lim_up = trans_hist_c1p2_ts_fimg->GetBinContent(20+i);
        Double_t trans_bin_lim_low = trans_hist_c1p2_ts_fimg->GetBinContent(20+i) - bin_err_low_fimg[i-1];
        Double_t trans_bin_content = (trans_bin_lim_low + trans_bin_lim_up)/2;
        Double_t trans_bin_error = (trans_bin_lim_up - trans_bin_lim_low)/2;
        trans_bkgd_hist_fimg->SetBinContent(i, trans_bin_content);
        trans_bkgd_hist_fimg->SetBinError(i, trans_bin_error);

        Double_t xsec_bin_lim_low = 0;
        Double_t xsec_bin_lim_up = 0;
        if (trans_bin_lim_up != 0)
        {
            xsec_bin_lim_low = - std::log(trans_bin_lim_up) / n_C_1p2cm;
        }
        if (trans_bin_lim_low != 0)
        {
            xsec_bin_lim_up = - std::log(trans_bin_lim_low) / n_C_1p2cm;
        }
        Double_t xsec_bin_content = (xsec_bin_lim_up + xsec_bin_lim_low)/2;
        Double_t xsec_bin_err = (xsec_bin_lim_up - xsec_bin_lim_low)/2;
        xsec_bkgd_hist_fimg->SetBinContent(i, xsec_bin_content);
        xsec_bkgd_hist_fimg->SetBinError(i, xsec_bin_err);
    }
}

void plot_subtraction_plots(TH1D* hist_1, const char* x_label, const char* legend_1, TH1D* hist_2, const char* legend_2, TH1D* sys_hist, const char* legend_sys_hist, TH1D* subtraction_plot, TH1D* subtraction_plot_sys, const char* plot_title, const char* output_file_name, Double_t legend_coord[4], Double_t max_e){

    TLine* zero_line = new TLine(1e-1, 0., max_e, 0.);
    zero_line->SetLineWidth(2);
    zero_line->SetLineColor(1);
    zero_line->SetLineStyle(2);

    //Plotting
    SetMArEXStyle();
    
    gStyle->SetStatX(0.27);
    gStyle->SetStatY(0.9);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.17);

    gStyle->SetCanvasDefW(800); //600
    gStyle->SetCanvasDefH(400); //500
    gStyle->SetPadRightMargin(0.05);

    TCanvas *canv = new TCanvas(Form("c%d", plot_index)," ");
    canv->cd();
    canv->Draw();

    TPad *p_upper = new TPad(Form("p_upper_%d", plot_index), Form("p_upper_%d", plot_index), 0., 0.4, 1., 1.);
    p_upper->SetFillColor(kWhite);
    p_upper->SetBottomMargin(0.00001);
    p_upper->SetBorderMode(0);
    p_upper->Draw();
    p_upper->cd();

    cout << "plotting sys_hist " << endl;
    sys_hist->GetYaxis()->SetTitle(x_label);
    sys_hist->GetYaxis()->SetLabelSize(0.08);
    sys_hist->GetYaxis()->SetTitleSize(0.08);
    sys_hist->GetYaxis()->SetTitleOffset(0.5);
    // X Axis
    sys_hist->GetXaxis()->SetLabelOffset(999);
    sys_hist->GetXaxis()->SetLabelSize(0);
    sys_hist->SetTitle(plot_title);
    sys_hist->SetLineWidth(1);
    sys_hist->SetFillColor(kCyan);
    sys_hist->Draw("e2");
    sys_hist->SetStats(0);
    gPad->SetGrid();
    gPad->SetLogx();
    // gPad->SetLogy();

    cout << "plotting hist_1 " << endl;
    hist_1->SetTitle("");
    hist_1->Draw("SAME");

    cout << "plotting hist_2 " << endl;
    hist_2->SetLineWidth(2);
    hist_2->SetLineColor(2);
    // hist_2->GetXaxis()->SetRangeUser(1e-1, max_e);
    hist_2->Draw("][SAME");

    TLegend *sub_legend = new TLegend(legend_coord[0], legend_coord[1], legend_coord[2], legend_coord[3]);
    sub_legend->AddEntry(sys_hist, legend_sys_hist, "f");
    sub_legend->AddEntry(hist_1, legend_1, "l");
    sub_legend->AddEntry(hist_2, legend_2, "l");
    sub_legend->Draw();

    canv->cd(0);
    TPad *p_lower = new TPad(Form("p_upper_%d", plot_index), Form("p_upper_%d", plot_index), 0., 0., 1., 0.4);
    p_lower->SetFillColor(kWhite);
    p_lower->SetTopMargin(0.00001);
    p_lower->SetBottomMargin(0.3);
    p_lower->SetBorderMode(0);
    p_lower->Draw();
    p_lower->cd();

    cout << "plotting subtraction_plot_sys " << endl;
    subtraction_plot_sys->SetTitle("");
    // X Axis
    subtraction_plot_sys->GetXaxis()->SetTitle("Energy (in eV)");
    subtraction_plot_sys->GetXaxis()->SetTitleSize(0.12);
    subtraction_plot_sys->GetXaxis()->SetLabelOffset(0.01);
    subtraction_plot_sys->GetXaxis()->SetTitleOffset(1.13); //decrease to move up
    subtraction_plot_sys->GetXaxis()->SetLabelSize(0.12);
    // Y Axis
    subtraction_plot_sys->GetYaxis()->SetTitle("Difference");
    subtraction_plot_sys->GetYaxis()->SetTitleSize(0.12);
    subtraction_plot_sys->GetYaxis()->SetTitleOffset(0.32);
    subtraction_plot_sys->GetYaxis()->SetLabelSize(0.12);
    subtraction_plot_sys->GetYaxis()->SetNdivisions(5);
    // TGaxis::SetExponentOffset(-0.06, -0.8, "y"); // X and Y offset for Y axis
    // subtraction_plot_sys->GetXaxis()->SetRangeUser(1e-1, max_e);
    // gPad->SetGrid();
    subtraction_plot_sys->SetFillColor(kCyan);
    subtraction_plot_sys->Draw("e2");
    gPad->SetLogx();

    cout << "plotting subtraction_plot " << endl;
    subtraction_plot->Draw("SAME");
    cout << "plotting zero_line " << endl;
    zero_line->Draw("SAME");
    cout << "printing canvas " << endl;
    // canv->Print(output_file_name);

    plot_index++;
    cout << "increased plot index " << endl;
    return;

}

void c1p2_ana(){

    trans_hist_c1p2_ts_ptbc = Get_1D_hist("../rootFiles/trans_xsec_hists_c1p2_ts_20bpd.root", "transmission_hist_e_PTBC");
    trans_hist_c1p2_ts_fimg = Get_1D_hist("../rootFiles/trans_xsec_hists_c1p2_ts_20bpd.root", "transmission_hist_e_FIMG");
    xsec_hist_c1p2_ts_ptbc = Get_1D_hist("../rootFiles/trans_xsec_hists_c1p2_ts_20bpd.root", "cross_section_hist_e_PTBC");
    xsec_hist_c1p2_ts_fimg = Get_1D_hist("../rootFiles/trans_xsec_hists_c1p2_ts_20bpd.root", "cross_section_hist_e_FIMG");
    endf_trans_hist = Get_1D_hist("../rootFiles/trans_xsec_hists_c1p2_ts_20bpd.root", "endf_trans_hist");
    endf_xsec_hist = Get_1D_hist("../rootFiles/trans_xsec_hists_c1p2_ts_20bpd.root", "endf_xsec_hist");

    fill_bkgd_sys();

    TH1D* trans_subtraction_hist_ptbc = subtractHists_theoryEval(trans_hist_c1p2_ts_ptbc, endf_trans_hist);
    TH1D* trans_subtraction_hist_fimg = subtractHists_theoryEval(trans_hist_c1p2_ts_fimg, endf_trans_hist);
    TH1D* xsec_subtraction_hist_ptbc = subtractHists_theoryEval(xsec_hist_c1p2_ts_ptbc, endf_xsec_hist);
    TH1D* xsec_subtraction_hist_fimg = subtractHists_theoryEval(xsec_hist_c1p2_ts_fimg, endf_xsec_hist);

    Int_t num_bins_ptbc = trans_bkgd_hist_ptbc->GetNbinsX();
    Int_t num_bins_fimg = trans_bkgd_hist_fimg->GetNbinsX();

    TH1D* trans_subtraction_hist_ptbc_sys = (TH1D*)trans_bkgd_hist_ptbc->Clone("trans_subtraction_hist_ptbc_sys");
    TH1D* trans_subtraction_hist_fimg_sys = (TH1D*)trans_bkgd_hist_fimg->Clone("trans_subtraction_hist_fimg_sys");
    TH1D* xsec_subtraction_hist_ptbc_sys = (TH1D*)xsec_bkgd_hist_ptbc->Clone("xsec_subtraction_hist_ptbc_sys");
    TH1D* xsec_subtraction_hist_fimg_sys = (TH1D*)xsec_bkgd_hist_fimg->Clone("xsec_subtraction_hist_fimg_sys");

    for (Int_t i = 1; i < num_bins_ptbc+1; i++)
    {
        Double_t trans_val = trans_hist_c1p2_ts_ptbc->GetBinContent(20+i);
        Double_t trans_bkgd_err_up = trans_bkgd_hist_ptbc->GetBinError(i) + trans_bkgd_hist_ptbc->GetBinContent(i) - trans_val;
        Double_t trans_bkgd_err_low = trans_bkgd_hist_ptbc->GetBinError(i) - trans_bkgd_hist_ptbc->GetBinContent(i) + trans_val;

        Double_t xsec_val = xsec_hist_c1p2_ts_ptbc->GetBinContent(20+i);
        Double_t xsec_bkgd_err_up = xsec_bkgd_hist_ptbc->GetBinError(i) + xsec_bkgd_hist_ptbc->GetBinContent(i) - xsec_val;
        Double_t xsec_bkgd_err_low = xsec_bkgd_hist_ptbc->GetBinError(i) - xsec_bkgd_hist_ptbc->GetBinContent(i) + xsec_val;

        trans_subtraction_hist_ptbc_sys->SetBinContent(i, trans_subtraction_hist_ptbc->GetBinContent(20+i) + (trans_bkgd_err_up - trans_bkgd_err_low)/2);
        trans_subtraction_hist_ptbc_sys->SetBinError(i, (trans_bkgd_err_up + trans_bkgd_err_low)/2);

        xsec_subtraction_hist_ptbc_sys->SetBinContent(i, xsec_subtraction_hist_ptbc->GetBinContent(20+i) + (xsec_bkgd_err_up - xsec_bkgd_err_low)/2);
        xsec_subtraction_hist_ptbc_sys->SetBinError(i, (xsec_bkgd_err_up + xsec_bkgd_err_low)/2); 
    }

    for (Int_t i = 1; i < num_bins_fimg+1; i++)
    {
        Double_t trans_val = trans_hist_c1p2_ts_fimg->GetBinContent(20+i);
        Double_t trans_bkgd_err_up = trans_bkgd_hist_fimg->GetBinError(i) + trans_bkgd_hist_fimg->GetBinContent(i) - trans_val;
        Double_t trans_bkgd_err_low = trans_bkgd_hist_fimg->GetBinError(i) - trans_bkgd_hist_fimg->GetBinContent(i) + trans_val;

        Double_t xsec_val = xsec_hist_c1p2_ts_fimg->GetBinContent(20+i);
        Double_t xsec_bkgd_err_up = xsec_bkgd_hist_fimg->GetBinError(i) + xsec_bkgd_hist_fimg->GetBinContent(i) - xsec_val;
        Double_t xsec_bkgd_err_low = xsec_bkgd_hist_fimg->GetBinError(i) - xsec_bkgd_hist_fimg->GetBinContent(i) + xsec_val;

        trans_subtraction_hist_fimg_sys->SetBinContent(i, trans_subtraction_hist_fimg->GetBinContent(20+i) + (trans_bkgd_err_up - trans_bkgd_err_low)/2);
        trans_subtraction_hist_fimg_sys->SetBinError(i, (trans_bkgd_err_up + trans_bkgd_err_low)/2);

        xsec_subtraction_hist_fimg_sys->SetBinContent(i, xsec_subtraction_hist_fimg->GetBinContent(20+i) + (xsec_bkgd_err_up - xsec_bkgd_err_low)/2);
        xsec_subtraction_hist_fimg_sys->SetBinError(i, (xsec_bkgd_err_up + xsec_bkgd_err_low)/2); 
    }

    trans_bkgd_hist_fimg->GetXaxis()->SetRangeUser(1e-1, 1e6);
    trans_bkgd_hist_ptbc->GetXaxis()->SetRangeUser(1e-1, 1e8);
    xsec_bkgd_hist_fimg->GetXaxis()->SetRangeUser(1e-1, 1e6);
    xsec_bkgd_hist_ptbc->GetXaxis()->SetRangeUser(1e-1, 1e8);

    trans_bkgd_hist_ptbc->GetYaxis()->SetRangeUser(0.42,1.05);
    trans_bkgd_hist_fimg->GetYaxis()->SetRangeUser(0.42,0.95);
    // xsec_bkgd_hist_ptbc->GetYaxis()->SetRangeUser(-0.6,7.1);
    xsec_bkgd_hist_fimg->GetYaxis()->SetRangeUser(0.1,7.1);

    trans_hist_c1p2_ts_fimg->GetXaxis()->SetRangeUser(1e-1, 1e6);
    trans_hist_c1p2_ts_ptbc->GetXaxis()->SetRangeUser(1e-1, 1e8);
    xsec_hist_c1p2_ts_fimg->GetXaxis()->SetRangeUser(1e-1, 1e6);
    xsec_hist_c1p2_ts_ptbc->GetXaxis()->SetRangeUser(1e-1, 1e8);

    trans_hist_c1p2_ts_ptbc->GetYaxis()->SetRangeUser(0.42,1.05);
    trans_hist_c1p2_ts_fimg->GetYaxis()->SetRangeUser(0.42,0.95);
    xsec_hist_c1p2_ts_ptbc->GetYaxis()->SetRangeUser(-0.6,7.1);
    xsec_hist_c1p2_ts_fimg->GetYaxis()->SetRangeUser(0.1,7.1);

    trans_subtraction_hist_ptbc_sys->GetXaxis()->SetRangeUser(1e-1, 1e8);
    trans_subtraction_hist_fimg_sys->GetXaxis()->SetRangeUser(1e-1, 1e6);
    xsec_subtraction_hist_ptbc_sys->GetXaxis()->SetRangeUser(1e-1, 1e8);
    xsec_subtraction_hist_fimg_sys->GetXaxis()->SetRangeUser(1e-1, 1e6);

    trans_subtraction_hist_ptbc_sys->GetYaxis()->SetRangeUser(-0.19,0.19);
    trans_subtraction_hist_fimg_sys->GetYaxis()->SetRangeUser(-0.19,0.19);
    xsec_subtraction_hist_ptbc_sys->GetYaxis()->SetRangeUser(-1.9,1.9);
    xsec_subtraction_hist_fimg_sys->GetYaxis()->SetRangeUser(-1.9,1.9);

    trans_subtraction_hist_ptbc->GetXaxis()->SetRangeUser(1e-1, 1e8);
    trans_subtraction_hist_fimg->GetXaxis()->SetRangeUser(1e-1, 1e6);
    xsec_subtraction_hist_ptbc->GetXaxis()->SetRangeUser(1e-1, 1e8);
    xsec_subtraction_hist_fimg->GetXaxis()->SetRangeUser(1e-1, 1e6);

    trans_subtraction_hist_ptbc->GetYaxis()->SetRangeUser(-0.19,0.19);
    trans_subtraction_hist_fimg->GetYaxis()->SetRangeUser(-0.19,0.19);
    xsec_subtraction_hist_ptbc->GetYaxis()->SetRangeUser(-1.9,1.9);
    xsec_subtraction_hist_fimg->GetYaxis()->SetRangeUser(-1.9,1.9);

    //calculating subtraction plot sys


    Double_t xsec_legend_coord[4] = {0.15,0.1,0.45,0.35};
    Double_t trans_legend_coord[4] = {0.15,0.55,0.45,0.8};

    plot_subtraction_plots(trans_hist_c1p2_ts_fimg, "Transmission", "1.2cm C", endf_trans_hist, "ENDF-VIII", trans_bkgd_hist_fimg, "Systematic Error", trans_subtraction_hist_fimg, trans_subtraction_hist_fimg_sys, "Transmission - C 1.2cm - Micromegas", "../plots/results_plots/trans_c1p2_ts_fimg_20bpd.png", trans_legend_coord, 1e6);

    plot_subtraction_plots(trans_hist_c1p2_ts_ptbc, "Transmission", "1.2cm C", endf_trans_hist, "ENDF-VIII", trans_bkgd_hist_ptbc, "Systematic Error", trans_subtraction_hist_ptbc, trans_subtraction_hist_ptbc_sys, "Transmission - C 1.2cm - Fission Chamber", "../plots/results_plots/trans_c1p2_ts_ptbc_20bpd.png", trans_legend_coord, 1e8);

    plot_subtraction_plots(xsec_hist_c1p2_ts_fimg, "Cross Section", "1.2cm C", endf_xsec_hist, "ENDF-VIII", xsec_bkgd_hist_fimg, "Systematic Error", xsec_subtraction_hist_fimg, xsec_subtraction_hist_fimg_sys, "Cross Section - C 1.2cm - Micromegas", "../plots/results_plots/xsec_c1p2_ts_fimg_20bpd.png", xsec_legend_coord, 1e6);

    plot_subtraction_plots(xsec_hist_c1p2_ts_ptbc, "Cross Section", "1.2cm C", endf_xsec_hist, "ENDF-VIII", xsec_bkgd_hist_ptbc, "Systematic Error",  xsec_subtraction_hist_ptbc, xsec_subtraction_hist_ptbc_sys, "Cross Section - C 1.2cm - Fission Chamber", "../plots/results_plots/xsec_c1p2_ts_ptbc_20bpd.png", xsec_legend_coord, 1e8);

    // cout << xsec_subtraction_hist_ptbc_sys->GetBinContent(100) << endl;

    // TCanvas *canv = new TCanvas(Form("c%d", plot_index)," ");
    // canv->cd();
    // canv->Draw();
    // xsec_bkgd_hist_ptbc->Draw();
    // // subtraction_plot_sys->SetFillColor(kCyan);
    // // subtraction_plot_sys->Draw("e2");
    // gPad->SetLogx();
}