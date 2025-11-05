import numpy as np
import matplotlib.pyplot as plt
# from Luminosity_function import Phi_bol, log_PhiStar, MStar,gamma_1, gamma_2,log_LStar

def bolometric_correction(logLbol):
    # Bolometric Luminosity correction for UV (1450A) band
    # from Shen et al. 2020 (with Hopkins.)
    # L in unit of L_sun.
    c1, c2 = 1.862, 4.876
    k1, k2 = -0.361, -0.0063
    L = 10 ** logLbol
    Lsun = 3.839 * 1e33
    LbolLband = c1 * (L / (1e10 * Lsun)) ** k1 + c2 * (L / (1e10 * Lsun)) ** k2
    return np.log10(LbolLband)

def dispersion(logLbol):
    # dispersion on Bolometric Correction.
    sig1, sig2, logL0, sig3 = -0.372, 0.405, 42.31, 2.310
    import scipy.special as ss
    err_func = ss.erf((logLbol - logL0) / (np.sqrt(2) * sig3))
    sigma = sig2 + sig1 * ((1/2) + (1/2) * err_func)
    return sigma

# logLbol = 48
# gaus_dist = np.random.normal(loc=bolometric_correction(logLbol),scale=dispersion(logLbol),size=1000)
# hist, bin_edges = np.histogram(gaus_dist,bins=50)
# delta_bin = np.diff(bin_edges)[0]
# print(hist * delta_bin)
# area = np.sum(hist * delta_bin)
# print(area)
# plt.hist(gaus_dist,bins=50,density=True)
# plt.show()

# reproducing Fig.2 in Shen+2020.

# logLbol = np.linspace(38,48,100)
# logLbol_Lband = bolometric_correction(logLbol)
# sigma = dispersion(logLbol)
# Lbol_top = logLbol_Lband + sigma
# Lbol_bottom = logLbol_Lband - sigma
#
# # plt.plot(logLbol,logLbol_Lband,color='black',label='UV')
# # plt.fill_between(logLbol,Lbol_top,Lbol_bottom,alpha=0.2,color='gray',label=r'1$\sigma$')
# plt.plot(logLbol,sigma,label=r'$\sigma$')
# plt.xlim(38,48)
# # plt.ylim(0.5,3)
# plt.ylim(0.0,0.5)
#
# plt.xlabel('log($L_{bol}$[erg/s])')
# # plt.ylabel(r'log($L_{bol}$/$L_{1450\AA}$)')
# plt.ylabel(r'$\sigma_{corr}$')
# plt.legend()
# # sigma = dispersion(logLbol,sig1, sig2, logL0, sig3)
# # plt.plot(logLbol,sigma)
# # plt.scatter(logLbol,sigma)
#
# # plt.savefig('../Figures/L1450_correction.pdf', bbox_inches='tight')
#
# plt.show()

def L1450_wStdev(logLbol):
    # L1450 obtained from Bolometric correction
    # with dispersion (Shen+2020 eq.4 and 5).
    C = bolometric_correction(logLbol)
    logL1450 = logLbol - C
    mean_corr = logL1450
    stdev = dispersion(logLbol)
    corr_upp, corr_low = logLbol + stdev, logLbol - stdev
    return mean_corr, corr_upp, corr_low

def M1450(logLbol):
    # convert L1450 to M1450.
    Lsun = 3.839 * 1e33
    Msun = 4.83
    L_1450 = 10 ** L1450_wStdev(logLbol)[0]
    # print(L_1450)
    return -2.5 * np.log10(L_1450 / Lsun) + Msun

# logLbol = np.linspace(43,47,20)
# M_1450 = M1450(logLbol)
# print(M_1450)

from Kcorr_K2019 import m
# Kulkarni+2019's conversion from M1450 to i-mag for different redshift.

# i_mag = m(M_1450,0.5)
# print(i_mag)