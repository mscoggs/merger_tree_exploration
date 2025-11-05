import numpy as np
import matplotlib.pyplot as plt
from BolometricCorrection import M1450,L1450_wStdev,bolometric_correction,dispersion
from Kcorr_K2019 import m
# Luminosity functions from Shen et al. 2020 (with Hopkins.)
# From Shen+2020 eq(14):
def T_0(x):
    return 1

def T_1(x):
    return x

def T_2(x):
    return 2 * (x ** 2) - 1

def gamma_1(z):
    a0, a1, a2 = 0.8396, -0.2519, 0.0198
    return a0 * T_0(1+z) + a1 * T_1(1+z) + a2 * T_2(1+z)

def gamma_2(z):
    b0, b1, b2 = 2.5432, -1.0528, 1.1284
    z_ref = 2.
    numerator = 2 * b0
    denominator = ((1+z) / (1 + z_ref)) ** b1 + ((1+z) / (1 + z_ref)) ** b2
    return numerator / denominator

def log_LStar(z):
    c0, c1, c2 = 13.0124, -0.5777, 0.4545
    z_ref = 2.
    numerator = 2 * c0
    denominator = ((1+z) / (1 + z_ref)) ** c1 + ((1+z) / (1 + z_ref)) ** c2
    return numerator / denominator

def log_PhiStar(z):
    d0, d1 = -3.5148, -0.4045
    return d0 * T_0(1+z) + d1 * T_1(1+z)

def MStar(z,c0,c1,c2,z_ref):
    Msun = 4.83
    Mstar = -2.5 * log_LStar(z,c0,c1,c2,z_ref) + Msun
    return Mstar

def Phi_Lbol(logLbol, z):
    Lbol = 10 ** logLbol
    # Best-fit 11 global parameters from Table 4 in Shen2020.
    Lsun = 3.839 * 1e33
    # --------------- dn_dlog(Lbol) Calculation -------------- #
    log_phiStar = log_PhiStar(z)
    # print(log_phiStar)
    LStar = Lsun * 10 ** log_LStar(z)

    gamma1 = gamma_1(z)
    gamma2 = gamma_2(z)
    # print(gamma1, gamma2,log_phiStar)

    alpha = - (gamma1 + 1)
    beta = - (gamma2 + 1)

    term_LStar = ((Lbol / LStar) ** gamma1 + (Lbol / LStar) ** gamma2)
    log_dndLbol = - np.log10(term_LStar) + log_phiStar

    return log_dndLbol


# print(Phi_Lbol(46, 1))
# L_bol = np.linspace(42, 49, 200)
# M_1450 = M1450(L_bol)
# print(M_1450)

# z = [0.5,1.5,2.5,3.5,4.5,5.5,6.5]
z = [0.31,0.50,0.72,0.91,1.10,1.30,1.50,1.71,1.98,2.30,2.45,2.55,2.65,
            2.75,2.85,2.95,3.05,3.15,3.25,3.34,3.44,3.88,4.35,4.92,6.00]
z_25 = np.linspace(0.1,6,25)

# for i in range(25):
#     redshift = z[i]
#     imag = m(M_1450,redshift)
#     log_phi = Phi_Lbol(L_bol, redshift)
#     np.savetxt('../Data/LF_imag_Shen2020/lf_imag_z{}.txt'.format('%.2f' % redshift), np.c_[imag,log_phi],fmt='% .3f % .6f')

fig, axs = plt.subplots(nrows=5,ncols=5,figsize=(9,7))

# fig, axs = plt.subplots(nrows=2,ncols=4,figsize=(9,4.5))
# fig.subplots_adjust(wspace=0.02,hspace=0.02)
# fig.delaxes(axs[1,3])
# # fig.delaxes(axs[2,2])

for i, ax in enumerate(fig.axes):

    ######## Shen-Hopkins 2020 LF #####
    # shen = np.loadtxt('../Data/LF_imag_Shen2020/lf_imag_z{}.txt'.format('%.2f'%z[i]))
    # i_mag = shen[:,0]
    # lf_s = shen[:,1]
    k19_intp = np.loadtxt('../Data/K19-DoublePL/LF_new/lf-z{}.txt'.format('%.2f' % z[i]))
    i_mag = k19_intp[:, 1]
    lf_s = k19_intp[:, 2]

    # shen = np.loadtxt('../Data/LF_M1450_Shen2020/lf_m1450_z{}.txt'.format('%.1f' % z[i]))
    # m = shen[:,0]
    # lf_s = shen[:, 1]

    ###### Kulkarni 2019 LF ######
    kul = np.loadtxt('../Data/LF_imag_kulkarni2019/bestfit_lf_{}.txt'.format('%.2f'%z[i]))
    mag, lf = kul[:,0], kul[:,1]
    # lf = np.loadtxt('../Data/LF_imag_kulkarni2019/bestfit-lsst-bins/bestfit_lf_{}.txt'.format('%.1f'%z[i]))
    # mag = np.loadtxt('../Data/LF_imag_kulkarni2019/i-mag/i-mag_z{}.txt'.format('%.1f'%z[i]))

    # mag = np.loadtxt('../Data/LF_imag_kulkarni2019/i-mag/i-mag_z{}.txt'.format('%.1f'%z[i]))
    # lf = np.loadtxt('../Data/LF_imag_kulkarni2019/bestfit-lsst-bins/bestfit_lf_{}.txt'.format('%.1f'%z[i]))
    m1450 = np.linspace(-32,-16,200)
    #
    # print(m1450)
    # ax.plot(i_mag,lf_s,color='red',label='z=Shen+2020')
    ax.plot(i_mag, lf_s, color='blue', linestyle='--',label='K19 (Interp)')
    ax.plot(mag,lf,color='blue',label='Kulkarni+2019')
    # ax.plot(m, lf_s, color='red', label='z=Shen+2020')
    # ax.plot(m1450,lf,color='blue',label='Kulkarni+2019')
    ax.text(23,-10,'<z>={}'.format('%.2f' % (z[i])),fontsize=8)
    ax.set_xlabel('i-mag')
    ax.set_xlim(15,25)
    ax.set_ylim(-12,-4)
    ax.invert_xaxis()
    ax.grid(color='gray', linestyle='--', linewidth=0.5, zorder=0)
    # ax.legend(loc=3,markerscale=2, prop={'size': 6},edgecolor='none')

fig.text(0.06,0.5,'log($\phi$[mag$^{-1}$cMpc$^{-3}$])',rotation='vertical')
plt.legend(prop={'size': 9},edgecolor='none',bbox_to_anchor=(0.5,6.5))

# plt.savefig('../Figures/lf-s20-k19_convolv(i-mag)-fine.pdf', bbox_inches='tight')

# lf = np.loadtxt('../Data/LF_imag_kulkarni2019/bestfit-lsst-bins/bestfit_lf_0.5.txt')
# m1450 = np.linspace(-32, -16, 200)
# plt.plot(m1450, lf, color='blue', label='Kulkarni+2019')
# plt.xlim(-31,-17)
# plt.ylim(-12,-4)
# plt.gca().invert_xaxis()
plt.show()

