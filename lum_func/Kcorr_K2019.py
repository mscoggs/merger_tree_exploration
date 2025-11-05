# K-correction used for g-, i-, and z-band in Kulkarni+2019 (Fig.2)
# Use this file to make i-band correction, with K2019 eq. (1)
# cosmological parameters: H0 = 70, Omega_m = 0.3, Omega_lambda = 0.7
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.cosmology import WMAP9 as cosmo
import scipy.interpolate as intp

kcorri_l15 = np.loadtxt('../Data/kcorr_imag_kulkarni2019/kcorri_l15.txt')
kcorri_t02 = np.loadtxt('../Data/kcorr_imag_kulkarni2019/kcorri_t02.txt')
kcorri_v01 = np.loadtxt('../Data/kcorr_imag_kulkarni2019/kcorri_v01.txt')
z_1, kcorri_1 = kcorri_l15[:,1], kcorri_l15[:,2]
z_2, kcorri_2 = kcorri_t02[:,1], kcorri_t02[:,2]
z_3, kcorri_3 = kcorri_v01[:,1], kcorri_v01[:,2]

# print(z_1)
# print(kcorri_1)

# plot Fig. 2.
plt.plot(z_1,kcorri_1,label='l15')
plt.plot(z_2,kcorri_2,label='t02')
plt.plot(z_3,kcorri_3,label='v01')
plt.legend()
plt.show()

naptgrg

# Define eq.(1) in two ways:
# 1. Given apparent (observed) mag, m, find absolute mag, M1450
# 2. Given absolute mag, M1450, find apparent mag, m.
# Cosmological parameters:
H0 = 70 # km s-1 Mpc-1
Omega_m, Omega_lambda = 0.3, 0.7

c = 2.998 * 1e10
def Mbh(i_mag, q_s):
    # compute Mbh at z=2, from Haiman et al. 2009
    bh_mass = 3 * 1e6 * 10 ** ((24-i_mag) / 2.5)
    # print(np.log10(bh_mass))
    logM1 = bh_mass / (1 + q_s)
    return np.log10(logM1) # in unit of Msun

def imag(logM1, q_s,z):
    f_edd = 0.01
    d_z = float(cosmo.luminosity_distance(z) / (u.parsec))
    d_2 = float(cosmo.luminosity_distance(2) / (u.parsec))
    return 26 + 2.5 * np.log10((f_edd / 0.01) ** (-1) * (10 ** logM1 / (3 * 1e7)) ** (-1) * (d_z / d_2) ** 2)

ratio = 0.1

# mag = np.linspace(14,29,200)
# m1 = Mbh(mag,ratio)
# # print(m1)
#
# log_m1 = np.linspace(4.5,10,200)
# i_mag_fine = imag(log_m1,ratio)
# print(i_mag_fine)

# i_mag2 = np.arange(16,28,1)
# print(Mbh(i_mag2,q))

def M1450(m,z):
    # print(z)
    # m is apparent i-band magnitude, Kcorr is i-band K-correction to m at 1450 A.
    # From z, find luminosity distance in Mpc.
    D = cosmo.luminosity_distance(z)
    D_L = float(D / (1e6*u.parsec))  # luminosity distance in Mpc.
    # print(D_L)
    if z <= 4.7:
        ind = np.where(z_1 == z)[0]
        Kcorr = kcorri_1[ind][0][-1]

    if z > 4.7:
        Kcorr = np.interp(z,z_1,kcorri_1)
        # print(Kcorr)

    M1450 = m - 5 * np.log10(D_L) - 25 - Kcorr
    # print(M1450)
    return M1450

# print(M1450(26.5,2))

def m(M1450,z,linear):
    # M1450 is absolute magnitude, Kcorr is i-band K-correction to m at 1450 A.
    # From z, find luminosity distance in Mpc.
    D = cosmo.luminosity_distance(z)
    D_L = float(D / (u.parsec))  # luminosity distance in Mpc.
    # print(D_L)
    if z <= 4.700:
        ind = np.where(z_1 == z)[0]
        # print(ind)
        Kcorr = kcorri_1[ind]
        # print(Kcorr)
    if z > 4.700:
        if linear == True:
            # Kcorr = -2.2
            f = np.poly1d(np.polyfit(z_1,kcorri_1, 1))
            Kcorr = f(z)

            # f = intp.interp1d(z_1,kcorri_1,fill_value='extrapolate')
            # Kcorr = f(z)
            # Kcorr = np.interp(z, z_1, kcorri_1)
            # print(Kcorr)
        else:
            Kcorr = -2.2
    # return K/corr
    m = M1450 + 5 * np.log10(D_L / 1e6) + 25 + Kcorr
    return m

# redshift = np.linspace(0.1,7.5,60)
# M_1450 = np.linspace(-32,-8,200)

# K = []
# K_flat = []
# for z in redshift:
#     K.append(float(m(1,round(z,2),linear=True)))
#     K_flat.append(float(m(1,round(z,2),linear=False)))
#
#
# # print(redshift)
# # print(K)
# plt.plot(redshift,K,zorder=1,linestyle='solid',label='linear extrapolation')
# plt.plot(redshift,K_flat,zorder=0,linestyle='--',label='constant extrapolation')
# plt.xlabel('z')
# plt.ylabel('K-correction')
# plt.legend()
# # plt.savefig('../Figures/K-corr.pdf')
# plt.show()


# print(m(-32,2))
# i = m(-31.040,5.5)[0]
# print(i)

#
# lf = np.loadtxt('../Data/bestfit_lf.txt')
# z_eff = 0.5
# z_0 = [0.31,0.50,0.72,0.91,1.10,1.30,1.50,1.71,1.98,2.30,2.45,2.55,2.65,
#             2.75,2.85,2.95,3.05,3.15,3.25,3.34,3.44,3.88,4.35,4.92,6.00]
# print(len(z_0))

# z_0 = [0.5,1.5,2.5,3.5,4.5,5.5,6.5]
# print(z_0)
# M_i = []
# for q_z in z_0:
#     print(q_z)
#     print(round(M1450(25,q_z),2))
# print(np.round(M_i,1))

redshift = np.linspace(0.1,6,200)
# redshift = np.linspace(0.9,1.1,10)
# redshift = np.linspace(0.1,7,60)
M_1450 = np.linspace(-32,-8,200)

# i_mag = []
# for q in M_1450:
#     mag = m(q,1,linear=False)
#     print(mag[0])
#     i_mag.append(mag[0])
# np.savetxt('../../../LSST-LISA/data/m-z1.txt',np.c_[M_1450,i_mag],fmt='%.4f %.4f')
#
# plt.plot(i_mag,M_1450)
# plt.title('z=0.3')
# plt.show()

for z in redshift:
    i_mag = []
    for q in M_1450:
        mag = m(q,round(z,2),linear=False)
        print(mag)
        i_mag.append(mag)
    # np.savetxt('../Data/K19-DoublePL/imag_q{0}_k19/m-z{1}.txt'.format('%.1f' % ratio,'%.2f' % z),np.c_[M_1450,i_mag],fmt='%.4f %.4f')

    plt.plot(i_mag,M_1450)
    plt.title('z={}'.format(z))
    plt.show()

# lf_mid = (lf[1:] + lf[:-1]) / 2
# delta = -np.diff(M_1450)
#
# Nqso = np.sum(lf_mid * delta)
# print(Nqso)

# shen = np.loadtxt('../Data/LF_M1450_shen2020/lf_m1450_z6.5.txt')
# m1450 = shen[:,0]
# lf = shen[:,1]
# z_eff = 6.5
# print(len(m1450))

# i_mag = []
# for q in range(200):
#     i = m(M_1450[q],z_eff)[0]
#     i_mag.append(i)
#
# print(i_mag)
# np.savetxt('../Data/LF_imag_Shen2020/lf_imag_z6.5.txt',np.c_[i_mag,lf],fmt='%.6f %.6f')
# frac = i_mag / M_1450

# lf_i = lf / frac
# np.savetxt('../Data/lf_i_z1-1.2.txt',np.c_[i_mag,lf_i],fmt='% .4f % .4f')

# # plt.plot(i_mag,M_1450)
# plt.plot(i_mag,lf)
# # plt.plot(M_1450,lf)
# plt.gca().invert_xaxis()
# plt.show()

# def dV_dzdOmega(z):
#     pre_fac = (c / H0) * 1e-5
#     D = cosmo.luminosity_distance(z)
#     d_L = float(D / (u.parsec)) # in Mpc
#     denom = ((1+z) ** 2) * ((Omega_m * ((1+z) ** 3) + Omega_lambda) ** 0.5)
#     result = pre_fac * (d_L ** 2) / denom
#     return result # in Mpc^3
#
# def dV_dz(z,A):
#     # A is survey area in deg^2.
#     result = dV_dzdOmega(z) * A * 4 * np.pi / 41253
#     return result

def volume(z, area):
    omega = (area/41253.0)*4.0*np.pi # str
    vol_perstr = cosmo.differential_comoving_volume(z) # cMpc^3 str^-1 dz^-1
    volperstr = float(vol_perstr * u.sr * (1e6*u.parsec)**(-3))
    return omega*volperstr # cMpc^3 dz^-1







