import numpy as np
import math
import sys
import glob
import pandas as pd
import random
from scipy.interpolate import griddata
from scipy import integrate

from scipy import interpolate




#SWITCHES 
sys.setrecursionlimit(9000)
ACH_TEMP = 4e3 #1e4 is standard
ACH_MAX_TEMP = 1e4
COOLING_FACTOR = 20 #1, >1 increase the number of dcbh candidates



OMEGA_B = 0.0486
OMEGA_M = 0.315
OMEGA_L = 1-OMEGA_M
SIGMA_8 = 0.811
delta_c = 1.686
h=0.67
rho_crit = 1.5*1e11 #Msol/mpc^3
T_HUBBLE = 4.55e17 #s

K = 1.380649e-16 #ergs/K
SEC_PER_MYR = 3.15576e+13
m_p = 1.6e-27
msol = 2e30
mpc = 3.086e+24
n = rho_crit*msol/mpc**3*OMEGA_B/m_p #cm^-3
N_GAS = 0.2
u=1.22 #for neutral primordial gas
MPC_TO_CM = 3.086e+24




def get_redshift(file):
    f = open(file, "r")
    z=[]
    for line in f:
        a=float(line.strip("\n"))
        z.append(1/a - 1)
    return z

def f_shield(n_h2,n,T):
    x = n_h2/5e14
    b = 3 #km/s
    b5 = b#/1e5 #cm/s
    alph = alpha(n,T)
    a = 0.965/np.power(1+x/b5, alph)
    b = 0.035/np.power(1+x, 0.5)
    c = np.exp(-8.5 *1e-4*np.power(1+x, 0.5))
    f = np.array(a + b*c)
    f[np.where(f>= 1)] = 1
    return f

def alpha(n,T):
    c1,c2,c3,c4,c5 = 0.2856,0.8711, 1.928, 0.9639, 3.892
    A1 = c2*np.log10(T) - c3
    A2 = -c4*np.log10(T) + c5
    return A1 * np.exp(-c1 *np.log10(n))  + A2



def r_vir(m_h,z):
    #returns mpc
    #h=0.7

    OMEGA_K = 1-OMEGA_M-OMEGA_L
    OMEGA_M_Z=OMEGA_M*np.float_power(1+z,3)/ (OMEGA_M*np.float_power(1+z,3) + OMEGA_L + OMEGA_K*np.float_power(1+z,2))
    d = OMEGA_M_Z-1
    delta_c = 18*np.pi*np.pi + 82*d -39*d*d

    return 1e-3*0.784/h* np.float_power(m_h*h/1e8, 1.0/3)*np.float_power(OMEGA_M/OMEGA_M_Z * (delta_c/(18*np.pi*np.pi)), -1.0/3) * (10/(1+z))


def calc_N_h2_col(m,z,J,T):
    dist = r_vir(m,z)*MPC_TO_CM
    n_h2 = calc_n_h2(J,T)
    return dist*n_h2


def calc_J_shield(m,z,J):
    T = T_vir(m,z)
    n_h2_col = calc_N_h2_col(m,z,J,T)
    n = calc_n_H(T)/(0.76*u)
    f = f_shield(n_h2_col,n,T)
    print(f)
    return J*f

def T_vir(m,z):
    #h=0.7
    #eq 26 from https://arxiv.org/pdf/astro-ph/0010468.pdf
    #returns virial temp in units of Kelvin
    OMEGA_K = 1-OMEGA_M-OMEGA_L
    OMEGA_M_Z=OMEGA_M*np.float_power(1+z,3)/ (OMEGA_M*np.float_power(1+z,3) + OMEGA_L + OMEGA_K*np.float_power(1+z,2))
    d = OMEGA_M_Z-1

    delta_c = 18*np.pi*np.pi + 82*d -39*d*d
    return 1.98e4*(u/0.6)*np.float_power(m*h/1e8, 2.0/3)*np.float_power(OMEGA_M/OMEGA_M_Z * delta_c/(18*np.pi*np.pi),1.0/3) *(1+z)/10.0

def calc_k_LW(J_lw):
    return 1.39e-12*J_lw

def z_to_time(z):
        a=2./(3.*np.sqrt(OMEGA_L))
        b = np.sqrt(OMEGA_L/OMEGA_M)*np.float_power(1+z, -1.5)
        c = 71*1e5*3.154e+7/3.086e18
        return a/c*np.log(b+np.sqrt(1+b*b))

def calc_n_H(T):
    u = 1.22
    return N_GAS*0.76*u*6*np.power(T/1000, 1.5) #cm^-3

def calc_n_e(T):
    return 1.2e-5*np.sqrt(OMEGA_M)/(h*OMEGA_B)*calc_n_H(T)

def calc_k9(T):
    a = np.power(T, 0.928)
    b = np.exp(-T/16200)
    k9 = 1.4e-18*a*b
    return k9

def calc_n_h2(J_lw,T):
    k9 = calc_k9(T)
    n_e = calc_n_e(T)
    n_H = calc_n_H(T)
    k_LW = calc_k_LW(J_lw)
    return k9*n_H*n_e/k_LW

def U(n_h, temps):  #ergs/cm^3
    #V = 4/3*np.pi*np.power(r,3)
    u=1.22
    n = n_h/(0.76*u)
    return 1.5*n*K*temps#*V

#
# def calc_t_cool(z,ms, j_lw):
#     temps = T_vir(ms,z)
#     n_h = calc_n_H(temps)
#     n_h2 = calc_n_h2(j_lw, temps)
#     u = U(n_h, temps)
#     c = total_cooling_rate(temps,n_h)*n_h*n_h2
#     t_cool = 1/(5/3 - 1)*u/(c)
#     return t_cool


def cooling_LTE(temps, n_h):
    T = temps/1.0e3
    a = 9.5e-22*np.power(T, 3.76)/(1+0.12*np.power(T, 2.1))*np.exp(-2.197e-3/np.power(T,3))
    b = 3e-24*np.exp(-0.51/T) + 6.7e-19*np.exp(-5.86/T)
    c = 1.6e-18*np.exp(-11.7/T)
    return (a+b+c)/n_h

def cooling_GP(T):
    T=np.log10(T)
    return np.power(10,-103.0 + 97.59*T -48.05*T**2 + 10.80*T**3-0.9032*T**4)

def total_cooling_rate(temps, n_h):
    a = cooling_LTE(temps, n_h)
    b = cooling_GP(temps)
    return a/(1+a/b)

def dynamical_heating(mass_list,dm_list, dt_list,z_list, temp_list):

    dmdt = dm_list/dt_list
    dedt =  temp_list/mass_list * K/(5.0/3.0-1.0) * dmdt/SEC_PER_MYR
    return dedt

def atomic_cooling_halos(temps):
    no_index = np.where(temps < ACH_TEMP)
    arr = np.full(len(temps), True)
    arr[no_index] = False
    return arr

def calc_t_cool(mass, dm,dt,zs, temps, j_lw):
    n_h = calc_n_H(temps)
    n = n_h*13.0/12.0
    n_h2 = calc_n_h2(j_lw, temps)

    u = U(n, temps)
    c = total_cooling_rate(temps,n_h)*n_h2*n_h
    h = dynamical_heating(mass, dm, dt,zs, temps)*n_h
    ch = c-h

    ch[np.where(ch==0)] = 1e-50
    t_cool = u/(ch)
    t_cool = np.array(t_cool)
    t_cool[np.where(np.isnan(t_cool) == True)] = 1e20
    t_cool[np.where(np.isinf(t_cool) == True)] = 1e20
    injection = np.where(t_cool <= 0)
    t_cool[injection] = 1e20
    return t_cool

def cooling_condition(t_cool,prog_counts, times):
     no_index = np.where(t_cool < times*SEC_PER_MYR/COOLING_FACTOR)
     arr = np.full(len(t_cool), True)
     arr[no_index] = False
     yes_index = np.where(prog_counts == 0)
     arr[yes_index] = True
     return arr









#
#
# def mcrit_mihir(j_21, v_bc,z):
#     #eq 2 onward from https://arxiv.org/pdf/2010.04169.pdf
#
#     m0 = 1.96e5
#     b1,b2,b3 = 0.8,1.83, -0.06
#     a0,y1,y2,y3 = 1.64, 0.36, -0.62,.13
#     j0,v0,jv0 = 1,30,3
#     a = a0*np.power(1+j_21/j0,y1)*np.power(1+v_bc/v0, y2)*np.power(1+j_21*v_bc/jv0, y3)
#     m_z20 = m0*np.power(1+j_21/j0,b1)*np.power(1+v_bc/v0, b2)*np.power(1+j_21*v_bc/jv0, b3)
#     m_crit = m_z20*np.power((1+z)/20, -a)
#     return m_crit

'''
def cooling_condition(mass, z, j_lw, v_bc, prog_counts,dm,dt):
    dcbh_bools = np.full(len(mass), False)


    mcrit = mcrit_mihir(j_lw, v_bc, z)
    yes_index = np.where((mass < mcrit))
    dcbh_bools[yes_index] = True

    thub = z_to_time(z)*SEC_PER_MYR/np.sqrt(200)
    no_index = np.where(dcbh_bools == False)[0]
    for nogo in no_index:
        mc,mh,t,z1,= mcrit[nogo], mass[nogo], thub[nogo], z[nogo]
        T_crit,T_h = T_vir(mc,z1), T_vir(mh,z1)
        T_scaling = T_h/T_crit
        U_crit = 1.5*K*T_crit  #omitting n_h, cancels out
        Lam_N_h2 = U_crit/t
        U_h    = 1.5*K*T_h
        heating = dynamical_heating(mh, dm[nogo], dt[nogo],z1, T_h)
        n_halo, n_crit = calc_n_H(T_h), calc_n_H(T_crit)
        Lam_halo, Lam_crit = total_cooling_rate(T_h,n_halo),total_cooling_rate(T_crit,n_crit)
        Lam_N_h2_halo = Lam_N_h2 * (n_halo*Lam_halo)/ (n_crit*Lam_crit)  #Assuming lam cooling is more complicated than our simple model, scaling mihir's value using our simple prescription
        t_h = U_h/(Lam_N_h2_halo - heating)
        if(t_h > thub[nogo]):
            dcbh_bools[nogo] = True
    return dcbh_bools
'''


J_MEAN_BASE = 100
def draw_jlw_base():
    upperbound =7
    final = integrate.quad(calc_N, 0, upperbound, args=(np.log10(J_MEAN_BASE)))[0]
    x_vals, cdf_vals = [],[]
    for x in np.arange(0,upperbound,0.01):
        cdf_vals.append(integrate.quad(calc_N, 0, x, args=(np.log10(J_MEAN_BASE))) [0]/final)
        x_vals.append(x)
    cdf_vals.append(1),x_vals.append(upperbound)
    return x_vals, cdf_vals
#
# def draw_jlw(j_mean, x_vals, cdf_vals, N_draws = 1):
#
#     draw =np.random.uniform(0,1.0,N_draws)
#     new_xvals = x_vals+np.log10(j_mean/J_MEAN_BASE)
#     f1 = interpolate.interp1d(cdf_vals, new_xvals,kind = 'linear')(draw)
#     return np.power(10.0,f1)

def draw_jlw_ratio(x_vals, cdf_vals, N_draws = 1):

    draw =np.random.uniform(0,1.0,N_draws)
    f1 = interpolate.interp1d(cdf_vals, x_vals,kind = 'linear')(draw)
    return np.power(10.0,f1)/J_MEAN_BASE

def calc_N(x, j_mean_log):
    if(x<j_mean_log): N = np.power(10, np.log10(500) - (j_mean_log-x)*2)
    else: N= np.power(10,np.log10(500)+ (j_mean_log-x)*2)
    return N



def main():
    print(sys.getrecursionlimit())
    block = int(sys.argv[1])

    for tree_num in range(block*100, block*100+100):

        print("reading tree chunk ", str(tree_num))

        #ds =  pd.read_csv(f,delim_whitespace=True,skiprows = start_index,nrows=n_rows, names=["ind","snap_num","desc","prog_count","first_prog","next_prog","mass/1e10"])


        ds =  pd.read_csv("/scratch/08288/tg875874/trees/tree_"+str(tree_num),delim_whitespace=True,skiprows = 1,names=["ind","snap_num","desc","prog_count","first_prog","next_prog","mass/1e10"])


        inds = np.array(ds["ind"])
        z = get_redshift("Snapshot_alist_mod")

        snaps = np.array(ds["snap_num"])
        zs = np.array(z)[snaps]
        mass = np.array(ds["mass/1e10"]*1e10)/h
        temps = T_vir(mass,zs)
        inds = np.array(ds["ind"])


        times=z_to_time(zs)
        z_eval = zs#[0:1000000]
        m_eval = mass#[0:1000000]
        df_j = pd.read_csv("J_lw.txt", delim_whitespace=True)
        z_jlw,m_jlw,jlw = df_j["z"], df_j["mass"], df_j["j"]
        grid_jlw = np.array(griddata(np.array([z_jlw,m_jlw]).T, jlw, (z_eval, m_eval), method='nearest'))




        achs = atomic_cooling_halos(temps)

        desc = np.array(ds["desc"])
        prog = np.array(ds["first_prog"])
        prog2 = np.array(ds["next_prog"])
        no_prog = np.where(prog == -1)[0]
        prog[no_prog] = no_prog
        dm = np.array(mass)-np.array(mass)[prog]
        dt = np.array(times)-np.array(times)[prog]
        dm[no_prog] = 0
        dt[no_prog] = -1
        prog = np.array(ds["first_prog"])

        x_base, cdf_base = draw_jlw_base()
        jlw_ratios = draw_jlw_ratio(x_base, cdf_base, N_draws= np.size(grid_jlw))



        ####PROPOGATING JRATIO TO PROGS
        end_ids = np.where((np.array(ds["prog_count"]) == 0)&(achs == False))[0]
        all_indices_for_jratio = []
        all_jratios = []
        for end_id in end_ids:

            indice_subset = []
            halo_id = end_id
            desc_id = desc[halo_id]
            indice_subset.append(halo_id)

            while(True):

                halo_id = desc_id
                desc_id = desc[halo_id]
                indice_subset.append(halo_id)

                #if(achs[halo_id]==True):
                if(temps[halo_id]>=ACH_MAX_TEMP): #Changing this to propogate the jratio all the way to the max temp rather than the lowest
                    ratio = jlw_ratios[halo_id]
                    all_jratios.extend([ratio for count in range(len(indice_subset))])
                    all_indices_for_jratio.extend(indice_subset)
                    break


        for index in all_indices_for_jratio:
            jlw_ratios[index] = all_jratios[index]

        drawn_jlw = grid_jlw*jlw_ratios

        j_lw = calc_J_shield(mass,zs,drawn_jlw)
        ds["j_lw"] = drawn_jlw




        t_cools = calc_t_cool(mass, dm,dt,zs, temps, j_lw)
        ds["t_cools"] = t_cools
        ds["tcool/thub"] = t_cools/(times*SEC_PER_MYR)
        #ds["t_hubs"] = times
        #dcbh_sites = cooling_condition(t_cools,np.array(ds["prog_count"]), times)

        dcbh_sites = cooling_condition(t_cools,ds["prog_count"], times)
        dcbh_pure_branch = np.full(len(dcbh_sites), 1)

        end_ids = np.where((np.array(ds["prog_count"]) == 0)&(achs == False))[0]
        total = np.size(end_ids)

        count = 0
        print("TOTAL NUMBER OF ENDS", total)

        for end_id in end_ids:

		    #3 cases:
		    #1 the halo doesn't have a non-ach prog
		    #2 the halo has non-ach progs and they all meet the conditions
		    #3 the halo has been polluted
			#starting from a dcbh-cand, propogate the failure up the tree

            desc_id = desc[end_id]
            halo_id = end_id


            if(count%1000 ==0):print(desc_id, halo_id)
            current_state = 2
            while(True):

                dcbh_site = dcbh_sites[halo_id]
                if(dcbh_site == False): current_state = 3
                old_state = dcbh_pure_branch[halo_id]
                if(old_state==3): break

                if(old_state < current_state): dcbh_pure_branch[halo_id] = current_state

                if(desc_id == -1): break

                halo_id = desc_id
                desc_id = desc[desc_id]
                if(halo_id == 0 and desc_id == 0):
                    print("ERROR")
                    return

            count+=1

        prog_achs = achs[prog]
        prog_purity = dcbh_pure_branch[prog]
        #b = np.where((achs == True) & (prog_achs == False) & (prog_purity==2) & (snaps>2))[0]
        #MODIFYING THIS FOR LEO PROJECT
        b = np.where((achs == True) & (prog_achs == False) & (prog_purity==2) & (snaps>2))[0]
        locations = []
        total_dcbh_sites = np.size(b)
        for x in range(np.size(b)):
            locations.append(b[x])
            desc_id = desc[b[x]]
            while(desc_id != -1):
                locations.append(desc_id)
                desc_id = desc[desc_id]

            #getting all progenitors of these dcbh candidates
            #prog_indices, secondary = get_all_progs(prog, prog2, b[x])
            #locations.extend(prog_indices)
            locations.extend(get_all_progs(prog, prog2, b[x]))

            #locations.extend(get_all_progs(prog_ids1))
            #progenitor_indices = get_all_progenitor_indices(desc,inds, b[x])
            #print("prog_indices: ", progenitor_indices)
            #print("and their index ", inds[np.array(progenitor_indices)])
            #print("\n\n")
            #add these to location#
            #locations.extend(progenitor_indices)




        locations = np.unique(np.array(locations))
        ds_write = ds.iloc[locations]

        ds_write.to_csv('dcbh_files/tree'+str(tree_num)+'_sitenum_'+str(total_dcbh_sites)+'.csv',index=True)


def get_all_progs(prog, prog2, i):
    '''
    prog_id = prog[i]
    prog_id2 = prog2[i]
    locations = []
    secondary = []
    while(prog_id != -1):
        print(prog_id)
        locations.append(prog_id)
        if(prog_id2 != -1): locations.append(prog_id2)
        prog_id = prog[prog_id]
    return locations, secondary   

    '''
    prog_id1 = prog[i]
    prog_id2 = prog2[i]
    list_ = [i]
    if(prog_id1 == -1): return list_
    if(prog_id2 != -1): list_.extend(get_all_progs(prog, prog2, prog_id2))
    list_.extend(get_all_progs(prog, prog2, prog_id1))
    return list_
    



    
def get_all_progenitor_indices(desc_list, inds_list, target_indice):
    target_ind = inds_list[target_indice]
    progenitor_indices = np.where(desc_list == target_ind)[0]
    #print("pi ", progenitor_indices)
    #print(target_ind)
    #print(target_ind in desc_list)
    #print(np.size(progenitor_indices))
    


    final_list = [target_indice]
    if (np.size(progenitor_indices) == 0): return final_list
    for pi in progenitor_indices:
        final_list.extend(get_all_progenitor_indices(desc_list, inds_list, pi))
    return final_list



if __name__ == "__main__":
    main()
