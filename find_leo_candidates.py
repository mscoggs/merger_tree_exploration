import numpy as np
import math
import sys
from astropy import units as u
from astropy import constants as const
import pandas as pd
import astropy
from astropy import units as u
from astropy.constants import G
import glob

rho_crit = 1.5*1e11 #Msol/mpc^3
T_HUBBLE = 4.55e17 #s
OMEGA_M = 0.315
OMEGA_B = 0.0486
OMEGA_M = 0.307
OMEGA_L = 1-OMEGA_M
SIGMA_8 = 0.811
delta_c = 1.686
h=0.677
rho_crit = 1.5*1e11 #Msol/mpc^3
T_HUBBLE = 4.55e17 #s
MPC_TO_CM = 3.086e+24
K = 1.380649e-16 #ergs/K
SEC_PER_MYR = 3.15576e+13
m_p = 1.6e-27
msol = 2e30
mpc = 3.086e+24
n = rho_crit*msol/mpc**3*OMEGA_B/m_p #cm^-3
N_GAS=0.2


r_leo = 250*u.kpc #kpc, not radius of leo but distance between leo and mw
mass_MW = 0.9e12*u.Msun
mass_Leo = 7e8*u.Msun
mass_tolerance = 1e8*u.Msun
mt = mass_tolerance*3.5

tff = np.pi/2.0 * (np.power(r_leo,1.5)/np.sqrt(2*G*(mass_MW +mass_Leo)))
tff = tff.to("Myr")




def get_redshift(file):
    f = open(file, "r")
    z=[]
    for line in f:
        a=float(line.strip("\n"))
        z.append(1/a - 1)
    return z

def z_to_time(z):
        a=2./(3.*np.sqrt(OMEGA_L))
        b = np.sqrt(OMEGA_L/OMEGA_M)*np.float_power(1+z, -1.5)
        c = 71*1e5*3.154e+7/3.086e18
        return a/c*np.log(b+np.sqrt(1+b*b))

def T_vir(m,z):
    #h=0.7

    u = 1.22
    #eq 26 from https://arxiv.org/pdf/astro-ph/0010468.pdf
    #returns virial temp in units of Kelvin
    OMEGA_K = 1-OMEGA_M-OMEGA_L
    OMEGA_M_Z=OMEGA_M*np.float_power(1+z,3)/ (OMEGA_M*np.float_power(1+z,3) + OMEGA_L + OMEGA_K*np.float_power(1+z,2))
    d = OMEGA_M_Z-1

    delta_c = 18*np.pi*np.pi + 82*d -39*d*d
    return 1.98e4*(u/0.6)*np.float_power(m*h/1e8, 2.0/3)*np.float_power(OMEGA_M/OMEGA_M_Z * delta_c/(18*np.pi*np.pi),1.0/3) *(1+z)/10.0


def T_vir_inverted(T,z):
    #h=0.7
    u = 1.22
    #eq 26 from https://arxiv.org/pdf/astro-ph/0010468.pdf
    #returns virial temp in units of Kelvin
    OMEGA_K = 1-OMEGA_M-OMEGA_L
    OMEGA_M_Z=OMEGA_M*np.float_power(1+z,3)/ (OMEGA_M*np.float_power(1+z,3) + OMEGA_L + OMEGA_K*np.float_power(1+z,2))
    d = OMEGA_M_Z-1

    delta_c = 18*np.pi*np.pi + 82*d -39*d*d
    m = np.float_power(T/(1.98e4*(u/0.6)) / (np.float_power(OMEGA_M/OMEGA_M_Z * delta_c/(18*np.pi*np.pi),1.0/3) *(1+z)/10.0),3.0/2)*1e8/h
    return m

def get_all_progs(inds, prog, prog2, i):
    prog_id1 = prog[i]
    prog_id2 = prog2[i]
    list_ = [i]
#     print(prog_id2 in inds)
#     print(prog_id1, prog_id2)
    if((prog_id1 == -1) or (prog_id1 not in inds)): return list_
    if(prog_id2 != -1):
        if(prog_id2 in inds):
            prog_id2 =np.where(inds ==prog_id2)[0][0]
            list_.extend(get_all_progs(inds, prog, prog2, prog_id2))
    prog_id1 =np.where(inds ==prog_id1)[0][0]
    list_.extend(get_all_progs(inds,prog, prog2, prog_id1))
    return list_


def save_data(data, filenum):
    columns = ['merger_ratio', 'target_temp', 'min_t_ratio', 'count']
    df = pd.DataFrame(data, columns=columns)
    df.to_csv("leo_candidacy_files/leo_candidacy_"+str(filenum),index=False)

def get_grid_ind(x,y,z):
    return TEMP_N*TCOOL_N* x + TCOOL_N*y + z



def make_grid():
    m_list, temp_list, tratio_list = [],[],[]
    count_list = []
    for x in range(MERGER_N):
        for y in range(TEMP_N):
            for z in range(TCOOL_N):
                m_list.append(MERGER_RANGE[x])
                temp_list.append(TEMP_RANGE[y])
                tratio_list.append(TCOOL_RANGE[z])
                count_list.append(0)
    m_list = np.array(m_list)
    temp_list = np.array(temp_list)
    tratio_list = np.array(tratio_list)
    count_list = np.array(count_list)
    grid = np.array([m_list, temp_list, tratio_list, count_list]).T
    return grid
    
def check_all_progs(inds, prog, prog2, tvirs):
    bool_list = np.full(np.size(inds), True)
    tvir_filter1 = np.where((tvirs <= np.max(TEMP_RANGE)))[0]
    tvir_filter2 = np.where((tvirs <= np.max(TEMP_RANGE)) & (tvirs >= np.min(TEMP_RANGE)))[0]
    sub_inds = inds[tvir_filter1]
    sub_prog = prog[tvir_filter2]
    sub_prog2 = prog2[tvir_filter2]
    true = []

    size = np.size(sub_prog)
    for p1,p2,num in zip(sub_prog, sub_prog2, range(size)):
        #print("working on number", num, "out of", size)

        
        if(p1 == -1): true.append(True)
        elif((p1 in sub_inds) and ((p2 == -1) or (p2 in sub_inds))): true.append(True)
        else: true.append(False)

    for i in range(len(tvir_filter2)):
        bool_list[tvir_filter2[i]] = true[i]
    return bool_list


def find_leo_dcbhs(start=0, stop=1000):
    dcbh_files = np.sort(list(glob.iglob("dcbh_files/*.csv")))
    dcbh_files_number = np.array([int(b.split("/tree")[-1].split("_")[0]) for b in dcbh_files])
    dcbh_files_order = np.argsort(dcbh_files_number)
    
    dcbh_files = dcbh_files[dcbh_files_order]
    dcbh_files_number = dcbh_files_number[dcbh_files_order]

    dcbh_files = dcbh_files[start:stop]
    dcbh_files_number = dcbh_files_number[start:stop]
    inds_keep = []
    print(dcbh_files_number)


    for filenum, f in zip(dcbh_files_number, dcbh_files):#range(Ntrees):

        grid = make_grid()

        print("\n\n\nworking on file", f, str(filenum)+"/"+str(len(dcbh_files)))
        print("\n\n")

        ds =  pd.read_csv(f,delimiter=",",skiprows = 1, names=["filler","ind","snap_num","desc","prog_count","first_prog","next_prog","mass/1e10", "jlw", "t_cools", "tcool/thub"])

        inds = np.array(ds["ind"])
        z = get_redshift("Snapshot_alist_mod")
        snaps = np.array(ds["snap_num"])
        zs = np.array(z)[snaps]
        mass = np.array(ds["mass/1e10"]*1e10)/h
        t_ratio = np.array(ds["tcool/thub"])

        inds = np.array(ds["ind"])
        desc = np.array(ds["desc"])
        prog = np.array(ds["first_prog"])
        prog2 = np.array(ds["next_prog"])
        times=z_to_time(zs)
        tvirs = T_vir(mass,zs)
        print("file info read, checking progs")



        has_all_progs = check_all_progs(inds, prog, prog2, tvirs)

        print("filtering for leo candidates")


        selection_filter = np.where((times >= (np.max(times)-tff.value)) & (mass <= (mass_Leo+mt).value) & (mass >= (mass_Leo-mt).value))
        filtered_mass = mass[selection_filter]
        filtered_progs = prog[selection_filter]
        filtered_inds = inds[selection_filter]

        print("removing duplicate branches")
        print(np.size(selection_filter))

        if(np.size(selection_filter) == 0):
            blank_grid = make_grid()
            save_data(blank_grid, filenum)
            print("no leo candidates, skipping fnum", filenum)
            continue
        #remove duplicate candidates (different snaps of the same halo)
        branch_ids = []
        print(len(filtered_progs))
        for x,current_prog in zip(range(len(filtered_progs)), filtered_progs):
            print(current_prog)

            if any(current_prog in sub_branch for sub_branch in branch_ids): 
                print("HIT CONTINUE SECTION", current_prog)
                continue

            this_progs_ids = [filtered_inds[x]]
            prog_i = filtered_progs[x]

            while(True):
                index = np.where(inds == prog_i)[0]
                if(index.size == 0): break
                this_progs_ids.append(prog_i)

                prog_i = prog[index[0]]

            branch_ids.append(this_progs_ids)


        #find the number of leo candidates
        max_mergers = 3
        data = []
        print("iterating sub branches")
        print(branch_ids)
        #print(np.size(branch_ids))
        print("thats the size of the branch ids")

        leo_indices = []

        for sub_branch in branch_ids:
            indices = []
            for ind in sub_branch:
                a=np.where(inds ==ind)[0][0]
                indices.append(a)

            count_for_branch = 0
            for y in range(TEMP_N):
                branch_mass = mass[indices]
                branch_z = zs[indices]
                branch_tvirs= tvirs[indices]
                branch_has_all_progs = has_all_progs[indices]


                target_temp = TEMP_RANGE[y]
                tvir_filter = np.where(branch_tvirs >= target_temp)[0]
                branch_mass = branch_mass[tvir_filter]
                branch_z = branch_z[tvir_filter]
                branch_has_all_progs = branch_has_all_progs[np.where(branch_tvirs <= target_temp)[0]]
                if(np.any(branch_has_all_progs == False)):
                    print("failed progs for temp:", target_temp)
                    continue


                index_of_ach = (np.array(indices)[tvir_filter])[-1]
                
                print(index_of_ach)
                print(prog[index_of_ach])
                #print(np.where(inds ==prog[index_of_ach]))
                #print(prog[index_of_ach] in inds)
                #print(prog[index_of_ach] in indices)
                #print(tvirs[index_of_ach], zs[index_of_ach], mass[index_of_ach])
                index_of_ach = np.where(inds ==prog[index_of_ach])
                if(np.size(index_of_ach) == 0):
                    print("failed to find index of ach, skipping")
                    continue

                index_of_ach = index_of_ach[0][0]
                #print(tvirs[index_of_ach], zs[index_of_ach], mass[index_of_ach])

                for x in range(MERGER_N):
                    for z in range(TCOOL_N):
                        grid_ind = get_grid_ind(x,y,z)
                        merger_ratio, target_temp, min_t_ratio, count = grid[grid_ind]
                        print("working on mr, tt, mt:", merger_ratio, target_temp, min_t_ratio)

                        large_merger_count = 0
                        for xi in range(len(branch_mass) -1):
                            if(branch_mass[xi] >= merger_ratio* branch_mass[xi+1]):
                                large_merger_count += 1
                        lmc = large_merger_count
                        if((lmc <= max_mergers) and (lmc > 0)):
                            #check the cooling times
                            #get all the progs
                            locations = []
                            locations.extend(get_all_progs(inds,prog, prog2, index_of_ach))
                            ach_t_ratio = t_ratio[locations]
                            if(np.all(ach_t_ratio > min_t_ratio)):
                                count += 1

                        count_for_branch += count
                        grid[grid_ind][-1] = count

            if(count_for_branch > 0): leo_indices.extend(indices)


                
        ds_write = ds.iloc[leo_indices]
        ds_write.to_csv("leo_files/tree"+str(filenum)+".csv", index=True)
        save_data(grid, filenum)





MERGER_N = 20
TEMP_N = 20
TCOOL_N = 10
MERGER_RANGE = np.linspace(1.1, 2, MERGER_N)
TEMP_RANGE = np.linspace(4e3, 1e4, TEMP_N)
#tvir_filter = np.where(branch_tvirs >= target_temp)[0]
TCOOL_RANGE = np.linspace(0.05, 1.0, TCOOL_N)

args = sys.argv
print(args)
a, b = 0,1000
if(len(args) > 1):
    a = int(args[1])
    b = int(args[2])

print("this is start and stop index for tree files:", a, b)

find_leo_dcbhs(start = a, stop = b)
