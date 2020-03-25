# IMPORT
import re
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import os
import math
from math import factorial
from scipy.optimize import curve_fit
from scipy import stats
from collections import OrderedDict



def calculate_correction_barrier(reorganization, free_energy, couplings):
    if couplings > reorganization/2:
        return reorganization/4
    else:
        if free_energy == 0.0:
            return couplings - couplings**2/reorganization
        else:
            print "Free energy != O not implemented"
            raise SystemExit

def calculate_factor(coupling, reorganization, free_energy, temperature, frequency):
    exposent = (np.pi**1.5 * coupling**2) / ( 2*np.pi * np.sqrt(reorganization*temperature) * frequency)
    Plz = 1 - np.exp(-exposent)
    #print "Hab", coupling
    #print "exposent", exposent
    #print "plz", Plz
    #print 2*Plz / (1 + Plz)
    if free_energy >= - reorganization:
        return 2*Plz / (1 + Plz)

def calculate_rate(coupling, reorganization, free_energy, temperature, method, frequency = 1):
    #Convert everything in atomic units
    coupling = coupling  / 27211.399
    reorganization = reorganization  / 27211.399
    free_energy = free_energy  /     27211.399
    temperature = temperature / 315777.09
    frequency = frequency * 0.0000046 / (2*np.pi)

    barrier = np.square( ( reorganization + free_energy ) ) / ( 4 * reorganization)
    if method == 'NA':
        factor = (2 * np.pi) * np.square(coupling) / np.sqrt( 4 * np.pi * reorganization * temperature)
    elif method == 'AD':
        barrier += - calculate_correction_barrier(reorganization, free_energy, coupling)
        factor = frequency
    elif 'TOT' in method:
        barrier += - calculate_correction_barrier(reorganization, free_energy, coupling)
        factor = frequency * calculate_factor(coupling, reorganization, free_energy, temperature, frequency)
        #print frequency
        #print factor
    else:
        print "method should be NA or AD"
   
    if barrier < 0:
        barrier = 0.0
    #print barrier
    expo = np.exp( - barrier / temperature)
    #print "expo", expo
    rate = expo*factor / 0.024188  # 1au = 0.02418884 fs
    return rate # in fs-1

import math
import numpy as np


def calculate_rate_MLJ(coupling, lambS, lambI, free_energy, W0, temperature):
    """ This calculated the quantized Marcus-Levitch-Jortner rate according to Cupellini 2017 paper
        Inputs:
          - Coupling [meV]
          - lambS [meV] : Classical reorganization energy (called also solvent reorganization)
          - lambI [meV]: Quantum reorganization energy (that could be derived from the 4 point scheme)
          - free_energy [meV]
          - W0 [s-1]: which is the angular frequency of the quantized mode
          - temperature [K]

        Outputs:
          - rate in s-1

    """

    hbar  = 6.58211928E-16 # eV*s
    kboltz = 8.617333262145E-5  # eV

    #Convert everything in eV
    coupling = coupling*1e-3
    lambS = lambS*1e-3
    lambI = lambI*1e-3
    free_energy = free_energy*1e-3
    
    # Huang-Rhys factor for the quantum mode
    S = lambI/(hbar*W0)
    KT = temperature*kboltz

    FCS = 0.0
    control = 0.0
    j = 0
    while True:
        FCS += math.exp(-S)*S**j/(math.factorial(j))*math.exp(-(free_energy +lambS+j*hbar*W0)**2/(4*lambS*KT))
        control += math.exp(-S)*S**j/(math.factorial(j))
        #print j,control,FCS
        if 1.0 - control < 1e-14: break
        j += 1
    
    
    J = math.sqrt(1/(4*math.pi*lambS*KT))*FCS
    rate_s   = 2*math.pi*(J*coupling**2)/hbar # s^-1
    #print 'k(M.L.J.) =', rate_s, 's^-1'
    
    rate = rate_s*1e-15 # convert in fs-1
    return rate


def get_eigen(matrix):
    # we use np.linalg.eigh to deal with simmetric matrix
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    idx = eigenvalues.argsort()
    eigenvalues = np.array(eigenvalues[idx])
    eigenvectors = np.array(eigenvectors[:, idx])
    #print "all_eigen", eigenvectors
    return eigenvalues, eigenvectors # eigenvectors[0, :], eigenvectors[:, 0], eigenvectors[:, 1]

def time_evolution(n, Kinetic_matrix, tstep, nstep, starting_pos):
    #find eigenvalues and eigenvectors
    L, U = get_eigen(Kinetic_matrix)
    print "DIAGONALIZATION FINISHED"

    #Define initial condition
    P0 = np.zeros(n)
    P0[starting_pos] = 1.0
    
    #Define time step serie
    t = np.arange(nstep+1)*tstep

    Uinv = np.linalg.inv(U)

    # Find y(0)
    y0 = np.dot(Uinv,P0)
    # Solving P(t) = P(O)*U*exp(Lt)*U-1
    exps = np.exp(np.einsum('i,j->ij',L,t))
    Pt  = np.dot(U,np.einsum('ij,i->ij',exps,y0))
    return t, Pt

def read_coordinates(filename):
    com = []
    with open(filename, "r") as f:
        for line in f.readlines():
            com.append([float(i) for i in line.split()])
    return com
    
def read_connect(filename):
    f = open(filename, "r")
    file_ = f.readlines()
    natom = int(file_[0])
    connectivity_list = []
    for line in file_[1:]:
        list_line = [int(i) for i in line.split()]
        connectivity_list.append(list_line)
    return natom, connectivity_list

def isclose(a, b, rel_tol=1e-09, abs_tol=0.1):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def check_connect_from_distance(msd_length, com):
    connectivity = []
    #start looping over one mol
    count1 = 0
    for point1 in com:
        count1 += 1
        count2 = count1
        for point2 in com[count1:]:
            count2 += 1
            distance = abs(np.array(point2)-np.array(point1))
            distance =  np.linalg.norm(distance)
            # loop over list of interaction to do the check: if there are two equivalent interactions we need to do something else
            for iter_3, check_dist in enumerate(msd_length):
                #print "check_dist", check_dist
                diff = isclose(distance, check_dist)
                if diff:
                    connectivity.append([count1, count2, iter_3+1])   
    return connectivity


def build_kinetic_matrix(natom, connectivity_list, interaction_dict):
    #build kinetic matrix able to deal with both PBC and no PBC depending if their are present or not in connectivity.dat 
    k_matrix = np.zeros((natom,natom))
    
    #loop over the connectivity instead of a double loop over the kinetic matrix is faster
    for il in connectivity_list:
        #the interactions are 1-based but python counts from zero
        k = il[0]-1
        l = il[1]-1
        type_inter = str(il[2])
        mean = interaction_dict[type_inter]
        # deal with PBC in connectivity file
        if k_matrix[k,l] == 0.0:
            k_matrix[k,l] = mean
            k_matrix[l,k] = k_matrix[k,l]
        else:
            k_matrix[k,l] += mean
            k_matrix[l,k] = k_matrix[k,l]

    # construct the diagonal
    for diag in range(natom):
        # sum the column elements to keep detailed balance
        k_matrix[diag,diag] = -sum(k_matrix[diag,:])
    #print "\n"
    #print "Kmatrix \n"
    #print('\n'.join([''.join(['{:20}'.format(item) for item in row]) 
    #  for row in k_matrix]))
    return k_matrix


#def build_kinetic_matrix_NOPBC(natom, connectivity_list, interaction_dict):
#    #build kinetic matrix no PBC 
#    k_matrix = np.zeros((natom,natom))
#    
#    #loop over the connectivity instead of a double loop over the kinetic matrix is faster
#    for il in connectivity_list:
#        #the interactions are 1-based but python counts from zero
#        k = il[0]-1
#        l = il[1]-1
#        type_inter = str(il[2])
#        mean = interaction_dict[type_inter]
#    
#        k_matrix[k,l] = mean
#        k_matrix[l,k] = k_matrix[k,l]
#    
#    # construct the diagonal
#    for diag in range(natom):
#        # sum the column elements to keep detailed balance
#        k_matrix[diag,diag] = -sum(k_matrix[diag,:])
#    #print "Kmatrix", k_matrix
#    print('\n'.join([''.join(['{:4}'.format(item) for item in row]) 
#      for row in k_matrix]))
#    return k_matrix

#<L^2> = \sum_i Pi*(iL)^2 #correct expression
def _calc_3D_msd3(Pt, com_chain):
    populations = Pt.T
    com = np.array(com_chain)
    results = []
    vect0 = np.array([0.0, 0.0, 0.0])
    for (x,y) in zip(populations[0], com):
        vect0 +=  x *y

    for serie in populations:
        msd = np.zeros(9)
        for alpha in range(3):
            for beta in range(3):
                for (x, y) in zip(serie, com):
                    index = beta * 3 + alpha
                    msd[index] += x * (y[alpha] - vect0[alpha]) * (y[beta] - vect0[beta])
        results.append(list(msd))
    return results

#check if the matrix is symmetric and use this info for the diagonalization
def check_symmetric(a, rtol=1e-05, atol=1e-08):
    return np.allclose(a, a.T, rtol=rtol, atol=atol)





############################## MAIN ######################################################################

if __name__ == '__main__':

    """This code can solve the Master equations for a charge travelling in a given plane or chain depending 
      on the topology. The available rates are Non-adiabatic Marcus rate, Adiabatic Marcus rate (see Blumberger 2015) """



    general_path = "./"
    
    select_list_dir_2_plot = {
        
        'RUB_AOM':{    
            'npairs'  : 6, # provide number of pairs
             'Coupling' : [101.8, 10, 20, 30, 40, 50], #meV ---- provide list of couplings for all pairs
             'msd_length' : [8.75, ],  # distances Angstrom for each pair of molecules considered, not necessary is PBC is used 
                                       # this is only needed when one wants to reconstruct connectivity directly from pairs distance
             'Reorg energy' : [152.0, 152.0, 152.0, 152.0, 152.0, 152.0, ], #meV ----provide list for all pairs
             'free_energy_bias' : [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], #meV ----provide list for all pairs  
             'frequency' : [1601, 1601,1601,1601,1601,1601,],  # frequency cm-1 -----provide list for all pairs
             'Temperature' : 300, #Kelvin
            'connectivity' : general_path + "connectivity.dat",     
            'coordinates' : general_path + "coordinates.dat",     
             'connectivity_create' : general_path + "connectivity_created.dat" ,
            'start_position' : 0, # 0-based stating position
            'number of steps' : 1000, #tot number of steps, 
            'range linear fit' : slice(500, 1000), #range for the linear fit depending on the number of total steps
            'Time step' : 0.1 , #fs
            'Read_connectivity' : True, # if true you mast supply connectivity file in which pbc are or are not included
             'Method rate calc' : 'NA', #this can be Marcus = "NA", Jacob mod = "TOT", or MLJ_rate for quantized rate
            # if MLJ_rate is provided 
              "lambda_classical"  : [152.0, 152.0, 152.0, 152.0, 152.0, 152.0, ], #meV
              "lambda_quantum"    : [152.0, 152.0, 152.0, 152.0, 152.0, 152.0, ], # meV 
              'frequency_quantum' : [1601, 1601,1601,1601,1601,1601,],  # angualar freq in s-1
        },
        
    }
    
    # Multiple systems can be done one after the other by providing the appropriate parameters
    systems_list = [ 'RUB_AOM']
    
    #############################################   CODE ############################################


    params = {
       'axes.labelsize': 30,
       'font.size': 30,
       'legend.fontsize': 30,
       'xtick.labelsize': 30,
       'ytick.labelsize': 30,
       'text.usetex': False,
       'figure.figsize': [7,7],
        'axes.linewidth' : 3
       }
    rcParams.update(params)
    
    #hardcoded
    #free_energy = 0.0
    
    for system in systems_list:
        print "SYSTEM:", system
        
        n_pairs = select_list_dir_2_plot[system]['npairs']
        msd_length = select_list_dir_2_plot[system]['msd_length']
        hab = select_list_dir_2_plot[system]['Coupling']
        lambda_ = select_list_dir_2_plot[system]['Reorg energy']
        frequency = select_list_dir_2_plot[system]['frequency']
        temperature = select_list_dir_2_plot[system]['Temperature']
        nstep = select_list_dir_2_plot[system]['number of steps']
        range_ = select_list_dir_2_plot[system]['range linear fit']
        tstep = select_list_dir_2_plot[system]['Time step']
        pbc =  select_list_dir_2_plot[system]['Read_connectivity']
        created_conn =  select_list_dir_2_plot[system]['connectivity_create']
        method =  select_list_dir_2_plot[system]['Method rate calc']
        free_energy = select_list_dir_2_plot[system]['free_energy_bias'] 
    
        coord_file = select_list_dir_2_plot[system]['coordinates']
        starting_pos = select_list_dir_2_plot[system]['start_position']

        if (method=="MLJ_rate"):
            lambS = select_list_dir_2_plot[system]["lambda_classical"]
            lambI =  select_list_dir_2_plot[system]["lambda_quantum"]
            W0 =    select_list_dir_2_plot[system]['frequency_quantum']
        
        # check if you want to use PBC or not
        if pbc:
            print "CONNECTIVITY IS READ FROM EXTERNAL FILE"
            connectivity = select_list_dir_2_plot[system]['connectivity']
            nmol, connectivity_list = read_connect(connectivity)
            nmol = int(nmol)
            print "number of interactions :", len(connectivity_list)
        else:
            print "CONNECTIVITY IS RECONSTRUCTED FROM DISTANCES AND PBC NOT USED"
            coord_ = read_coordinates(coord_file)
            nmol = int(len(coord_))
            #chech for duplicate 
            msd_length = list(OrderedDict.fromkeys(msd_length))
            ##### here we check if there are duplicated interaction that we cannot handle using connect from distance for now!!!
            print "adjusted number of interactions :", msd_length
            print "NMOL:", nmol
            connectivity_list = check_connect_from_distance(msd_length, coord_)
            # write constructed connect to an external file
            with open(created_conn, 'w') as f:
                f.write("%s \n" %nmol)
                for item in connectivity_list:
                    f.write("%s  %s  %s\n" % (item[0],item[1],item[2]) )
        
        # write interaction list
        interaction_dict = {}
        for i in range(n_pairs):
            #v = hrr[i] # get coupling
            #below is the line actually implemented in the matlab code from where I took this
            #w =(v*v/6.5821192569654e-16)*((3.1415926/(e1[i]*0.026))**0.5)*np.exp(-e1[i]/(4*temperature*8.617333*1e-5))    
            #you can just give more line as input my rate matches the published
            if (method=="MLJ_rate"):
                # calculate the quantum MLJ rate
                rate = calculate_rate_MLJ(coupling=hab[i], lambS=lambS[i], lambI=lambI[i], 
                                                         free_energy=free_energy[i], W0=W0[i], temperature=temperature)
            else:
                # evaluate rate as specified in the input 
                rate = calculate_rate(coupling=hab[i], reorganization=lambda_[i], 
                                                              free_energy=free_energy[i], temperature=temperature,
                               method=method, frequency=frequency[i])###*1e15 #convertion to second-1
            print "rate %s: %s fs-1" %(str(i+1), rate)
            
            interaction_dict.update({str(i+1) : rate})
    
            
        # for now the number of interactions is just determined by the connectivity file
        k_matrix = build_kinetic_matrix(nmol, connectivity_list, interaction_dict)
        #print k_matrix
        simmetric = check_symmetric(k_matrix, rtol=1e-07, atol=1e-10)
        print "SIMMETRY", simmetric
    
        # calculate time evolution
        t, Pt = time_evolution(nmol, k_matrix, tstep, nstep, starting_pos)
        
        # plot population evolution for n states
        #for i in range(5):
        #    plt.plot(t, Pt[i])
        #plt.show()
        
        # read coordinates of the center of masses    
        com_chain = read_coordinates(coord_file)
        #print "COM", com_chain
        
        # calculate msd list (direction) of list(snaph): [[xx,xy,xz,yx,yy,yz, zx,zy,zz],[xx,xy,xz,yx,yy...]] 
        msd_list = _calc_3D_msd3(Pt, com_chain)
        # transpose to have time snapshots for each direction
        msd_snap = np.array(msd_list).T
        
        
        #fitting msd <L^2> xx
        fit = np.polyfit( t[range_], msd_snap[0][range_], 1 )
        #plot msd <L^2> xx
        fit_fn = np.poly1d(fit)
        plt.plot(t[range_], fit_fn(t[range_]), color = "k")
        plt.plot(t, msd_snap[0], label = "<L^2> dir. XX", linewidth=5, color="r")
        plt.plot(t[range_], fit_fn(t[range_]), color = "k")
        #mobility calculation
        mobility = fit[0] * 10**15 * 10**(-16) / (2 * 0.0000861728 * float(temperature) )
        print "Linear mob. from <L^2> in xx = \sum_i Pi*(iL)^2 is: %.5f cm2.(V.s)-1" % mobility
        
        #fitting msd <L^2> xy
        fit = np.polyfit( t[range_], msd_snap[1][range_], 1 )
        #plot msd <L^2> xxy
        fit_fn = np.poly1d(fit)
        plt.plot(t[range_], fit_fn(t[range_]), color = "k")
        plt.plot(t, msd_snap[1], label = "<L^2> dir. XY", linewidth=5, color="y")
        plt.plot(t[range_], fit_fn(t[range_]), color = "k")
        #mobility calculation
        mobility = fit[0] * 10**15 * 10**(-16) / (2 * 0.0000861728 * float(temperature) )
        print "Linear mob. from <L^2> in xy = \sum_i Pi*(iL)^2 is: %.5f cm2.(V.s)-1" % mobility
        
        #fitting msd <L^2> xz
        fit = np.polyfit( t[range_], msd_snap[2][range_], 1 )
        #plot msd <L^2> xz
        fit_fn = np.poly1d(fit)
        plt.plot(t[range_], fit_fn(t[range_]), color = "k")
        plt.plot(t, msd_snap[2], label = "<L^2> dir. XZ", linewidth=5, color="k")
        plt.plot(t[range_], fit_fn(t[range_]), color = "k")
        #mobility calculation
        mobility = fit[0] * 10**15 * 10**(-16) / (2 * 0.0000861728 * float(temperature) )
        print "Linear mob. from <L^2> in xz = \sum_i Pi*(iL)^2 is: %.5f cm2.(V.s)-1" % mobility
        
        #fitting msd <L^2> yy
        fit = np.polyfit( t[range_], msd_snap[4][range_], 1 )
        #plot msd <L^2> yy
        fit_fn = np.poly1d(fit)
        plt.plot(t[range_], fit_fn(t[range_]), color = "k")
        plt.plot(t, msd_snap[4], label = "<L^2> dir. YY", linewidth=5, color="g")
        plt.plot(t[range_], fit_fn(t[range_]), color = "k")
        #mobility calculation
        mobility = fit[0] * 10**15 * 10**(-16) / (2 * 0.0000861728 * float(temperature) )
        print "Linear mob. from <L^2> in yy = \sum_i Pi*(iL)^2 is: %.5f cm2.(V.s)-1" % mobility
        
    
        #fitting msd <L^2> yz
        fit = np.polyfit( t[range_], msd_snap[5][range_], 1 )
        #plot msd <L^2> yz
        fit_fn = np.poly1d(fit)
        plt.plot(t[range_], fit_fn(t[range_]), color = "k")
        plt.plot(t, msd_snap[5], label = "<L^2> dir. YZ", linewidth=5, color="m")
        plt.plot(t[range_], fit_fn(t[range_]), color = "k")
        #mobility calculation
        mobility = fit[0] * 10**15 * 10**(-16) / (2 * 0.0000861728 * float(temperature) )
        print "Linear mob. from <L^2> in yz = \sum_i Pi*(iL)^2 is: %.5f cm2.(V.s)-1" % mobility
    
    
        #fitting msd <L^2> zz
        fit = np.polyfit( t[range_], msd_snap[8][range_], 1 )
        #plot msd <L^2> zz
        fit_fn = np.poly1d(fit)
        plt.plot(t[range_], fit_fn(t[range_]), color = "k")
        plt.plot(t, msd_snap[8], label = "<L^2> dir. ZZ", linewidth=5, color="c")
        plt.plot(t[range_], fit_fn(t[range_]), color = "k")
        #mobility calculation
        mobility = fit[0] * 10**15 * 10**(-16) / (2 * 0.0000861728 * float(temperature) )
        print "Linear mob. from <L^2> in zz = \sum_i Pi*(iL)^2 is: %.5f cm2.(V.s)-1" % mobility
    
        plt.title(system)
        plt.xlabel('Time (fs)')
        plt.ylabel( r"MSD ($\AA^2$)")
        lgd = plt.legend(frameon=False, bbox_to_anchor=(1.1, 1.0))
        plt.tight_layout()
        plt.savefig('MSD-%s' % system, bbox_inches='tight', bbox_extra_artists=(lgd,))
        #plt.show()
        
        

            
