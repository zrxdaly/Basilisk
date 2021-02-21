#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 10:10:01 2020

OBJ: compare the buoyancy field between the no operation and different operation

@author: dai
"""

import matplotlib.pyplot as plt
import numpy as np
import glob as glob

# plt.rc('font', family='serif')
plt.rc('font',**{'family':'serif','serif':['Times']})
plt.rc('xtick', labelsize='15')
plt.rc('ytick', labelsize='15')
plt.rc('text', usetex=True)
lw = 2
Color_R2 = "#000000"         
Color_R3 = '#fcf18f'          
Color_R4 = '#377eb8'         
Color_R5 = '#008000'   
Color_R6 = '#084594'             # 6.25m
Color_R7 = '#ff7f00'             # 8m
Color_R8 = '#808080'             # 10m     

# Color_R1 = '#0000a0'             # 1m
# Color_R2 = "#000000"             # 2m
# Color_R3 = '#e41a1c'             # 3.125m
# Color_R4 = '#377eb8'             # 4m
# Color_R5 = '#008000'             # 5m

# Color_R9 = '#808080'             # 12.5m

#%%
dir1 = '/home/dai/software/basilisk/Exp_field_Basilisk/3DGoal_vrlab/'
Goal_dir = sorted(glob.glob(dir1 + 'goal_fan*/'))

FR_file = sorted(glob.glob(Goal_dir[0] + 't=*'))   # FR: Full rotation
HLR_file = sorted(glob.glob(Goal_dir[1] + 't=*'))  # HRR: half left rotation
HRR_file = sorted(glob.glob(Goal_dir[2] + 't=*'))  # HLR: Half right rotation
NR_file = sorted(glob.glob(Goal_dir[3] + 't=*'))   # NR: no rotation

# calculation of the max warming and reference warming
# in 2D
G = 9.81
T_ref = 273
KCR = 273.15
INV = 0.3
T_hub = G/T_ref*INV*10.5

def K2C(T_hub):
    return(T_hub*T_ref/G + T_ref - KCR)

TC_hub = K2C(T_hub)

Dt = 1    # s  the time interval
T0 = 30     # s  time length for no operation
T1 = 990    # s  time length for fan operation

#%% average the data at different height over time
def D2_N(num, Tstart, Tend):
    # 1234 represent different case
    FR_05 = sorted(glob.glob(Goal_dir[num] + 't=*y=005'))
    FR_10 = sorted(glob.glob(Goal_dir[num] + 't=*y=010'))
    FR_15 = sorted(glob.glob(Goal_dir[num] + 't=*y=015'))
    FR_20 = sorted(glob.glob(Goal_dir[num] + 't=*y=020'))
    FR_25 = sorted(glob.glob(Goal_dir[num] + 't=*y=025'))
    FR_30 = sorted(glob.glob(Goal_dir[num] + 't=*y=030'))
    FR_35 = sorted(glob.glob(Goal_dir[num] + 't=*y=035'))
    FR_40 = sorted(glob.glob(Goal_dir[num] + 't=*y=040'))
    FR_45 = sorted(glob.glob(Goal_dir[num] + 't=*y=045'))
    FR_50 = sorted(glob.glob(Goal_dir[num] + 't=*y=050'))
    FR_series_05 = []
    FR_series_10 = []
    FR_series_15 = []
    FR_series_20 = []
    FR_series_25 = []
    FR_series_30 = []
    FR_series_35 = []
    FR_series_40 = []
    FR_series_45 = []
    FR_series_50 = []
    
    FR_series_05_append = FR_series_05.append
    FR_series_10_append = FR_series_10.append
    FR_series_15_append = FR_series_15.append
    FR_series_20_append = FR_series_20.append
    FR_series_25_append = FR_series_25.append
    FR_series_30_append = FR_series_30.append
    FR_series_35_append = FR_series_35.append
    FR_series_40_append = FR_series_40.append
    FR_series_45_append = FR_series_45.append
    FR_series_50_append = FR_series_50.append

    SUM_FR_N_H = []
    for i in np.arange(Tstart, Tend):
        FR_Bu05N = np.loadtxt(FR_05[i], dtype='f',skiprows=2)
        FR_Bu10N = np.loadtxt(FR_10[i], dtype='f',skiprows=2)
        FR_Bu15N = np.loadtxt(FR_15[i], dtype='f',skiprows=2)
        FR_Bu20N = np.loadtxt(FR_20[i], dtype='f',skiprows=2)
        FR_Bu25N = np.loadtxt(FR_25[i], dtype='f',skiprows=2)
        FR_Bu30N = np.loadtxt(FR_30[i], dtype='f',skiprows=2)
        FR_Bu35N = np.loadtxt(FR_35[i], dtype='f',skiprows=2)
        FR_Bu40N = np.loadtxt(FR_40[i], dtype='f',skiprows=2)
        FR_Bu45N = np.loadtxt(FR_45[i], dtype='f',skiprows=2)
        FR_Bu50N = np.loadtxt(FR_50[i], dtype='f',skiprows=2)

        FR_series_05_append(K2C(FR_Bu05N.mean()))
        FR_series_10_append(K2C(FR_Bu10N.mean()))
        FR_series_15_append(K2C(FR_Bu15N.mean()))
        FR_series_20_append(K2C(FR_Bu20N.mean()))
        FR_series_25_append(K2C(FR_Bu25N.mean()))
        FR_series_30_append(K2C(FR_Bu30N.mean()))
        FR_series_35_append(K2C(FR_Bu35N.mean()))
        FR_series_40_append(K2C(FR_Bu40N.mean()))
        FR_series_45_append(K2C(FR_Bu45N.mean()))
        FR_series_50_append(K2C(FR_Bu50N.mean()))
    
    SUM_FR_N_H = [np.average(FR_series_05),
                  np.average(FR_series_10),
                  np.average(FR_series_15),
                  np.average(FR_series_20),
                  np.average(FR_series_25),
                  np.average(FR_series_30),
                  np.average(FR_series_35),
                  np.average(FR_series_40),
                  np.average(FR_series_45),
                  np.average(FR_series_50)]
    
    return(SUM_FR_N_H)

SUM_FR_N_H = D2_N(0, 0, T0)
SUM_FR_N = np.mean(SUM_FR_N_H)

SUM_HLR_N_H = D2_N(1, 0, T0)
SUM_HLR_N = np.mean(SUM_HLR_N_H)

SUM_HRR_N_H = D2_N(2, 0, T0)
SUM_HRR_N = np.mean(SUM_HRR_N_H)

SUM_NR_N_H = D2_N(3, 0, T0)
SUM_NR_N = np.mean(SUM_NR_N_H)

# more cases added after the more data coming in

#%%

def D2(num, Tstart, Tend, SUM_FR_N_H):
    # 1234 represent different case
    FR_05 = sorted(glob.glob(Goal_dir[num] + 't=*y=005'))
    FR_10 = sorted(glob.glob(Goal_dir[num] + 't=*y=010'))
    FR_15 = sorted(glob.glob(Goal_dir[num] + 't=*y=015'))
    FR_20 = sorted(glob.glob(Goal_dir[num] + 't=*y=020'))
    FR_25 = sorted(glob.glob(Goal_dir[num] + 't=*y=025'))
    FR_30 = sorted(glob.glob(Goal_dir[num] + 't=*y=030'))
    FR_35 = sorted(glob.glob(Goal_dir[num] + 't=*y=035'))
    FR_40 = sorted(glob.glob(Goal_dir[num] + 't=*y=040'))
    FR_45 = sorted(glob.glob(Goal_dir[num] + 't=*y=045'))
    FR_50 = sorted(glob.glob(Goal_dir[num] + 't=*y=050'))
    FR_series_05 = []
    FR_series_10 = []
    FR_series_15 = []
    FR_series_20 = []
    FR_series_25 = []
    FR_series_30 = []
    FR_series_35 = []
    FR_series_40 = []
    FR_series_45 = []
    FR_series_50 = []

    FR_series_05_append = FR_series_05.append
    FR_series_10_append = FR_series_10.append
    FR_series_15_append = FR_series_15.append
    FR_series_20_append = FR_series_20.append
    FR_series_25_append = FR_series_25.append
    FR_series_30_append = FR_series_30.append
    FR_series_35_append = FR_series_35.append
    FR_series_40_append = FR_series_40.append
    FR_series_45_append = FR_series_45.append
    FR_series_50_append = FR_series_50.append

    D2 = []
    for i in np.arange(Tstart, Tend):
        FR_Bu05N = np.loadtxt(FR_05[i], dtype='f',skiprows=2)
        FR_Bu10N = np.loadtxt(FR_10[i], dtype='f',skiprows=2)
        FR_Bu15N = np.loadtxt(FR_15[i], dtype='f',skiprows=2)
        FR_Bu20N = np.loadtxt(FR_20[i], dtype='f',skiprows=2)
        FR_Bu25N = np.loadtxt(FR_25[i], dtype='f',skiprows=2)
        FR_Bu30N = np.loadtxt(FR_30[i], dtype='f',skiprows=2)
        FR_Bu35N = np.loadtxt(FR_35[i], dtype='f',skiprows=2)
        FR_Bu40N = np.loadtxt(FR_40[i], dtype='f',skiprows=2)
        FR_Bu45N = np.loadtxt(FR_45[i], dtype='f',skiprows=2)
        FR_Bu50N = np.loadtxt(FR_50[i], dtype='f',skiprows=2)
        # substraction
        FR_series_05_append(K2C(FR_Bu05N.mean()) - SUM_FR_N_H[0])
        FR_series_10_append(K2C(FR_Bu10N.mean()) - SUM_FR_N_H[1])
        FR_series_15_append(K2C(FR_Bu15N.mean()) - SUM_FR_N_H[2])
        FR_series_20_append(K2C(FR_Bu20N.mean()) - SUM_FR_N_H[3])
        FR_series_25_append(K2C(FR_Bu25N.mean()) - SUM_FR_N_H[4])
        FR_series_30_append(K2C(FR_Bu30N.mean()) - SUM_FR_N_H[5])
        FR_series_35_append(K2C(FR_Bu35N.mean()) - SUM_FR_N_H[6])
        FR_series_40_append(K2C(FR_Bu40N.mean()) - SUM_FR_N_H[7])
        FR_series_45_append(K2C(FR_Bu45N.mean()) - SUM_FR_N_H[8])
        FR_series_50_append(K2C(FR_Bu50N.mean()) - SUM_FR_N_H[9])
    
    D2 = np.array([FR_series_05, FR_series_10, FR_series_15, FR_series_20,
          FR_series_25, FR_series_30, FR_series_35, FR_series_40,
          FR_series_45, FR_series_50])
    
    return(D2)

#%%
FR_2D = D2(0, T0, T1, SUM_FR_N_H)
HLR_2D = D2(1, T0, T1, SUM_HLR_N_H)
HRR_2D = D2(2, T0, T1, SUM_HRR_N_H)
NR_2D = D2(3, T0, T1, SUM_NR_N_H)
# this require some time
#%% plot the statistics

FR_DS = FR_2D.mean(axis = 0)
HLR_DS = HLR_2D.mean(axis = 0)
HRR_DS = HRR_2D.mean(axis = 0)
NR_DS = NR_2D.mean(axis = 0)

#%% the efficiency then can be expressed as 

EF_FR = FR_DS/(TC_hub - SUM_FR_N)
EF_HRR = HRR_DS/(TC_hub - SUM_HRR_N)
EF_HLR = HLR_DS/(TC_hub - SUM_HLR_N)
EF_NR = NR_DS/(TC_hub - SUM_NR_N)

# EF_FR_H = np.array(FR_2D)/(TC_hub - np.array(SUM_FR_N_H))
# EF_HRR_H = np.array(HRR_2D)/(TC_hub - np.array(SUM_HRR_N_H))
# EF_HLR_H = np.array(HLR_2D)/(TC_hub - np.array(SUM_HLR_N_H))
# EF_NR_H = np.array(NR_2D)/(TC_hub - np.array(SUM_NR_N_H))


#%%
T = np.arange(T0,T1, Dt)

fig1 = plt.figure(figsize=(6.4, 4.8))
ax = fig1.add_subplot(1,1,1)

h1 = plt.plot(T, EF_FR, "-", label="full",linewidth = lw,color = Color_R2)  
h2 = plt.plot(T, EF_HLR, "-", label="half left",linewidth = lw,color = Color_R3) 
h3 = plt.plot(T, EF_HRR, "-", label="half right",linewidth = lw,color = Color_R4) 
h4 = plt.plot(T, EF_NR, "-", label="fixed right",linewidth = lw,color = Color_R5)  
# ax.set_xlim([230,275])
# ax.set_ylim([0,300])
ax.set_xlabel('Time [s]',fontsize=18)
ax.set_ylabel(r"$\Delta T/\Delta T_{max}$",fontsize=18)
plt.title("whole domain, time serie comparision",fontsize=18)
legend = plt.legend(loc='lower right', frameon=False,fontsize = 18,ncol = 1)
plt.tight_layout()
plt.savefig("domain3D.pdf")

#%% plot the height statistics

def list2Series(FR_2D,SUM_FR_N_H):
    FR_2D = np.array(FR_2D)
    FR_H1 = (FR_2D[0]+FR_2D[1])/2
    FR_H2 = (FR_2D[2]+FR_2D[3])/2
    FR_H3 = (FR_2D[4]+FR_2D[5])/2
    FR_H4 = (FR_2D[6]+FR_2D[7])/2
    FR_H5 = (FR_2D[8]+FR_2D[9])/2
    
    FR_NH1 = (SUM_FR_N_H[0]+SUM_FR_N_H[1])/2
    FR_NH2 = (SUM_FR_N_H[2]+SUM_FR_N_H[3])/2
    FR_NH3 = (SUM_FR_N_H[4]+SUM_FR_N_H[5])/2
    FR_NH4 = (SUM_FR_N_H[6]+SUM_FR_N_H[7])/2
    FR_NH5 = (SUM_FR_N_H[8]+SUM_FR_N_H[9])/2
    
    EF_FR_H1 = FR_H1/(TC_hub - FR_NH1)
    EF_FR_H2 = FR_H2/(TC_hub - FR_NH2)
    EF_FR_H3 = FR_H3/(TC_hub - FR_NH3)
    EF_FR_H4 = FR_H4/(TC_hub - FR_NH4)
    EF_FR_H5 = FR_H5/(TC_hub - FR_NH5)
    return(EF_FR_H1,EF_FR_H2,EF_FR_H3,EF_FR_H4,EF_FR_H5)


FR_H1,FR_H2,FR_H3,FR_H4,FR_H5 = list2Series(FR_2D,SUM_FR_N_H)
HLR_H1,HLR_H2,HLR_H3,HLR_H4,HLR_H5 = list2Series(HLR_2D,SUM_HLR_N_H)
HRR_H1,HRR_H2,HRR_H3,HRR_H4,HRR_H5 = list2Series(HRR_2D,SUM_HRR_N_H)
NR_H1,NR_H2,NR_H3,NR_H4,NR_H5 = list2Series(NR_2D,SUM_NR_N_H)


#%%
fig2 = plt.figure(figsize=(6.4, 4.8))
ax = fig2.add_subplot(1,1,1)
h1 = plt.plot(T, FR_H1, "-", label="0.75m",linewidth = lw,color = Color_R2)  
h2 = plt.plot(T, FR_H2, "-", label="1.75m",linewidth = lw,color = Color_R3) 
h3 = plt.plot(T, FR_H3, "-", label="2.75m",linewidth = lw,color = Color_R4) 
h4 = plt.plot(T, FR_H4, "-", label="3.75m",linewidth = lw,color = Color_R5)  
h5 = plt.plot(T, FR_H5, "-", label="4.75m",linewidth = lw,color = Color_R6)  
ax.set_ylim([0,0.013])
ax.set_xlabel('Time [s]',fontsize=18)
ax.set_ylabel(r"$\Delta T/\Delta T_{max}$",fontsize=18)
plt.title("height, time serie, full rotation",fontsize=18)
legend = plt.legend(loc='lower center', frameon=False,fontsize = 18,ncol = 2)
plt.tight_layout()
plt.savefig("Full_3D.pdf")


fig2 = plt.figure(figsize=(6.4, 4.8))
ax = fig2.add_subplot(1,1,1)
h1 = plt.plot(T, HRR_H1, "-", label="0.75m",linewidth = lw,color = Color_R2)  
h2 = plt.plot(T, HRR_H2, "-", label="1.75m",linewidth = lw,color = Color_R3) 
h3 = plt.plot(T, HRR_H3, "-", label="2.75m",linewidth = lw,color = Color_R4) 
h4 = plt.plot(T, HRR_H4, "-", label="3.75m",linewidth = lw,color = Color_R5)  
h5 = plt.plot(T, HRR_H5, "-", label="4.75m",linewidth = lw,color = Color_R6)  
ax.set_ylim([0,0.013])
ax.set_xlabel('Time [s]',fontsize=18)
ax.set_ylabel(r"$\Delta T/\Delta T_{max}$",fontsize=18)
plt.title("height, time serie, right half",fontsize=18)
legend = plt.legend(loc='lower center', frameon=False,fontsize = 18,ncol = 2)
plt.tight_layout()
plt.savefig("HR_3D.pdf")


fig2 = plt.figure(figsize=(6.4, 4.8))
ax = fig2.add_subplot(1,1,1)
h1 = plt.plot(T, HLR_H1, "-", label="0.75m",linewidth = lw,color = Color_R2)  
h2 = plt.plot(T, HLR_H2, "-", label="1.75m",linewidth = lw,color = Color_R3) 
h3 = plt.plot(T, HLR_H3, "-", label="2.75m",linewidth = lw,color = Color_R4) 
h4 = plt.plot(T, HLR_H4, "-", label="3.75m",linewidth = lw,color = Color_R5)  
h5 = plt.plot(T, HLR_H5, "-", label="4.75m",linewidth = lw,color = Color_R6)  
ax.set_ylim([0,0.013])
ax.set_xlabel('Time [s]',fontsize=18)
ax.set_ylabel(r"$\Delta T/\Delta T_{max}$",fontsize=18)
plt.title("height, time serie, left half",fontsize=18)
legend = plt.legend(loc='upper left', frameon=False,fontsize = 18,ncol = 2)
plt.tight_layout()
plt.savefig("HL_3D.pdf")


fig2 = plt.figure(figsize=(6.4, 4.8))
ax = fig2.add_subplot(1,1,1)
h1 = plt.plot(T, NR_H1, "-", label="0.75m",linewidth = lw,color = Color_R2)  
h2 = plt.plot(T, NR_H2, "-", label="1.75m",linewidth = lw,color = Color_R3) 
h3 = plt.plot(T, NR_H3, "-", label="2.75m",linewidth = lw,color = Color_R4) 
h4 = plt.plot(T, NR_H4, "-", label="3.75m",linewidth = lw,color = Color_R5)  
h5 = plt.plot(T, NR_H5, "-", label="4.75m",linewidth = lw,color = Color_R6)  
ax.set_ylim([0,0.013])
ax.set_xlabel('Time [s]',fontsize=18)
ax.set_ylabel(r"$\Delta T/\Delta T_{max}$",fontsize=18)
plt.title("height, time serie, fixed right",fontsize=18)
legend = plt.legend(loc='upper left', frameon=False,fontsize = 18,ncol = 2)
plt.tight_layout()
plt.savefig("NR_3D.pdf")

#%%
fig2 = plt.figure(figsize=(6.4, 4.8))
ax = fig2.add_subplot(1,1,1)
h1 = plt.plot(T, FR_H1, "-", label="full",linewidth = lw,color = Color_R2)
h2 = plt.plot(T, HRR_H1, "-", label="right",linewidth = lw,color = Color_R3)
h3 = plt.plot(T, HLR_H1, "-", label="left",linewidth = lw,color = Color_R4)
h4 = plt.plot(T, NR_H1, "-", label="fix",linewidth = lw,color = Color_R5)
ax.set_ylim([0,0.012])
ax.set_xlabel('Time [s]',fontsize=18)
ax.set_ylabel(r"$\Delta T/\Delta T_{max}$",fontsize=18)
plt.title("height, time serie, height 0.75m",fontsize=18)
legend = plt.legend(loc='best', frameon=False,fontsize = 18,ncol = 2)
plt.tight_layout()
plt.savefig("RO_H075.pdf")


fig2 = plt.figure(figsize=(6.4, 4.8))
ax = fig2.add_subplot(1,1,1)
h1 = plt.plot(T, FR_H2, "-", label="full",linewidth = lw,color = Color_R2)
h2 = plt.plot(T, HRR_H2, "-", label="right",linewidth = lw,color = Color_R3)
h3 = plt.plot(T, HLR_H2, "-", label="left",linewidth = lw,color = Color_R4)
h4 = plt.plot(T, NR_H2, "-", label="fix",linewidth = lw,color = Color_R5)
ax.set_ylim([0,0.012])
ax.set_xlabel('Time [s]',fontsize=18)
ax.set_ylabel(r"$\Delta T/\Delta T_{max}$",fontsize=18)
plt.title("height, time serie, height 1.75m",fontsize=18)
legend = plt.legend(loc='best', frameon=False,fontsize = 18,ncol = 2)
plt.tight_layout()
plt.savefig("RO_H175.pdf")


fig2 = plt.figure(figsize=(6.4, 4.8))
ax = fig2.add_subplot(1,1,1)
h1 = plt.plot(T, FR_H3, "-", label="full",linewidth = lw,color = Color_R2)
h2 = plt.plot(T, HRR_H3, "-", label="right",linewidth = lw,color = Color_R3)
h3 = plt.plot(T, HLR_H3, "-", label="left",linewidth = lw,color = Color_R4)
h4 = plt.plot(T, NR_H3, "-", label="fix",linewidth = lw,color = Color_R5)
ax.set_ylim([0,0.012])
ax.set_xlabel('Time [s]',fontsize=18)
ax.set_ylabel(r"$\Delta T/\Delta T_{max}$",fontsize=18)
plt.title("height, time serie, height 2.75m ",fontsize=18)
legend = plt.legend(loc='best', frameon=False,fontsize = 18,ncol = 2)
plt.tight_layout()
plt.savefig("RO_H275.pdf")


fig2 = plt.figure(figsize=(6.4, 4.8))
ax = fig2.add_subplot(1,1,1)
h1 = plt.plot(T, FR_H4, "-", label="full",linewidth = lw,color = Color_R2)
h2 = plt.plot(T, HRR_H4, "-", label="right",linewidth = lw,color = Color_R3)
h3 = plt.plot(T, HLR_H4, "-", label="left",linewidth = lw,color = Color_R4)
h4 = plt.plot(T, NR_H4, "-", label="fix",linewidth = lw,color = Color_R5)
ax.set_ylim([0,0.012])
ax.set_xlabel('Time [s]',fontsize=18)
ax.set_ylabel(r"$\Delta T/\Delta T_{max}$",fontsize=18)
plt.title("height, time serie, height 3.75m",fontsize=18)
legend = plt.legend(loc='best', frameon=False,fontsize = 18,ncol = 2)
plt.tight_layout()
plt.savefig("RO_H375.pdf")

fig2 = plt.figure(figsize=(6.4, 4.8))
ax = fig2.add_subplot(1,1,1)
h1 = plt.plot(T, FR_H5, "-", label="full",linewidth = lw,color = Color_R2)
h2 = plt.plot(T, HRR_H5, "-", label="right",linewidth = lw,color = Color_R3)
h3 = plt.plot(T, HLR_H5, "-", label="left",linewidth = lw,color = Color_R4)
h4 = plt.plot(T, NR_H5, "-", label="fix",linewidth = lw,color = Color_R5)
ax.set_ylim([0,0.012])
ax.set_xlabel('Time [s]',fontsize=18)
ax.set_ylabel(r"$\Delta T/\Delta T_{max}$",fontsize=18)
plt.title("height, time serie, height 4.75m",fontsize=18)
legend = plt.legend(loc='best', frameon=False,fontsize = 18,ncol = 2)
plt.tight_layout()
plt.savefig("RO_H475.pdf")













