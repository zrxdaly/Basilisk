#!/usr/bin/env python3
# -*)/2- coding: utf-8 -*-
"""
Created on Thu Nov 19 10:10:01 2020

OBJ: compare the buoyancy field between the no operation and different operation

@author: dai
"""

import matplotlib.pyplot as plt
import numpy as np
import glob as glob

plt.rc('font',**{'family':'serif','serif':['Times']})
plt.rc('xtick', labelsize='15')
plt.rc('ytick', labelsize='15')
plt.rc('text', usetex=True)
lw = 2
Color_R2 = "#000000"
Color_R3 = '#fcf18f'
Color_R4 = '#377eb8'
Color_R5 = '#008000'

#%%
dir1 = '/home/dai/software/basilisk/Exp_field_Basilisk/'
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

#%% 
Dt = 1    # s  the time interval
T0 = 60     # s  time length for no operation
T1 = 1020-60    # s  time length for fan operation

def f2sum(FR_file):
    FR_Buo_N = np.loadtxt(FR_file[0], dtype='f',skiprows=1)
    # average over no operation
    for i in np.arange(int(T0/Dt-1)):
        FR_Buo_N = FR_Buo_N + np.loadtxt(FR_file[i+1], dtype='f',skiprows=1)
    
    # no operation 2D average
    FR_Buo_N_AV = FR_Buo_N/(T0/Dt)
    # no operation sum
    SUM_FR_N = FR_Buo_N_AV.mean()
    SUM_FR_N = K2C(SUM_FR_N)
    SUM_FR_N_H = FR_Buo_N_AV.mean(axis = 1)
    SUM_FR_N_H = K2C(SUM_FR_N_H)   # to celsius
    return(SUM_FR_N,SUM_FR_N_H)

SUM_FR_N,SUM_FR_N_H = f2sum(FR_file)
SUM_HLR_N,SUM_HLR_N_H = f2sum(HLR_file)
SUM_HRR_N,SUM_HRR_N_H = f2sum(HRR_file)
# SUM_NR_N,SUM_NR_N_H = f2sum(NR_file)

#%%
# average over whole domain
def do_A(FR_file, SUM_FR_N):
    SUM_FR_Y = []
    for i in np.arange(int(T0/Dt),int((T0+T1)/Dt)):
        FR_Buo_Y = np.loadtxt(FR_file[i], dtype='f',skiprows=1)
        SUM_FR_Y.append(K2C(FR_Buo_Y.mean())-SUM_FR_N)  # mean over domain and transfer to celsius
    return(SUM_FR_Y)

# average over different height
def do_Height_A(FR_file,SUM_FR_N_H):
    SUM_FR_YH = []
    for i in np.arange(int(T0/Dt),int((T0+T1)/Dt)):
        FR_Buo_Y = np.loadtxt(FR_file[i], dtype='f',skiprows=1)
        SUM_FR_YH.append(K2C(FR_Buo_Y.mean(axis = 1))-SUM_FR_N_H)
    return(SUM_FR_YH)

SUM_FR_Y = do_A(FR_file, SUM_FR_N)
SUM_HRR_Y = do_A(HRR_file, SUM_HRR_N)
SUM_HLR_Y = do_A(HLR_file, SUM_HLR_N)
# SUM_NR_Y = do_A(NR_file, SUM_NR_N)

SUM_FR_YH = do_Height_A(FR_file,SUM_FR_N_H)
SUM_HRR_YH = do_Height_A(HRR_file,SUM_HRR_N_H)
SUM_HLR_YH = do_Height_A(HLR_file,SUM_HLR_N_H)
# SUM_NR_YH = do_Height_A(NR_file,SUM_NR_N_H)

#%% the efficiency then can be expressed as 
EF_FR = SUM_FR_Y/(TC_hub - SUM_FR_N)
EF_HRR = SUM_HRR_Y/(TC_hub - SUM_HRR_N)
EF_HLR = SUM_HLR_Y/(TC_hub - SUM_HLR_N)
# EF_NR = SUM_NR_Y/(TC_hub - SUM_NR_N)

# over different height
EF_FR_H = SUM_FR_YH/(TC_hub - SUM_FR_N_H)
EF_HRR_H = SUM_HRR_YH/(TC_hub - SUM_HRR_N_H)
EF_HLR_H = SUM_HLR_YH/(TC_hub - SUM_HLR_N_H)
# EF_NR_H = SUM_NR_YH/(TC_hub - SUM_NR_N_H)

#%% plot the statistics
# time series
T = np.arange(T0,T1+T0, Dt)

fig1 = plt.figure(figsize=(6.4, 4.8))
ax = fig1.add_subplot(1,1,1)

h1 = plt.plot(T, EF_FR, "-", label="full",linewidth = lw,color = Color_R2)  
h2 = plt.plot(T, EF_HRR, "-", label="half left",linewidth = lw,color = Color_R3) 
h3 = plt.plot(T, EF_HLR, "-", label="half right",linewidth = lw,color = Color_R4) 
# h4 = plt.plot(T, EF_NR, "-", label="fixed right",linewidth = lw,color = Color_R5)  
# ax.set_xlim([230,275])
ax.set_ylim([0,1])
ax.set_xlabel('Time [s]',fontsize=18)
ax.set_ylabel(r"$\Delta T/\Delta T_{max}$",fontsize=18)
plt.title("whole domain, time serie comparision",fontsize=18)
legend = plt.legend(loc='upper left', frameon=False,fontsize = 18,ncol = 1)
plt.tight_layout()
plt.savefig("domain.png")

#%% plot the height statistics
def list2Series(EF_FR_H):
    FR_H1 = np.zeros(len(EF_FR_H))
    FR_H2 = np.zeros(len(EF_FR_H))
    FR_H3 = np.zeros(len(EF_FR_H))
    FR_H4 = np.zeros(len(EF_FR_H))
    FR_H5 = np.zeros(len(EF_FR_H))
    for i in np.arange(len(EF_FR_H)):
        FR_H1[i] = (EF_FR_H[i][0]+EF_FR_H[i][1])/2
        FR_H2[i] = (EF_FR_H[i][2]+EF_FR_H[i][3])/2
        FR_H3[i] = (EF_FR_H[i][4]+EF_FR_H[i][5])/2
        FR_H4[i] = (EF_FR_H[i][6]+EF_FR_H[i][7])/2
        FR_H5[i] = (EF_FR_H[i][8]+EF_FR_H[i][9])/2
    
    return(FR_H1,FR_H2,FR_H3,FR_H4,FR_H5)

FR_H1,FR_H2,FR_H3,FR_H4,FR_H5 = list2Series(EF_FR_H)
HRR_H1,HRR_H2,HRR_H3,HRR_H4,HRR_H5 = list2Series(EF_HRR_H)
HLR_H1,HLR_H2,HLR_H3,HLR_H4,HLR_H5 = list2Series(EF_HLR_H)
# NR_H1,NR_H2,NR_H3,NR_H4,NR_H5 = list2Series(EF_NR_H)

#%%
fig2 = plt.figure(figsize=(6.4, 4.8))
ax = fig2.add_subplot(1,1,1)
h1 = plt.plot(T, FR_H1, "-", label="0.75m",linewidth = lw,color = Color_R2)  
h2 = plt.plot(T, FR_H2, "-", label="1.75m",linewidth = lw,color = Color_R3) 
h3 = plt.plot(T, FR_H3, "-", label="2.75m",linewidth = lw,color = Color_R4) 
h4 = plt.plot(T, FR_H4, "-", label="3.75m",linewidth = lw,color = Color_R5)  
h5 = plt.plot(T, FR_H5, "-", label="4.75m",linewidth = lw,color = Color_R5)  
ax.set_ylim([0,1])
ax.set_xlabel('Time [s]',fontsize=18)
ax.set_ylabel(r"$\Delta T/\Delta T_{max}$",fontsize=18)
plt.title("height, time serie, full rotation",fontsize=18)
legend = plt.legend(loc='upper left', frameon=False,fontsize = 18,ncol = 1)
plt.tight_layout()
plt.savefig("Full.png")


fig2 = plt.figure(figsize=(6.4, 4.8))
ax = fig2.add_subplot(1,1,1)
h1 = plt.plot(T, HRR_H1, "-", label="0.75m",linewidth = lw,color = Color_R2)  
h2 = plt.plot(T, HRR_H2, "-", label="1.75m",linewidth = lw,color = Color_R3) 
h3 = plt.plot(T, HRR_H3, "-", label="2.75m",linewidth = lw,color = Color_R4) 
h4 = plt.plot(T, HRR_H4, "-", label="3.75m",linewidth = lw,color = Color_R5)  
h5 = plt.plot(T, HRR_H5, "-", label="4.75m",linewidth = lw,color = Color_R5)  
ax.set_ylim([0,1])
ax.set_xlabel('Time [s]',fontsize=18)
ax.set_ylabel(r"$\Delta T/\Delta T_{max}$",fontsize=18)
plt.title("height, time serie, left half",fontsize=18)
legend = plt.legend(loc='upper left', frameon=False,fontsize = 18,ncol = 1)
plt.tight_layout()
plt.savefig("HL.png")


fig2 = plt.figure(figsize=(6.4, 4.8))
ax = fig2.add_subplot(1,1,1)
h1 = plt.plot(T, HLR_H1, "-", label="0.75m",linewidth = lw,color = Color_R2)  
h2 = plt.plot(T, HLR_H2, "-", label="1.75m",linewidth = lw,color = Color_R3) 
h3 = plt.plot(T, HLR_H3, "-", label="2.75m",linewidth = lw,color = Color_R4) 
h4 = plt.plot(T, HLR_H4, "-", label="3.75m",linewidth = lw,color = Color_R5)  
h5 = plt.plot(T, HLR_H5, "-", label="4.75m",linewidth = lw,color = Color_R5)  
ax.set_ylim([0,1])
ax.set_xlabel('Time [s]',fontsize=18)
ax.set_ylabel(r"$\Delta T/\Delta T_{max}$",fontsize=18)
plt.title("height, time serie, right half",fontsize=18)
legend = plt.legend(loc='upper left', frameon=False,fontsize = 18,ncol = 1)
plt.tight_layout()
plt.savefig("HR.png")


# fig2 = plt.figure(figsize=(6.4, 4.8))
# ax = fig2.add_subplot(1,1,1)
# h1 = plt.plot(T, NR_H1, "-", label="0.75m",linewidth = lw,color = Color_R2)  
# h2 = plt.plot(T, NR_H2, "-", label="1.75m",linewidth = lw,color = Color_R3) 
# h3 = plt.plot(T, NR_H3, "-", label="2.75m",linewidth = lw,color = Color_R4) 
# h4 = plt.plot(T, NR_H4, "-", label="3.75m",linewidth = lw,color = Color_R5)  
# h5 = plt.plot(T, NR_H5, "-", label="4.75m",linewidth = lw,color = Color_R5)  
# # ax.set_ylim([0,110])
# ax.set_xlabel('Time [s]',fontsize=18)
# ax.set_ylabel(r"$\Delta T/\Delta T_{max}$",fontsize=18)
# plt.title("height, time serie, fixed right",fontsize=18)
# legend = plt.legend(loc='upper left', frameon=False,fontsize = 18,ncol = 1)
# plt.tight_layout()
# plt.savefig("NR.png")














