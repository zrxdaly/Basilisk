#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 11:42:26 2020

OBJ: The evaluation of the rotation parameter

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

#%% read the folder
dir1 = "/home/dai/software/basilisk/Exp_field_Basilisk/2DGoal_FACTOR/RO_Period/fan_goal_180R/"
Goal_dir = sorted(glob.glob(dir1 + 'R*/goal_fan/'))

R140_file = sorted(glob.glob(Goal_dir[0] + 't=*'))
R160_file = sorted(glob.glob(Goal_dir[1] + 't=*'))
R180_file = sorted(glob.glob(Goal_dir[2] + 't=*'))
R240_file = sorted(glob.glob(Goal_dir[3] + 't=*'))
R288_file = sorted(glob.glob(Goal_dir[4] + 't=*'))
R360_file = sorted(glob.glob(Goal_dir[5] + 't=*'))

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
T0 = 60     # s  time length for no operation
T1 = 1500-60    # s  time length for fan operation

#%% no operation func: read files to sum up the field
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

#%% 
R140_N,R140_N_H = f2sum(R140_file)
R160_N,R160_N_H = f2sum(R160_file)
R180_N,R180_N_H = f2sum(R180_file)
R240_N,R240_N_H = f2sum(R240_file)
R288_N,R288_N_H = f2sum(R288_file)
R360_N,R360_N_H = f2sum(R360_file)

#%% average over whole domain
def do_A(FR_file, SUM_FR_N):
    SUM_FR_Y = []
    for i in np.arange(int(T0/Dt),int((T0+T1)/Dt)):
        FR_Buo_Y = np.loadtxt(FR_file[i], dtype='f',skiprows=1)
        SUM_FR_Y.append(K2C(FR_Buo_Y.mean())-SUM_FR_N)
    return(SUM_FR_Y)

# average over different height
def do_Height_A(FR_file,SUM_FR_N_H):
    SUM_FR_YH = []
    for i in np.arange(int(T0/Dt),int((T0+T1)/Dt)):
        FR_Buo_Y = np.loadtxt(FR_file[i], dtype='f',skiprows=1)
        SUM_FR_YH.append(K2C(FR_Buo_Y.mean(axis = 1))-SUM_FR_N_H)
    return(SUM_FR_YH)

#%% the relative effect of warning 
SUM_R140_Y = do_A(R140_file, R140_N)
SUM_R160_Y = do_A(R160_file, R160_N)
SUM_R180_Y = do_A(R180_file, R180_N)
SUM_R240_Y = do_A(R240_file, R240_N)
SUM_R288_Y = do_A(R288_file, R288_N)
SUM_R360_Y = do_A(R360_file, R360_N)

SUM_R140_YH = do_Height_A(R140_file,R140_N_H)
SUM_R160_YH = do_Height_A(R160_file,R160_N_H)
SUM_R180_YH = do_Height_A(R180_file,R180_N_H)
SUM_R240_YH = do_Height_A(R240_file,R240_N_H)
SUM_R288_YH = do_Height_A(R288_file,R288_N_H)
SUM_R360_YH = do_Height_A(R360_file,R360_N_H)

#%% the efficiency then can be expressed as 
EF140 = SUM_R140_Y/(TC_hub - R140_N)
EF160 = SUM_R160_Y/(TC_hub - R160_N)
EF180 = SUM_R180_Y/(TC_hub - R180_N)
EF240 = SUM_R240_Y/(TC_hub - R240_N)
EF288 = SUM_R288_Y/(TC_hub - R288_N)
EF360 = SUM_R360_Y/(TC_hub - R360_N)

EF_R140_H = SUM_R140_YH/(TC_hub - R140_N_H)
EF_R160_H = SUM_R160_YH/(TC_hub - R160_N_H)
EF_R180_H = SUM_R180_YH/(TC_hub - R180_N_H)
EF_R240_H = SUM_R240_YH/(TC_hub - R240_N_H)
EF_R288_H = SUM_R288_YH/(TC_hub - R288_N_H)
EF_R360_H = SUM_R360_YH/(TC_hub - R360_N_H)

#%% plot the statistics
# time series
T = np.arange(T0,T1+T0, Dt)
fig1 = plt.figure(figsize=(6.4, 4.8))
ax = fig1.add_subplot(1,1,1)
h1 = plt.plot(T, EF140, "-", label="R140",linewidth = lw,color = Color_R2)  
h2 = plt.plot(T, EF160, "-", label="R160",linewidth = lw,color = Color_R3) 
h3 = plt.plot(T, EF180, "-", label="R180",linewidth = lw,color = Color_R4) 
h4 = plt.plot(T, EF240, "-", label="R240",linewidth = lw,color = Color_R5)  
h4 = plt.plot(T, EF288, "-", label="R288",linewidth = lw,color = Color_R6)  
h4 = plt.plot(T, EF360, "-", label="R360",linewidth = lw,color = Color_R7)  
ax.set_xlabel('Time [s]',fontsize=18)
ax.set_ylabel(r"$\Delta T/\Delta T_{max}$",fontsize=18)
plt.title("whole domain, time serie comparision",fontsize=18)
legend = plt.legend(loc='upper left', frameon=False,fontsize = 18,ncol = 1)
plt.tight_layout()
plt.savefig("RO_domain.pdf")

#%% plot the height statistics
# add the data of two height together
def list2Series(SUM_FR_YH):
    FR_H1 = np.zeros(len(SUM_FR_YH))
    FR_H2 = np.zeros(len(SUM_FR_YH))
    FR_H3 = np.zeros(len(SUM_FR_YH))
    FR_H4 = np.zeros(len(SUM_FR_YH))
    FR_H5 = np.zeros(len(SUM_FR_YH))
    for i in np.arange(len(SUM_FR_YH)):
        FR_H1[i] = (SUM_FR_YH[i][0]+SUM_FR_YH[i][1])/2
        FR_H2[i] = (SUM_FR_YH[i][2]+SUM_FR_YH[i][3])/2
        FR_H3[i] = (SUM_FR_YH[i][4]+SUM_FR_YH[i][5])/2
        FR_H4[i] = (SUM_FR_YH[i][6]+SUM_FR_YH[i][7])/2
        FR_H5[i] = (SUM_FR_YH[i][8]+SUM_FR_YH[i][9])/2
    
    return(FR_H1,FR_H2,FR_H3,FR_H4,FR_H5)

#%% 
R030_H1,R030_H2,R030_H3,R030_H4,R030_H5 = list2Series(EF_R140_H)
R060_H1,R060_H2,R060_H3,R060_H4,R060_H5 = list2Series(EF_R160_H)
R090_H1,R090_H2,R090_H3,R090_H4,R090_H5 = list2Series(EF_R180_H)
R120_H1,R120_H2,R120_H3,R120_H4,R120_H5 = list2Series(EF_R240_H)
R150_H1,R150_H2,R150_H3,R150_H4,R150_H5 = list2Series(EF_R288_H)
R180_H1,R180_H2,R180_H3,R180_H4,R180_H5 = list2Series(EF_R360_H)


#%%
fig2 = plt.figure(figsize=(6.4, 4.8))
ax = fig2.add_subplot(1,1,1)
h1 = plt.plot(T, R030_H1, "-", label="0.75m",linewidth = lw,color = Color_R2)  
h2 = plt.plot(T, R030_H2, "-", label="1.75m",linewidth = lw,color = Color_R3) 
h3 = plt.plot(T, R030_H3, "-", label="2.75m",linewidth = lw,color = Color_R4) 
h4 = plt.plot(T, R030_H4, "-", label="3.75m",linewidth = lw,color = Color_R5)  
h5 = plt.plot(T, R030_H5, "-", label="4.75m",linewidth = lw,color = Color_R6)  
ax.set_ylim([-0.01,0.42])
ax.set_xlabel('Time [s]',fontsize=18)
ax.set_ylabel(r"$\Delta T/\Delta T_{max}$",fontsize=18)
plt.title("height, time serie, R140",fontsize=18)
legend = plt.legend(loc='upper left', frameon=False,fontsize = 18,ncol = 1)
plt.tight_layout()
plt.savefig("RS_R140.pdf")


fig2 = plt.figure(figsize=(6.4, 4.8))
ax = fig2.add_subplot(1,1,1)
h1 = plt.plot(T, R060_H1, "-", label="0.75m",linewidth = lw,color = Color_R2)  
h2 = plt.plot(T, R060_H2, "-", label="1.75m",linewidth = lw,color = Color_R3) 
h3 = plt.plot(T, R060_H3, "-", label="2.75m",linewidth = lw,color = Color_R4) 
h4 = plt.plot(T, R060_H4, "-", label="3.75m",linewidth = lw,color = Color_R5)  
h5 = plt.plot(T, R060_H5, "-", label="4.75m",linewidth = lw,color = Color_R6)  
ax.set_ylim([-0.01,0.42])
ax.set_xlabel('Time [s]',fontsize=18)
ax.set_ylabel(r"$\Delta T/\Delta T_{max}$",fontsize=18)
plt.title("height, time serie, R160",fontsize=18)
legend = plt.legend(loc='upper left', frameon=False,fontsize = 18,ncol = 1)
plt.tight_layout()
plt.savefig("RS_R160.pdf")


fig2 = plt.figure(figsize=(6.4, 4.8))
ax = fig2.add_subplot(1,1,1)
h1 = plt.plot(T, R090_H1, "-", label="0.75m",linewidth = lw,color = Color_R2)  
h2 = plt.plot(T, R090_H2, "-", label="1.75m",linewidth = lw,color = Color_R3) 
h3 = plt.plot(T, R090_H3, "-", label="2.75m",linewidth = lw,color = Color_R4) 
h4 = plt.plot(T, R090_H4, "-", label="3.75m",linewidth = lw,color = Color_R5)  
h5 = plt.plot(T, R090_H5, "-", label="4.75m",linewidth = lw,color = Color_R6)  
ax.set_ylim([-0.01,0.42])
ax.set_xlabel('Time [s]',fontsize=18)
ax.set_ylabel(r"$\Delta T/\Delta T_{max}$",fontsize=18)
plt.title("height, time serie, R180",fontsize=18)
legend = plt.legend(loc='upper left', frameon=False,fontsize = 18,ncol = 1)
plt.tight_layout()
plt.savefig("RS_R180.pdf")


fig2 = plt.figure(figsize=(6.4, 4.8))
ax = fig2.add_subplot(1,1,1)
h1 = plt.plot(T, R120_H1, "-", label="0.75m",linewidth = lw,color = Color_R2)  
h2 = plt.plot(T, R120_H2, "-", label="1.75m",linewidth = lw,color = Color_R3) 
h3 = plt.plot(T, R120_H3, "-", label="2.75m",linewidth = lw,color = Color_R4) 
h4 = plt.plot(T, R120_H4, "-", label="3.75m",linewidth = lw,color = Color_R5)  
h5 = plt.plot(T, R120_H5, "-", label="4.75m",linewidth = lw,color = Color_R6)  
ax.set_ylim([-0.01,0.42])
ax.set_xlabel('Time [s]',fontsize=18)
ax.set_ylabel(r"$\Delta T/\Delta T_{max}$",fontsize=18)
plt.title("height, time serie, R240",fontsize=18)
legend = plt.legend(loc='upper left', frameon=False,fontsize = 18,ncol = 1)
plt.tight_layout()
plt.savefig("RS_R240.pdf")


fig2 = plt.figure(figsize=(6.4, 4.8))
ax = fig2.add_subplot(1,1,1)
h1 = plt.plot(T, R150_H1, "-", label="0.75m",linewidth = lw,color = Color_R2)  
h2 = plt.plot(T, R150_H2, "-", label="1.75m",linewidth = lw,color = Color_R3) 
h3 = plt.plot(T, R150_H3, "-", label="2.75m",linewidth = lw,color = Color_R4) 
h4 = plt.plot(T, R150_H4, "-", label="3.75m",linewidth = lw,color = Color_R5)  
h5 = plt.plot(T, R150_H5, "-", label="4.75m",linewidth = lw,color = Color_R6)  
ax.set_ylim([-0.01,0.42])
ax.set_xlabel('Time [s]',fontsize=18)
ax.set_ylabel(r"$\Delta T/\Delta T_{max}$",fontsize=18)
plt.title("height, time serie, R288",fontsize=18)
legend = plt.legend(loc='upper left', frameon=False,fontsize = 18,ncol = 1)
plt.tight_layout()
plt.savefig("RS_R288.pdf")


fig2 = plt.figure(figsize=(6.4, 4.8))
ax = fig2.add_subplot(1,1,1)
h1 = plt.plot(T, R180_H1, "-", label="0.75m",linewidth = lw,color = Color_R2)  
h2 = plt.plot(T, R180_H2, "-", label="1.75m",linewidth = lw,color = Color_R3) 
h3 = plt.plot(T, R180_H3, "-", label="2.75m",linewidth = lw,color = Color_R4) 
h4 = plt.plot(T, R180_H4, "-", label="3.75m",linewidth = lw,color = Color_R5)  
h5 = plt.plot(T, R180_H5, "-", label="4.75m",linewidth = lw,color = Color_R6)  
ax.set_ylim([-0.01,0.42])
ax.set_xlabel('Time [s]',fontsize=18)
ax.set_ylabel(r"$\Delta T/\Delta T_{max}$",fontsize=18)
plt.title("height, time serie, R360",fontsize=18)
legend = plt.legend(loc='upper left', frameon=False,fontsize = 18,ncol = 1)
plt.tight_layout()
plt.savefig("RS_R360.pdf")












