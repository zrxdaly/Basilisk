#%%
from os import sep
import matplotlib.pyplot as plt
import numpy as np
import glob 
import xarray as xr

from matplotlib.ticker import NullFormatter

from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# plt.rc('font', family='serif')
plt.rc('font',**{'family':'serif','serif':['Times']})
plt.rc('xtick', labelsize='15')
plt.rc('ytick', labelsize='15')
# plt.rc('text', usetex=True)
lw = 2
Color_R2 = '#000000'
Color_R3 = '#fcf18f'
Color_R4 = '#377eb8'
Color_R5 = '#008000'
Color_R6 = '#084594'
Color_R7 = '#ff7f00'
Color_R8 = '#808080'
#%%
dir1 = "/home/leaf/software/basilisk/basilisk/Exp_field_Basilisk/ventilator/resultslice/"
dir_buo = dir1 + "buo/"
dir_vel_u = dir1 + "vel_u/"
dir_vel_v = dir1 + "vel_v/"
dir_vel_w = dir1 + "vel_w/"

G = 9.81
T_ref = 273
KCR = 273.15
INV = 0.3
T_hub = G/T_ref*INV*10.5

def K2C(T_hub):
    return(T_hub*T_ref/G + T_ref - KCR)

TC_hub = K2C(T_hub)

Dt = 1    # s  the time interval
T0 = 1200     # s  time length for no operation
T1 = 1020    # s  time length for fan operation


#%% get the reference state of 3D data at each height
def D2_N(filedir, Tstart, Tend):
    # 1234 represent different case
    FR_05 = sorted(glob.glob(filedir + 't=*y=350'))

    # FR_Bu05N = np.loadtxt(FR_05[Tstart], dtype='f',skiprows=2)
    FR_Bu05N = np.loadtxt(FR_05[Tstart], dtype='f')

    for i in np.arange(Tstart+1, Tend):
        # FR_Bu05N = FR_Bu05N + np.loadtxt(FR_05[i], dtype='f',skiprows=2)
        FR_Bu05N = FR_Bu05N + np.loadtxt(FR_05[i], dtype='f')
        print(i)
    # FR_Bu05N_TA = K2C(FR_Bu05N/(Tend - Tstart))
    FR_Bu05N_TA = FR_Bu05N/(Tend - Tstart)
    
    return(FR_Bu05N_TA)

#%%
# reference state of the buoyancy data (the temp differences over height)
buo_R2 = D2_N(dir_vel_u, 100, 900)
#%% 

Fsim_buo = np.transpose(Fsim_buo, (1, 0))

# %% get the measurement data (temp)
in_folder = '/home/leaf/Documents/Data_Krabbendijke_2021/Exp_Data_cali/05-07-2020'

# %%
C1_ff_files = in_folder + "/calibrated_WEST"
C3_ff_files = in_folder + "/calibrated_EAST"
nc1_file = sorted(glob.glob(C1_ff_files+'/*.nc'))
nc3_file = sorted(glob.glob(C3_ff_files+'/*.nc'))

# %%  # index of 2000 is around H23T10
West_file = xr.open_mfdataset(nc1_file, concat_dim="time",
            data_vars=["tmpw"], coords='all', compat='override')

East_file = xr.open_mfdataset(nc3_file, concat_dim="time",
            data_vars=["tmpw"], coords='all', compat='override')

Int_sec = slice("2021-05-07T20:00:00","2021-05-08T01:00")

West_file_int = West_file.sel(time = Int_sec)
East_file_int = East_file.sel(time = Int_sec)

W_tmpad_15S = West_file_int.tmpw.resample(time="10.77S").mean()
E_tmpad_15S = East_file_int.tmpw.resample(time="10.77S").mean()

# %% get the cable section
X_W06 = np.array([4228, 3963])
X_W12 = np.array([3625, 3358])
X_W18 = np.array([3018, 2757])
X_W24 = np.array([2411, 2144])
X_W31 = np.array([1788, 1526])
X_W39 = np.array([1164, 888])

X_E00 = np.array([4243, 3967])
X_E06 = np.array([3623, 3350])
X_E12 = np.array([2933, 2661])
X_E18 = np.array([2312, 2051])
X_E24 = np.array([1713, 1443])
X_E30 = np.array([1136, 893])


def find_near(array,value):
    idx = (np.abs(array - value)).argmin()
    return(array[idx])

def get_T_list(X_W12, dsW_01, WE):
    # this equation is based on the georeference [x_cable, x_plot] = [39499, 623.9], [730, 3853.7]
    if WE == "W":
        X_W12 = X_W12 * (-1.003355) + 4586.149
    elif WE == "E":
        X_W12 = X_W12 * (-1.003586) + 4567.56325

    idx_W12 = np.zeros_like(X_W12)

    for index,value in enumerate(X_W12):
        idx_W12[index] = find_near(dsW_01["x"],value)

    tmp_W12 = dsW_01.sel(x=slice(idx_W12[0], idx_W12[1]))
    return(tmp_W12)

tmp_W06 = get_T_list(X_W06, W_tmpad_15S, "W")
tmp_W12 = get_T_list(X_W12, W_tmpad_15S, "W")
tmp_W18 = get_T_list(X_W18, W_tmpad_15S, "W")
tmp_W24 = get_T_list(X_W24, W_tmpad_15S, "W")
tmp_W31 = get_T_list(X_W31, W_tmpad_15S, "W")
tmp_W39 = get_T_list(X_W39, W_tmpad_15S, "W")

tmp_E00 = get_T_list(X_E00, E_tmpad_15S, "E")
tmp_E06 = get_T_list(X_E06, E_tmpad_15S, "E")
tmp_E12 = get_T_list(X_E12, E_tmpad_15S, "E")
tmp_E18 = get_T_list(X_E18, E_tmpad_15S, "E")
tmp_E24 = get_T_list(X_E24, E_tmpad_15S, "E")
tmp_E30 = get_T_list(X_E30, E_tmpad_15S, "E")

x_len = 134
def intep(tmp_W06, WE):
    # len_x = np.size(tmp_W06["x"])
    x_inp = np.linspace(tmp_W06["x"][0], tmp_W06["x"][-1], x_len)
    tmp_W06_inp = tmp_W06.interp(x = x_inp)
    tmp_W06_inp["x"] = np.linspace(0,270,x_len)
    if WE == "W":
        return(tmp_W06_inp)
    elif WE == "E":
        return(tmp_W06_inp)

tmp_W06_inp = intep(tmp_W06, "W")
tmp_W12_inp = intep(tmp_W12, "W")
tmp_W18_inp = intep(tmp_W18, "W")
tmp_W24_inp = intep(tmp_W24, "W")
tmp_W31_inp = intep(tmp_W31, "W")
tmp_W39_inp = intep(tmp_W39, "W")

tmp_E00_inp = intep(tmp_E00, "E")
tmp_E06_inp = intep(tmp_E06, "E")
tmp_E12_inp = intep(tmp_E12, "E")
tmp_E18_inp = intep(tmp_E18, "E")
tmp_E24_inp = intep(tmp_E24, "E")
tmp_E30_inp = intep(tmp_E30, "E")

#%%
WE_tmp = xr.concat([tmp_W39_inp, tmp_W31_inp, tmp_W24_inp, tmp_W18_inp, tmp_W12_inp, 
                   tmp_W06_inp, tmp_E00_inp, tmp_E06_inp, tmp_E12_inp, tmp_E18_inp, tmp_E24_inp, tmp_E30_inp], 
                  "y").T
#%%
OP_T_section = slice("2021-05-07T23:40", "2021-05-08T01:00:00")
WU_T_section = slice("2021-05-07T23:00", "2021-05-07T23:40:00")
RE2_T_section = slice("2021-05-07T21:41", "2021-05-07T23:00:00")
RE1_T_section = slice("2021-05-07T20:22", "2021-05-07T21:41:00")

WE_tmp_RE1 = WE_tmp.sel(time = RE1_T_section)
WE_tmp_RE2 = WE_tmp.sel(time = RE2_T_section)
WE_tmp_WU = WE_tmp.sel(time = WU_T_section)
WE_tmp_OP = WE_tmp.sel(time = OP_T_section)

# %%
F_rota = WE_tmp_OP[0:112,:].mean("time")

Y = np.array([-100, -80, -60, -40, -20, 0,20,40,60,80,110,150])[::-1] * (-1)

F_rota = F_rota.assign_coords({"y": Y})

#%%
F_rota_inp = F_rota.interp(y=np.linspace(-150, 100, 125))

n_n = 57
s_n = 77
xx = np.linspace(0, 270, 134)
yy = np.linspace(-150, 100, 125)
# Fsim_buo[0, 175-77:175+57, 175-75:175+50]

#%%
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize = (6, 11))

tp1 = ax1.contourf(yy, xx, F_rota_inp.values, levels = np.arange(0, 8,0.1),cmap = "magma")
tp2 = ax2.contourf(yy, xx, Fsim_buo[175-77:175+57, 175-75:175+50], levels = np.arange(0, 8,0.1),cmap = "magma")
tp3 = ax3.contourf(yy, xx, F_rota_inp.values - Fsim_buo[175-77:175+57, 175-75:175+50], levels = np.arange(-1.5, 3,0.1),cmap = "magma")

for c in tp1.collections:
    c.set_edgecolor("face")
for c in tp2.collections:
    c.set_edgecolor("face")
for c in tp3.collections:
    c.set_edgecolor("face")
cb = fig.colorbar(tp1, ax=(ax1,ax2))
cb2 = fig.colorbar(tp3, ax=(ax3))
cb.set_label(r"temp")
ax1.set_ylabel("exp",fontsize=18)
ax2.set_ylabel("simu",fontsize=18)
ax3.set_ylabel("exp-simu",fontsize=18)








# %%  load the wind data and 
dir_sonic = "/home/leaf/Documents/Data_Krabbendijke_2021/Raw_data/Krabbendijke_experiment_2021_0507/sonic_data/"
file_sonic_COM6 = sorted(glob.glob(dir_sonic + "data_COM6*.dat"))
file_sonic_COM7 = sorted(glob.glob(dir_sonic + "data_COM7*.dat"))
# COM7 seems to be W31
# %%
import pandas as pd

def get_pdF(num):
    C6_H23 = np.genfromtxt(file_sonic_COM6[num], usecols=(1, 2, 3, 4))
    C7_H23 = np.genfromtxt(file_sonic_COM7[num], usecols=(1, 2, 3, 4))
    C6_date23 = np.genfromtxt(file_sonic_COM6[num],delimiter=",", usecols=0, dtype="str")
    C7_date23 = np.genfromtxt(file_sonic_COM7[num],delimiter=",", usecols=0, dtype="str")
    C6_H23 = pd.DataFrame(data = C6_H23, index=C6_date23, columns=["u", "v", "w", "T"])
    C7_H23 = pd.DataFrame(data = C7_H23, index=C7_date23, columns=["u", "v", "w", "T"])
    return(C6_H23, C7_H23)
#%%
C6, C7 = get_pdF(0)

for i in range(len(file_sonic_COM6)-1):
    C6 = pd.concat([C6, get_pdF(i+1)[0]])
    C7 = pd.concat([C7, get_pdF(i+1)[1]])
    print(i)

#%% get all the data and save to netcdf

C6.index = pd.to_datetime(C6.index)
C7.index = pd.to_datetime(C7.index)
#%%
C6["u"] = -C6["u"]
C6["v"] = -C6["v"]
#%%
C77 = C7.copy()
th = 9.85/180 * np.pi
C77["u"] = C7["u"] * np.cos(th) - C7["v"] * np.sin(th)
C77["v"] = C7["u"] * np.sin(th) + C7["v"] * np.cos(th)
C7 = C77
# C6 = C6.to_xarray()
# C7 = C77.to_xarray()

#%%
# C6 = C6.load()
# C6.to_netcdf(
#     path="C6_sonic.nc"
# )
# C6.close()

# C7 = C7.load()
# C7.to_netcdf(
#     path="C7_sonic.nc"
# )
# C7.close()

#%% some visualization of the sonic data
C6.u.resample("2T").mean().plot()
C6.v.resample("2T").mean().plot()
C6.w.resample("2T").mean().plot()
#%%
C7.u.resample("2T").mean().plot()
C7.v.resample("2T").mean().plot()
C7.w.resample("2T").mean().plot()
#%% calcualtion of the angle
C6["L"] = (C6["u"] ** + C6["v"] ** 2)**0.5
mask1 = (C6["u"]<=0) & (C6["v"]<=0)
mask2 = (C6["u"]<=0) & (C6["v"]>=0)
mask3 = (C6["u"]>=0) & (C6["v"]>=0)
mask4 = (C6["u"]>=0) & (C6["v"]<=0)

value1 = np.arcsin(-C6["u"]/C6["L"])/np.pi * 180
value2 = (np.arcsin(C6["v"]/C6["L"]) + np.pi/2)/np.pi * 180
value3 = (np.arcsin(C6["u"]/C6["L"]) + np.pi)/np.pi * 180
value4 = (np.arcsin(-C6["v"]/C6["L"]) + np.pi * 3/2)/np.pi * 180

C6['ang'] = np.where(mask1, value1, 
            np.where(mask2, value2, 
            np.where(mask3, value3,
            np.where(mask4, value4, 0))))
#%%
C7["L"] = (C7["u"] ** + C7["v"] ** 2)**0.5
mask1 = (C7["u"]<=0) & (C7["v"]<=0)
mask2 = (C7["u"]<=0) & (C7["v"]>=0)
mask3 = (C7["u"]>=0) & (C7["v"]>=0)
mask4 = (C7["u"]>=0) & (C7["v"]<=0)

value1 = np.arcsin(-C7["u"]/C7["L"])/np.pi * 180
value2 = (np.arcsin(C7["v"]/C7["L"]) + np.pi/2)/np.pi * 180
value3 = (np.arcsin(C7["u"]/C7["L"]) + np.pi)/np.pi * 180
value4 = (np.arcsin(-C7["v"]/C7["L"]) + np.pi * 3/2)/np.pi * 180

C7['ang'] = np.where(mask1, value1, 
            np.where(mask2, value2, 
            np.where(mask3, value3,
            np.where(mask4, value4, 0))))
#%%
# C6.drop('L', axis=1)
# C7.drop('L', axis=1)

# C6 = C6.to_xarray()
# C7 = C7.to_xarray()

# C6 = C6.load()
# C6.to_netcdf(
#     path="C6_sonic.nc"
# )
# C6.close()

# C7 = C7.load()
# C7.to_netcdf(
#     path="C7_sonic.nc"
# )
# C7.close()
#%%
C6["KE"] = C6["u"]**2 + C6["v"]**2 + C6["w"]**2
C7["KE"] = C7["u"]**2 + C7["v"]**2 + C7["w"]**2
# %% plot the KE plot 
C6.KE.resample("1S").mean().plot()
C7.KE.resample("1S").mean().plot()

# C6.KE.plot()
# C7.KE.plot()

#%% get the simualtion data 
#%% get the reference state of 3D data at each height
def get_vel(filedir, Tstart, Tend):
    # 1234 represent different case
    FR_15 = sorted(glob.glob(filedir + 't=*y=003'))

    W12 = []
    W31 = []

    for i in np.arange(Tstart, Tend):
        W12 = np.append(W12, np.loadtxt(FR_15[i], dtype='f',skiprows=2)[175, 155])
        W31 = np.append(W31, np.loadtxt(FR_15[i], dtype='f',skiprows=2)[175, 125])
        print(i)
    W12 = K2C(W12)
    W31 = K2C(W31)
    
    return(W12, W31)

#%%
# reference state of the buoyancy data (the temp differences over height)
W12_U, W31_U = get_vel(dir_vel_u, 0, T0)
W12_V, W31_V = get_vel(dir_vel_v, 0, T0)
W12_W, W31_W = get_vel(dir_vel_w, 0, T0)
#%%
C6_slice = C6.loc['2021-05-07 23:40:00':'2021-05-07 23:59:59'].resample("1S").mean()
C6_slice["sim_KE"] = W12_U**2 + W12_V**2 + W12_W**2

C7_slice = C7.loc['2021-05-07 23:40:00':'2021-05-07 23:59:59'].resample("1S").mean()
C7_slice["sim_KE"] = W31_U**2 + W31_V**2 + W31_W**2

# %% comparison between sim and exp
fig1 = plt.figure(figsize=(6.4, 4.8))
C6_slice.KE.plot(legend = "exp_KE")
C6_slice.sim_KE.plot(legend = "sim U")
plt.title("W12")
# %%
fig2 = plt.figure(figsize=(6.4, 4.8))
C7_slice.KE.plot(legend = "exp_KE")
C7_slice.sim_KE.plot(legend = "sim U")
plt.title("W31")
# %%
