from __future__ import print_function

import pyPamtra
import shutil
import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import imp
from netCDF4 import Dataset

#_________LMDZ data___________________
path = "/home/grzegorc/AWACA/LMDZ/OUT_golden_case_v2"
nc_file = path + "/TEST-amip-ERA5-LAM.01_20250212_20250218_INS_histinsD17.nc"
nc_data = Dataset(nc_file, "r")

lon = nc_data.variables['lon'][:]
lat = nc_data.variables['lat'][:]
Alt = nc_data.variables['zfull'][0,:,0,0] #if 'Alt' in nc_data.variables else None
T = nc_data.variables['temp'][:]       # shape: (time, level, lat, lon) — check yours
RH = nc_data.variables['rhl'][:]
p = nc_data.variables['pres'][:]

print('max rh',np.shape(T))

p[p<1.]=1. #minimum value accepted by PAMTRA
z = nc_data.variables['zfull'][:]
time = nc_data.variables['time_counter'][:]

#_________cosp data subcolumns________
path="/home/grzegorc/AWACA/COSP/COSPv2.0_lmdz_hillman/driver/run"
nc_file = path+"/hydro_output_golden_case_v2.nc"
nc_data = Dataset(nc_file, "r")

Qi=nc_data['I_LSCICE'][:][::-1,:,:]
Ql=nc_data['I_LSCLIQ'][:][::-1,:,:]
Qr=nc_data['I_LSRAIN'][:][::-1,:,:]
Qs=nc_data['I_LSSNOW'][:][::-1,:,:]

ncol=np.shape(Qs)[1]

Qi=np.transpose(Qi, (2, 1, 0))
Qs=np.transpose(Qs, (2, 1, 0))
Qr=np.transpose(Qr, (2, 1, 0))
Ql=np.transpose(Ql, (2, 1, 0))

#_________format the shape of data_____

T=np.repeat(T[:, np.newaxis, :,0,0], ncol, axis=1)
RH=np.repeat(RH[:, np.newaxis, :,0,0], ncol, axis=1)
p=np.repeat(p[:, np.newaxis, :,0,0], ncol, axis=1)
z=np.repeat(z[:, np.newaxis, :,0,0], ncol, axis=1)

#air density
Rd=287.
Rho_air=p/(Rd*T)

lon=np.full(np.shape(T)[:-1],lon)
lat=np.full(np.shape(T)[:-1],lat)

q_hydro=np.zeros(np.shape(T))
q_hydro=np.repeat(q_hydro[:,:,:,np.newaxis],5, axis=3)

#separate ice of LS scheme and cirrus parameterization
Qi_ice_only=np.zeros(np.shape(Qi))
Qi_cirrus=np.zeros(np.shape(Qi))

#to consider cirrus cloud
#if cirrus==True:
Qi_ice_only[np.where(T>235.15)]=Qi[np.where(T>235.15)]
Qi_cirrus[np.where(T<=235.15)]=Qi[np.where(T<=235.15)]
#else:

#______Definition hydrometeors_____
id_liq=0
id_rain=1
id_snow=2
id_cir=3
id_ice=4 #id larger than 4 are of ice particles

#binned distribution of ice crystals for id_ice*>3
N_bin=12

for i in range(0,N_bin):
    globals()["id_ice_"+str(i)]=id_ice+i
    print('id_ice', globals()["id_ice_"+str(i)])

## Define array of q_hydro
q_hydro=np.zeros(np.shape(T))
q_hydro=np.repeat(q_hydro[:,:,:,np.newaxis],4+N_bin, axis=3)

q_hydro[:,:,:,id_liq]=Ql
q_hydro[:,:,:,id_rain]=Qr
q_hydro[:,:,:,id_snow]=Qs
q_hydro[:,:,:,id_cir]=Qi_cirrus
# ice is defined later

#________load PAMTRA______________
imp.reload(pyPamtra)

pam = pyPamtra.pyPamtra()

#_______Create profile_______

pamData = dict()
jsel=420
isel=340

pamData["lon"] = lon[isel:jsel,:]
pamData["lat"] = lat[isel:jsel,:]
pamData["temp"] = T[isel:jsel,:,:]
pamData["relhum"] = RH[isel:jsel,:,:]
pamData["hgt"] = z[isel:jsel,:,:]
pamData["press"] = p[isel:jsel,:,:]
pamData["hydro_q"] = q_hydro[isel:jsel,:,:]

#________hydrometeor input________
##___Liq_properties

r_liq=2e-5
Rho_liq=1000.

N_liq=(q_hydro[isel:jsel,:,:,id_liq]*Rho_air[isel:jsel,:,:])/(Rho_liq*4/3*np.pi*r_liq**3)
pam.df.addHydrometeor(("liq", -99., 1, Rho_liq, -99., -99., -99., -99. ,3,   1, "mono", -99., -99., -99., -99.,2*r_liq, -99.,"mie-sphere", "khvorostyanov01_drops", -99.))

##___Ice_properties

r_ice=5e-5
Rho_ice=917.
AR_ice=1.
# Variables for the ice exponential distribution
RTT = 273.16
eps=0. #TKE dissip (don't know why it is here)
tempvig1 = -21.06 + RTT
tempvig2 = -30.35 + RTT
naero5=0.5
RPI=np.pi
gamma_snwretro=0. #snow retroaction: don't know why it exist

# size doubling bins definition
#OLD BIN DEFINITION Bin_center=5*2**(np.arange(0,N_bin+1,1))*1e-6 #center of each bin

#New bins definition
end = 2e-3     

bin_0 = end / (2**(N_bin-1))

Bin_center = bin_0 * 2**np.arange(N_bin+1)

Bin_edges=[0.,]
for i in range(0,N_bin):
    Bin_edges.append(np.exp((np.log(Bin_center[i]) + np.log(Bin_center[i+1])) / 2)) #bin boundaries
Bin_edges[0] = Bin_edges[1]/2. #bottom boundary of the first bin

Bin_center=Bin_center[:-1]
Delta_bin=np.array(Bin_edges[1:])-np.array(Bin_edges[:-1])

# DeMott
if (naero5<0):
    if(T>tempvig1):
        nb_crystals = 1e3 * 10**(-0.14*(T-tempvig1) - 2.88)
    elif(T>tempvig2):
        nb_crystals = 1e3 * 10**(-0.31*(T-tempvig1) - 2.88)
    else:
        nb_crystals = 1e3 * 10**(0.)
else:
    nb_crystals = 1e3 * 5.94e-5 * ( RTT - T )**3.33 * naero5**(0.0264*(RTT-T)+0.0033)

qiceini_incl  = Qi_ice_only + gamma_snwretro * Qs+ 1e-30 #1e-30 avoid 0 values
nb_crystals[qiceini_incl==1e-30]=0.
#qiceini_incl[qiceini_incl<eps]=eps

lambda_PSD  = ((RPI*Rho_ice*nb_crystals) / (Rho_air * qiceini_incl )) ** (1./3.)
lambda_PSD[np.isnan(lambda_PSD)] = 0
N0_PSD  = nb_crystals * lambda_PSD

#number of ice crystal in each bin
nb_crystals_bin = np.zeros(np.shape(lambda_PSD))
nb_crystals_bin = np.expand_dims(nb_crystals_bin, axis=-1)  # shape becomes (72,100,95,1)
nb_crystals_bin = np.tile(nb_crystals_bin, (1,1,1,N_bin))   # shape becomes (72,100,95,N_bin)

for i in range(0,N_bin):
    nb_crystals_bin[:,:,:,i]=N0_PSD*np.exp(-lambda_PSD*Bin_center[i])*Delta_bin[i]
#   nb_crystals_bin[:,:,:,i]=nb_crystals_bin[:,:,:,i]*nb_crystals/np.sum(nb_crystals_bin,3) #concentration of ice crystals per bin


#add each bin as one hydrometeor
for i in range(0,N_bin):
    q_hydro[:,:,:,globals()["id_ice_"+str(i)]]=nb_crystals_bin[:,:,:,i]/Rho_air*Rho_ice*4/3*np.pi*(Bin_center[i]/2)**3
    pam.df.addHydrometeor(("ice_bin_"+str(i), -99., -1 , Rho_ice,  -99., -99. ,-99., -99. , 3 ,1, "mono", -99., -99., -99., -99., Bin_center[i], -99., "mie-sphere", "heymsfield10_particles",0.0))


##___Rain_properties___

r_rain=0.0005
rain_fallspeed=4.
Rho_rain=1000.

pam.df.addHydrometeor(("rain",-99.,  1 , Rho_rain , -99., -99.,-99. , -99., 3 ,1,"mono",-99.0, -99.0, -99.0, -99.0,2*r_rain,-99.0,"mie-sphere",rain_fallspeed,0.0))

##___Snow_properties___

r_snow=0.001
AR_snow=1.
snow_fallspeed=1.
Rho_snow = 1.e3 * 0.178 * ( r_snow * 2 * 1000. )**(-0.922)
N_snow=(q_hydro[isel:jsel,:,:,id_snow]*Rho_air[isel:jsel,:,:])/(Rho_snow*4/3*np.pi*r_snow**3)

pam.df.addHydrometeor(("snow",-99., -1 , Rho_snow, -99., -99.,-99. , -99., 3 ,1,"mono",-99.0, -99.0, -99.0, -99.0,2*r_snow,-99.0,"mie-sphere",snow_fallspeed,0.0))

##___Cirrus_properties___

C_cirrus=0.5
Rho_cir=920.
AR_cir=0.5
N_cir=np.zeros(np.shape(q_hydro[isel:jsel, :, :, id_cir]))
N_cir[q_hydro[isel:jsel, :, :, id_cir]>0]= 0.3e6 #0.3 particle per cm-3 to m-3
r_cir=((q_hydro[isel:jsel,:,:,id_cir]*Rho_air[isel:jsel,:,:])/(N_cir*Rho_cir*4/3*np.pi+1e-30))**(1/3)

pam.df.addHydrometeor(("cirrus", C_cirrus, -1 , Rho_ice,  130., 3.0 ,0.684, 2.  , 3 ,1, "mono_cosmo_ice", -99., -99., -99., -99., 2*r_ice, -99., "mie-sphere", "heymsfield10_particles",0.0))

pam.createProfile(**pamData)


## ___ PAMTRA FULL SPECTRA ___
FULL_SPECTRA=True
if FULL_SPECTRA==True:
    pam.nmlSet["hydro_fullspec"] = True
    pam.df.addFullSpectra()
    #print("fullspec",pam.df.dataFullSpec.keys())
    #fullspec dict_keys(['rho_ds', 'd_ds', 'd_bound_ds', 'n_ds', 'mass_ds', 'area_ds', 'as_ratio', 'canting', 'fallvelocity', 'rg_kappa_ds', 'rg_beta_ds', 'rg_gamma_ds', 'rg_zeta_ds'])

    ##___LIQ:
    pam.df.dataFullSpec["rho_ds"][:,:,:,id_liq,:]=1000.
    pam.df.dataFullSpec["d_ds"][:,:,:,id_liq,:]=2*r_liq
    pam.df.dataFullSpec["d_bound_ds"][:,:,:,id_liq,:]=[r_liq,4*r_liq]
    pam.df.dataFullSpec["n_ds"][:,:,:,id_liq,0]=N_liq
    pam.df.dataFullSpec["mass_ds"][:,:,:,id_liq,0]=q_hydro[isel:jsel,:,:,id_liq]
    pam.df.dataFullSpec["area_ds"][:,:,:,id_liq,:]=np.pi*r_liq**2
    pam.df.dataFullSpec["as_ratio"][:,:,:,id_liq,:]=1.
    pam.df.dataFullSpec["canting"][:,:,:,id_liq,:]=0.
    pam.df.dataFullSpec["fallvelocity"][:,:,:,id_liq,:]=0.#liq_fallspeed
    pam.df.dataFullSpec["rg_beta_ds"][:,:,:,id_liq,:]=-99.
    pam.df.dataFullSpec["rg_kappa_ds"][:,:,:,id_liq,:]=-99.
    pam.df.dataFullSpec["rg_gamma_ds"][:,:,:,id_liq,:]=-99.
    pam.df.dataFullSpec["rg_zeta_ds"][:,:,:,id_liq,:]=-99.

    ##___RAIN:
    pam.df.dataFullSpec["rho_ds"][:,:,:,id_rain,:]=1000.
    pam.df.dataFullSpec["d_ds"][:,:,:,id_rain,:]=2*r_rain
    pam.df.dataFullSpec["d_bound_ds"][:,:,:,id_rain,:]=[r_rain,4*r_rain]
    pam.df.dataFullSpec["n_ds"][:,:,:,id_rain,0]=q_hydro[isel:jsel,:,:,id_rain]/(1000*4/3*np.pi*r_rain**3)
    pam.df.dataFullSpec["mass_ds"][:,:,:,id_rain,0]=q_hydro[isel:jsel,:,:,id_rain]
    pam.df.dataFullSpec["area_ds"][:,:,:,id_rain,:]=np.pi*r_liq**2
    pam.df.dataFullSpec["as_ratio"][:,:,:,id_rain,:]=1.
    pam.df.dataFullSpec["canting"][:,:,:,id_rain,:]=0.
    pam.df.dataFullSpec["fallvelocity"][:,:,:,id_rain,:]=rain_fallspeed
    pam.df.dataFullSpec["rg_beta_ds"][:,:,:,id_rain,:]=-99.
    pam.df.dataFullSpec["rg_kappa_ds"][:,:,:,id_rain,:]=-99.
    pam.df.dataFullSpec["rg_gamma_ds"][:,:,:,id_rain,:]=-99.
    pam.df.dataFullSpec["rg_zeta_ds"][:,:,:,id_rain,:]=-99.

    ##___SNOW:
    pam.df.dataFullSpec["rho_ds"][:,:,:,id_snow,:]=Rho_snow
    pam.df.dataFullSpec["d_ds"][:,:,:,id_snow,:]=2*r_snow
    pam.df.dataFullSpec["d_bound_ds"][:,:,:,id_snow,:]=[r_snow,4*r_snow]
    pam.df.dataFullSpec["n_ds"][:,:,:,id_snow,0]=N_snow
    pam.df.dataFullSpec["mass_ds"][:,:,:,id_snow,0]=q_hydro[isel:jsel,:,:,id_snow]
    pam.df.dataFullSpec["area_ds"][:,:,:,id_snow,:]=-99 #np.pi*r_snow**2
    pam.df.dataFullSpec["as_ratio"][:,:,:,id_snow,:]=AR_snow
    pam.df.dataFullSpec["canting"][:,:,:,id_snow,:]=0.
    pam.df.dataFullSpec["fallvelocity"][:,:,:,id_snow,:]=snow_fallspeed
    pam.df.dataFullSpec["rg_beta_ds"][:,:,:,id_snow,:]=-99.
    pam.df.dataFullSpec["rg_kappa_ds"][:,:,:,id_snow,:]=-99.
    pam.df.dataFullSpec["rg_gamma_ds"][:,:,:,id_snow,:]=-99.
    pam.df.dataFullSpec["rg_zeta_ds"][:,:,:,id_snow,:]=-99.

    ##___CIRRUS
    pam.df.dataFullSpec["rho_ds"][:,:,:,id_cir,:]=Rho_cir
    pam.df.dataFullSpec["d_ds"][:,:,:,id_cir,0]=2*r_cir
    pam.df.dataFullSpec["d_bound_ds"][:,:,:,id_cir,:]=[0,4]#4*np.max(r_cir)]
    pam.df.dataFullSpec["n_ds"][:,:,:,id_cir,0]=N_cir
    pam.df.dataFullSpec["mass_ds"][:,:,:,id_cir,0]=q_hydro[isel:jsel,:,:,id_cir]
    pam.df.dataFullSpec["area_ds"][:,:,:,id_cir,:]=-99 #np.pi*r_cir**2
    pam.df.dataFullSpec["as_ratio"][:,:,:,id_cir,:]=AR_cir
    pam.df.dataFullSpec["canting"][:,:,:,id_cir,:]=0.
    pam.df.dataFullSpec["fallvelocity"][:,:,:,id_cir,:]=0.#snow_fallspeed
    pam.df.dataFullSpec["rg_beta_ds"][:,:,:,id_cir,:]=-99.
    pam.df.dataFullSpec["rg_kappa_ds"][:,:,:,id_cir,:]=-99.
    pam.df.dataFullSpec["rg_gamma_ds"][:,:,:,id_cir,:]=-99.
    pam.df.dataFullSpec["rg_zeta_ds"][:,:,:,id_cir,:]=-99.

    ##___ICE
    for i in range(0,N_bin):
        print('fullspec ice',i,globals()["id_ice_"+str(i)])
        pam.df.dataFullSpec["rho_ds"][:,:,:,globals()["id_ice_"+str(i)],:]=Rho_ice
        pam.df.dataFullSpec["d_ds"][:,:,:,globals()["id_ice_"+str(i)],0]=Bin_center[i]
        pam.df.dataFullSpec["d_bound_ds"][:,:,:,globals()["id_ice_"+str(i)],:]=[Bin_edges[i],Bin_edges[i+1]]#4*np.max(r_cir)]
        pam.df.dataFullSpec["n_ds"][:,:,:,globals()["id_ice_"+str(i)],0]=nb_crystals_bin[isel:jsel,:,:,i]
        pam.df.dataFullSpec["mass_ds"][:,:,:,globals()["id_ice_"+str(i)],0]=q_hydro[isel:jsel,:,:,globals()["id_ice_"+str(i)]]
        pam.df.dataFullSpec["area_ds"][:,:,:,globals()["id_ice_"+str(i)],:]=-99 #np.pi*r_cir**2
        pam.df.dataFullSpec["as_ratio"][:,:,:,globals()["id_ice_"+str(i)],:]=AR_ice
        pam.df.dataFullSpec["canting"][:,:,:,globals()["id_ice_"+str(i)],:]=0.
        pam.df.dataFullSpec["fallvelocity"][:,:,:,globals()["id_ice_"+str(i)],:]=0.#snow_fallspeed
        pam.df.dataFullSpec["rg_beta_ds"][:,:,:,globals()["id_ice_"+str(i)],:]=-99.
        pam.df.dataFullSpec["rg_kappa_ds"][:,:,:,globals()["id_ice_"+str(i)],:]=-99.
        pam.df.dataFullSpec["rg_gamma_ds"][:,:,:,globals()["id_ice_"+str(i)],:]=-99.
        pam.df.dataFullSpec["rg_zeta_ds"][:,:,:,globals()["id_ice_"+str(i)],:]=-99.

#__________namelist___________

pam.nmlSet['tmatrix_db'] = 'file'
pam.nmlSet['tmatrix_db_path'] = 'example_db/'
pam.nmlSet["passive"] = False

#__________scattering method____________
pam.df.data["scat_name"][:] = "tmatrix"
#pam.df.data["as_ratio"][:] = 1.0

pam.set["pyVerbose"] = 1
Run_pamtra=False
Run_pamtra=True#False
if Run_pamtra==True:
    print("Start to run")
    pam.runParallelPamtra(35.0,
                      pp_deltaX=1,    # 2 profiles in X per worker
                      pp_deltaY=1,    # 1 profile in Y per worker
                      pp_deltaF=1,    # 1 frequency per worker
                      pp_local_workers="auto")  # detect CPU cores
    print("SHAPE",np.shape(pam.r["Ze"]))

    plt.figure('hydro')
    plt.imshow(np.sum(Qi,1))
    plt.colorbar()
    print("Run end")
    plt.show()


    plt.figure('test')
    plt.imshow(pam.r["Ze"][:,0,:,0,0,0],aspect='auto',vmin=-60,vmax=20)
    plt.colorbar()
    print("Run end")
    plt.show()

# Output NetCDF file path
#output_file = "/home/grzegorc/AWACA/PAMTRA/pamtra/Ouput_cosp_final_no_fullspec.nc"
output_file = "/home/grzegorc/AWACA/PAMTRA/pamtra/Ouput_golden_case_D17_v2_new_bin_wateronly.nc"
output_file = "/home/grzegorc/AWACA/PAMTRA/pamtra/test.nc"#Ouput_golden_case_D17_v2_new_bin_wateronly.nc"
with Dataset(output_file, "w", format="NETCDF4") as nc_out:

    # Dimensions
    nc_out.createDimension("time", len(time[isel:jsel]))
    nc_out.createDimension("col", ncol)
    nc_out.createDimension("level", T.shape[2])
    nc_out.createDimension("hydro", 4+N_bin)

    # Variables simples
    nc_out.createVariable("time", "f8", ("time",))[:] = time[isel:jsel]
    nc_out.createVariable("col", "i4", ("col",))[:] = np.arange(ncol)
    nc_out.createVariable("level", "i4", ("level",))[:] = Alt
    nc_out.createVariable("hydro", "i4", ("hydro",))[:] = np.arange(4+N_bin)

    # Écriture des champs
    nc_out.createVariable("lon", "f4", ("time", "col"))[:] = pamData["lon"]
    nc_out.createVariable("lat", "f4", ("time", "col"))[:] = pamData["lat"]
    nc_out.createVariable("temp", "f4", ("time", "col", "level"))[:] = pamData["temp"]
    nc_out.createVariable("relhum", "f4", ("time", "col", "level"))[:] = pamData["relhum"]
    nc_out.createVariable("hgt", "f4", ("time", "col", "level"))[:] = pamData["hgt"]
    nc_out.createVariable("press", "f4", ("time", "col", "level"))[:] = pamData["press"]
    nc_out.createVariable("hydro_q", "f4", ("time", "col", "level", "hydro"))[:] = pamData["hydro_q"]
    nc_out.createVariable("Ze", "f4", ("time", "col", "level"))[:] = pam.r["Ze"][:,:,:,0,0,0]
    print("PAMTRA output saved as NetCDF file  ",output_file)
