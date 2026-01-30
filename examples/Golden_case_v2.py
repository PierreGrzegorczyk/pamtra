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
 
#______Definition hydrometeors_____
id_liq=0
id_rain=1
id_snow=2
id_ice=3

## Define array of q_hydro
q_hydro=np.zeros(np.shape(T))
q_hydro=np.repeat(q_hydro[:,:,:,np.newaxis],4, axis=3)

q_hydro[:,:,:,id_liq]=0.#Ql
q_hydro[:,:,:,id_rain]=0.#Qr
q_hydro[:,:,:,id_snow]=Qs
q_hydro[:,:,:,id_ice]=0.#Qi
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

r_ice=50e-6
Rho_ice=917.
AR_ice=1.
N_ice=(q_hydro[isel:jsel,:,:,id_ice]*Rho_air[isel:jsel,:,:])/(Rho_ice*4/3*np.pi*r_ice**3)
pam.df.addHydrometeor(("ice", -99., -1, Rho_ice,  -99., -99. ,-99., -99. , 3 ,1, "mono", -99., -99., -99., -99., 2*r_ice, -99., "mie-sphere", "heymsfield10_particles",0.0))

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
    ##___ICE
    pam.df.dataFullSpec["rho_ds"][:,:,:,id_ice,:]=Rho_ice
    pam.df.dataFullSpec["d_ds"][:,:,:,id_ice,:]=2*r_ice
    pam.df.dataFullSpec["d_bound_ds"][:,:,:,id_ice,:]=[r_ice,4*r_ice]
    pam.df.dataFullSpec["n_ds"][:,:,:,id_ice,0]=N_ice
    pam.df.dataFullSpec["mass_ds"][:,:,:,id_ice,0]=q_hydro[isel:jsel,:,:,id_ice]
    pam.df.dataFullSpec["area_ds"][:,:,:,id_ice,:]=-99 #np.pi*r_cir**2
    pam.df.dataFullSpec["as_ratio"][:,:,:,id_ice,:]=AR_ice
    pam.df.dataFullSpec["canting"][:,:,:,id_ice,:]=0.
    pam.df.dataFullSpec["fallvelocity"][:,:,:,id_ice,:]=0.#snow_fallspeed
    pam.df.dataFullSpec["rg_beta_ds"][:,:,:,id_ice,:]=-99.
    pam.df.dataFullSpec["rg_kappa_ds"][:,:,:,id_ice,:]=-99.
    pam.df.dataFullSpec["rg_gamma_ds"][:,:,:,id_ice,:]=-99.
    pam.df.dataFullSpec["rg_zeta_ds"][:,:,:,id_ice,:]=-99.


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
    plt.imshow(np.sum(Qi,1).T,aspect='auto')
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
output_file = "/home/grzegorc/AWACA/PAMTRA/pamtra/Ouput_golden_case_D17_v2_new_bin_monoice.nc"
with Dataset(output_file, "w", format="NETCDF4") as nc_out:

    # Dimensions
    nc_out.createDimension("time", len(time[isel:jsel]))
    nc_out.createDimension("col", ncol)
    nc_out.createDimension("level", T.shape[2])
    nc_out.createDimension("hydro", 4)

    # Variables simples
    nc_out.createVariable("time", "f8", ("time",))[:] = time[isel:jsel]
    nc_out.createVariable("col", "i4", ("col",))[:] = np.arange(ncol)
    nc_out.createVariable("level", "i4", ("level",))[:] = Alt
    nc_out.createVariable("hydro", "i4", ("hydro",))[:] = np.arange(4)

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
