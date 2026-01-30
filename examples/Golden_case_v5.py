from __future__ import print_function

import pyPamtra
import shutil
import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import imp
from netCDF4 import Dataset

#_________LMDZ data___________________
path = "/home/grzegorc/AWACA/LMDZ/OUT_golden_case_v5"
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
path="/home/grzegorc/AWACA/COSP/COSPv2.0_lmdz_hillman_precip/driver/run"
nc_file = path+"/hydro_output_golden_case_v5.nc"
nc_data = Dataset(nc_file, "r")

Qi=nc_data['I_LSCICE'][:][::-1,:,:]
Ql=nc_data['I_LSCLIQ'][:][::-1,:,:]
Qr=nc_data['I_LSRAIN'][:][::-1,:,:]
Qs=nc_data['I_LSSNOW'][:][::-1,:,:]
print("number of extreme values", np.sum(Qs>0.05))
plt.figure('Qi')
plt.imshow(np.sum(Qi,1),aspect='auto')
plt.colorbar()
print("Show Qi")

plt.figure('Qs')
plt.imshow(np.sum(Qs,1),aspect='auto')
plt.colorbar()
print("Show Qs")
plt.show()




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
q_hydro=np.repeat(q_hydro[:,:,:,np.newaxis],4, axis=3)


#______Definition hydrometeors_____
id_liq=0
id_rain=1
id_snow=2
id_ice=3 



## Define array of q_hydro
q_hydro=np.zeros(np.shape(T))
q_hydro=np.repeat(q_hydro[:,:,:,np.newaxis],4, axis=3)

#q_hydro[:,:,:,id_liq]=Ql
q_hydro[:,:,:,id_rain]=Qr
q_hydro[:,:,:,id_snow]=Qs
q_hydro[:,:,:,id_ice]=Qi
# ice is defined later

#________load PAMTRA______________
imp.reload(pyPamtra)

pam = pyPamtra.pyPamtra()

#_______Create profile_______

pamData = dict()
jsel=420
isel=340

jsel=426
isel=425

#isel=0#3*24*3
#jsel=-1#6*24*3

pamData["lon"] = lon[isel:jsel,:]
pamData["lat"] = lat[isel:jsel,:]
pamData["temp"] = T[isel:jsel,:,:]
pamData["relhum"] = RH[isel:jsel,:,:]
pamData["hgt"] = z[isel:jsel,:,:]
pamData["press"] = p[isel:jsel,:,:]
pamData["hydro_q"] = q_hydro[isel:jsel,:,:]

#________hydrometeor input________
##___Liq_properties
nbin_liq=100
r_liq=1e-5
Rho_liq=1000.

N_liq=(q_hydro[isel:jsel,:,:,id_liq]*Rho_air[isel:jsel,:,:])/(Rho_liq*4/3*np.pi*r_liq**3)
pam.df.addHydrometeor(("liq", -99., 1, Rho_liq, -99., -99., -99., -99. ,3,   nbin_liq, "mono", -99., -99., -99., -99.,2*r_liq, -99.,"mie-sphere", "khvorostyanov01_drops", -99.))

##___Rain_properties___
nbin_rain=100
r_rain=0.0005
rain_fallspeed=4.
Rho_rain=1000.

N_rain=q_hydro[isel:jsel,:,:,id_rain]/(Rho_rain*4/3*np.pi*r_rain**3)
pam.df.addHydrometeor(("rain",-99.,  1 , Rho_rain , -99., -99.,-99. , -99., 3 ,nbin_rain,"mono",-99.0, -99.0, -99.0, -99.0,2*r_rain,-99.0,"mie-sphere",rain_fallspeed,0.0))

##___Snow_properties___
nbin_snow=100
r_snow=0.001
AR_snow=1.
snow_fallspeed=1.
Rho_snow = 1.e3 * 0.178 * ( r_snow * 2 * 1000. )**(-0.922)
N_snow=(q_hydro[isel:jsel,:,:,id_snow]*Rho_air[isel:jsel,:,:])/(Rho_snow*4/3*np.pi*r_snow**3)

pam.df.addHydrometeor(("snow",-99., -1 , Rho_snow, -99., -99.,-99. , -99., 3 ,nbin_snow,"mono",-99.0, -99.0, -99.0, -99.0,2*r_snow,-99.0,"mie-sphere","heymsfield10_particles",0.0))

##___Ice_properties___
nbin_ice=100
C_ice=1.
Rho_ice=920.
AR_ice=1.#
r_ice=50e-6#((q_hydro[isel:jsel,:,:,id_cir]*Rho_air[isel:jsel,:,:])/(N_cir*Rho_cir*4/3*np.pi+1e-30))**(1/3)
N_ice=(q_hydro[isel:jsel,:,:,id_ice]*Rho_air[isel:jsel,:,:])/(Rho_ice*4/3*np.pi*r_ice**3)

pam.df.addHydrometeor(("ice", C_ice, -1 , Rho_ice,  130., 3.0 ,0.684, 2.  , 3 ,nbin_ice, "mono_cosmo_ice", -99., -99., -99., -99., 2*r_ice, -99., "mie-sphere", "heymsfield10_particles",0.0))






d_bound_hydro=np.linspace(10e-6,r_snow*2+10e-6,nbin_snow+1)
d_hydro=np.linspace(2*r_liq,r_snow*2,nbin_snow)



print("d_bound_hydro",d_bound_hydro)
print("d_hydro",d_hydro)
pam.createProfile(**pamData)
print('check',pamData["hydro_q"].shape)
## ___ PAMTRA FULL SPECTRA ___
FULL_SPECTRA=True
if FULL_SPECTRA==True:
    pam.nmlSet["hydro_fullspec"] = True
    pam.df.addFullSpectra()

#print("fullspec",pam.df.dataFullSpec.keys())
    #fullspec dict_keys(['rho_ds', 'd_ds', 'd_bound_ds', 'n_ds', 'mass_ds', 'area_ds', 'as_ratio', 'canting', 'fallvelocity', 'rg_kappa_ds', 'rg_beta_ds', 'rg_gamma_ds', 'rg_zeta_ds'])

    pam.nmlSet["radar_mode"] = "spectrum"
    pam.nmlSet["radar_min_v"] = -15.
    pam.nmlSet["radar_max_v"] =  15.
    pam.nmlSet["radar_nfft"]  = 256
    pam.nmlSet["passive"] = False
    pam.nmlSet["randomseed"] = 10
    pam.nmlSet["radar_mode"] = "spectrum"
    pam.nmlSet["radar_aliasing_nyquist_interv"] = 3
    pam.nmlSet["hydro_adaptive_grid"] = False
    pam.nmlSet["conserve_mass_rescale_dsd"] = False
    pam.nmlSet["radar_use_hildebrand"] = True
    pam.nmlSet["radar_noise_distance_factor"] = -2

    pam.nmlSet["save_psd"] = True
    #pam.nmlSet['tmatrix_db'] = 'file'
    #pam.nmlSet['tmatrix_db_path'] = 'example_db/'
    pam.nmlSet["passive"] = False 


    print("shape of fullspectra", np.shape(pam.df.dataFullSpec["d_bound_ds"][:,:,:,:,:]))
##___LIQ:
    pam.df.dataFullSpec["d_bound_ds"][:,:,:,id_liq,:]=d_bound_hydro
    pam.df.dataFullSpec["d_ds"][:,:,:,id_liq,:]=d_hydro

    ibin=np.where(d_hydro==2*r_liq)[0][0]
    print("IBIN",ibin)
    pam.df.dataFullSpec["rho_ds"][:,:,:,id_liq,:]=1000.
    pam.df.dataFullSpec["n_ds"][:,:,:,id_liq,ibin]=N_liq
    pam.df.dataFullSpec["mass_ds"][:,:,:,id_liq,ibin]=q_hydro[isel:jsel,:,:,id_liq]
    pam.df.dataFullSpec["area_ds"][:,:,:,id_liq,:]=np.pi/4. *  pam.df.dataFullSpec["d_ds"][:,:,:,id_liq,:] ** 2
    pam.df.dataFullSpec["as_ratio"][:,:,:,id_liq,:]=1.
    #pam.df.dataFullSpec["canting"][:,:,:,id_liq,:]=0.
    #pam.df.dataFullSpec["fallvelocity"][:,:,:,id_liq,:]=rain_fallspeed#0.#liq_fallspeed
    #pam.df.dataFullSpec["rg_beta_ds"][:,:,:,id_liq,:]=-99.
    #pam.df.dataFullSpec["rg_kappa_ds"][:,:,:,id_liq,:]=-99.
    #pam.df.dataFullSpec["rg_gamma_ds"][:,:,:,id_liq,:]=-99.
    #pam.df.dataFullSpec["rg_zeta_ds"][:,:,:,id_liq,:]=-99.

    ##___RAIN:
    ibin=np.where(d_hydro==2*r_rain)[0][0]
    print("IBIN1",ibin)

    pam.df.dataFullSpec["d_bound_ds"][:,:,:,id_rain,:]=d_bound_hydro
    pam.df.dataFullSpec["d_ds"][:,:,:,id_rain,:]=d_hydro
    pam.df.dataFullSpec["rho_ds"][:,:,:,id_rain,:]=1000.
    pam.df.dataFullSpec["n_ds"][:,:,:,id_rain,ibin]=N_rain#q_hydro[isel:jsel,:,:,id_rain]/(1000*4/3*np.pi*r_rain**3)
    pam.df.dataFullSpec["mass_ds"][:,:,:,id_rain,ibin]=q_hydro[isel:jsel,:,:,id_rain]
    pam.df.dataFullSpec["area_ds"][:,:,:,id_rain,:]=np.pi/4. *  pam.df.dataFullSpec["d_ds"][:,:,:,id_rain,:] ** 2
    pam.df.dataFullSpec["as_ratio"][:,:,:,id_rain,:]=1.
    #pam.df.dataFullSpec["canting"][:,:,:,id_rain,:]=0.
    #pam.df.dataFullSpec["fallvelocity"][:,:,:,id_rain,:]=rain_fallspeed
    #pam.df.dataFullSpec["rg_beta_ds"][:,:,:,id_rain,:]=-99.
    #pam.df.dataFullSpec["rg_kappa_ds"][:,:,:,id_rain,:]=-99.
    #pam.df.dataFullSpec["rg_gamma_ds"][:,:,:,id_rain,:]=-99.
    #pam.df.dataFullSpec["rg_zeta_ds"][:,:,:,id_rain,:]=-99.

    ##___SNOW:

    ibin=np.where(d_hydro==2*r_snow)[0][0]
    print("IBIN2",ibin)
    pam.df.dataFullSpec["d_bound_ds"][:,:,:,id_snow,:]=d_bound_hydro
    pam.df.dataFullSpec["d_ds"][:,:,:,id_snow,:]=d_hydro
    pam.df.dataFullSpec["rho_ds"][:,:,:,id_snow,:]=Rho_snow
    pam.df.dataFullSpec["n_ds"][:,:,:,id_snow,ibin]=0.#N_snow
    pam.df.dataFullSpec["mass_ds"][:,:,:,id_snow,ibin]=q_hydro[isel:jsel,:,:,id_snow]
    pam.df.dataFullSpec["as_ratio"][:,:,:,id_snow,:]=AR_snow
    pam.df.dataFullSpec["area_ds"][:,:,:,id_snow,:]=np.pi/4.*pam.df.dataFullSpec["d_ds"][:,:,:,id_snow,:] ** 2
    #pam.df.dataFullSpec["canting"][:,:,:,id_snow,:]=0.
    #pam.df.dataFullSpec["fallvelocity"][:,:,:,id_snow,:]=snow_fallspeed
    #pam.df.dataFullSpec["rg_beta_ds"][:,:,:,id_snow,:]=-99.
    #pam.df.dataFullSpec["rg_kappa_ds"][:,:,:,id_snow,:]=-99.
    #pam.df.dataFullSpec["rg_gamma_ds"][:,:,:,id_snow,:]=-99.
    #pam.df.dataFullSpec["rg_zeta_ds"][:,:,:,id_snow,:]=-99.

    ##___ICE
    #print("rice",2*r_ice,d_hydro[4])
    #ibin=np.where(d_hydro==2*r_ice)
    #print("IBIN3",ibin)
    ibin=4#ibin[0][0]

    pam.df.dataFullSpec["d_bound_ds"][:,:,:,id_ice,:]=d_bound_hydro
    pam.df.dataFullSpec["d_ds"][:,:,:,id_ice,:]=d_hydro
    pam.df.dataFullSpec["rho_ds"][:,:,:,id_ice,:]=Rho_ice
    #pam.df.dataFullSpec["n_ds"][:,:,:,id_ice,ibin-3]=0.25*N_ice
    #pam.df.dataFullSpec["n_ds"][:,:,:,id_ice,ibin-2]=0.25*N_ice
    #pam.df.dataFullSpec["n_ds"][:,:,:,id_ice,ibin-1]=0.25*N_ice
    pam.df.dataFullSpec["n_ds"][:,:,:,id_ice,ibin]=0.#25*N_ice
    #pam.df.dataFullSpec["mass_ds"][:,:,:,id_ice,ibin-3]=0.25*q_hydro[isel:jsel,:,:,id_ice]
    #pam.df.dataFullSpec["mass_ds"][:,:,:,id_ice,ibin-2]=0.25*q_hydro[isel:jsel,:,:,id_ice]
    #pam.df.dataFullSpec["mass_ds"][:,:,:,id_ice,ibin-1]=0.25*q_hydro[isel:jsel,:,:,id_ice]
    pam.df.dataFullSpec["mass_ds"][:,:,:,id_ice,ibin]=0.#25*q_hydro[isel:jsel,:,:,id_ice]
    pam.df.dataFullSpec["area_ds"][:,:,:,id_ice,:]=np.pi/4. *  pam.df.dataFullSpec["d_ds"][:,:,:,id_ice,:] ** 2
    pam.df.dataFullSpec["as_ratio"][:,:,:,id_ice,:]=AR_ice
    #pam.df.dataFullSpec["canting"][:,:,:,id_ice,:]=0.
    #pam.df.dataFullSpec["fallvelocity"][:,:,:,id_ice,:]=snow_fallspeed#snow_fallspeed
    #pam.df.dataFullSpec["rg_beta_ds"][:,:,:,id_ice,:]=-99.
    #pam.df.dataFullSpec["rg_kappa_ds"][:,:,:,id_ice,:]=-99.
    #pam.df.dataFullSpec["rg_gamma_ds"][:,:,:,id_ice,:]=-99.
    #pam.df.dataFullSpec["rg_zeta_ds"][:,:,:,id_ice,:]=-99.
    #plt.imshow(N_ice[:,0,:]/1000,aspect='auto')
    plt.show()
    pam.df.dataFullSpec["fallvelocity"][:,:,:,id_liq,:]  = rain_fallspeed
    pam.df.dataFullSpec["fallvelocity"][:,:,:,id_rain,:] = rain_fallspeed
    pam.df.dataFullSpec["fallvelocity"][:,:,:,id_snow,:] = snow_fallspeed
    pam.df.dataFullSpec["fallvelocity"][:,:,:,id_ice,:]  = snow_fallspeed  # or ice fall speed

#__________PAMTRA namelist___________

#pam.nmlSet['tmatrix_db'] = 'file'
#pam.nmlSet['tmatrix_db_path'] = 'example_db/'
pam.nmlSet["passive"] = False
pam.nmlSet["radar_mode"] = "spectrum"#to deactivate the doppler spectra: simple


#__________scattering method____________
#pam.df.data["scat_name"][:] = "tmatrix"
#pam.df.data["as_ratio"][:] = 1.0

pam.set["pyVerbose"] = 1
Run_pamtra=False
Run_pamtra=True#False
if Run_pamtra==True:
    print("Start to run")
    pam.runParallelPamtra(35.0,
                      pp_deltaX=2,    # 2 profiles in X per worker
                      pp_deltaY=2,    # 1 profile in Y per worker
                      pp_deltaF=1,    # 1 frequency per worker
                      pp_local_workers="auto")  # detect CPU cores


#    pam.runPamtra([35.0], checkData=False)
    print("Run end")
    print("variables",pam.r.keys())

    print("Ze",np.shape(pam.r["Ze"]))#.keys())
    plt.figure('test')
    plt.imshow(pam.r["Ze"][:,0,:,0,0,0].T,aspect='auto',vmin=-30,vmax=30,cmap="plasma")
    plt.colorbar()
    plt.show()



    Ze_profile = pam.r['radar_moments'][0, :, :, 0, 0, 0, 0]

    plt.figure('test2')
    plt.imshow(Ze_profile.T,aspect='auto',cmap="plasma")
    plt.colorbar()
    plt.show()


    print("radar_spectra",np.shape(pam.r["radar_spectra"]))#.keys())
    plt.figure('doppler test')
    plt.imshow(np.sum(pam.r["radar_spectra"][0,:,:,0,0,:],1).T,aspect='auto')#,vmin=-30,vmax=30,cmap="plasma")
    plt.show()


    print("radar_moments",np.shape(pam.r["radar_moments"]))#.keys())
    plt.figure('radar_moments')
    plt.imshow(pam.r["radar_moments"][0,:,:,0,0,0,1].T,aspect='auto')#,vmin=-30,vmax=30,cmap="plasma")
    plt.show()


# Output NetCDF file path
#output_file = "/home/grzegorc/AWACA/PAMTRA/pamtra/Ouput_cosp_final_no_fullspec.nc"
output_file = "/home/grzegorc/AWACA/PAMTRA/pamtra/Ouput_golden_case_D17_v5.nc"
#output_file = "/home/grzegorc/AWACA/PAMTRA/pamtra/test.nc"#Ouput_golden_case_D17_v2_new_bin_wateronly.nc"
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
    nc_out.createVariable("N_ice", "f4", ("time", "col", "level"))[:] = N_ice
    nc_out.createVariable("N_liq", "f4", ("time", "col", "level"))[:] = N_liq
    nc_out.createVariable("N_snow", "f4", ("time", "col", "level"))[:] = N_snow
    nc_out.createVariable("N_rain", "f4", ("time", "col", "level"))[:] = N_rain
    nc_out.createVariable("Ze", "f4", ("time", "col", "level"))[:] = pam.r["Ze"][:,:,:,0,0,0]
    print("PAMTRA output saved as NetCDF file  ",output_file)
