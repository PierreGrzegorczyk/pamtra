from __future__ import print_function

import pyPamtra
import shutil
import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import imp
from netCDF4 import Dataset

#________Switch on different options


Run_pamtra=False
Run_pamtra=True


Run_parall=False
Run_parall=True

Run_spectra=True
Run_spectra=False

Write_output=True
Write_output=False




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
tke = nc_data.variables['tke'][:]
u = nc_data.variables['vitu'][:]
v = nc_data.variables['vitv'][:]
w = nc_data.variables['vitw'][:]

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
print("extreme", np.sum(Qs>0.05))


##_________quicklook data input_________

#plt.figure('Qi')
#plt.imshow(np.sum(Qi,1),aspect='auto')
#plt.colorbar()
#print("Show Qi")
#plt.show()

#plt.figure('Qs')
#plt.imshow(np.sum(Qs,1),aspect='auto')
#plt.colorbar()
#print("Show Qs")
#plt.show()


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
tke=np.repeat(tke[:, np.newaxis, :,0,0], ncol, axis=1)
u=np.repeat(u[:, np.newaxis, :,0,0], ncol, axis=1)
v=np.repeat(v[:, np.newaxis, :,0,0], ncol, axis=1)
w=np.repeat(w[:, np.newaxis, :,0,0], ncol, axis=1)

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

q_hydro[:,:,:,id_liq]=Ql
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

isel=0#3*24*3
jsel=-1#6*24*3

jsel=491
isel=490

pamData["lon"] = lon[isel:jsel,:]
pamData["lat"] = lat[isel:jsel,:]
pamData["temp"] = T[isel:jsel,:,:]
pamData["relhum"] = RH[isel:jsel,:,:]
pamData["hgt"] = z[isel:jsel,:,:]
pamData["press"] = p[isel:jsel,:,:]
pamData["hydro_q"] = q_hydro[isel:jsel,:,:]
pamData["airturb"] = T[isel:jsel,:,:]/T[isel:jsel,:,:]*0.01
tke[tke<0.012]=0.012
#tke[tke>5]=5
pamData["airturb"] = tke[isel:jsel,:,:]#/T[isel:jsel,:,:]*0.011
pamData["wind_w"] = T[isel:jsel,:,:]/T[isel:jsel,:,:]*0#0.011
print('keys pamdata',pamData.keys())
pamData["wind_w"] =-w[isel:jsel,:,:]#/T[isel:jsel,:,:]*0.1

plt.figure('w')
plt.plot(w[isel:jsel,0,:][0],z[isel:jsel,0,:][0])

plt.figure('tke')
plt.plot(tke[isel:jsel,0,:][0],z[isel:jsel,0,:][0])
plt.xlabel('tke')
plt.xscale('log')


#________hydrometeor input________
##___Liq_properties

r_liq=2e-5
Rho_liq=1000.

N_liq=(q_hydro[isel:jsel,:,:,id_liq]*Rho_air[isel:jsel,:,:])/(Rho_liq*4/3*np.pi*r_liq**3)
pam.df.addHydrometeor(("liq", -99., 1, Rho_liq, -99., -99., -99., -99. ,3,   1, "mono", -99., -99., -99., -99.,2*r_liq, -99.,"mie-sphere", "khvorostyanov01_drops", -99.))

##___Rain_properties___

r_rain=0.0005
rain_fallspeed=4.
Rho_rain=1000.
N_rain=q_hydro[isel:jsel,:,:,id_rain]/(Rho_rain*4/3*np.pi*r_rain**3)
pam.df.addHydrometeor(("rain",-99.,  1 , Rho_rain , -99., -99.,-99. , -99., 3 ,1,"mono",-99.0, -99.0, -99.0, -99.0,2*r_rain,-99.0,"mie-sphere",rain_fallspeed,0.0))

##___Snow_properties___
C_snow=1.
r_snow=0.001
AR_snow=1.
snow_fallspeed=1.
Rho_snow = 1.e3 * 0.178 * ( r_snow * 2 * 1000. )**(-0.922)
N_snow=(q_hydro[isel:jsel,:,:,id_snow]*Rho_air[isel:jsel,:,:])/(Rho_snow*4/3*np.pi*r_snow**3)

pam.df.addHydrometeor(("snow",C_snow, -1 , Rho_snow, -99., -99., np.pi/4., 2. ,  3 ,1,"mono",-99.0, -99.0, -99.0, -99.0,2*r_snow,-99.0,"ss-rayleigh-gans","lmdz_snow",0.0))
#pam.df.addHydrometeor(("snow",C_snow, -1 , Rho_snow, 130., 3.0 ,0.684, 2. ,  3 ,1,"mono_cosmo_ice",-99.0, -99.0, -99.0, -99.0,2*r_snow,-99.0,"ss-rayleigh-gans","heymsfield10_particles",0.0)) #gives strange things
#pam.df.addHydrometeor(("snow",C_snow, -1 , Rho_snow, 130., 3.0 ,0.684, 2. ,  3 ,1,"mono",-99.0, -99.0, -99.0, -99.0,2*r_snow*1e-3,-99.0,"mie-sphere","heymsfield10_particles",0.0))
print(pam.df)
##___Cirrus_properties___

C_ice=1.
Rho_ice=917.
AR_ice=1.#
r_ice=50e-6#((q_hydro[isel:jsel,:,:,id_cir]*Rho_air[isel:jsel,:,:])/(N_cir*Rho_cir*4/3*np.pi+1e-30))**(1/3)
N_ice=(q_hydro[isel:jsel,:,:,id_ice]*Rho_air[isel:jsel,:,:])/(Rho_ice*4/3*np.pi*r_ice**3)


#pam.df.addHydrometeor(("ice", C_ice, -1 , Rho_ice,  130., 3.0 ,0.684, 2.  , 3 ,1, "mono_cosmo_ice", -99., -99., -99., -99., 2*r_ice, -99., "ss-rayleigh-gans", "heymsfield10_particles",0.0))
pam.df.addHydrometeor(("ice", C_ice, -1 , Rho_ice,  -99,-99 ,np.pi/4, 2.  , 3 ,1, "mono", -99., -99., -99., -99., 2*r_ice, -99., "ss-rayleigh-gans", "heymsfield10_particles",0.0))
pam.createProfile(**pamData)

## ___ PAMTRA FULL SPECTRA ___
FULL_SPECTRA=False
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

    ##___RAIN:
    pam.df.dataFullSpec["rho_ds"][:,:,:,id_rain,:]=1000.
    pam.df.dataFullSpec["d_ds"][:,:,:,id_rain,:]=2*r_rain
    pam.df.dataFullSpec["d_bound_ds"][:,:,:,id_rain,:]=[r_rain,4*r_rain]
    pam.df.dataFullSpec["n_ds"][:,:,:,id_rain,0]=N_rain#q_hydro[isel:jsel,:,:,id_rain]/(1000*4/3*np.pi*r_rain**3)
    pam.df.dataFullSpec["mass_ds"][:,:,:,id_rain,0]=q_hydro[isel:jsel,:,:,id_rain]
    pam.df.dataFullSpec["area_ds"][:,:,:,id_rain,:]=np.pi*r_liq**2
    pam.df.dataFullSpec["as_ratio"][:,:,:,id_rain,:]=1.
    pam.df.dataFullSpec["canting"][:,:,:,id_rain,:]=0.
    pam.df.dataFullSpec["fallvelocity"][:,:,:,id_rain,:]=rain_fallspeed

    ##___SNOW:
    pam.df.dataFullSpec["rho_ds"][:,:,:,id_snow,:]=Rho_snow
    pam.df.dataFullSpec["d_ds"][:,:,:,id_snow,:]=2*r_snow
    pam.df.dataFullSpec["d_bound_ds"][:,:,:,id_snow,:]=[r_snow,3*r_snow]
    pam.df.dataFullSpec["n_ds"][:,:,:,id_snow,0]=N_snow
    pam.df.dataFullSpec["mass_ds"][:,:,:,id_snow,0]=q_hydro[isel:jsel,:,:,id_snow]
    pam.df.dataFullSpec["area_ds"][:,:,:,id_snow,:]=np.pi*r_snow**2
    pam.df.dataFullSpec["as_ratio"][:,:,:,id_snow,:]=AR_snow
    pam.df.dataFullSpec["canting"][:,:,:,id_snow,:]=0.
    pam.df.dataFullSpec["fallvelocity"][:,:,:,id_snow,:]=snow_fallspeed
#    pam.df.dataFullSpec["rg_beta_ds"][:,:,:,id_snow,:]=-99.
#    pam.df.dataFullSpec["rg_kappa_ds"][:,:,:,id_snow,:]=-99.
#    pam.df.dataFullSpec["rg_gamma_ds"][:,:,:,id_snow,:]=-99.
#    pam.df.dataFullSpec["rg_zeta_ds"][:,:,:,id_snow,:]=-99.
#    RG_BETA  = 0.22
#    RG_KAPPA = 2.52
#    RG_GAMMA = 2.36
#    RG_ZETA  = 0.049

#"   pam.df.dataFullSpec["rg_beta_ds"][:,:,:,id_snow,:]  = RG_BETA
#   pam.df.dataFullSpec["rg_kappa_ds"][:,:,:,id_snow,:] = RG_KAPPA
#    pam.df.dataFullSpec["rg_gamma_ds"][:,:,:,id_snow,:] = RG_GAMMA
#    pam.df.dataFullSpec["rg_zeta_ds"][:,:,:,id_snow,:]  = RG_ZETA



##___ICE
    pam.df.dataFullSpec["rho_ds"][:,:,:,id_ice,:]=Rho_ice
    pam.df.dataFullSpec["d_ds"][:,:,:,id_ice,0]=2*r_ice
    pam.df.dataFullSpec["d_bound_ds"][:,:,:,id_ice,:]=[0,4]#4*np.max(r_cir)]
    pam.df.dataFullSpec["n_ds"][:,:,:,id_ice,0]=N_ice
    pam.df.dataFullSpec["mass_ds"][:,:,:,id_ice,0]=q_hydro[isel:jsel,:,:,id_ice]
    pam.df.dataFullSpec["area_ds"][:,:,:,id_ice,:]=np.pi*r_ice**2
    pam.df.dataFullSpec["as_ratio"][:,:,:,id_ice,:]=AR_ice
    pam.df.dataFullSpec["canting"][:,:,:,id_ice,:]=0.
    pam.df.dataFullSpec["fallvelocity"][:,:,:,id_ice,:]=0.#snow_fallspeed
    plt.imshow(N_ice[:,0,:]/1000,aspect='auto')
    plt.show()

#__________namelist___________

#pam.nmlSet['tmatrix_db'] = 'file'
#pam.nmlSet['tmatrix_db_path'] = 'example_db/'
pam.nmlSet["randomseed"] = 10
pam.nmlSet["passive"] = False
pam.nmlSet['radar_airmotion']=True
#pam.nmlSet['radar_allow_negative_dD_dU'] = True

#new for doppler
pam.nmlSet["radar_mode"] = "spectrum"
pam.nmlSet["radar_aliasing_nyquist_interv"] = 10
pam.nmlSet["hydro_adaptive_grid"] = False
pam.nmlSet["conserve_mass_rescale_dsd"] = False
pam.nmlSet["radar_use_hildebrand"] = True
pam.nmlSet["radar_noise_distance_factor"] = -2
#pam.nmlSet["radar_noise_distance_factor"] = 0#-6#0.5
pam.nmlSet["radar_save_noise_corrected_spectra"]=   False
pam.nmlSet["radar_nfft"]=256
pam.nmlSet["radar_use_wider_peak"]=True
#pam.nmlSet['radar_nPeaks']=1.
#pam.nmlSet['radar_peak_min_bins']=-2.
#pam.nmlSet['radar_pnoise0']= -38.23
           #pam.nmlSet["radar_smooth_spectrum"]=True
#pam.nmlSet["radar_use_wider_peak"] = True
#pam.nmlSet["radar_noise_distance_factor"] = 0.#-2
#pam.nmlSet['radar_max_v']= 18.
#pam.nmlSet['radar_min_v']= -18.
#pam.nmlSet["save_psd"] = True

#__________scattering method____________
#pam.df.data["scat_name"][:] = "tmatrix"
#pam.df.data["scat_name"][:] = "ssrga"
#pam.df.data["as_ratio"][:] = 1.0

pam.set["pyVerbose"] = 2

if Run_pamtra==True:
    print("Start to run")

    if Run_parall==True:
        pam.runParallelPamtra(35.0,
                      pp_deltaX=1,    # 2 profiles in X per worker
                      pp_deltaY=1,    # 1 profile in Y per worker
                      pp_deltaF=1,    # 1 frequency per worker
                      pp_local_workers="auto")  # detect CPU cores
    else:
        pam.runPamtra(35.0,checkData=False)

    print("Run end")

    print("Pam.r keys",pam.r.keys())

    print("SHAPE",np.shape(pam.r["Ze"]))
    print("MAX",np.max(pam.r["Ze"]))


    print('shape',np.shape(pam.r["radar_vel"]),np.shape(pam.r["radar_spectra"]))
    print("max min vel",np.min(pam.r["radar_vel"][:]),np.max(pam.r["radar_vel"][:]))
    print("max dBZ",np.max(pam.r["radar_spectra"][:]))


    #pam.r["radar_spectra"][pam.r["radar_spectra"]<-35]=0.
    plt.figure('Quicklook reflectivity')
    plt.imshow(pam.r["Ze"][:,0,:,0,0,0].T,aspect='auto',vmin=-30,vmax=30,cmap="jet")
    #plt.pcolormesh(np.array([1]),Alt/1000,pam.r["Ze"][:,0,:,0,0,0].T,vmin=-50,vmax=30,cmap="jet")
    plt.colorbar() 


    plt.figure('reflectivity profile')
    plt.plot(pam.r["Ze"][0,0,:,0,0,0],Alt/1000)
    plt.xlim(-30,30) 
    plt.ylim(0,10) 


    plt.figure('Quicklook doppler spectra subcol 1')
    plt.pcolormesh(pam.r["radar_vel"][:],Alt/1000,pam.r["radar_spectra"][0,0,:,0,0,:],vmin=-50,vmax=30,cmap="plasma")
    plt.ylim(0,10)
    #plt.imshow(pam.r["radar_spectra"][0,0,:,0,0,:],vmin=-50,vmax=30,cmap="plasma",aspect='auto')
    plt.colorbar()



    plt.figure('Quicklook doppler spectra all subcol')
    plt.pcolormesh(pam.r["radar_vel"][:],Alt/1000,np.mean(pam.r["radar_spectra"][0,:,:,0,0,:],0),vmin=-50,vmax=30,cmap="plasma")
    plt.ylim(0,10)
    #plt.imshow(pam.r["radar_spectra"][0,0,:,0,0,:],vmin=-50,vmax=30,cmap="plasma",aspect='auto')
    plt.colorbar()
    plt.show()







# Output NetCDF file path
#output_file = "/home/grzegorc/AWACA/PAMTRA/pamtra/Ouput_cosp_final_no_fullspec.nc"
output_file = "/home/grzegorc/AWACA/PAMTRA/pamtra/Ouput_golden_case_D17_v5_ssrga.nc"
#output_file = "/home/grzegorc/AWACA/PAMTRA/pamtra/test.nc"#Ouput_golden_case_D17_v2_new_bin_wateronly.nc"

if Write_output==True:
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
