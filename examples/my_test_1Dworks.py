from __future__ import print_function

import pyPamtra
import shutil
import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import imp
from netCDF4 import Dataset

#_________LMDZ data___________________
path = "/home/grzegorc/AWACA/COSP/output_lmdz"
nc_file = path + "/histhf.nc"

nc_data = Dataset(nc_file, "r")

lon = nc_data.variables['lon'][:]
lat = nc_data.variables['lat'][:]
Alt = nc_data.variables['Alt'][:] #if 'Alt' in nc_data.variables else None
T = nc_data.variables['temp'][:]       # shape: (time, level, lat, lon) — check yours
RH = nc_data.variables['rhum'][:]*100
p = nc_data.variables['pres'][:]
p[p<1.]=1. #minimum value accepted by PAMTRA
z = nc_data.variables['zfull'][:]
time = nc_data.variables['time_counter'][:]


#_________cosp data subcolumns________
path="/home/grzegorc/AWACA/COSP/COSPv2.0_lmdz_hillman/driver/run"
nc_file = path+"/hydro_output.nc"
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

lon=np.full(np.shape(T)[:-1],lon)
lat=np.full(np.shape(T)[:-1],lat)

q_hydro=np.zeros(np.shape(T))
q_hydro=np.repeat(q_hydro[:,:,:,np.newaxis],5, axis=3)

q_hydro[:,:,:,0]=Ql
q_hydro[:,:,:,1]=Qi
q_hydro[:,:,:,2]=Qr
q_hydro[:,:,:,3]=Qs



#________load PAMTRA______________
imp.reload(pyPamtra)
pam = pyPamtra.pyPamtra()

jsel=66 #select time
#________hydrometeor input________

descriptorFile = "../descriptorfiles/descriptor_file_COSMO_1mom.txt"
descriptorFile = "../descriptorfiles/LMDZ2.txt"

#pam.df.readFile(descriptorFile)

#Hydrometeor____properties

#SNOW
r_snow=0.001
Rho_snow = 1.e3 * 0.178 * ( r_snow * 2 * 1000. )**(-0.922)
snow_fallspeed=1.
print('rho_snow',Rho_snow)
#snow are assumed to be spherical

#RAIN
r_rain=0.0005
rain_fallspeed=4.

#ICE
Rho_ice=917.
r_ice=5e-5

pam.df.addHydrometeor(("liq", -99., 1, 1000, -99., -99., -99., -99. ,3,   1, "mono", -99., -99., -99., -99.,2e-5, -99.,"mie-sphere", "khvorostyanov01_drops", -99.))
pam.df.addHydrometeor(("ice", -99., -1 , Rho_ice,  130., 3.0 ,0.684, 2.  , 3 ,1, "mono_cosmo_ice", -99., -99., -99., -99., 2*r_ice, -99., "mie-sphere", "heymsfield10_particles",0.0))
pam.df.addHydrometeor(("rain",-99.,  1 , 1000 , -99., -99.,-99. , -99., 3 ,1,"mono",-99.0, -99.0, -99.0, -99.0,2*r_rain,-99.0,"mie-sphere",rain_fallspeed,0.0))
pam.df.addHydrometeor(("snow",-99., -1 , Rho_snow, -99., -99.,-99. , -99., 3 ,1,"mono",-99.0, -99.0, -99.0, -99.0,2*r_snow,-99.0,"mie-sphere",snow_fallspeed,0.0))
pam.df.addHydrometeor(("cirrus", -99., -1 , 920.,  130., 3.0 ,0.684, 2.  , 3 ,1, "mono_cosmo_ice", -99., -99., -99., -99., 2*r_ice, -99., "mie-sphere", "heymsfield10_particles",0.0))

#pam.df.addHydrometeor(("cirrus", -99., 1 , 1000,  -99., -99.,-99. ,-99.  , 3 , 1, "mono_cosmo_ice", -99., -99., -99., -20e-6, -99., -99., "mie-sphere", "heymsfield10_particles",0.0))



print('pam.df.nhydro',pam.df.nhydro)
print('ice', pam.df)



#________input data_____________
#pam.readPamtraProfile("../profile/example_input.lay")
#pam.filterProfiles(np.array([[True, False], [False, False]]))
#for k in pam.p.keys():
#    print(k, np.shape(pam.p[k]))

pamData = dict()
pamData["lon"] = lon[jsel,0]
pamData["lat"] = lat[jsel,0]
pamData["temp"] = T[jsel,0,:]
pamData["relhum"] = RH[jsel,0,:]
pamData["hgt"] = z[jsel,0,:]
pamData["press"] = p[jsel,0,:]
pamData["hydro_q"] = np.max(q_hydro[:,0,:],0)

print("z shape",np.shape(pamData["lon"]),np.shape(pamData["relhum"]))
pam.createProfile(**pamData)

#print('pam.p.keys()',pam.p.keys())
#print(pam.p)
#print(pam.p['hydro_d'])


#data figure
old_figure=False
if old_figure == True: 
    fig = plt.figure("meteo profile",figsize=[12,4])
    fig.add_axes([0.1,0.1,0.25,0.8])
    plt.plot(pam.p['press'][0,0,:]/100.,pam.p['hgt'][0,0,:])
    plt.ylabel('height [m]')
    plt.xlabel('pressure [hPa]')
    fig.add_axes([0.425,0.1,0.25,0.8])
    plt.plot(pam.p['temp'][0,0,:],pam.p['hgt'][0,0,:])
    plt.xlabel('temperature [K]')
    fig.add_axes([0.75,0.1,0.25,0.8])
    plt.plot(pam.p['relhum'][0,0,:],pam.p['hgt'][0,0,:])
    plt.xlabel('rel. humidity [%]')
    plt.show()


#__________namelist___________

pam.nmlSet['tmatrix_db'] = 'file'
pam.nmlSet['tmatrix_db_path'] = 'example_db/'
pam.nmlSet["passive"] = False
#print('namelist',pam.nmlSet)


#__________scattering method____________
pam.df.data["scat_name"][:] = "tmatrix"
pam.df.data["as_ratio"][:] = 1.0

print("tmatrix or not", pam.df.data["scat_name"])
#pam.set["pyVerbose"] = 10
if 0==0: 
    print("Start to run")
    pam.runPamtra(94.0)
    #print("End to run")
    print(np.shape(pam.r["Ze"]))
    plt.figure("reflectivity")
    plt.subplot(121)
    plt.plot(pam.r['Ze'][0,0,:,0,0,0],pam.p['hgt'][0,0,:]/1000)
    plt.ylabel('height [m]')
    plt.xlim(-30,30)
    plt.xlabel('reflectivity [dBZ]')

    plt.subplot(122)
    plt.plot(pam.p['hydro_q'][0,0,:]*1000,pam.p['hgt'][0,0,:]/1000)
    plt.ylabel('height [m]')
    plt.xlim(0,1)
    plt.xlabel('qi')

    plt.show()
    print(pam.r)
    print("Run end")

    #shutil.rmtree('example_db')
    #print('cleaning up done')
