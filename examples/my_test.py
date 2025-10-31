from __future__ import print_function

import pyPamtra
import shutil
import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import imp
from netCDF4 import Dataset

#_________LMDZ data___________________
path = "/mnt/c/Users/grzegorczyk/AWACA/COSP/output_lmdz"
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
path="/mnt/c/Users/grzegorczyk/AWACA/COSP/COSPv2.0_lmdz_hillman/driver/run"
nc_file = path+"/hydro_output.nc"
nc_data = Dataset(nc_file, "r")

Qi=nc_data['I_LSCICE'][:]
Ql=nc_data['I_LSCLIQ'][:]
Qr=nc_data['I_LSRAIN'][:]
Qs=nc_data['I_LSSNOW'][:]

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
#q_hydro[:,:,:,0]=Ql
#q_hydro[:,:,:,1]=Qi
#q_hydro[:,:,:,2]=Qr
#q_hydro[:,:,:,3]=Qs

#________load PAMTRA______________
imp.reload(pyPamtra)


pam = pyPamtra.pyPamtra()
#________hydrometeor input________
descriptorFile = "../descriptorfiles/LMDZ.txt"
pam.df.readFile(descriptorFile)
print('pam.df.nhydro',pam.df.nhydro)

#________input data_____________
#pam.readPamtraProfile("../profile/example_input.lay")
#pam.filterProfiles(np.array([[True, False], [False, False]]))
#for k in pam.p.keys():
#    print(k, np.shape(pam.p[k]))

pamData = dict()
pamData["lon"] = lon
pamData["lat"] = lat
pamData["temp"] = T
pamData["relhum"] = RH
pamData["hgt"] = z
#print("z shape",np.shape(z))


pamData["press"] = p
pamData["hydro_q"] = q_hydro
pam.createProfile(**pamData)

#print('pam.p.keys()',pam.p.keys())
#print(pam.p)
#print(pam.p['hydro_q'])

#data figure
old_figure=False
if old_figure == True: 
    fig = plt.figure(figsize=[12,4])
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


#pam.set["pyVerbose"] = 10
if 0==0: 
    print("Start to run")
    pam.runPamtra(89.0)
    #print("End to run")

    print(np.shape(pam.r["Ze"]))

    print("##########################")

    print('cleaning up')
    #shutil.rmtree('example_db')
    #print('cleaning up done')
