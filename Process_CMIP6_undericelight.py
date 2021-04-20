#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 11:29:04 2021

@author: stroeve
"""

import sys
sys.path.insert(0,'/anaconda2/pkgs')

import pandas as pd
import numpy as np
import numpy.ma as ma
import os


import xarray as xr
import netCDF4 as nc
from netCDF4 import Dataset, num2date

#these are used for the regridding tool
import pyproj as proj
from scipy.interpolate import griddata

#these are needed for plotting
import matplotlib as mpl
# from mpl_toolkits.basemap import Basemap
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib.cm as cm2

import cftime
from cftime import DatetimeNoLeap

import cartopy
import cartopy.crs as ccrs

#from cftime import num2date, date2num

from skimage import data, color
from skimage.transform import rescale, resize, downscale_local_mean
from pylab import *
# from PIL import Image


# import gzip

#from scipy import interpolate as sp
#from scipy.ndimage.interpolation import zoom
#from scipy.ndimage.interpolation import rotate
from datetime import datetime
#from scipy.stats.mstats import pearsonr
#from scipy.stats import linregress
#import scipy.stats as ss
#from matplotlib.ticker import FuncFormatter
#import matplotlib.mlab as mlab
from scipy.stats import norm
from matplotlib import ticker
from scipy import optimize


#%% LOAD Function
import datetime as dt
def output2Netcdf(data_slice,modeltime,fn,itype,model_name):
    try: ncfile.close()
    except: pass
#we will output all the under-ice par and transmittance in EASE-grid    
    xdim=361
    ydim=361
    
    ncfile=Dataset(fn,'a',format='NETCDF4')
    ncfile.title='Model data '+model_name
    
    #first create dimensions
    time_dim = ncfile.createDimension("time", None)  # can be appended to
    lat_dim=ncfile.createDimension("lat_EASE",xdim)
    lon_dim=ncfile.createDimension("lon_EASE",ydim)
    
    #create variables
    times=ncfile.createVariable("time","f4",("time",))
    times.units='year since '+str(yearstart)
    times.long_name='time'
    latitudes = ncfile.createVariable("lat", "f4", ("lat",))
    lat.units='degrees_north'
    lat.long_name='latitude'
    longitudes = ncfile.createVariable("lon", "f4", ("lon",))
    lon.units='degrees_ease'
    lon.long_name='longitude'
    
    #define a 3D variable to hold the data
    if itype == 'trans':
        value=ncfile.createVariable("transmissivity","f4",("time","lat","lon"),zlib='True')
        value.units="None"
        value.standard_name='Transmittance under sea ice'
    elif itype == 'par':
        value=ncfile.createVariable("PAR","f4",("time","lat","lon"),zlib='True')
        value.units='micromoles per square meter per second'
        value.standard_name='Under-ice PAR'
    
    value[:,:,:]=data_slice
    times=modeltime
    ncfile.close(); print(fn+' data set is closed')
    
    # ds1=xr.Dataset( data_vars={'PAR':(['x','y'],Fsw_TR_NEW),
    #                         'dXr':(['x','y'],dXr),
    #                         'dYr':(['x','y'],dYr)},

    #              coords =  {'lon':(['x','y'],ease_lons),
    #                         'lat':(['x','y'],ease_lats)})
    # ds2=xr.Dataset( data_vars={'Transmittance':(['x','y'],T_snow),
    #                         'dXr':(['x','y'],dXr),
    #                         'dYr':(['x','y'],dYr)},

    #              coords =  {'lon':(['x','y'],ease_lons),
    #                         'lat':(['x','y'],ease_lats)})


#%% LOAD Function
def get_under_ice_light(sic,sit,snd,alb,swd,latsy,lonsx,modelname,yearstart):
# #**********function to compute under-ice light*************
    #create the output netCDF data sets
    nyears=sic.shape[0]
    print('nyears ',nyears)
    year=[]
    big_value=1.e20
    xdim=sic.shape[1]
    ydim=sic.shape[2]

    transmittance=np.empty((nyears,361,361))
    par=np.empty((nyears,361,361))

#since we want to put all data in EASE we need to load those lat/lons
    # Get EASE grid
    EASE_grid = Dataset('/Users/stroeve/Documents/seaice/grid.nc')
    lons_ease = np.array(EASE_grid['lon'])
    lats_ease = np.array(EASE_grid['lat'])
    
    print('xdim and ydim ',xdim,ydim,nyears)
# surface transmission parameter
    i_0_snw  =  0.3;     # for snow 
    i_0_bi   =  0.3;     # for ice  (may change) 
    i_0_pond =  0.58;    # for melt pond
    i_0_ow   =  0.58;    # for open water
# surface transmission layer thickness (m)
    h_0_snw  =  0.03;     # for snow 
    h_0_bi   =  0.1;    # for bare ice ice  (may change !) 
    alb_pond = 0.27;  # for ponded ice 
# attenuation coeff ice    
    K_i  = 1;
# ITD pdf    
    hpdf = [0.0646, 0.1415, 0.173, 0.1272, 0.1114, 0.0824, 0.0665, \
        0.0541, 0.0429, 0.0347, 0.0287, 0.024,0.0194, 0.016, 0.0136]
    hpdf=np.asarray(hpdf)
    for iyears in range(nyears):        ## big loop on years (could also have been on months as Gaelle did)
        year.append(iyears+yearstart)
        print('processing year ',iyears+yearstart)
        SIC=sic[iyears,:,:]
        f_bi=reshape(SIC,xdim*ydim)      ## re-shape sea ice concentration as 1D array for loop
        # print('inside years loop for undericelight ',range(nyears),iyears,alb.shape)
    

    #### Snow depth ####
        snod=snd[iyears,:,:]
        h_s=reshape(snod,xdim*ydim)   #re-shape to 1D array for loop 
    
    #### Sea ice thickness ####
        SIT=sit[iyears,:,:]
        hi=reshape(SIT,xdim*ydim)
        
    #### albedo and incoming solar #### 
        ALB=alb[iyears,:,:]
        alb_new=reshape(ALB,xdim*ydim)

    #### Incoming solar irradiance ####
        FSW0=swd[iyears,:,:]
        Fsw0=reshape(FSW0,xdim*ydim)
        

## ITD 15 classes ##
        h_cutoff=3.   # maximum 3 times mean ice thickness       
        hi15=np.zeros([xdim*ydim,15]) ## ice thickness distribution
        for ii in range (1,16):
            factor2=(2*ii - 1)/15.    ## factor for distribution
            hi15[:,ii-1]=hi*factor2     
            hi15[:,ii-1]=(h_cutoff/2.)*hi15[:,ii-1] ## h_cutoff vary with the SIC distributed
        #done with initializing the 15 SIT distribution classes


# initialisation arrays transmittance, irradiance    
        T_ow=np.zeros(xdim*ydim)
        T_snow=np.zeros(xdim*ydim)
        Fsw_tr_new=np.zeros([xdim*ydim,16])
# hs array between 0 and 2 times mean snow depth
        hs=np.arange(0.,2.*nanmax(h_s),0.05)
# initialisation snow attenuation coeff
        K_s=np.zeros([len(hs),xdim*ydim])

# initialisation intermediate steps
        f_att_snow=np.zeros([xdim*ydim,len(hs),15])
        i_s_hom=np.zeros([xdim*ydim,16])
        t_s_hom=np.zeros([xdim*ydim,16])

    #loops on ITD and assign K_s according to wet/dry snow
        for jj in range(0,15): 
            for ii in range (0,xdim*ydim):
                for kk in range(0,len(hs)):
                    # if temp[ii] <= 0:   # dry snow
                    K_s[kk,ii] = 10;
                    f_att_snow[ii,kk,jj]  = exp(-K_s[kk,ii]*(hs[kk]-h_0_snw))* exp(-K_i*hi15[ii,jj]) ;
                    # elif hs[kk]<=0.03:  # if snow depth < SSL thickness
                    if hs[kk]<=0.03:  # if snow depth < SSL thickness
                        K_s[kk,ii] = 40; 
                        f_att_snow[ii,kk,jj]  = exp(-K_s[kk,ii]*(hs[kk]))* exp(-K_i*hi15[ii,jj]) ;
                        # else:
                        #     K_s[kk,ii]=7;   # wet snow
                        #     f_att_snow[ii,kk,jj]  = exp(-K_s[kk,ii]*(hs[kk]-h_0_snw))* exp(-K_i*hi15[ii,jj]) ;
        
                #end loop for ii in range (0,xdim*ydim) and kk in range (0,len(hs))
            index=where(f_att_snow>1)
            f_att_snow[index]=1
        
        
            for ii in range (0,xdim*ydim):
                # distribution snow depth 
                g_hs = zeros([len(hs)]) ; 
                g_hs[where(hs  < 2.0 * h_s[ii])]  = f_bi[ii] / ( 2.0 * h_s[ii] );
                i_s_hom[ii,jj] = sum(g_hs * 0.01 * f_att_snow[ii,:,jj]); # integrated transmission
                for kk in range (0,len(hs)):
                    # computing transmittance
                    if K_s[kk,ii]==40:
                        t_s_hom[ii,jj] = ((1-alb_new[ii]) * i_s_hom[ii,jj]) / f_bi[ii] ;
                    else: 
                        t_s_hom[ii,jj] = ((1-alb_new[ii]) * i_0_snw * i_s_hom[ii,jj]) / f_bi[ii] ;     
                
                T_ow[ii] = (1 - alb_new[ii]);  ## transmittance open water
                # under-ice irradiance and PAR calculation
                Fsw_tr_new[ii,jj] = Fsw0[ii]* ((t_s_hom[ii,jj] * f_bi[ii]*3.51) + (T_ow[ii] * (1-f_bi[ii]))*2.30); 
               #end loop for ii  in range (0,xdim*ydim)
        
# sum ITD 15 classes and apply pdf 
        #end loop for jj in range(0,15)
        for i in range (0,xdim*ydim): 
            t_s_hom[i,15]=sum(t_s_hom[i,0:15]*hpdf[0:15])
            Fsw_tr_new[i,15]=sum(Fsw_tr_new[i,0:15]*hpdf[0:15])

        #end loop for SIT classse
#       print('Mean under-ice PAR = ')
        # print(nanmean(Fsw_tr_new[:,15]))
        sys.stdout.flush()    

# change shape to 2D array for map plotting
    
        T_snow=np.zeros([xdim,ydim])
        Fsw_TR_NEW=np.zeros([xdim,ydim])  #under-ice PAR
        c=-1
        for ii in range (0,xdim):
            for jj in range (0,ydim):
                c=c+1
                T_snow[ii,jj]=t_s_hom[c,15]
                Fsw_TR_NEW[ii,jj]=Fsw_tr_new[c,15]

            T_snow=flipud(T_snow)
            Fsw_TR_NEW=flipud(Fsw_TR_NEW)

    
        c=-1
        H_S=np.zeros([xdim,ydim])
        for ii in range (0,xdim):
            for jj in range (0,ydim):
                c=c+1
                H_S[ii,jj]=h_s[c]
    
#   use mask on transmittance and under ice PAR using snow depth
            # Fsw_TR_NEW=np.ma.array(Fsw_TR_NEW,mask=(isnan(flipud(H_S))==True))
            # T_snow=np.ma.array(T_snow,mask=(isnan(flipud(H_S))==True))

            
        # print('testing the size of the under-ice PAR and what year we are in ', Fsw_TR_NEW.shape, iyears, str(year[iyears]))
        plt.imshow(Fsw_TR_NEW) 
        #regrid to EASE
        Fsw_TR_EASE=regrid(Fsw_TR_NEW,lonsx,latsy,lons_ease,lats_ease)
        T_snow_EASE=regrid(T_snow,lonsx,latsy,lons_ease,lats_ease)
        par[iyears,:,:]=Fsw_TR_EASE
        transmittance[iyears,:,:]=T_snow_EASE
        
        plt.imshow(Fsw_TR_EASE)
        # plot(Fsw_TR_EASE,lats_ease,lons_ease,'PAR')
        # fname='/Volumes/Lacie/CMIP6/Ecolight/UnderIcePAR_'+modelname+'_April_'+str(year[iyears])
        # print('filename ',fname)
        # plt.savefig(fname,dpi=300)
        # plt.close()
        
    #end the loop on all years
    
    return par,transmittance
    
    
#function to regrid data to the same grid (this is a problem for the incoming and outgoing solar radiation fields)
args = proj.Proj(proj="aeqd", lat_0=90, lon_0=0, datum="WGS84", units="m")
crs_wgs = proj.Proj(init='epsg:4326')  # assuming you're using WGS84 geographic

#%% LOAD Function
def regrid(data_in,
           lon_in,
           lat_in,
           lon_out,
           lat_out,
           method='nearest'):

    xout, yout = proj.transform(crs_wgs, args, np.array(lon_out),np.array(lat_out))

    xin, yin = proj.transform(crs_wgs, args, np.array(lon_in),np.array(lat_in))

    output = griddata((xin.ravel(),yin.ravel()),
                    np.array(data_in).ravel(),
                    (xout,yout),
                    method=method)
    
    return(output)

#%% LOAD Function
def plot(data,latsy,lonsx,string):
    #minv=data.min()
    #maxv=data.max()
    
    # cmap1=plt.cm.jet
    # cmap1.set_bad('white')
    # cmap2=plt.cm.bwr
    # cmap2.set_bad('white')
    
#     cmap2=plt.cm.bwr
#     cmap2.set_bad('white')
#    print('min and max of data values ',minv,maxv)
    if string == 'siflswutop' or string == 'siflswdtop':
        label_unit='W/m2'
        minv=50.
        maxv=400.
        
    elif string == 'siconc':
        label_unit=' %'
        minv=0.
        maxv=1.
    elif string == 'sithick':
        label_unit='m'
        minv=0.0
        maxv=5.0
    elif string == 'sisnthick':
        data=data*100.
        label_unit='cm'
        minv=0.0
        maxv=50.0    
    elif string == 'PAR':
        label_unit=r'$\mu mol.m^{-2}.s^{-1}$'
        minv=0
        maxv=12
    elif string == 'alb':
        label_unit='albedo'
        minv=0
        maxv=1.0
        
#switch to using cartopy for plotting
    fig=plt.figure(figsize=(6.5,6.5))
    ax=plt.axes(projection=ccrs.NorthPolarStereo())
    ax.set_extent([-180,180,90,66],ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, edgecolor='black',zorder=1)
    
    bg=ax.pcolormesh(lonsx,latsy,data,vmin=minv,vmax=maxv,transform=ccrs.PlateCarree(),cmap='plasma')
    cb=fig.colorbar(gb,orientation='vertical',shrink=1)
    cb.set_label(label_unit,fontsize='x-large')
        
# #    print('inside plot routine ',latsy.shape,lonsx.shape,data.shape)
#     lat_0 = 90 # Latitude of Origin
#     lon_0 = -45 # Central Meridian
# #bbox = [45,-180.,90,180.] # Bounding Box
#     bbox = [45,-45,45,35]
#     fig=plt.figure(figsize=(6.5,6.5))
#     ax = fig.add_subplot() 
#     mp = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # choice of lat, lon boundariessic
#     mp.fillcontinents(color='white',lake_color='white')
#     mp.drawcoastlines()
#     mp.drawparallels(np.arange(-80.,81.,20.))
#     mp.drawmeridians(np.arange(-180.,181.,20.))
#     mp.drawmapboundary(fill_color='white') 
# #    print('min and max values and step ',minv,maxv,(maxv-minv)/10.)
#     cf1 = mp.contourf(lonsx, latsy, data, np.arange(minv,maxv,(maxv-minv)/10.), latlon=True, cmap=plt.cm.viridis, extend='both')
# # Add an annotation for each sub-plot
#     box = dict(boxstyle='square,pad=0.3',fc='white',ec='white',lw=0.5) # white rectangle with white outline
#     ax.annotate(string, xy=(0.02,0.96), xycoords='axes fraction', va='top', ha= 'left', size=8, bbox=box)
    
# # # The above adds text plotted on top of the white box located in the upper-left corner of each map
#     hc=plt.colorbar(cf1,);  
#     hc.set_label(label_unit,size=12);

#this program will run through an emission scenario to find the models to compute under-ice light following Stroeve et al. 2021
#the light processing component is calling a function based off of Gaelle's work

#inialise location of the CMIP6 data and the variable names

filepath='/Volumes/LaCie/CMIP6/'
minlat=45 # in degrees sets the minimum latitude limit

#variables needed for each model name
ncvarlist=['siconc','sithick','sisnthick','siflswutop','siflswdtop']
experiment=['historical','ssp126','ssp585']
months=['April','May','June','July','August','September','October']

#get model names from .csv file
filenames=pd.read_csv(filepath+'CMIP6_modelnames.csv',header='infer')

#%% LOAD Function
def find_files(datapath,variable):
    listing = sorted(os.listdir(f'{datapath}'))
    good_files=[]
    for file in listing:
        if variable in file:
            good_files.append(file)       
    return(good_files)

#we can loop through the models but we only proceed if we find files for each model name for the 5 variables

#search for each file
datapath=filepath+experiment[1]+'/'

sic_files=find_files(datapath,ncvarlist[0])
sit_files=find_files(datapath,ncvarlist[1])
snd_files=find_files(datapath,ncvarlist[2])
swu_files=find_files(datapath,ncvarlist[3])
swd_files=find_files(datapath,ncvarlist[4])
print(len(sic_files),len(sit_files),len(snd_files),len(swu_files),len(swd_files))

#now we need to remove the variable name from each in order to find the subsets
subset_sic=[]
for f in sic_files:
    s=f[7::]
    subset_sic.append(s)
subset_sit=[]
for f in sit_files:
    s=f[8::]
    subset_sit.append(s)
subset_snd=[]
for f in snd_files:
    s=f[10::]
    subset_snd.append(s)
subset_swu=[]
for f in swu_files:
    s=f[11::]
    subset_swu.append(s)
subset_swd=[]
for f in swd_files:
    s=f[11::]
    subset_swd.append(s)

 #%% LOAD Function   
def get(inpath,string):
    if string == 'lon':
       ncf=nc.Dataset(inpath)
       try:
          lons=ncf.variables['lon'][:].data 
       except:
          lons=ncf.variables['longitude'][:].data
       ncf.close()
       return(lons)
               
    elif string == 'lat':
        ncf=nc.Dataset(inpath)
        try:
          lats=ncf.variables['lat'][:].data 
        except:
          lats=ncf.variables['latitude'][:].data
        ncf.close()
        return(lats)
    elif string == 'time':
        ncf=nc.Dataset(inpath)
        timeslen = ncf['time'].shape[0]
        ncf.close()
        return(time)
    
#****************************************************
#main body of code that loops through all models
set_of_siconc_models = set(subset_sic)    
set_of_sithick_models = set(subset_sit)
set_of_snd_models = set(subset_snd)
set_of_swu_models = set(subset_swu)
set_of_swd_models = set(subset_swd) 

#this finds a set of all models/ensemble runs that are there for each variable
models_with_all_variables = (set_of_siconc_models & set_of_sithick_models & set_of_snd_models & set_of_swu_models & set_of_swd_models) # Find intersection of sets
#convert set to list
model_list=list(models_with_all_variables)
model_list.sort()

big_value=1.e18
    
# Get EASE grid
EASE_grid = Dataset('/Users/stroeve/Documents/seaice/grid.nc')
lons_ease = np.array(EASE_grid['lon'])
lats_ease = np.array(EASE_grid['lat'])

month=['January','February','March','April','May','June']

for f in model_list[0:1]:
# for f in sorted(models_with_all_variables):
    print('proccessing model ',f)
#get the latitude/longitude using the first file from the siconc to set the latitude/longitude for all variables
    inpath=datapath+ncvarlist[0]+'_'+f   
    lats=get(inpath,'lat')
    lons=get(inpath,'lon')
            
#do this if the lats and lons are 1-d arrays instead of 2-d
    if lats.ndim == 1:
        lons, lats = np.meshgrid(lons,lats)
    lats[(lats > 90) | (lats < -90)] = np.nan
    lons[(lons > 360) | (lons < -360)] = np.nan

#convert lons to -180 to 180
    lons = np.where(lons > 180, lons-360, lons) # This replaces any lon greater than 180 with lons-360 (so -180 to 0)

#find all indices for NH so as not to run the code for all pixels
    irows = np.unique(np.where((np.isfinite(lats) == 1) & (lats >= minlat))[0])
    lonsx=lons[irows,:]
    latsy=lats[irows,:]
    
#Establish Time Units                
    # timeslen = get(inpath,'time')
    # print('Establish time units ',timeslen)
    
#for some reason this num2date doesn't work though it should according to the guide docs
    # time=ncf.variables['time']
    # time_convert=num2date(time[:],time.units, time.calendar)
    #instead use xarray to get the start year 
    f1=xr.open_dataset(inpath) #can just read the first file name for this
        #extract out the year 
    try:
        modeltime=f1.time #this gives an array of all the times
        timestart=modeltime.values[0]
        yearstart=pd.to_datetime(timestart).year     
    except:
        modeltime=f1.indexes['time'].to_datetimeindex() 
        timestart=modeltime.values[0]
        yearstart=pd.to_datetime(timestart).year
    
    output=[]  #declare array that all variables will be stored in to extract from later
# now loop over all the variables per model/ensemble
    for ncvar in ncvarlist: #loop through each variable name
        files2open=datapath+ncvar+'_'+f
        print(files2open)
        ncparts=[]
        ncf=nc.Dataset(files2open)
        
        #testing with xarray
        # test=xr.open_dataset(files2open)
        # test_var=test[ncvar].values
    
        print('processing variable ',ncvar, ' ', ncf[ncvar].shape)
        one_file=ncf[ncvar]  #this is an array of time x lon x lat
        
        #we first just check the shape in order to regrid if neccesary
        one_year=one_file[0,:,:]                     
        if (one_year.shape != lats.shape): #this seems to happen to the incoming solar/outgoing solar radiation
            print('in test for lat/lon ',files2open)
            newlats=get(files2open,'lat')
            newlons=get(files2open,'lon')             
        #do this if the lats are 1-d arrays instead of 2-d
            if newlats.ndim == 1:
                print('needing to convert to 2D array')
                newlons, newlats = np.meshgrid(newlons,newlats)
                newlats[(newlats > 90) | (newlats < -90)] = np.nan
                newlons[(newlons > 360) | (newlons < -360)] = np.nan
                #convert lons to -180 to 180
                newlons = np.where(newlons > 180, newlons-360, newlons) # This replaces any lon greater than 180 with lons-360 (so -180 to 0)
            
            i=0  #now we will have to fill the array if the data has to be regridded 
            t=one_file.shape[0]
            x=lats.shape[0]
            y=lats.shape[1]
            array_to_fill=np.empty((t,x,y))  #we fill an array as a function of time, lon,lat
            
            while i < len(one_file):  #loop over each month/year
                one_year=one_file[i,:,:]  #extract out 2D array for each year                   
                one_year=np.absolute(one_year) #for some reason MPI has negative sisflswutop and siflswdtop
                if (one_year.shape != lats.shape): #now replace the one_file with the regridded data if we need to convert swu and swd
                    newlons.shape
                    newlats.shape
                    new_data=regrid(one_year,newlons,newlats,lons,lats)  #regrid the data for incoming and outgoing solar to the sea ice fields
                    array_to_fill[i]=new_data #now reassign the regridded one_year data back to the one_file array for that variable
                    one_file=array_to_fill   

                i += 1
                       
        # plot(one_year,lats,lons,ncvar)
        ncstart=int(ncf['time'].units.strip('days since ')[0:4])  #note this gives 1850 
                     
        #load each part of the data for this variable
        ncparts.append(one_file[:,irows])
        output.append( np.concatenate(ncparts) ) #append the entire 3-D array of this variable to the main empty list
 
    #end loop for all variables

#separate list into components
    sic,sit,snd,swu,swd=output[0], output[1], output[2], output[3], output[4] 
    
    #convert sic to fraction
    sic=sic/100.
    #set big values to NaN
    sic[sic>big_value]=np.nan
    sit[sit>big_value]=np.nan
    snd[snd>big_value]=np.nan
    swd[swd>big_value]=np.nan
    swu[swu>big_value]=np.nan
    #if snow depth is in cm we have to convert to meters
    if (snd.max() > 10.0 ) :
        print('need to convert snow depth to m ',snd.max())
        snd=snd/100.
    
#first let's test this for the month of April
    imon=3
    sic_one=sic[imon::12,:,:] #start at month index of 3/4 (april/May) and skip every 12th value
    sit_one=sit[imon::12,:,:]
    snd_one=snd[imon::12,:,:]
    swu_one=swu[imon::12,:,:]
    swd_one=swd[imon::12,:,:]
#    swd_one[swd_one>big_value]=np.nan #set big values to Nan
    # swu_one[swu_one>big_value]=0.
    # sit_one[sit_one>big_value]=np.nan
    # snd_one[snd_one>big_value]=np.nan
    # sic_one[sic_one>big_value]=np.nan
    alb_one=swu_one/swd_one 
    times_one=modeltime[imon::12]
    # plot(alb_april[0,:,:],latsy,lonsx,'alb')
    model_name=f[6:35]

#this is getting the under-ice light for all years for that particular month
    par_model,trans_model = get_under_ice_light(sic_one,sit_one,snd_one,alb_one,swd_one,latsy,lonsx,model_name,yearstart)
#    fn_par='/Users/stroeve/Documents/IPCC/CMIP6/Ecolight/'+experiment[1]+'/'+model_name+'_PAR.nc'
    fn_par='/Volumes/LaCie/CMIP6/Ecolight/'+experiment[1]+'/'+model_name+'_PAR_'+month[imon]+'_'+str(yearstart)+'.nc'
    fn_trans='/Volumes/LaCie/CMIP6/Ecolight/'+experiment[1]+'/'+model_name+'_transmittance_'+month[imon]+'_'+str(yearstart)+'.nc'
    ds1=xr.Dataset(data_vars={'par':(['time','x','y'],par_model)},coords =  {'lon':(['x','y'],lons_ease),'lat':(['x','y'],lats_ease),'time':(['time'],times_one)})
    ds2=xr.Dataset(data_vars={'transmittance':(['time','x','y'],trans_model)},coords =  {'lon':(['x','y'],lons_ease),'lat':(['x','y'],lats_ease),'time':(['time'],times_one)})
    ds1.to_netcdf(fn_par)
    ds2.to_netcdf(fn_trans)
  #   output2Netcdf(par_model,times_one,fn_par,'par',model_name)
    #break
    
    
#    break
    # imon=3
    # for imon in range(len(months)):
    #     sic_mon=sic[imon::12,:,:]/100.
    #     sit_mon=sit[imon::12,:,:]
    #     snd_mon=snd[imon::12,:,:]
    #     swu_mon=swu[imon::12,:,:]
    #     swd_mon=swd[imon::12,:,:]
    #     swd_mon[swd_mon>big_value]=np.nan #set big values to Nan
    #     swu_mon[swu_mon>big_value]=0.
    #     sit_mon[sit_mon>big_value]=np.nan
    #     snd_mon[snd_mon>big_value]=np.nan
    #     sic_mon[sic_mon>big_value]=np.nan
    #     alb_mon=swu_mon/swd_mon 
    #     get_under_ice_light(sic_mon,sit_mon,snd_mon,alb_mon)
    #     imon +=1
        




# #   mapping under ice PAR
    
#     cmap1=plt.cm.jet
#     cmap1.set_bad('white')
    
#     cmap2=plt.cm.bwr
#     cmap2.set_bad('white')
    
# #create a figure, determining the size in inches    
#   #   fig=plt.figure(figsize=(6.5,6.5))

# # Create an axis on which to plot. If we only want one map, we say we want a 1 by 1 grid and
# # then we plot to plot in the 1st location. This seems silly now, but it will be helpful later.
#   #   ax = fig.add_subplot(1,1,1)
# # Create the Basemap Object. You need to supply a projection, the bounding box, 
# # and (when applicable) standard lines like a standard parallel or central meridian.
# # For Arctic work, I prefer to use Lambert's Azimuthal Equal-Area Grid (i.e. EASE2 or laea)
# #     lat_0 = 90 # Latitude of Origin
# #     lon_0 = 0 # Central Meridian
# #     bbox = [30,-45,30,135] # order here is lower-left lat, lower-left lon, upper-right lat, upper-right lon    
# #     mp = Basemap(projection='laea',lat_0=lat_0,lon_0=lon_0,llcrnrlat=bbox[0], llcrnrlon=bbox[1], urcrnrlat=bbox[2], urcrnrlon=bbox[3], resolution='c')
# # # The "resolution" here is 'c' for 'coarse'. It only impacts the resolution of boundaries that are built in to basemap (like coastlines and country outlines)

# # # Add other accessories   
# # # Add meridians every 30° and parallels every 10°. The color of '0.4' makes them a medium gray.
# #     mp.drawmeridians( np.arange(0,360,30), color="0.4", linewidth=0.8)
# #     mp.drawparallels( np.arange(0,90,10), color="0.4", linewidth=0.8)
# # # Add coastlines. I sometimes add ESRI shapefiles instead because I find the estuaries too generous.
# #     mp.drawcoastlines()

#     fig= plt.figure(figsize=(6.5,6.5))
#     mp = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l') # choice of lat, lon boundaries
#     mp.fillcontinents(color='white',lake_color='white')
#     mp.drawcoastlines()
#     mp.drawparallels(np.arange(-80.,81.,20.))
#     mp.drawmeridians(np.arange(-180.,181.,20.))
#     mp.drawmapboundary(fill_color='white')   
#     cf1 = mp.contourf(lonsx, latsy, Fsw_TR_NEW, np.arange(0,12,0.1), latlon=True, cmap=plt.cm.viridis, extend='both')
#     # Add an annotation for each sub-plot
#     box = dict(boxstyle='square,pad=0.3',fc='white',ec='white',lw=0.5) # white rectangle with white outline
#     ax.annotate(title[0], xy=(0.02,0.96), xycoords='axes fraction', va='top', ha= 'left', size=8, bbox=box)
#     # # The above adds text plotted on top of the white box located in the upper-left corner of each map
#     hc=plt.colorbar(cf1,);  
#     hc.set_label(r'$\mu mol.m^{-2}.s^{-1}$',size=12);
    
#  #   plt.savefig('/Users/stroeve/Documents/EcoLight/UnderIcePAR_{value2}{value}_SnowDistrib_Icedistrib_DiffITD.png'.format(value=year[iyears],value2=month[2]))
    
#  # We have 2 rows and 2 columns in our plot grid
