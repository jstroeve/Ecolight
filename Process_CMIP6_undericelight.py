#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 11:29:04 2021

@author: stroeve
"""

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date
import netCDF4 as nc
import numpy as np
import csv
import glob
from pathlib import Path
from io import StringIO
import os
import sys
sys.path.insert(0,'/anaconda2/pkgs')
import string

import numpy.ma as ma

import matplotlib.pyplot as plt
from skimage import data, color
from skimage.transform import rescale, resize, downscale_local_mean
from pylab import *
from PIL import Image
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap
import gzip
import matplotlib.cm as cm2
from scipy import interpolate as sp
from scipy.ndimage.interpolation import zoom
from scipy.ndimage.interpolation import rotate
from datetime import datetime
from scipy.stats.mstats import pearsonr
from scipy.stats import linregress
import scipy.stats as ss
from matplotlib.ticker import FuncFormatter
import matplotlib.mlab as mlab
from scipy.stats import norm
from matplotlib import ticker
from scipy import optimize




#this program will run through an emission scenario to find the models to compute under-ice light following Stroeve et al. 2021
#the light processing component is calling a function based off of Gaelle's work

#inialise location of the CMIP6 data and the variable names

filepath='/Volumes/LaCie/CMIP6/'
minlat=45 # in degrees sets the minimum latitude limit

#variables needed for each model name
ncvarlist=['siconc','sithick','sisnthick','siflswutop','siflswdtop']
experiment=['historical','ssp126','ssp585']

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
          lons=ncf.variables['lat'][:].data 
        except:
          lons=ncf.variables['latitude'][:].data
        ncf.close()
        return(lons)
    elif string == 'time':
        ncf=nc.Dataset(inpath)
        timeslen = ncf['time'].shape[0]
        ncf.close()
        return(time)

set_of_siconc_models = set(subset_sic)    
set_of_sithick_models = set(subset_sit)
set_of_swu_models = set(subset_swu)
set_of_swd_models = set(subset_swd)
set_of_snd_models = set(subset_snd)

#this finds a set of all models/ensemble runs that are there for each variable
models_with_all_variables = (set_of_siconc_models & set_of_sithick_models & set_of_snd_models & set_of_swu_models & set_of_swd_models) # Find intersection of sets

big_value=1.e18 

for f in models_with_all_variables:
    
#get the latitude/longitude using the first file from the siconc to set the latitude/longitude
    inpath=datapath+ncvarlist[0]+'_'+f
    lats=get(inpath,'lat')
    lons=get(inpath,'lon')
            
#do this if the lats are 1-d arrays instead of 2-d
    if lats.ndim == 1:
        lons, lats = np.meshgrid(lons,lats)
    lats[(lats > 90) | (lats < -90)] = np.nan
    lons[(lons > 360) | (lons < -360)] = np.nan

#convert lons to -180 to 180
    lons = np.where(lons > 180, lons-360, lons) # This replaces any lon greater than 180 with lons-360 (so -180 to 0)

    irows = np.unique(np.where((np.isfinite(lats) == 1) & (lats >= minlat))[0])
    lonsx=lons[irows,:]
    latsy=lats[irows,:]
       
#Establish Time Units                
    timeslen = get(inpath,'time')
     
    
    
#   file_sic_open.append(ncvarlist[0]+'_'+f)
#    print(file_sic_open)
    
    
# now loop over all the variables per model/ensemble
    output=[]
    for ncvar in ncvarlist: #loop through each variable name
        files2open=datapath+ncvar+'_'+f
        ncparts=[]
        for good_files in files2open: #for each file that matches for all variables, load the array
            ncf=nc.Dataset(good_files)
            ncstart=int(ncf['time'].units.strip('days since ')[0:4])
            times=ncf['time'][:].data
        #load each part of the data for this variable
        ncparts.append(ncf[ncvar][:,irows].data)
    output.append( np.concatenate(ncparts) ) #append the entire 3-D array of this variable to the main empty list

# #separate list into components
# sic,sit,snd,swu,swd=output[0]/100., output[1], output[2], output[3], output[4]
# #first let's test this for the month of April
# sic_april=sic[3::12,:,:] #start at month index of 3 (april) and skip every 12th value
# sit_april=sit[3::12,:,:]
# snd_april=snd[3::12,:,:]
# swu_april=swu[3::12,:,:]
# swd_april=swd[3::12,:,:]
# swd_april[swd_april>big_value]=np.nan #set big values to Nan
# swu_april[swu_april>big_value]=0.
# sit_april[sit_april>big_value]=np.nan
# snd_april[snd_april>big_value]=np.nan
# sic_april[sic_april>big_value]=np.nan
# alb_april=swu_april/swd_april 



# #**********function to compute under-ice light*************
# #************DONE GETTING ALL THE APRIL MODEL VARIABLES
# nyears=sic_april.shape[0]
# year=[]
# big_value=1.e20

# # surface transmission parameter
# i_0_snw  =  0.3;     # for snow 
# i_0_bi   =  0.3;     # for ice  (may change) 
# i_0_pond =  0.58;    # for melt pond
# i_0_ow   =  0.58;    # for open water
# # surface transmission layer thickness (m)
# h_0_snw  =  0.03;     # for snow 
# h_0_bi   =  0.1;    # for bare ice ice  (may change !) 
# alb_pond = 0.27;  # for ponded ice 
# # attenuation coeff ice    
# K_i  = 1;
# # ITD pdf    
# hpdf = [0.0646, 0.1415, 0.173, 0.1272, 0.1114, 0.0824, 0.0665, \
#         0.0541, 0.0429, 0.0347, 0.0287, 0.024,0.0194, 0.016, 0.0136]
# hpdf=np.asarray(hpdf)
# for iyears in range(164,nyears):        ## big loop on years (could also have been on months as Gaelle did)
#     year.append(iyears+1850)
#     sic=sic_april[iyears,:,:]
#     f_bi=reshape(sic,88*320)      ## re-shape sea ice concentration as 1D array for loop
    

#     #### Snow depth ####
#     snod=snd_april[iyears,:,:]
#     h_s=reshape(snod,88*320)   #re-shape to 1D array for loop 
    
#     #### Sea ice thickness ####
#     SIT=sit_april[iyears,:,:]
#     hi=reshape(SIT,88*320)
    
#     # for zz in range (0,88*320):
#     #     if f_bi[zz]>(big_value):
#     #         f_bi[zz]=np.float('nan')
#     #     if h_s[zz]>(big_value):
#     #         h_s[zz]=np.float('nan')
#     #     if hi[zz]>(big_value):
#     #         hi[zz]=np.float('nan')
      

# ## ITD 15 classes ##
#     h_cutoff=3.   # maximum 3 times mean ice thickness       
#     hi15=np.zeros([88*320,15]) ## ice thickness distributed
#     for ii in range (1,16):
#         factor2=(2*ii - 1)/15.    ## factor for distribution
#         hi15[:,ii-1]=hi*factor2     
#         hi15[:,ii-1]=(h_cutoff/2.)*hi15[:,ii-1] ## h_cutoff vary with the SIC distributed

#     #### albedo and incoming solar #### 
#     ALB=alb_april[iyears,:,:]
#     alb=reshape(ALB,88*320)

#     #### Incoming solar irradiance ####
#     FSW0=swd_april[iyears,:,:]
#     Fsw0=reshape(FSW0,88*320)


# # initialisation arrays transmittance, irradiance    
#     T_ow=np.zeros(88*320)
#     T_snow=np.zeros(88*320)
#     Fsw_tr_new=np.zeros([88*320,16])
# # hs array between 0 and 2 times mean snow depth
#     hs=np.arange(0.,2.*nanmax(h_s),0.05)
# # initialisation snow attenuation coeff
#     K_s=np.zeros([len(hs),88*320])

# # initialisation intermediate steps
#     f_att_snow=np.zeros([88*320,len(hs),15])
#     i_s_hom=np.zeros([88*320,16])
#     t_s_hom=np.zeros([88*320,16])

#     #loops on ITD and assign K_s according to wet/dry snow
#     for jj in range(0,15): 
#         for ii in range (0,88*320):
#             for kk in range(0,len(hs)):
# #                if temp[ii] <= 0:   # dry snow
#                 K_s[kk,ii] = 10;
#                 f_att_snow[ii,kk,jj]  = exp(-K_s[kk,ii]*(hs[kk]-h_0_snw))* exp(-K_i*hi15[ii,jj]) ;
#                 # elif hs[kk]<=0.03:  # if snow depth < SSL thickness
#                 if hs[kk]<=0.03:  # if snow depth < SSL thickness
#                     K_s[kk,ii] = 40; 
#                     f_att_snow[ii,kk,jj]  = exp(-K_s[kk,ii]*(hs[kk]))* exp(-K_i*hi15[ii,jj]) ;
#                 # else:
#                 #     K_s[kk,ii]=7;   # wet snow
#                 #     f_att_snow[ii,kk,jj]  = exp(-K_s[kk,ii]*(hs[kk]-h_0_snw))* exp(-K_i*hi15[ii,jj]) ;
        
#         index=where(f_att_snow>1)
#         f_att_snow[index]=1
        
        
#         for ii in range (0,88*320):
#                 # distribution snow depth 
#                 g_hs = zeros([len(hs)]) ; 
#                 g_hs[where(hs  < 2.0 * h_s[ii])]  = f_bi[ii] / ( 2.0 * h_s[ii] );
#                 i_s_hom[ii,jj] = sum(g_hs * 0.01 * f_att_snow[ii,:,jj]); # integrated transmission
#                 for kk in range (0,len(hs)):
#                     # computing transmittance
#                     if K_s[kk,ii]==40:
#                         t_s_hom[ii,jj] = ((1-alb[ii]) * i_s_hom[ii,jj]) / f_bi[ii] ;
#                     else: 
#                         t_s_hom[ii,jj] = ((1-alb[ii]) * i_0_snw * i_s_hom[ii,jj]) / f_bi[ii] ;     
                
#                 T_ow[ii]     = (1 - alb[ii]);  ## transmittance open water
#                 # under-ice irradiance and PAR calculation
#                 Fsw_tr_new[ii,jj] = Fsw0[ii]* ((t_s_hom[ii,jj] * f_bi[ii]*3.51) + (T_ow[ii] * (1-f_bi[ii]))*2.30); 
        
# # sum ITD 15 classes and apply pdf   
    
#     for i in range (0,88*320): 
#         t_s_hom[i,15]=sum(t_s_hom[i,0:15]*hpdf[0:15])
#         Fsw_tr_new[i,15]=sum(Fsw_tr_new[i,0:15]*hpdf[0:15])

#     print('Mean under-ice PAR = ')
#     print(nanmean(Fsw_tr_new[:,15]))
#     sys.stdout.flush()    

# # change shape to 2D array for map plotting
    
#     T_snow=np.zeros([88,320])
#     Fsw_TR_NEW=np.zeros([88,320])  #under-ice PAR
#     c=-1
#     for ii in range (0,88):
#         for jj in range (0,320):
#             c=c+1
#             T_snow[ii,jj]=t_s_hom[c,15]
#             Fsw_TR_NEW[ii,jj]=Fsw_tr_new[c,15]

#     #T_snow=flipud(T_snow)
#     #Fsw_TR_NEW=flipud(Fsw_TR_NEW)

    
#     c=-1
#     H_S=np.zeros([88,320])
#     for ii in range (0,88):
#       for jj in range (0,320):
#         c=c+1
#         H_S[ii,jj]=h_s[c]
    
# # #   use mask on transmittance and under ice PAR using snow depth
# #     Fsw_TR_NEW=np.ma.array(Fsw_TR_NEW,mask=(isnan(flipud(H_S))==True))
# #     T_snow=np.ma.array(T_snow,mask=(isnan(flipud(H_S))==True))


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
