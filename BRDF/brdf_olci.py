#!/usr/bin/env python
# coding: utf-8
"""
Created on Wed Jan 23 18:21:01 2019
Modified by Javier on July 15 2019

@author: mario
mario.benincasa@artov.isac.cnr.it
mbeninca@gmail.com
"""

"""
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

# brdf_mario.py
# Bidirectional Reflectance Distribution Function
# computation of the bdrf correction for OLCI sensor

# Wang2006.nc content
# netcdf Wang2006 {
# dimensions:
#     Coefficient = 4 ;
#     Sigma = 5 ;
#     Wavelength = 6 ;
# variables:
#     float Coefficient(Coefficient) ;
#     float Sigma(Sigma) ;
#     float Wavelength(Wavelength) ;
#     float Wang2006(Sigma, Coefficient, Wavelength) ;
# data:

#  Wavelength = 412, 443, 490, 510, 555, 670 ;
# }

# morel_fq.nc content 
# netcdf morel_fq {
# dimensions:
#     n_phi = 13 ;
#     n_senz = 17 ;
#     n_chl = 6 ;
#     n_solz = 6 ;
#     n_wave = 7 ;
# variables:
#     float phi(n_phi) ;
#         phi:units = "degrees" ;
#     float senz(n_senz) ;
#         senz:units = "degrees" ;
#     float chl(n_chl) ;
#         chl:units = "mg m^-3" ;
#     float solz(n_solz) ;
#         solz:units = "degrees" ;
#     float wave(n_wave) ;
#         wave:units = "um" ;
#     float foq(n_wave, n_solz, n_chl, n_senz, n_phi) ;

# // global attributes:
#         :description = "F/Q" ;
#         :history = "Created Tue Jul 21 13:43:50 2015" ;
#         :source = "cre_morel_ncfile.py" ;
#         :author = "Rick Healy richard.healy@nasa.gov" ;
# data:

#  wave = 412.5, 442.5, 490, 510, 560, 620, 660 ;
# }


"""
import with something like:
  sys.path.insert(0,'/path/to/dir/containing/us')
  import brdf_mario
"""

import os #, multiprocessing
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator, interp1d #, griddata
from multiprocessing import cpu_count

#from joblib import Parallel, delayed

#number of cores
n_cores=cpu_count()
#where we reside, where to look for LUT files Wang2006.nc and morel_fq.nc
sdir=os.path.dirname(os.path.realpath(__file__))

# constants
f0 = 0.9795218   # what is f0?
riw = 1.334      # refractive_index_water
sigmac= 0.0731   # constant to multiply square root of windspeed (?)
powers= np.array([1., 2., 3., 4.], dtype=np.float32)

# build interpolators
# Wang 2006
wang2006_file=Dataset(os.path.join(sdir,'Wang2006.nc'),'r')
wang2006_Coefficient=wang2006_file.variables['Coefficient'][:]
wang2006_Sigma=wang2006_file.variables['Sigma'][:]
wang2006_Wavelength=wang2006_file.variables['Wavelength'][:]
wang2006_nbands=wang2006_file.dimensions['Wavelength'].size
wang2006_table=wang2006_file.variables['Wang2006'][:]
wang2006_file.close()
# wang2006_table=wang2006_table(sigma, coefficent, wavelength)
#swap first and last dimensions, to obtain wang2006_table(wavelength,sigma,coefficent)
wang2006_table=np.moveaxis(wang2006_table,2,0)

# f over Q Morel Gentili
foq_file=Dataset(os.path.join(sdir,'morel_fq.nc'),'r')
foq_rel_azimuths=foq_file.variables['phi'][:]
foq_sensor_zeniths=foq_file.variables['senz'][:]
foq_lchl_levels=np.log10(foq_file.variables['chl'][:])
foq_solar_zeniths=foq_file.variables['solz'][:]
foq_bands=foq_file.variables['wave'][:]
foq_nbands=foq_file.dimensions['n_wave'].size
foq_table=foq_file.variables['foq'][:]
foq_file.close()
foqint_morel=dict()
for windex,w in np.ndenumerate(foq_bands):
    widx=windex[0]  #windex is a tuple of 1 element (because array is 1 column)
    foqint_morel[widx]=RegularGridInterpolator(
                        (foq_solar_zeniths, foq_lchl_levels, foq_sensor_zeniths, foq_rel_azimuths),
                        foq_table[widx,:,:,:,:],method='linear', 
                        bounds_error=False, fill_value=None
                       )
# /interpolators

def fresnel_sensor(angle):
    """
    Compute the effects of the air-sea transmittance for sensor view zenith angle
    Angle (in degrees) can be single float or array.
    Returns a scalar for every sensor zenith angle.
    """
    mu=np.cos(np.deg2rad(angle))
    mu2=mu**2
    sq=np.sqrt(riw**2 - 1 + mu2)
    musq=mu*sq
    r2=((mu-sq)/(mu+sq))**2
    q1=(1-mu2-musq)/(1-mu2+musq)
    fres=r2*(q1**2+1)/2
    tf=1-fres
    return f0/tf

def fresnel_solar(windspeed,angle):

    """
    Compute the effects of the air-sea transmittance for solar path 
    for the windspeed module and the solar zenith angle.
    Windspeed and angle can be single scalars or (same dim) arrays
    Returns the xxx factor on 7 bands (instead of six) to directly compare with foq
    """
    from joblib import Parallel, delayed
    sigma=sigmac*np.sqrt(windspeed)
    s=np.asarray(angle,dtype=np.float32)    # to handle the same way single float and array
    s[s>80]=80
    s=np.log(np.cos(np.deg2rad(s)))    
    if s.shape:
        log_cosine_vector=np.ndarray((s.shape[0],wang2006_Coefficient.size),dtype=np.float32)
        for i,p in enumerate(wang2006_Coefficient):
            log_cosine_vector[:,i]=s**powers[i]
        brdf_values=np.zeros((wang2006_nbands,wang2006_Sigma.size,s.shape[0]),dtype=np.float32)
        brdf_values=1+np.tensordot(wang2006_table,log_cosine_vector.transpose(), axes=1)
        result=np.zeros((s.shape[0],wang2006_nbands+1),dtype=np.float32)
        def myf(p):
            brdf_interf=interp1d(wang2006_Sigma,brdf_values[:,:,p],fill_value='extrapolate')
            #extrapolate necessario!
            a=brdf_interf(sigma[p])
            return a
        #(n_jobs=n_cores, prefer="threads")
        #investigation ... 20190330
        result[:,0:wang2006_nbands]= Parallel(n_jobs=n_cores)(delayed(myf)(i) for i in range(s.shape[0]))
#        for i in range(s.shape[0]):
#            result[i,0:wang2006_nbands]=myf(i)
        result[:,wang2006_nbands]=result[:,wang2006_nbands-1]
        """
        for p in range(s.shape[0]):
            brdf_interpolator=interp1d(wang2006_Sigma,brdf_values[:,:,p])
            result[p,0:wang2006_nbands]=brdf_interpolator(sigma[p])
            result[p,wang2006_nbands]=result[p,wang2006_nbands-1]
            #or, slower (almost twice the time)
            #result[p,0:wang2006_nbands]=griddata(wang2006_Sigma,brdf_values[:,:,p].transpose(),sigma[p])
            #result[p,wang2006_nbands]=result[p,wang2006_nbands-1]
        """
    else:
        log_cosine_vector = s**powers
        brdf_values=np.zeros((wang2006_nbands,wang2006_Sigma.size),dtype=np.float32)
        brdf_values=1+np.dot(wang2006_table,log_cosine_vector)        
        result=np.zeros((wang2006_nbands+1),dtype=np.float32)
        brdf_interpolator=interp1d(wang2006_Sigma,brdf_values)
        result[0:wang2006_nbands]=brdf_interpolator(sigma)
        result[wang2006_nbands]=result[wang2006_nbands-1]
        #or 
        #result[0:wang2006_nbands]=griddata(wang2006_Sigma,brdf_values.transpose(),sigma)
        #result[wang2006_nbands]=result[wang2006_nbands-1]
    return result

def brdf(ws0, ws1, chl, sza, saa, vza, vaa):
    """
    Compute the BRDF correction for OLCI, according to Morel 2002 and Wang 2006
    arguments:
        wind speed components (ws0,ws1)
        chlorophil concentrations (chl) (NOT log)
        solar zenith angle (deg) (sza)
        solar azimuth angle
        view (sensor) zenith angle (vza)
        view (sensor) azimuth angle (vaa)
    can be scalar or 1-D arrays (same dim, of course)
    Returns a vector of 7 (1 per band) brdf correction coefficients per pixel
    """
    #number of pixels we were called on
    if np.asarray(ws0).shape:
        npixels=ws0.shape[0] 
    else:
        npixels=1
    ws = np.sqrt(ws0**2 + ws1**2)
    #chl log10, limit lchl to allowed max in interpolation table
    lchl=np.asarray(np.log10(chl))
    lchl[lchl<foq_lchl_levels.min()]=foq_lchl_levels.min()
    lchl[lchl>foq_lchl_levels.max()]=foq_lchl_levels.max()
    
    specBRDF=np.ones((npixels,foq_nbands), dtype=np.float32)
    for widx in range(foq_nbands):
        specBRDF[:,widx]*=fresnel_sensor(vza)
    specBRDF*=fresnel_solar(ws,sza)
    #relative azimuth angle
    raa=np.asarray(np.abs(saa-vaa-180)) #### why this -180???
    raa[raa>180] -= 360
    raa=np.abs(raa)
    #wsza = np.rad2deg(np.arcsin(np.sin(np.deg2rad(sza))/riw))  # refractive_index_water
    #wsza=np.asarray(wsza)
    #wsza[wsza>foq_solar_zeniths.max()]=foq_solar_zeniths.max()
    vza=np.asarray(vza)
    vza[vza<foq_sensor_zeniths.min()]=foq_sensor_zeniths.min()
    vza[vza>foq_sensor_zeniths.max()]=foq_sensor_zeniths.max()
    wvza = np.rad2deg(np.arcsin(np.sin(np.deg2rad(vza))/riw))
    wvza=np.asarray(wvza)
    sza=np.asarray(sza)
    sza[sza<foq_solar_zeniths.min()]=foq_solar_zeniths.min()
    sza[sza>foq_solar_zeniths.max()]=foq_solar_zeniths.max()
    if npixels==1:
        arg0=np.array([0,lchl,0,0],dtype=np.float32)
        #arg =np.array([wsza,lchl,vza,raa],dtype=np.float32) 
        arg =np.array([sza,lchl,wvza,raa],dtype=np.float32) 
    else:
        zerov=np.zeros(npixels,dtype=np.float32)
        arg0=np.stack((zerov,lchl,zerov,zerov),axis=1)
        #arg=np.stack((wsza,lchl,vza,raa),axis=1)
        arg=np.stack((sza,lchl,wvza,raa),axis=1)
    for widx in range(foq_nbands):
        specBRDF[:,widx]*=foqint_morel[widx](arg0)
        specBRDF[:,widx]/=foqint_morel[widx](arg)
    return specBRDF

