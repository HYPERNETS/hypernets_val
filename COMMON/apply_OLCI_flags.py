#!/usr/bin/env python
# coding: utf-8
"""
Created on Wed Jan 23 18:21:01 2019
Modified by Javier on July 15 2019

@author: Marco Bracaglia
Modified by Javier A. Concha July 19 2019
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
import numpy as np
import numpy.ma as ma
#%% Flags values
'''
INVALID = Invalid flag: instrument data missing or invalid
WATER = Water (marine) with clear sky conditions, i.e. no clouds
CLOUD = Cloudy pixel
CLOUD_AMIBUOUS = Possibly a cloudy pixel, the flag removes semi-transparent clouds and other ambiguous cloud signatures
CLOUD_MARGIN = Cloud edge pixel, the flag provides an a-priori margin on the ‘CLOUD or CLOUD_AMBIGUOUS’ flag of 2 pixels at RR and 4 pixels at FR
SNOW_ICE = Possible sea-ice or snow contamination
INLAND_WATER = Fresh inland waters flag (from L1B); these pixels will also be flagged as LAND rather than WATER.
TIDAL = Pixel is in a tidal zone (from L1B)
COSMETIC = Cosmetic flag (from L1B)
SUSPECT = Suspect flag (from L1B)
HISOLZEN = High solar zenith: SZA > 70°
SATURATED = Saturation flag: saturated within any band from 400 to 754 nm or in bands 779, 865, 885 and 1020 nm
MEGLINT = Flag for pixels corrected for sun glint
HIGHGLINT = Flag for when the sun glint correction is not reliable
WHITECAPS = Flag for when the sea surface is rough enough for there to be whitecaps, which cause a brightening of the water-leaving reflectance
ADJAC = reserved for future use for an adjacency correction, so always set to false
WV_FAIL = IWV retrieval algorithm failed
AC_FAIL = BAC atmospheric correction is suspect
OC4ME_FAIL = OC4Me algorithm failed 
OCNN_FAIL = NN algorithm failed
KDM_ FAIL = KD490 algorithm failed

BPAC_ON = Bright Pixel Correction converged and a NIR signal was determined
WHITE_SCATT = "White" scatterer within the water e.g. coccoliths
LOWRW = Water-leaving reflectance at 560 nm is less than a defined threshold or HIINLD_F raised (flag for low pressure water i.e., high altitude inland waters)
HIGHRW = High water-leaving reflectance at 560 nm or the TSM retrieved as part of the BPAC is above a threshold
ANNOT = Annotation flags for the quality of the atmospheric correction, including:
      ANNOT_ANGSTROM (Ångström exponent cannot be computed);
      ANNOT_AERO_B (blue aerosols);
      ANNOT_ABSO_D (desert dust absorbing aerosols);
      ANNOT_ACLIM (aerosol model does not match aerosol climatology);
      ANNOT_ABSOA (absorbing aerosols);
      ANNOT_MIXR1 (aerosol mixing ratio is equal to 0 or 1);
      ANNOT_DROUT (minimum absolute value of the reflectance error at 510 nm is greater than a defined threshold);
      ANNOT_TAU06 (aerosol optical thickness is greater than a defined threshold)
RWNEG_O01 to RWNEG_O21 = Provides a "negative water-leaving reflectance" flag for each band’s water-leaving reflectance: the value below which pixels are flagged varies according to the band, with the threshold stored within a Look-Up Table
'''
INVALID	          =  1                     
WATER             =  2                   
LAND              =  4                   
CLOUD             =  8                   
CLOUD_AMBIGUOUS   =  8388608             
CLOUD_MARGIN      =  16777216            
SNOW_ICE          =  16                  
INLAND_WATER      =  32                  
TIDAL             =  64                  
COSMETIC          =  128                 
SUSPECT           =  256                 
HISOLZEN          =  512                 
SATURATED         =  1024                
MEGLINT           =  2048                
HIGHGLINT         =  4096                
WHITECAPS         =  8192                
ADJAC             =  16384               
WV_FAIL           =  32768               
PAR_FAIL          =  65536               
AC_FAIL           =  131072              
OC4ME_FAIL        =  262144              
OCNN_FAIL         =  524288              
KDM_FAIL          =  2097152             
BPAC_ON           =  33554432            
WHITE_SCATT       =  67108864            
LOWRW             =  134217728           
HIGHRW            =  268435456           
ANNOT_ANGSTROM    =  4294967296          
ANNOT_AERO_B      =  8589934592          
ANNOT_ABSO_D      =  17179869184         
ANNOT_ACLIM       =  34359738368         
ANNOT_ABSOA       =  68719476736         
ANNOT_MIXR1       =  137438953472        
ANNOT_DROUT       =  274877906944        
ANNOT_TAU06       =  549755813888        
RWNEG_O1          =  1099511627776       
RWNEG_O2          =  2199023255552       
RWNEG_O3          =  4398046511104       
RWNEG_O4          =  8796093022208       
RWNEG_O5          =  17592186044416      
RWNEG_O6          =  35184372088832      
RWNEG_O7          =  70368744177664      
RWNEG_O8          =  140737488355328     
RWNEG_O9          =  281474976710656     
RWNEG_O10         =  562949953421312     
RWNEG_O11         =  1125899906842624    
RWNEG_O12         =  2251799813685248    
RWNEG_O16         =  4503599627370496    
RWNEG_O17         =  9007199254740992    
RWNEG_O18         =  18014398509481984   
RWNEG_O21         =  36028797018963968     

#%%
def create_mask(flag_file):
    flag_file = np.int64(flag_file)
    flag_mask=ma.zeros(flag_file.shape)
    flag_mask[ma.where((flag_file & INVALID)\
		| (flag_file & LAND)\
		| (flag_file & CLOUD)\
		| (flag_file & CLOUD_AMBIGUOUS)\
		| (flag_file & CLOUD_MARGIN)\
		| (flag_file & SNOW_ICE)\
		| (flag_file & SUSPECT)\
#		| (flag_file & HISOLZEN)\
		| (flag_file & SATURATED)\
		| (flag_file & HIGHGLINT)\
		| (flag_file & WHITECAPS)\
		| (flag_file & AC_FAIL)\
		# | (flag_file & RWNEG_O2)\
		# | (flag_file & RWNEG_O3)\
		# | (flag_file & RWNEG_O4)\
		# | (flag_file & RWNEG_O5)\
		# | (flag_file & RWNEG_O6)\
		# | (flag_file & RWNEG_O7)\
		# | (flag_file & RWNEG_O8)\
		)]=1
    flag_mask = np.int64(flag_mask)
    return flag_mask 
#%%	
def main():
	print('This is the main() of the apply_flags_OLCI.py')
#%%
if __name__ == '__main__':
   main()  	