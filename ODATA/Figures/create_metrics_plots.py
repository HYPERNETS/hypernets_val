#!/usr/bin/env python3
# coding: utf-8
"""
Created on Thu Dec 10 21:31:16 2020

@author: javier.concha
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
import os
import pandas as pd
import sys
from matplotlib import pyplot as plt
import numpy as np
fs = 14
plt.rc('xtick',labelsize=fs)
plt.rc('ytick',labelsize=fs)

# User Defined Functions
code_home = os.path.abspath('../../')
sys.path.append(code_home)
import COMMON.common_functions as cfs
#%%
idir = '/Users/javier.concha/Desktop/Javier/2019_Roma/CNR_Research/HYPERNETS_D7p2/MDB_py/ODATA/Figures/'

S3A_7p00 = pd.read_csv(os.path.join(idir,'metrics_S3A_WFR_07p00.csv'))
S3B_7p00 = pd.read_csv(os.path.join(idir,'metrics_S3B_WFR_07p00.csv'))
S3A_6p13 = pd.read_csv(os.path.join(idir,'metrics_S3A_WFR_06p13.csv'))
S3B_6p13 = pd.read_csv(os.path.join(idir,'metrics_S3B_WFR_06p13.csv'))

C =[]
for wl in S3A_6p13['wl']:
    C.append(np.pi/cfs.get_F0(wl))

C = np.array(C)

# RMSE
plt.figure(figsize=(8,5))
plt.plot(S3A_6p13['wl'][S3A_6p13['wl']<900],C[S3A_6p13['wl']<900]*S3A_6p13['rmse'][S3A_6p13['wl']<900],'c.--')
plt.plot(S3A_7p00['wl'][S3A_7p00['wl']<900],C[S3B_7p00['wl']<900]*S3A_7p00['rmse'][S3A_7p00['wl']<900],'co-')
plt.plot(S3B_6p13['wl'][S3B_6p13['wl']<900],C[S3A_6p13['wl']<900]*S3B_6p13['rmse'][S3B_6p13['wl']<900],'m.--')
plt.plot(S3B_7p00['wl'][S3B_7p00['wl']<900],C[S3B_7p00['wl']<900]*S3B_7p00['rmse'][S3B_7p00['wl']<900],'mo-')
plt.xlabel('Wavelength (nm)',fontsize=fs)
plt.ylabel('RMSD (-)',fontsize=fs)

plt.legend(['S3A WFR IPF-OL-2 version 06.13','S3A WFR IPF-OL-2 version 07.00',
            'S3B WFR IPF-OL-2 version 06.13','S3B WFR IPF-OL-2 version 07.00'])
plt.show()

# RMSE
plt.figure(figsize=(8,5))
plt.plot(S3A_6p13['wl'][S3A_6p13['wl']<800],S3A_6p13['rmse'][S3A_6p13['wl']<800],'c.--')
plt.plot(S3A_7p00['wl'][S3A_7p00['wl']<800],S3A_7p00['rmse'][S3A_7p00['wl']<800],'co-')
plt.plot(S3B_6p13['wl'][S3B_6p13['wl']<800],S3B_6p13['rmse'][S3B_6p13['wl']<800],'m.--')
plt.plot(S3B_7p00['wl'][S3B_7p00['wl']<800],S3B_7p00['rmse'][S3B_7p00['wl']<800],'mo-')
plt.xlabel('Wavelength (nm)',fontsize=fs)
plt.ylabel('RMSD',fontsize=fs)

plt.legend(['S3A WFR IPF-OL-2 version 06.13','S3A WFR IPF-OL-2 version 07.00',
            'S3B WFR IPF-OL-2 version 06.13','S3B WFR IPF-OL-2 version 07.00'])
plt.show()

# MAPD
plt.figure(figsize=(8,5))
plt.plot(S3A_6p13['wl'][S3A_6p13['wl']<800],S3A_6p13['MAPD'][S3A_6p13['wl']<800],'c.--')
plt.plot(S3A_7p00['wl'][S3A_7p00['wl']<800],S3A_7p00['MAPD'][S3A_7p00['wl']<800],'co-')
plt.plot(S3B_6p13['wl'][S3B_6p13['wl']<800],S3B_6p13['MAPD'][S3B_6p13['wl']<800],'m.--')
plt.plot(S3B_7p00['wl'][S3B_7p00['wl']<800],S3B_7p00['MAPD'][S3B_7p00['wl']<800],'mo-')
plt.xlabel('Wavelength (nm)',fontsize=fs)
plt.ylabel('MAPD (%)',fontsize=fs)

plt.legend(['S3A WFR IPF-OL-2 version 06.13','S3A WFR IPF-OL-2 version 07.00',
            'S3B WFR IPF-OL-2 version 06.13','S3B WFR IPF-OL-2 version 07.00'])
plt.show()

# MPD
plt.figure(figsize=(8,5))
plt.plot(S3A_6p13['wl'][S3A_6p13['wl']<800],S3A_6p13['MPD'][S3A_6p13['wl']<800],'c.--')
plt.plot(S3A_7p00['wl'][S3A_7p00['wl']<800],S3A_7p00['MPD'][S3A_7p00['wl']<800],'co-')
plt.plot(S3B_6p13['wl'][S3B_6p13['wl']<800],S3B_6p13['MPD'][S3B_6p13['wl']<800],'m.--')
plt.plot(S3B_7p00['wl'][S3B_7p00['wl']<800],S3B_7p00['MPD'][S3B_7p00['wl']<800],'mo-')
plt.xlabel('Wavelength (nm)',fontsize=fs)
plt.ylabel('MPD (%)',fontsize=fs)

plt.legend(['S3A WFR IPF-OL-2 version 06.13','S3A WFR IPF-OL-2 version 07.00',
            'S3B WFR IPF-OL-2 version 06.13','S3B WFR IPF-OL-2 version 07.00'])
plt.show()

# r^2
plt.figure(figsize=(8,5))
plt.plot(S3A_6p13['wl'][S3A_6p13['wl']<800],S3A_6p13['r^2'][S3A_6p13['wl']<800],'c.--')
plt.plot(S3A_7p00['wl'][S3A_7p00['wl']<800],S3A_7p00['r^2'][S3A_7p00['wl']<800],'co-')
plt.plot(S3B_6p13['wl'][S3B_6p13['wl']<800],S3B_6p13['r^2'][S3B_6p13['wl']<800],'m.--')
plt.plot(S3B_7p00['wl'][S3B_7p00['wl']<800],S3B_7p00['r^2'][S3B_7p00['wl']<800],'mo-')
plt.xlabel('Wavelength (nm)',fontsize=fs)
plt.ylabel(r'$r^2$',fontsize=fs)

plt.legend(['S3A WFR IPF-OL-2 version 06.13','S3A WFR IPF-OL-2 version 07.00',
            'S3B WFR IPF-OL-2 version 06.13','S3B WFR IPF-OL-2 version 07.00'])
plt.show()