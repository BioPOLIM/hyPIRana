# -*- coding: utf-8 -*-
"""
hyPirana is a program for analysis of PiFM hyPIR spectra acquired using VistaScan
Created on Fri Apr 24 08:06:26 2020

@author: ungersebastian

@contributions: TauDan, sinaravi

modified on Fri Oct 16 by Daniela Taeuber for application to the spectral range of one tuner only
modified by Mohammad Soltaninezhad for rescaling intensities and saving figures
last modified on Fri August 13 2021 by Daniela Taeuber: implemented export of loadings as text-files

hyPirana.py can do:
- read a hyperspectral data set from a Vistascope
- use AreaSelect for selecting a region of interest in the data via a GUI
- use a substrate spectrum for calibrating the hyperspectral data. The substrate text-file should have one line with column titles and consist of two entries (wavelengths and intensities) separated by tabs "\t"
- calculate mean spectra
- run a PCA on the data set
"""

#%% imports & parameters

import hyPIRana as ir

from os.path import join
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime

#%% generate absolute paths of ressource data

path_dir = os.path.dirname(os.path.abspath(__file__))

path_import = r'hyPIRana\resources\single-fibrillar F-Actin'
headerfile  = 'F-actinDy4900007.txt'

path_final = join(path_dir, path_import)

today = datetime.strftime(datetime.now(), "%Y%m%d")

save_path = join(path_dir, 'export')

#%% loads data and plots associated VistaScan parameter images

my_data = ir.IRAFM(path_final, headerfile) 

pos =  [my_file['Caption']=='hyPIRFwd' for my_file in my_data['files']]
hyPIRFwd = np.array(my_data['files'])[pos][0]
data = np.reshape(hyPIRFwd['data'], (hyPIRFwd['data'].shape[0]*hyPIRFwd['data'].shape[1], hyPIRFwd['data'].shape[2]))

my_data.plot_all()

#%%

my_wl  = my_data['wavelength']
my_sum = np.sum(data, axis = 1)
coord = np.arange(len(my_sum))
zeros = np.zeros(len(my_sum))
data = data[my_sum != 0]
coord = coord[my_sum != 0]
my_sum = my_sum[my_sum != 0]

#Calibration using CaF:
    

path_caf = join(path_dir, path_import,'CalibrationFile_single fribrillar.txt')
caf_file = pd.read_csv(path_caf, delimiter = '\t')
caf_spc = np.array(caf_file)[:,1]

caf_spc = caf_spc[830:931]

spc_norm = np.array([(spc/caf_spc)/s for spc, s in zip(data, my_sum)])

#%%
# PlotIt
mean_spc = np.mean(spc_norm, axis = 0)
std_spc = np.std(spc_norm, axis = 0)
my_fig = plt.figure()
ax = plt.subplot(111)
ax.fill_between(x = my_data['wavelength'], y1 = mean_spc+std_spc, y2 = mean_spc-std_spc, alpha = 0.6)
plt.gca().invert_xaxis() #inverts values of x-axis
ax.plot(my_data['wavelength'], mean_spc)
ax.set_xlabel('wavenumber (1/cm)')
ax.set_ylabel('intensity (normalized)')
#ax.set_yticklabels([])
plt.title('mean spectrum')
my_fig.savefig( 'mean spectrum.png' )
my_fig.tight_layout()
my_spc = my_data.return_spc()
my_wl  = my_data['wavelength']

np.savetxt(r'export\MEAN.txt', mean_spc)

#%% 
from sklearn.decomposition import PCA

ncomp = 2

model = PCA(n_components=ncomp)

transformed_data = model.fit(spc_norm-mean_spc).transform(spc_norm-mean_spc).T
loadings = model.components_

print("loadingtype",type(model))
print("mysum",my_sum.shape)
print("data",data.shape)
print("spc_norm",spc_norm.shape)
print("meanspc",mean_spc.shape)
print("std",std_spc.shape)
print("transform",transformed_data.shape)
print("loading",loadings.shape)

my_fig = plt.figure()
ax = plt.subplot(111)
plt.gca().invert_xaxis() #inverts values of x-axis
for icomp in range(ncomp):
    ax.plot(my_data['wavelength'], loadings[icomp], label='PC'+str(icomp+1) )
ax.set_xlabel('wavenumber (1/cm)')
ax.set_ylabel('intensity (normalized)')
#ax.set_yticklabels([])
ax.legend()
plt.title('PCA-Loadings')

my_fig.savefig( r'export\PCA-Loadings.png' )

my_fig = plt.figure()
ax = plt.subplot(111)
ax.plot(transformed_data[0], transformed_data[1], '.')
ax.set_xlim(np.quantile(transformed_data[0], 0.05),np.quantile(transformed_data[0], 0.95))
ax.set_ylim(np.quantile(transformed_data[1], 0.05),np.quantile(transformed_data[1], 0.95))
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
plt.title('scatterplot')

my_fig.savefig( r'export\scatterplot.png' )

my_fig.tight_layout()
vmax = 100
maps = [zeros.copy() for icomp in range(ncomp)]
for icomp in range(ncomp):
    maps[icomp][coord] = transformed_data[icomp]
    maps[icomp] = np.reshape(maps[icomp], (my_data['xPixel'], my_data['yPixel']) )
    
    Km = maps[icomp][1:32,:]
    
    my_fig= plt.figure()
    ax = plt.subplot(111)
    plt.colorbar( ax.imshow(maps[icomp], cmap = 'coolwarm', extent = my_data.extent(), vmin = -vmax, vmax = vmax) )
    
    ax.set_xlabel('x scan ['+my_data['XPhysUnit']+']')
    ax.set_ylabel('y scan ['+my_data['YPhysUnit']+']')
    plt.title('factors PC'+str(icomp+1))
    
    my_fig.savefig( r'export\factors PC.png' )

    my_fig.tight_layout()

plt.show()



