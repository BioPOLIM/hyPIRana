
"""" monIRana is a program for analysis of a set of spectral data provided in a tabulated textfile
Created on 4 July 2021
@author: Mohammad Soltaninezhad
@contributor: Maryam Ali (code editing, calibration + scatter plot clustering)
@supervisor: Daniela Täuber

monIRana can do:
- Calibration
- calculate mean spectra of the complete data set
- run a PCA on the data set and provide plots of results
- clustered presentation of Amide rich & glycan rich analyzed spectra
"""


from os.path import join
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from sklearn import preprocessing

path1 =r"MonIRana\resourses\control(13)Amide(32)Sugar(19)_MA.txt"
pathCali1 =r"MonIRana\resourses\CaF2 for BacVan30.txt"
pathCali2 =r"MonIRana\resourses\CaF2 for Control30.txt"
df1 = pd.read_csv(path1,skiprows=[0], delimiter="\t" , header = None)
dfCali1 = pd.read_csv(pathCali1,skiprows=[0], delimiter="\t" , header = None)
dfCali2 = pd.read_csv(pathCali2,skiprows=[0], delimiter="\t" , header = None)
df1.to_numpy()
dfCali1.to_numpy()
dfCali2.to_numpy()

#%% Read Data

D = df1.values.T
dataCali1 = dfCali1.values.T
dataCali2 = dfCali2.values.T

#Classifying Amide rich spectra (Data1) & Glycan rich spectra (Data2)
Data1 = D[15:65]
Data2 = D[1:14]

print('data1:', Data1.shape)
print('data2:', Data2.shape)

#%% Calibration
# This runs if the improted data is not calibrated

#Data1 = Data1 / dataCali1
#Data2 = Data2 / dataCali2


#%% Normalization (L1)

my_sum1 = np.sum(Data1 , axis=1 ) #sum of every axis
my_sum1 = my_sum1[my_sum1 != 0]

my_sum2 = np.sum(Data2 , axis=1 ) #sum of every axis
my_sum2 = my_sum2[my_sum2 != 0]

spc_norm1 = np.array( [spc / s for spc, s in zip( Data1, my_sum1 )] )
spc_norm2 = np.array( [spc / s for spc, s in zip( Data2, my_sum2 )] )


spc_norm = np.concatenate((spc_norm1, spc_norm2), axis = 0)

mean_spc = np.mean(spc_norm, axis = 0)
mean_spc1 = np.mean(spc_norm1, axis = 0)
mean_spc2 = np.mean(spc_norm2, axis = 0)

my_wl = D[0]  # wavelength

#%% plotting mean spectra

mean_fig = plt.figure()
plt.plot(my_wl, mean_spc1.T, c = 'limegreen', label = 'treated')
plt.plot(my_wl, mean_spc2.T, c = 'magenta', label = 'control')
plt.legend()
plt.show()

#%% PCA

ncomp=2

model = PCA( ncomp )
print('spc_norm', spc_norm.shape)
# For running PCA on treated data only use spc_norm1 & mean_spc1 
transformed_data = model.fit( spc_norm - mean_spc ).transform( spc_norm - mean_spc ).T
loadings = model.components_

#%% pca_var = PCA( ncomp )

pca_var = PCA( ncomp )
data_var = pca_var.fit( spc_norm - mean_spc ).transform( spc_norm - mean_spc ).T
pca_var.fit( data_var )
pca_data_var = pca_var.transform( data_var )
per_var = np.round( pca_var.explained_variance_ratio_ * 100, decimals=1 )
labels = ['pc' + str( x ) for x in range( 1, len( per_var ) + 1 )]
plt.bar( x=range( 1, len( per_var ) + 1 ), height=per_var, tick_label=labels )
plt.ylabel( "percantage of explained variance" )
plt.xlabel( "Principle Components" )
plt.title( "Bar plot" )
plt.show()
data_acc = []
i_old = 0
for i in per_var:
    i_old = i_old + i
    data_acc.append( i_old )
plt.bar( x=range( 1, len( data_acc ) + 1 ), height=data_acc, tick_label=labels )
plt.ylabel( "accumulate variance" )
plt.xlabel( "Principle Components" )
plt.title( "Bar plot" )
plt.show()

#%% Plotting

my_fig = plt.figure()
ax = plt.subplot( 111 )
plt.gca().invert_xaxis()  # inverts values of x-axis
x_min = np.amin(transformed_data[0])
y_min = np.amin(transformed_data[1])
x_max = np.amax(transformed_data[0])
y_max = np.amax(transformed_data[1])
print("PC1_min  :", x_min)
print("PC1_max :", x_max)
print("PC2_max :", y_max)
print("PC2_min :", y_min)

# PC Loadings
for icomp in range( ncomp ):
    ax.plot(my_wl , loadings[icomp], label='PC' + str( icomp + 1 ) )
ax.set_xlabel( 'wavenumber ' )
ax.set_ylabel( 'intensity (normalized)' )
ax.set_yticklabels( [] )
ax.legend()
plt.title( 'PCA-Loadings' )
#my_fig.savefig( 'PCA-Loadings.png' )
#np.savetxt('LoadingsPC1.txt', loadings[0], delimiter='\t')
#np.savetxt('LoadingsPC2.txt', loadings[1], delimiter='\t')

#scatter plot
my_fig = plt.figure()
ax = plt.subplot( 111 )
ax.plot( transformed_data[0], transformed_data[1], '.' )
ax.set_xlabel( 'PC1 scores x1000')
ax.set_ylabel( 'PC2 scores x1000')
my_fig.tight_layout()
plt.show()

# clustered scatter plot
my_fig2 = plt.figure()
ax2 = plt.subplot( 111 )
ax2.plot((transformed_data[0, 0:31])*1000, (transformed_data[1, 0:31])*1000, '.', c = 'limegreen', label = 'Amide', markersize = 45)
ax2.plot((transformed_data[0, 32:50])*1000, (transformed_data[1, 32:50])*1000, '.', c = 'magenta', label = 'Glycan', markersize = 45)
ax2.plot((transformed_data[0, 51:65])*1000, (transformed_data[1, 51:65])*1000, '.', c = 'b', label = 'Control', markersize = 40)
ax2.set_xlabel( 'PC1 scores x1000')
ax2.set_ylabel( 'PC2 scores x1000')
my_fig2.tight_layout()
plt.show()
