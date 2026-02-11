
"""
hyPirana is a program for  applying Principal Component Analysis (PCA) of PiFM hyPIR spectra acquired using VistaScan
Created on Fri Apr 24 08:06:26 2020

@author: ungersebastian
New features added by Maryam Ali: Ability to read single\multiple datasets, calibration, 
analyze with combined PCA for multiple datasets, and clustered scatter plots
"""
#%% imports & parameters

from os.path import join
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from sklearn.preprocessing import Normalizer
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import non_negative_factorization

from IRAFM import IRAFM as ir


path_import = r'hyPIRana\resources'

headerfile = ['BacVan30_0011.txt','BacVan60_0013.txt'] 
Calib_file = [(pd.read_csv('avCaF2_BacVan30.txt', delimiter = '\t')), 
              (pd.read_csv('avCaF2_BacVan60.txt', delimiter = '\t'))]

path_final = join(path_import, path_import)


#%% Choosing no. of data sets to be read

n_set = 2


# loads data and plots associated VistaScan parameter images
Data = []

for n in range(n_set):
    my_data = ir(path_final, headerfile[n])
    my_data.plot_all()
    Data.append(my_data)



#%% Reading data & calibration files
# Normalization by L1 vector normalization
# plotting the mean spectrum after smoothing by savitzky-golay filter (2-5-5) 
# can be changed by 'windows_length' value

NORM = [] 
SUM= []
Img_pix =[]
plt.figure()

for n in range(n_set):
    pos =  [my_file['Caption']=='hyPIRFwd' for my_file in Data[n]['files']]
    hyPIRFwd = np.array(Data[n]['files'])[pos][0]
    #Artifact removing (first lines, because of scanning)
    hyPIRFwd['data'] = hyPIRFwd['data'][1:31, 1:32, :]
    Img_pix.append(hyPIRFwd['data'][:, :, 0].shape)
    data = np.reshape(hyPIRFwd['data'], (hyPIRFwd['data'].shape[0]*hyPIRFwd['data'].shape[1], hyPIRFwd['data'].shape[2]))
    
    
    #Calibration File for each dataset
    Cali_data = np.array(Calib_file[n])
    Cali_Spc = Cali_data[:, 1]
    
    data = data / Cali_Spc
    my_sum = np.sum(data, axis = 1)
    SUM.append(len(my_sum))
    # zeros were excluded prior normalization and analysis
    data = data[my_sum != 0]
    my_sum = my_sum[my_sum != 0]
        
    spc_norm = np.array([(spc)/s for spc, s in zip(data, my_sum)])
    NORM.append(spc_norm)    
    mean_spc = np.mean(spc_norm, axis = 0)
    plt.plot(Data[n]['wavelength'], 
             signal.savgol_filter(mean_spc, window_length=11, polyorder=2, mode="nearest"),
             label = 'Data'+str(n+1))

plt.legend()
plt.xlabel('Wavenumber (cm\u207b1)')
plt.ylabel('intensity (normalized)')
plt.title('mean spectrum')
plt.tight_layout()
plt.gca().invert_xaxis()   
    
#%% Merging data
# no. of data sets to be included in PCA (ordered from the read data)

set_PCA = 2
NORMSpc = NORM[0]
shape = [NORMSpc.shape]

if set_PCA == 1:
    NORMSpc = NORMSpc
else:
    for a in range(set_PCA-1):
#        sh = NORMSpc.shape
#        shape.append(sh)
        NORMSpc = np.concatenate((NORMSpc, NORM[a+1]), axis = 0)
        sh = NORMSpc.shape
        shape.append(sh)
        print(NORMSpc.shape)

coord = np.arange(len(NORMSpc))
MEANSpc = np.mean(NORMSpc, axis=0)

#%%

# And Apply PCA
# ncomp refers to the desired numbers of principal components

from sklearn.decomposition import PCA
ncomp = 3
model = PCA(n_components=ncomp)

transformed_data = model.fit(NORMSpc - MEANSpc).transform(NORMSpc - MEANSpc).T

     


#%% Loadings:
# Plotted after smoothing by savitzky-golay filter (2-5-5) 
# can be changed by 'windows_length' value

loadings = model.components_

my_fig = plt.figure()
ax = plt.subplot(111)
plt.gca().invert_xaxis() #inverts values of x-axis
for icomp in range(ncomp):
    ax.plot(Data[0]['wavelength'], 
            signal.savgol_filter(loadings[icomp], window_length=11, polyorder=2, mode="nearest"),
            label='PC'+str(icomp+1) )
    ax.legend()
ax.set_xlabel('Wavenumber (cm\u207b1)')
ax.set_ylabel('intensity (normalized)')
plt.title('PCA-Loadings')


#%% Scatter plots (clustered):

limits = [0]

for j in range(set_PCA):    
    limit = shape[j][0]
    limits.append(limit)


for sc in range(ncomp-1):
    my_fig2 = plt.figure()
    ax2 = plt.subplot(111)
    for j in range(set_PCA):
        ax2.plot(transformed_data[sc, limits[j]:limits[j+1]], transformed_data[sc+1, limits[j]:limits[j+1]], '.', label = ('Data ' + str(j+1)))        
    ax2.legend()
    ax2.set_xlabel('PC'+str(sc+1))
    ax2.set_ylabel('PC'+str(sc+2))
    plt.title('scatterplot PC' + str(sc+1) +'& PC'+ str(sc+2))
    my_fig.tight_layout()
   
    if sc == 1:
        my_fig3 = plt.figure()
        ax3 = plt.subplot(111)
        for j in range(set_PCA):
            ax3.plot(transformed_data[sc-1, limits[j]:limits[j+1]], transformed_data[sc+1, limits[j]:limits[j+1]], '.', label = ('Data ' + str(j+1)))        
        ax3.set_xlabel('PC'+str(sc))
        ax3.set_ylabel('PC'+str(sc+2))
        ax3.legend()
        plt.title('scatterplot PC' + str(sc) +'& PC'+ str(sc+2))
        my_fig.tight_layout()


#%% Factor plots:
# Overview of all PC factor maps [for all data sets in one plot]

total_sum = 0
for t in range (set_PCA):
    total_sum = total_sum + SUM[t]
 

zeros = np.zeros((total_sum))    
maps = [zeros.copy() for icomp in range(ncomp)]


xpix = Img_pix[0][0]
ypix = Img_pix[0][1]
Img = []

if set_PCA > 1:
    for p in range(set_PCA-1):
        xpix = xpix + Img_pix[p][0]
        Img.append(Img_pix[p])
        
b = 0


for icomp in range(ncomp):   
    maps[icomp][coord] = transformed_data[icomp]
    maps[icomp] = np.reshape(maps[icomp],(xpix, ypix))
    
    my_fig = plt.figure()
    ax = plt.subplot(111)
    #vmin and vmax were chosen for better visualization of our data sets (can be changed)
    plt.colorbar( ax.imshow(maps[icomp], cmap = 'coolwarm', vmin = -0.005, vmax = 0.005))
    ax.set_xlabel('x scan ['+Data[0]['XPhysUnit']+']')
    ax.set_ylabel('y scan ['+Data[0]['XPhysUnit']+']')
    ax.legend_ = None
    plt.title('factors PC  '+str(icomp+1))
    my_fig.tight_layout()

#%%

#Individual plots per each dataset

pix = Img_pix[0][0]

if set_PCA > 1:   
    for i in range(set_PCA):
        for icomp in range(ncomp): 
            maps_ind = maps[icomp][b:pix, :]       
            my_fig = plt.figure()
            ax = plt.subplot(111)
            plt.colorbar( ax.imshow(maps_ind, cmap = 'coolwarm', vmin = -0.005, vmax = 0.005))
            ax.set_xlabel('x scan ['+Data[0]['XPhysUnit']+']')
            ax.set_ylabel('y scan ['+Data[0]['XPhysUnit']+']')
            ax.legend_ = None
            plt.title('Data '+ str(i+1)+'PC'+str(icomp+1))
            my_fig.tight_layout()


#%% Heirarchical Clustering Analysis HCA

# Create a dendrogram for the ccombined data
# by calculating the distances between the spectra using Eclidean distance & ward method 

clustering = linkage(NORMSpc, method = "ward", metric="euclidean")
plt.figure()
dendrogram(clustering, color_threshold=0.07, above_threshold_color='grey')
plt.title('Dendrogram of Clusters')
plt.xlabel('Spectra')
plt.ylabel('Distances')
# Optional: plotting a line at where the user can estimate number of desired clusters 
plt.axhline(y = 0.1, color = 'black', linestyle = 'solid'  )
plt.show()

#%% Agglomerative clustering
# Here, the same steps are done to calculate the distances and generate clusters
# but the algorithm stops at a chosen number of clusters

# based on the line the user plot abive, number of desired clusters should be chosen here
# supervised clustering

n_clusters = 3

clustering_model = AgglomerativeClustering(n_clusters= n_clusters, metric='euclidean')
clustering_model.fit_predict(NORMSpc)

data_labels = clustering_model.labels_

#to check number of spectra per cluster
for i in range(n_clusters):
    number_spc = NORMSpc[data_labels == i].shape
    print('Cluster ' + str(i+1) + '  has:  ' + str(number_spc[0]))
    
#%% plotting the mean spectra for each cluster


hca_means = []
hca_means_sm = [] 


plt.figure()
for i in range (n_clusters):
    cluster_mean = np.mean(NORMSpc[data_labels == i], axis = 0)
    cluster_mean_sm = signal.savgol_filter((cluster_mean), window_length=11, polyorder=2, mode="nearest")
    plt.plot(Data[0]['wavelength'], cluster_mean_sm.T, label = 'cluster' + str(i+1))
    plt.legend()
    hca_means.append(cluster_mean)
    hca_means_sm.append(cluster_mean_sm)
    
plt.gca().invert_xaxis()
plt.title(' HCA Mean Spectra')
plt.xlabel('Wavenumber (cm\u207b1)')
plt.ylabel('intensity (normalized)')
plt.show()

#Exporting mean spectra
#np.savetxt('HCA_means_smoothed.txt', hca_means_sm, delimiter='\t')

#%% cluster maps
# created using non negative factorization NNF

H, W, n =  non_negative_factorization(np.abs(NORMSpc), n_components = n_clusters, init='custom',
                                      update_H=False, H=np.abs(hca_means), max_iter = 10000)

pix = Img_pix[0][0]
b = 0
res = H.__copy__().T
maps_hca = []

for ic, r in enumerate(res):
    my_fig= plt.figure()
    ax = plt.subplot(111)
    plt.colorbar(ax.imshow(np.reshape(r, (xpix, ypix)), cmap = 'gray', vmin = 0, vmax = 1))  
    plt.title(' cluster ' + str(ic+1) + ' map')
    plt.yticks([])
    plt.xticks([])
    plt.tight_layout()
    maps_hca.append(np.reshape(r, (xpix, ypix)))



#Individual plots per each dataset

pix = Img_pix[0][0]
b = 0

if set_PCA > 1:   
    maps_hca = np.array(maps_hca)
    for d in range(set_PCA):
        for n_c in range(n_clusters):            
            my_fig = plt.figure()
            ax = plt.subplot(111)
            plt.colorbar(ax.imshow(maps_hca[n_c, b:pix, :], cmap = 'gray', vmin = 0, vmax = 1))
            ax.set_xlabel('x scan ['+Data[0]['XPhysUnit']+']')
            ax.set_ylabel('y scan ['+Data[0]['XPhysUnit']+']')
            ax.legend_ = None
            plt.title('Data '+ str(d+1)+' cluster ' + str(n_c+1))
            my_fig.tight_layout()
        
        
        b = b + Img_pix[d][0]
        pix = pix + Img_pix[d][0]
            

        
        b = b + Img_pix[i][0]
        pix = pix + Img_pix[i][0]

        
