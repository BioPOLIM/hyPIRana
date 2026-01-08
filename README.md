# hyPIRana
Analysis software development for analysis of hyperspectral data in mid-IR photo-induced force microscopy (mid-IR PiFM, PiF-IR) 
-- MonIRana analyzes data sets of single PiF-IR spectra --


# hyPIRana v2.1.0 release adding new features for cluster analysis to hyPIRana

version v2.1.0 date 2025-01-08

Contributions to this update of hypIRana have been made by Maryam Ali
For details see headers of the code.

The new release of hyPIRana includes an approach for Heirarchical Clustering Analysis (HCA), doing:
- Generating a HCA dendrogram between the hyperspectra forming clusters.
- Supervised Agglomerative clustering based on user's desired number of spectra (referring to created dendrogram) 
- Plot the mean spectra for each cluster.
- generate the cluster maps by non negative factorization algorithm.

# hyPIRana v2.0.0 release extending hyPIRana to the combined analysis of multiple hyperspectral data sets 
version v2.0.0 date 2025-05-14

Contributions to this update of hypIRana have been made by Maryam Ali, Sebastian Unger, and Daniela Täuber. 
For details see headers of the code.

This new release of hypIRana can do:
- Calibration 
- Calculate mean spectra of the individual hyperspectral data sets
- Run a combined PCA on the data sets and provide plots of results
- Clustering of data sets by color in combined PCA scatter plots

Two additional hyperspectral data sets are available in this update of hyPIRana.py:
- Bacillus Subtilis incubated with Vancomycin for 30 minutes
- Bacillus Subtilis incubated with Vancomycin for 60 minutes


# hyPIRana v1.1.0 release including newly added MonIRana
version v1.1.0 date 2025-04-08

MonIRana is a program for analysis of a set of spectral data provided in a tabulated textfile
This code has been developed to analyze single spectral dataset from mid-IR Photo-induced Force Microscopy (PiF-IR) Contributions in this code so far have been made by: Mohammad Soltaninezhad, Sebastian Unger, Maryam Ali, and Daniela Täuber.

monIRana can do:
 - Calibration
 - Calculate mean spectra of the complete data set
 - Run a PCA on the data set and provide plots of results
 - Clustered presentation of Amide rich & glycan rich analyzed spectra

The available raw data in this repository are for:
 - Calibrated PiF-IR spectra for treated Bacillus Subtilis (classified to Amide rich & glycan rich spectra)
 - Calibrated PiF-IR spectra for untreated Bacillus Subtilis (control) from the same experiment


# hyPIRana v1.0.1 initial release MIT License
added the MIT license in the release
version v1.0.1 date 2024-10-15

# hyPIRana v1.0.0 initial release

see also: [hyPIRana on ungersebastian](https://github.com/ungersebastian/hyPIRana)

This sofware is part of hyperspectral analysis methods used in the published article: "Nanoscale chemical characterization of secondary protein structure of F-Actin using mid-infrared photoinduced force microscopy (PiF-IR)" by Jesvin Joseph, Lukas Spantzel, Maryam Ali, Dijo Moonnukandathil Joseph, Sebastian Unger, Katharina Reglinski, Christoph Krafft, Anne-Dorothea Müller, Christian Eggeling, Rainer Heintzmann, Michael Börsch, Adrian T. Press, Daniela Täuber. Spectrochimica Acta part A: Molecular and Biomolecular Spectroscopy, 306, 123612, 2024. https://doi.org/10.1016/j.saa.2023.123612

Contributions in this code so far have been made by Sebastian Unger, Maryam Ali, René Lachmann, Rainer Heintzmann, Mohammad Soltaninezhad and Daniela Täuber. For details see headers of the code.

This code analyzes the hyperspectral data by:
- Plotting generated images: AFM Topography, integrated PiF-IR hyperspectral image, and AFM phase image.
- Plotting the mean spectrum on PiF-IR intensity.
- Principal Component Analysis PCA represented by loaded components plot, scatter plot, and individual component factor plots.

In the context of the article by Joseph et al. Spectrochimica Acta part A: Molecular and Biomolecular Sepctroscopy, 306, 123612, 2024, two PiF-IR hyperspectral data sets have been analysed using hyPIRana. The Raw data have been available in this repository:
- single-fibrillar F-Actin
- crosslinked F-Actin

version v1.0.0 date 2024-10-13
