# MATLAB-Research: 
Snippets of MATLAB codes written for my research projects.


## COMSOL_simulation 
This contains disorder_slab.m file that interact with COMSOL to generate complex slab structures at random as pictured below, the extreme left is the entire structure, the rest is the zoom-in view of the central area of the structure, the red line is where measurement should be taken.

<p align="center">
  <img width="721" height="434" src="https://github.com/luoqiaoen/MATLAB-Research/blob/master/COMSOL_simulation/simulated_structure.png">
</p>

TM.m will compute transmission matrix of the structure based on input and ouput fields, and RMT.m will analyze the transmission matrix with eigen-decomposition (PCA) following random matrix theory.

## speckle_corr_sample
This contains sample code to process intensity correlation data for image reconstruction, will need phase retreival algorithm.

The sample and reconstruction can be like below:
<p align="center">
  <img width="500" height="200" src="https://github.com/luoqiaoen/MATLAB-Research/blob/master/speckle_corr_sample/recon.png">
</p>
