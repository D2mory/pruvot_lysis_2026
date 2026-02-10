# README

### Overview
Code for: *Connecting models, networks and experiments: Revisiting the role of viruses in marine carbon cycling*
Authors: Maxine Pruvôt, Bart Haegeman, Fabien Joux, David Demory.

All codes are written in MATLAB (version R2011a or later). To run the analyses or generate figures, ensure your current working directory is the folder containing the code and data files.

### Folder descriptions
The codes are organised according to the sections of the article:
- **section4_model_comparison**: codes for *Common ground to compare modelling approaches* (Section 4).
- **section5_virus_removal**: codes for *Virus removal as a tool to quantify their role* (Section 5).

Each section folder contains the scripts used to generate the main text and supplementary figures.

### Section 4
**Data:**
- The data used in section 4 is available in the file `data_model_comparison.mat`, it contains a single variable `res` :
	- `res(1)`: WIL = data from Wilhelm & Suttle (1999)
	- `res(2)`: WEI = data from Weitz _et al._ (2015)
	- `res(3)`: FUH = data from Fuhrman (2 datasets, 500 copies each), in order: Fuhrman 1992 and Furhman 1999
	- `res(4)`: MOJ = data from Mojica & Brussaard (2015) (2 datasets, 500 copies each), in order: MOJ_S (southern region) and MOJ_N (northern region)
	- `res(5)`: XIE = data from Xie _et al._ (2022) (2 datasets, 500 copies each), in order: XIE_H (HOT sampling site) and XIE_A (Arabian Sea sampling site)
	
	- `res.flx` are the flux matrices used for the analysis
	- `res.fmn` are average flux matrices of the 1000 realisations for the models of WIL and WEI
	
- If you want to generate the data yourself, run `data_generation_model_comparison.m`. It will generate the data file presented above. This code extracts the flux matrix for each of the selected models, and compile and save them in `data_model_comparison.mat`.

**Figures:**
-  `fig_3_PCA.m`: PCA and plotting of the figure 3
-  `fig_S1.m`: plotting of the figure S1

**Model functions:** used to generate the data from each model presented in the article
- `FUH.m`: Fuhrman (1992, 1999). `FUH('B')` for Fuhrman (1992) and `FUH('PB')` for Fuhrman (1999)
- `MOJ.m`: Mojica (2015). `MOJ('S')` for the southern region, and `MOJ('N')` for the northern region
- `XIE.m`: Xie _et al._ (2022). `XIE('H')` for the HOT site, and `XIE('A')` for the Arabian Sea site
- `WIL.m`: Wilhelm & Suttle (1999)
- `WEI.m`: Weitz _et al._ (2015), with its associated functions starting with "wei_"
	- `wei_eqs.m`: function containing the equilibrium expressions for the model of Weitz _et al._ (2015)
	- `wei_rhs.m`: function computing the dynamical equations from Weitz _et al._ (2015), and generating the flux matrix of its model
	- `wei_smp.m`: function sampling uniformly parameters values within ranges

### Section 5
**Data:**
- The data used in section 5 is available in the file `data_virus_removal.mat`.
- If you want to generate the data yourself, run `data_generation_virus_removal.m`. It will generate the data file presented above. Running this script should take several minutes. This code: 
	- Generates 10000 parameters sets leading to realisations realisations from a simplified dynamical model
	- Removes viruses from the model and then simulates it without viruses
	- Extracts and saves the fluxes, parameters sets, and equilibrium of the simulation that reached equilibrium

**Functions:**
- `gauss_trf.m`: function of gaussianization
- `hexbin_plot.m`: function for visualisation of the figures 4, 5, and S2: hexagonal binning and plot
- `vr_chk.m`: function insuring that all variables steady states are within desired (realistic) ranges
- `vr_eqs.m`: function calculating the steady state values for the different variables
- `vr_flx_trf.m`: function that transforms the matrix into the common format
- `vr_rhsC.m`: function that gives the dynamical equations and flux matrix in carbon units
- `vr_smp.m`: function sampling uniformly parameters values within ranges
- `vr_trf.m`:  function transforming variables from abundance to carbon variables or the other way around

**Figures:**
- `fig_4.m`: figure 4
- `fig_5.m`: figure 5 
- `fig_S2.m`: figure S2
- `fig_S3.m`: figure S3

### Discussion
**Data:** The data presented here is used only in the figure S4, they are compiled viral lysis rates from Mojica & Brussaard (2026):
- `MDA_data.xlsx` : data from the Modified Dilution Assay
- `VP_data.xlsx` : data from the Virus Reduction Assay

**Figures:** 
- `fig_S4.m`: figure S4


### References
- Fuhrman, J.A. (1992). Bacterioplankton roles in cycling of organic matter: the microbial food web. In: _Primary Productivity and Biogeochemical Cycles in the Sea_. Springer, pp. 361–383. https://doi.org/10.1007/978-1-4899-0762-2_20
- Fuhrman, J.A. (1999). Marine viruses and their biogeochemical and ecological effects. _Nature_, 399, 541–548. https://doi.org/10.1038/21119.
- Mojica, K.D. (2015). The viral shunt in a stratified Northeast Atlantic Ocean. _Viral Lysis of Marine Microbes in Relation to Vertical Stratification_. PhD thesis, University of Amsterdam, chap. 7, pp. 207–222.
- Weitz, J.S., Stock, C.A., Wilhelm, S.W., Bourouiba, L., Coleman, M.L., Buchan, A., Follows, M.J., Fuhrman, J.A., Jover, L.F., Lennon, J.T. et al. (2015). A multitrophic model to quantify the effects of marine viruses on microbial food webs and ecosystem processes. _ISME Journal_, 9, 1352–1364. https://doi.org/10.1038/ismej.2014.220
- Wilhelm, S.W. & Suttle, C.A. (1999). Viruses and nutrient cycles in the sea: viruses play critical roles in the structure and function of aquatic food webs. _BioScience_, 49, 781–788. https://doi.org/10.2307/1313569.
- Xie, L., Zhang, R. & Luo, Y.W. (2022). Assessment of explicit representation of dynamic viral processes in regional marine ecological models. _Viruses_, 14, 1448. https://doi.org/10.3390/v14071448.

