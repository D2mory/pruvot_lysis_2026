# README

### Overview

Code for : *Connecting models, networks and experiments: Revisiting the role of viruses in marine carbon cycling*.  
Authors : M. Pruvôt, B. Haegeman, F. Joux, D. Demory.[![DOI](https://zenodo.org/badge/1154573909.svg)](https://doi.org/10.5281/zenodo.20813474)

The code is written in MATLAB. To run the analyses or generate the figures, ensure your current working directory is the folder containing code and data files.

The code is organised into two folders :

  -  **modelcomparison** : code for Figs. 3 and S1.
  -  **virusremoval** : code for Figs. 4, 5, S2, S3 and S4.


### Folder : modelcomparison

**Data:**

The data is available in the file `modelcomp_data.mat`. It contains a single variable `res` :

  - `res(1)`: data for models WIL, from Wilhelm & Suttle (1999).
  - `res(2)`: data for models WEI, from Weitz _et al._ (2015).
  - `res(3)`: data for models FUH, from Fuhrman (1992, 1999).
  - `res(4)`: data for models MOJ, from Mojica (2015) (two datasets : MOJ-S, southern region; MOJ-N, northern region).
  - `res(5)`: data for models XIE, from Xie _et al._ (2022) (two datasets : XIE-H, HOT sampling site; XIE-A, Arabian Sea sampling site).

The fields of the variable `res` :

  - `res.flx` : array of individual flux matrices (1000 realisations for models WIL and WEI).
  - `res.fmn` : average flux matrix (average over 1000 realisations for models WIL and WEI).
	
If you want to generate the _in silico_ data yourself, run `modelcomp_data.m`.

**Functions:**

  - `modelcomp_data.m` : generates the data for Figs. 3 and S1.
  - `modelcomp_fig3.m` : plots Fig. 3 of the paper.
  - `modelcomp_figS1.m` : plots Fig. S1 of the paper.
  - `model_WIL_flx.m` : model WIL; generates 1000 realisations.
  - `model_WEI_flx.m` : model WEI; generates 1000 realisations.
  - `model_WEI_eqs.m` : model WEI; computes the equilibrium state.
  - `model_WEI_rhs.m` : model WEI; evaluates the dynamical equations.
  - `model_WEI_smp.m` : model WEI; samples parameters values within specified ranges.
  - `model_FUH_flx.m` : models FUH; `model_FUH('B')` for Fuhrman (1992) and `model_FUH('BP')` for Fuhrman (1999).
  - `model_MOJ_flx.m` : models MOJ; `model_MOJ('S')` for the southern region and `model_MOJ('N')` for the northern region.
  - `model_XIE_flx.m` : models XIE; `model_XIE('H')` for the HOT site and `model_XIE('A')` for the Arabian Sea site.


### Folder : virusremoval

**Data:**

The data is available in the file `virusremoval_data.mat`. It contains 10000 parameterisations of the minimal marine ecosystem model, along with their corresponding equilibrium state and flux network representation. For most parameterisations, the model reaches an equilibrium state after virus removal. For these cases, the dataset also contains the equilibrium state and flux network representation after virus removal.

If you want to generate the _in silico_ data yourself, run `virusremoval_data.m`.

**Functions:**

  - `gauss_trf.m` : implements the gaussianisation transform.
  - `hexbin_plot.m` : plots a two-dimensional histogram using hexagonal binning.
  - `virusremoval_data.m` : generates the data for Figs. 4, 5, S2, S3 and S4.
  - `virusremoval_fig4.m` : plots Fig. 4 of the paper.
  - `virusremoval_fig5.m` : plots Fig. 5 of the paper.
  - `virusremoval_figS2.m` : plots Fig. S2 of the paper.
  - `virusremoval_figS3.m` : plots Fig. S3 of the paper.
  - `virusremoval_figS3.m` : plots Fig. S4 of the paper.
  - `virusremoval_model_chk.m` : checks whether a model parameterisation is consistent and realistic.
  - `virusremoval_model_eqs.m` : computes the equilibrium state.
  - `virusremoval_model_rhs.m` : evaluates the dynamical equations.
  - `virusremoval_model_rhsC.m` : evaluates the dynamical equations with variables expressed in carbon units.
  - `virusremoval_model_smp.m` : samples parameter values within specified ranges.
  - `virusremoval_model_std.m` : transforms flux matrices obtained from the model into standard format.


### References

  - Fuhrman, J.A. (1992). Bacterioplankton roles in cycling of organic matter: the microbial food web. In: _Primary Productivity and Biogeochemical Cycles in the Sea_. Springer, pp. 361–383. https://doi.org/10.1007/978-1-4899-0762-2_20
  - Fuhrman, J.A. (1999). Marine viruses and their biogeochemical and ecological effects. _Nature_, 399, 541–548. https://doi.org/10.1038/21119.
  - Mojica, K.D. (2015). The viral shunt in a stratified Northeast Atlantic Ocean. _Viral Lysis of Marine Microbes in Relation to Vertical Stratification_. PhD thesis, University of Amsterdam, chap. 7, pp. 207–222.
  - Weitz, J.S., Stock, C.A., Wilhelm, S.W., Bourouiba, L., Coleman, M.L., Buchan, A., Follows, M.J., Fuhrman, J.A., Jover, L.F., Lennon, J.T. et al. (2015). A multitrophic model to quantify the effects of marine viruses on microbial food webs and ecosystem processes. _ISME Journal_, 9, 1352–1364. https://doi.org/10.1038/ismej.2014.220
  - Wilhelm, S.W. & Suttle, C.A. (1999). Viruses and nutrient cycles in the sea: viruses play critical roles in the structure and function of aquatic food webs. _BioScience_, 49, 781–788. https://doi.org/10.2307/1313569.
  - Xie, L., Zhang, R. & Luo, Y.W. (2022). Assessment of explicit representation of dynamic viral processes in regional marine ecological models. _Viruses_, 14, 1448. https://doi.org/10.3390/v14071448.

